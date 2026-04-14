#!/usr/bin/env python3
"""Filter BLAST HitTable + FullSeq by numeric and RecName conditions, then export FASTA."""

from __future__ import annotations

import argparse
import operator
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

import yaml


OP_MAP = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
}


@dataclass
class Config:
    hit_table_files: List[str]
    fullseq_files: List[str]
    numeric_conditions: List[Dict[str, Any]]
    recname_filters: List[Dict[str, Any]]
    fasta_file: str


class ConfigError(Exception):
    pass


class ParseError(Exception):
    pass


def load_config(config_path: Path) -> Config:
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")

    try:
        data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"Failed to parse YAML: {config_path}") from exc

    if not isinstance(data, dict):
        raise ConfigError("YAML root must be a mapping")

    try:
        input_block = data["input"]
        filters_block = data["filters"]
        output_block = data["output"]
    except KeyError as exc:
        raise ConfigError(f"Missing required top-level key: {exc.args[0]}") from exc

    hit_table_files = input_block.get("hit_table_files", [])
    fullseq_files = input_block.get("fullseq_files", [])
    numeric_conditions = filters_block.get("numeric_conditions", [])
    recname_filters = data.get("recname_filters", [])
    fasta_file = output_block.get("fasta_file")

    if not hit_table_files:
        raise ConfigError("input.hit_table_files must contain at least one file")
    if not fullseq_files:
        raise ConfigError("input.fullseq_files must contain at least one file")
    if not fasta_file:
        raise ConfigError("output.fasta_file is required")

    return Config(
        hit_table_files=hit_table_files,
        fullseq_files=fullseq_files,
        numeric_conditions=numeric_conditions,
        recname_filters=recname_filters,
        fasta_file=fasta_file,
    )


def normalize_field_name(name: str) -> str:
    text = name.strip().lower()
    replacements = {
        "subject acc.ver": "subject",
        "% identity": "identity",
        "% positives": "positives",
        "query acc.ver": "query",
        "q. start": "q_start",
        "q. end": "q_end",
        "s. start": "s_start",
        "s. end": "s_end",
        "alignment length": "alignment_length",
        "gap opens": "gap_opens",
        "bit score": "bit_score",
    }
    if text in replacements:
        return replacements[text]
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text


def parse_hit_table_file(path: Path) -> tuple[Dict[str, List[Dict[str, Any]]], List[str]]:
    if not path.exists():
        raise ParseError(f"HitTable file not found: {path}")

    records_by_subject: Dict[str, List[Dict[str, Any]]] = {}
    columns: List[str] | None = None

    with path.open("r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                if line.startswith("# Fields:"):
                    raw_fields = [x.strip() for x in line.split(":", 1)[1].split(",")]
                    columns = [normalize_field_name(x) for x in raw_fields]
                continue

            parts = line.split("\t")
            if columns is None:
                raise ParseError(
                    f"HitTable header '# Fields:' not found before data in {path} (line {line_no})"
                )
            if len(parts) != len(columns):
                raise ParseError(
                    f"Column count mismatch in {path} line {line_no}: expected {len(columns)}, got {len(parts)}"
                )

            record = dict(zip(columns, parts))
            subject = record.get("subject")
            if not subject:
                raise ParseError(f"Missing 'subject' column in {path} line {line_no}")

            records_by_subject.setdefault(subject, []).append(record)

    if columns is None:
        raise ParseError(f"Could not find '# Fields:' header in {path}")

    required = {"subject", "identity", "positives"}
    missing = [x for x in required if x not in columns]
    if missing:
        raise ParseError(f"Required columns missing in {path}: {', '.join(missing)}")

    return records_by_subject, columns


def evaluate_numeric_conditions(
    records_by_subject: Dict[str, List[Dict[str, Any]]],
    conditions: List[Dict[str, Any]],
) -> tuple[Dict[str, List[Dict[str, Any]]], List[str]]:
    if not conditions:
        return records_by_subject, []

    errors: List[str] = []
    filtered: Dict[str, List[Dict[str, Any]]] = {}

    for subject, rows in records_by_subject.items():
        for row in rows:
            ok = True
            for cond in conditions:
                column = cond.get("column")
                op = cond.get("op")
                value = cond.get("value")

                if column is None or op is None or value is None:
                    raise ConfigError(f"Invalid numeric condition: {cond}")
                if op not in OP_MAP:
                    raise ConfigError(f"Unsupported operator '{op}' in condition: {cond}")
                if column not in row:
                    raise ParseError(f"Column '{column}' not found in HitTable record")

                try:
                    left = float(row[column])
                    right = float(value)
                except (TypeError, ValueError):
                    ok = False
                    errors.append(
                        f"Excluded row (subject={subject}): non-numeric value in column '{column}' -> '{row.get(column)}'"
                    )
                    break

                if not OP_MAP[op](left, right):
                    ok = False
                    break

            if ok:
                filtered.setdefault(subject, []).append(row)

    return filtered, errors


def parse_fasta_header(header: str) -> tuple[str, str]:
    header = header.strip()
    if not header:
        raise ParseError("Empty FASTA header")

    subject = header.split()[0]
    recname = ""

    m = re.search(r"RecName:\s*Full=([^;]+)", header)
    if m:
        recname = m.group(1).strip()
    else:
        rest = header[len(subject) :].strip()
        if rest:
            recname = re.sub(r"\s*\[[^\]]+\]\s*$", "", rest).strip()

    if not recname:
        recname = "UNKNOWN"

    return subject, recname


def parse_fullseq_file(path: Path) -> Dict[str, List[Dict[str, str]]]:
    if not path.exists():
        raise ParseError(f"FullSeq file not found: {path}")

    seq_dict: Dict[str, List[Dict[str, str]]] = {}
    current_header: str | None = None
    seq_chunks: List[str] = []

    def flush() -> None:
        nonlocal current_header, seq_chunks
        if current_header is None:
            return
        subject, recname = parse_fasta_header(current_header)
        sequence = "".join(seq_chunks).replace(" ", "")
        if not sequence:
            raise ParseError(f"No sequence found for FASTA header '{current_header}' in {path}")
        seq_dict.setdefault(subject, []).append(
            {"subject": subject, "RecName": recname, "Seq": sequence}
        )

    with path.open("r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_header = line[1:].strip()
                seq_chunks = []
            else:
                if current_header is None:
                    raise ParseError(f"Sequence line before FASTA header in {path} line {line_no}")
                seq_chunks.append(line)

    flush()
    return seq_dict


def match_recname(recname: str, cond: Dict[str, Any]) -> bool:
    keyword = cond.get("keyword")
    match_mode = cond.get("match")
    case_sensitive = bool(cond.get("case_sensitive", False))

    if keyword is None or match_mode is None:
        raise ConfigError(f"Invalid recname_filters entry: {cond}")

    if not case_sensitive:
        recname_cmp = recname.lower()
        keyword_cmp = str(keyword).lower()
    else:
        recname_cmp = recname
        keyword_cmp = str(keyword)

    if match_mode == "exact":
        return recname_cmp == keyword_cmp
    if match_mode == "partial":
        return keyword_cmp in recname_cmp
    if match_mode == "startswith":
        return recname_cmp.startswith(keyword_cmp)
    if match_mode == "regex":
        flags = 0 if case_sensitive else re.IGNORECASE
        return re.search(str(keyword), recname, flags=flags) is not None

    raise ConfigError(f"Unknown recname match mode: {match_mode}")


def filter_sequences_by_subject(
    filtered_hits: Dict[str, List[Dict[str, Any]]],
    seq_dict: Dict[str, List[Dict[str, str]]],
) -> Dict[str, List[Dict[str, str]]]:
    result: Dict[str, List[Dict[str, str]]] = {}
    for subject in filtered_hits:
        if subject in seq_dict:
            result[subject] = seq_dict[subject]
    return result


def filter_by_recname(
    subject_matched: Dict[str, List[Dict[str, str]]],
    recname_filters: List[Dict[str, Any]],
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for entries in subject_matched.values():
        for item in entries:
            if not recname_filters:
                rows.append(item)
                continue
            if any(match_recname(item["RecName"], cond) for cond in recname_filters):
                rows.append(item)
    return rows


def write_fasta(entries: List[Dict[str, str]], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8") as fw:
        for item in entries:
            fw.write(f">{item['subject']} | RecName={item['RecName']}\n")
            seq = item["Seq"]
            for i in range(0, len(seq), 80):
                fw.write(seq[i : i + 80] + "\n")


def merge_subject_dicts(dicts: List[Dict[str, List[Dict[str, Any]]]]) -> Dict[str, List[Dict[str, Any]]]:
    merged: Dict[str, List[Dict[str, Any]]] = {}
    for d in dicts:
        for k, v in d.items():
            merged.setdefault(k, []).extend(v)
    return merged


def run(config_path: Path) -> None:
    cfg = load_config(config_path)

    hit_dicts = [parse_hit_table_file(Path(p))[0] for p in cfg.hit_table_files]
    merged_hits = merge_subject_dicts(hit_dicts)
    filtered_hits, numeric_errors = evaluate_numeric_conditions(merged_hits, cfg.numeric_conditions)

    seq_dicts = [parse_fullseq_file(Path(p)) for p in cfg.fullseq_files]
    merged_seqs = merge_subject_dicts(seq_dicts)

    subject_matched = filter_sequences_by_subject(filtered_hits, merged_seqs)
    final_entries = filter_by_recname(subject_matched, cfg.recname_filters)

    write_fasta(final_entries, Path(cfg.fasta_file))

    for msg in numeric_errors:
        print(f"[WARN] {msg}", file=sys.stderr)

    print(f"[INFO] Hit subjects (all): {len(merged_hits)}")
    print(f"[INFO] Hit subjects (after numeric filter): {len(filtered_hits)}")
    print(f"[INFO] Subject-matched sequence groups: {len(subject_matched)}")
    print(f"[INFO] Output entries: {len(final_entries)}")
    print(f"[INFO] FASTA written: {cfg.fasta_file}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter HitTable + FullSeq and output FASTA")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml)",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    try:
        run(Path(args.config))
    except (ConfigError, ParseError, FileNotFoundError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"[ERROR] Unexpected failure: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
