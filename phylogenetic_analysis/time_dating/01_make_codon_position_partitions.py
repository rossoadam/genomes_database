#!/usr/bin/env python3
"""
Convert an IQ-TREE/RAxML-style CDS partition file into three codon-position
partitions while resetting codon phase at the start of every CDS.

Example input
-------------
DNA, gene1.fasta = 1-1734
DNA, gene2.fasta = 1735-3090

Example output
--------------
DNA, codon_pos1 = 1-1734\3,1735-3090\3
DNA, codon_pos2 = 2-1734\3,1736-3090\3
DNA, codon_pos3 = 3-1734\3,1737-3090\3

The script validates:
  * input syntax;
  * positive, ordered coordinates;
  * CDS lengths divisible by three;
  * overlaps;
  * gaps, optionally;
  * alignment length, when a PHYLIP alignment is supplied.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PARTITION_RE = re.compile(
    r"^\s*(?P<datatype>[^,]+?)\s*,\s*"
    r"(?P<name>.+?)\s*=\s*"
    r"(?P<start>\d+)\s*-\s*(?P<end>\d+)\s*$"
)


@dataclass(frozen=True)
class Partition:
    line_number: int
    datatype: str
    name: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert an IQ-TREE/RAxML-style CDS partition file into three "
            "locus-aware codon-position partitions."
        )
    )
    parser.add_argument(
        "input_partition",
        type=Path,
        help="Original IQ-TREE/RAxML-style partition file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "Output partition path. Default: "
            "<input_stem>.codon_positions.partitions"
        ),
    )
    parser.add_argument(
        "--format",
        choices=("raxml", "nexus"),
        default="raxml",
        help="Output syntax. Default: raxml.",
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        help=(
            "Optional PHYLIP alignment. Its header site count will be checked "
            "against the final partition coordinate."
        ),
    )
    parser.add_argument(
        "--report",
        type=Path,
        help=(
            "Optional TSV validation report. Default: "
            "<output>.validation.tsv"
        ),
    )
    parser.add_argument(
        "--allow-gaps",
        action="store_true",
        help=(
            "Allow uncovered coordinates between CDS partitions. By default, "
            "gaps are treated as an error."
        ),
    )
    parser.add_argument(
        "--allow-nontriplet",
        action="store_true",
        help=(
            "Allow CDS ranges whose lengths are not divisible by three. "
            "This is normally unsafe and is disabled by default."
        ),
    )
    parser.add_argument(
        "--datatype",
        default="DNA",
        help="Datatype label for RAxML-style output. Default: DNA.",
    )
    return parser.parse_args()


def default_output_path(input_path: Path, output_format: str) -> Path:
    suffix = ".codon_positions.nex" if output_format == "nexus" else \
             ".codon_positions.partitions"
    return input_path.with_name(input_path.name + suffix)


def read_partitions(path: Path) -> list[Partition]:
    if not path.is_file():
        raise ValueError(f"Input partition file not found: {path}")

    partitions: list[Partition] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            match = PARTITION_RE.match(line)
            if not match:
                raise ValueError(
                    f"Could not parse line {line_number}:\n{raw_line.rstrip()}\n"
                    "Expected syntax: DNA, partition_name = start-end"
                )

            start = int(match.group("start"))
            end = int(match.group("end"))

            if start < 1:
                raise ValueError(
                    f"Line {line_number}: start coordinate must be >= 1."
                )
            if end < start:
                raise ValueError(
                    f"Line {line_number}: end coordinate {end} is smaller "
                    f"than start coordinate {start}."
                )

            partitions.append(
                Partition(
                    line_number=line_number,
                    datatype=match.group("datatype").strip(),
                    name=match.group("name").strip(),
                    start=start,
                    end=end,
                )
            )

    if not partitions:
        raise ValueError(f"No partition records were found in {path}")

    return partitions


def validate_partitions(
    partitions: list[Partition],
    allow_gaps: bool,
    allow_nontriplet: bool,
) -> list[str]:
    warnings: list[str] = []

    previous: Partition | None = None

    for partition in partitions:
        if partition.length % 3 != 0:
            message = (
                f"Line {partition.line_number}, {partition.name}: "
                f"length {partition.length} is not divisible by 3."
            )
            if allow_nontriplet:
                warnings.append(message)
            else:
                raise ValueError(message)

        if previous is not None:
            if partition.start <= previous.end:
                raise ValueError(
                    f"Overlap detected: {previous.name} ends at {previous.end}, "
                    f"but {partition.name} starts at {partition.start}."
                )

            expected_start = previous.end + 1
            if partition.start != expected_start:
                message = (
                    f"Gap detected after {previous.name}: expected the next "
                    f"partition to start at {expected_start}, but "
                    f"{partition.name} starts at {partition.start}."
                )
                if allow_gaps:
                    warnings.append(message)
                else:
                    raise ValueError(message)

        previous = partition

    return warnings


def read_phylip_site_count(path: Path) -> int:
    if not path.is_file():
        raise ValueError(f"Alignment file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.strip().split()
            if not fields:
                continue
            if len(fields) < 2:
                break
            try:
                ntax = int(fields[0])
                nchar = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse PHYLIP header from {path}: "
                    f"{raw_line.rstrip()}"
                ) from exc

            if ntax < 1 or nchar < 1:
                raise ValueError(
                    f"Invalid PHYLIP dimensions in {path}: {ntax} {nchar}"
                )
            return nchar

    raise ValueError(f"No PHYLIP header was found in {path}")


def codon_ranges(
    partitions: Iterable[Partition],
    codon_position: int,
) -> list[str]:
    offset = codon_position - 1
    return [
        f"{partition.start + offset}-{partition.end}\\3"
        for partition in partitions
    ]


def write_raxml(
    output_path: Path,
    partitions: list[Partition],
    datatype: str,
) -> None:
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for position in (1, 2, 3):
            ranges = ",".join(codon_ranges(partitions, position))
            handle.write(
                f"{datatype}, codon_pos{position} = {ranges}\n"
            )


def write_wrapped_ranges(
    handle,
    ranges: list[str],
    indent: str = "        ",
    per_line: int = 8,
) -> None:
    for index in range(0, len(ranges), per_line):
        chunk = ranges[index:index + per_line]
        is_last = index + per_line >= len(ranges)
        ending = ";\n" if is_last else ",\n"
        handle.write(indent + ", ".join(chunk) + ending)


def write_nexus(output_path: Path, partitions: list[Partition]) -> None:
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#nexus\n")
        handle.write("begin sets;\n")

        for position in (1, 2, 3):
            handle.write(f"    charset codon_pos{position} =\n")
            write_wrapped_ranges(
                handle,
                codon_ranges(partitions, position),
            )

        handle.write("end;\n")


def write_report(
    report_path: Path,
    partitions: list[Partition],
) -> None:
    with report_path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "line_number",
                "datatype",
                "partition_name",
                "start",
                "end",
                "length",
                "length_mod_3",
                "codon_pos1_start",
                "codon_pos2_start",
                "codon_pos3_start",
            ]
        )

        for partition in partitions:
            writer.writerow(
                [
                    partition.line_number,
                    partition.datatype,
                    partition.name,
                    partition.start,
                    partition.end,
                    partition.length,
                    partition.length % 3,
                    partition.start,
                    partition.start + 1,
                    partition.start + 2,
                ]
            )


def main() -> int:
    args = parse_args()

    output_path = args.output or default_output_path(
        args.input_partition,
        args.format,
    )
    report_path = args.report or Path(
        str(output_path) + ".validation.tsv"
    )

    try:
        partitions = read_partitions(args.input_partition)
        warnings = validate_partitions(
            partitions,
            allow_gaps=args.allow_gaps,
            allow_nontriplet=args.allow_nontriplet,
        )

        if args.alignment is not None:
            alignment_sites = read_phylip_site_count(args.alignment)
            final_coordinate = partitions[-1].end

            if alignment_sites != final_coordinate:
                raise ValueError(
                    f"Alignment contains {alignment_sites} sites, but the "
                    f"last partition ends at {final_coordinate}."
                )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        if args.format == "raxml":
            write_raxml(output_path, partitions, args.datatype)
        else:
            write_nexus(output_path, partitions)

        write_report(report_path, partitions)

    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    total_sites = partitions[-1].end - partitions[0].start + 1

    print("Codon-position partition conversion completed.")
    print(f"Input partitions:   {len(partitions)}")
    print(
        f"Coordinate span:    "
        f"{partitions[0].start}-{partitions[-1].end}"
    )
    print(f"Coordinate count:   {total_sites}")
    print(f"Output format:      {args.format}")
    print(f"Partition output:   {output_path}")
    print(f"Validation report:  {report_path}")

    if warnings:
        print(f"Warnings:           {len(warnings)}")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("Validation:         passed")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
