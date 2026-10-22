#!/usr/bin/env python3
"""Safely normalize filenames created by 12_nhphyml_parallel_v2.sh.

Old: <ortholog_id>__<ortholog_id>.phylip_nhPhyml<suffix>
New: <ortholog_id>_nhPhyml<suffix>

Only regular files directly inside the supplied directory are candidates.
A dry run is the default. Use --execute to perform the renames.

Safety notes:
  * The complete rename plan is checked for name conflicts before execution.
  * Destination filenames are never overwritten, including in a race after
    preflight: an exclusive hard link is made before the old name is unlinked.
  * No file content is read or modified. The filesystem must support hard
    links for the renames to succeed.
  * A runtime error may leave a partially completed batch; rerunning is safe.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# Require *identical* ortholog IDs on both sides of the double underscore.
OLD_PATTERN = re.compile(
    r"(?P<ortholog_id>.+)__(?P=ortholog_id)\.phylip_nhPhyml(?P<suffix>.*)"
)
NEW_PATTERN = re.compile(r".+_nhPhyml.*")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rename existing v2 nhPhyml output files using concise ortholog "
            "IDs, without changing file contents. Defaults to --dry-run."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory containing nhPhyml outputs (top level only).",
    )
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--dry-run", action="store_true", help="Print the proposed changes only."
    )
    modes.add_argument(
        "--execute", action="store_true", help="Perform validated renames."
    )
    return parser.parse_args()


def print_summary(
    directory: Path,
    *,
    mode: str,
    inspected: int,
    matched: int,
    already_new: int,
    skipped: int,
    conflicts: int,
    renamed: int,
    errors: int,
) -> None:
    print("\n" + "=" * 66)
    print("nhPhyml output renaming summary")
    print("=" * 66)
    print(f"Directory scanned           : {directory}")
    print(f"Mode                        : {mode}")
    print(f"Regular files inspected     : {inspected}")
    print(f"Files matching old pattern  : {matched}")
    print(f"Files requiring renaming    : {matched}")
    print(f"Already new convention      : {already_new}")
    print(f"Skipped (other files)       : {skipped}")
    print(f"Filename conflicts          : {conflicts}")
    print(f"Successful renames          : {renamed}")
    print(f"Execution errors            : {errors}")
    print("=" * 66)


def main() -> int:
    args = parse_args()
    directory = args.output_dir.expanduser()
    mode = "EXECUTE" if args.execute else "DRY RUN"

    if not directory.is_dir():
        print(f"[ERROR] Not a directory: {directory}", file=sys.stderr)
        return 1

    plan: list[tuple[Path, Path]] = []
    inspected = 0
    already_new = 0
    skipped = 0

    try:
        # Sorting makes dry runs and logs reproducible.
        for source in sorted(directory.iterdir(), key=lambda p: p.name):
            # Do not follow symlinks (even symlinks pointing to regular files).
            if source.is_symlink() or not source.is_file():
                continue
            inspected += 1
            match = OLD_PATTERN.fullmatch(source.name)
            if match:
                ortholog = match.group("ortholog_id")
                suffix = match.group("suffix")
                destination = directory / f"{ortholog}_nhPhyml{suffix}"
                plan.append((source, destination))
            elif ".phylip_nhPhyml" not in source.name and NEW_PATTERN.fullmatch(source.name):
                already_new += 1
            else:
                skipped += 1
    except OSError as exc:
        print(f"[ERROR] Cannot inspect directory {directory}: {exc}", file=sys.stderr)
        return 1

    print(f"Mode: {mode}")
    print(f"Output directory: {directory}")
    print(f"Proposed renames: {len(plan)}")
    for source, destination in plan:
        print(f"  {source.name}  ->  {destination.name}")

    # Collect every collision before deciding whether anything may be renamed.
    collisions: dict[Path, list[str]] = defaultdict(list)
    destination_sources: dict[Path, list[Path]] = defaultdict(list)
    for source, destination in plan:
        destination_sources[destination].append(source)
        # lexists detects even dangling symlinks, preventing their replacement.
        if os.path.lexists(destination):
            collisions[destination].append("destination already exists")

    for destination, sources in destination_sources.items():
        if len(sources) > 1:
            collisions[destination].append(
                "multiple sources: " + ", ".join(path.name for path in sources)
            )

    if collisions:
        print("\n[ERROR] Conflicts detected; no files have been renamed:", file=sys.stderr)
        for destination, reasons in sorted(collisions.items()):
            print(f"  {destination.name}: {'; '.join(reasons)}", file=sys.stderr)
        print_summary(
            directory,
            mode=mode,
            inspected=inspected,
            matched=len(plan),
            already_new=already_new,
            skipped=skipped,
            conflicts=len(collisions),
            renamed=0,
            errors=0,
        )
        return 1

    if not args.execute:
        print("\n[DRY RUN] No files were changed. Use --execute to rename.")
        print_summary(
            directory,
            mode=mode,
            inspected=inspected,
            matched=len(plan),
            already_new=already_new,
            skipped=skipped,
            conflicts=0,
            renamed=0,
            errors=0,
        )
        return 0

    renamed = 0
    errors = 0
    for source, destination in plan:
        # Recheck that the source was not removed or replaced in the meantime.
        if source.is_symlink() or not source.is_file():
            print(
                f"[ERROR] Source disappeared or is no longer a regular file: "
                f"{source}",
                file=sys.stderr,
            )
            errors += 1
            continue

        # An ordinary Path.rename()/os.rename() could overwrite a destination
        # created AFTER preflight on POSIX. Creating a hard link first instead
        # guarantees the target path does not already exist, even in that race.
        try:
            os.link(source, destination, follow_symlinks=False)
        except FileExistsError:
            print(f"[ERROR] Destination appeared during execution: {destination}", file=sys.stderr)
            errors += 1
            continue
        except OSError as exc:
            print(f"[ERROR] Could not create {destination.name}: {exc}", file=sys.stderr)
            errors += 1
            continue

        try:
            source.unlink()
        except OSError as exc:
            print(
                f"[ERROR] Created {destination.name} but could not remove "
                f"{source.name}: {exc}",
                file=sys.stderr,
            )
            errors += 1
            # Best effort to restore the pre-rename state; do not touch source.
            try:
                destination.unlink()
            except OSError as rollback_exc:
                print(
                    f"[ERROR] Could not remove newly created name "
                    f"{destination.name}: {rollback_exc}. Both names may exist.",
                    file=sys.stderr,
                )
            continue

        renamed += 1
        print(f"[RENAMED] {source.name}  ->  {destination.name}")

    print_summary(
        directory,
        mode=mode,
        inspected=inspected,
        matched=len(plan),
        already_new=already_new,
        skipped=skipped,
        conflicts=0,
        renamed=renamed,
        errors=errors,
    )
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
