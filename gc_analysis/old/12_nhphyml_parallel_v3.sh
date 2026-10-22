#!/bin/bash

set -euo pipefail

# Run nhPhyml independently on each PHYLIP alignment, using xargs for parallelism.
#
# Usage:
#   bash 12_nhphyml_parallel_v3.sh \
#     /path/to/phylip_input_dir \
#     /path/to/rooted_tree.tre \
#     40 \
#     [/path/to/nhPhyml_binary] \
#     [/path/to/output_dir]
#
# Arguments (unchanged from v2):
#   1. phylip_input_dir   Directory containing *.phylip files.
#   2. tree_file          Rooted tree supplied to nhPhyml.
#   3. threads            Maximum number of simultaneous nhPhyml jobs.
#   4. nhphyml_binary     Optional binary path.
#   5. output_dir         Optional destination; default: 09_nhphyml_output
#                       alongside the input directory.
#
# Each job uses an isolated temporary working directory. Generated files named
# <ID>.phylip_nhPhyml<suffix> become <ID>_nhPhyml<suffix>, irrespective of suffix.
# Unexpected files are kept as <ID>__<original-name> with a warning.
#
# No existing destination (including a .status file) is overwritten. A lock on
# the output DIRECTORY coordinates v3 writers, and hard-link publication is
# no-clobber even if a noncooperating process creates a destination in the gap
# between checking and publishing it. A job's links are rolled back if the
# publication encounters an error. Files are first staged on the output FS so
# the final hard-link operation cannot fail due to a cross-filesystem link.
# Requires Bash, GNU find/xargs, and flock (Linux).

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <phylip_input_dir> <tree_file> <threads> [nhphyml_binary] [output_dir]"
    exit 1
fi

INPUT_DIR="$1"
TREE_FILE="$2"
THREADS="$3"
NHPHYML_BIN="${4:-/home/lepidodactylus/bin/nhPhyml/nhPhyml}"

if [ -n "${5:-}" ]; then
    OUTPUT_DIR="$5"
else
    PARENT_DIR=$(dirname "${INPUT_DIR%/}")
    OUTPUT_DIR="$PARENT_DIR/09_nhphyml_output"
fi

# Preserve the v2 input validation and CLI.
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR" >&2
    exit 1
fi

if [ ! -f "$TREE_FILE" ]; then
    echo "Error: Tree file not found: $TREE_FILE" >&2
    exit 1
fi

if [ ! -x "$NHPHYML_BIN" ]; then
    echo "Error: nhPhyml binary not executable: $NHPHYML_BIN" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer. Got: $THREADS" >&2
    exit 1
fi

if ! command -v flock >/dev/null 2>&1; then
    echo "Error: flock is required for safe concurrent output transfers." >&2
    exit 1
fi

mkdir -p -- "$OUTPUT_DIR"

echo "Input directory : $INPUT_DIR"
echo "Tree file       : $TREE_FILE"
echo "Threads         : $THREADS"
echo "nhPhyml binary  : $NHPHYML_BIN"
echo "Output directory: $OUTPUT_DIR"
echo "-----------------------------------"

# Count filenames safely, including unusual names with embedded newlines.
PHYLIP_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.phylip' -printf '.' | wc -c)
PHYLIP_COUNT=$(echo "$PHYLIP_COUNT" | tr -d '[:space:]')

if [ "$PHYLIP_COUNT" -eq 0 ]; then
    echo "No .phylip files found in: $INPUT_DIR"
    exit 1
fi

# Per-invocation records; existing .status files cannot affect the summary.
RUN_MARKERS_DIR=$(mktemp -d "${TMPDIR:-/tmp}/nhphyml_invocation_XXXXXXXX")
trap 'rm -rf -- "$RUN_MARKERS_DIR"' EXIT

echo "Found $PHYLIP_COUNT phylip files."
echo "Launching parallel nhPhyml jobs..."
echo "-----------------------------------"

# Called ONLY after a job has failed. Existing status files, including SUCCESS,
# are left untouched. With no existing status, create a new FAILED status using
# Bash noclobber; the directory lock coordinates concurrent v3 processes.
record_failed_status() (
    local base_name="$1"
    local status_file="$OUTPUT_DIR/${base_name}.status"
    local status_fd
    local old_status

    if ! exec {status_fd}<"$OUTPUT_DIR" || ! flock -x "$status_fd"; then
        echo "[ERROR] [$base_name] Cannot lock output directory to write FAILED status." >&2
        return 1
    fi

    if [[ -e "$status_file" || -L "$status_file" ]]; then
        old_status=$(head -n 1 -- "$status_file" 2>/dev/null || true)
        echo "[WARNING] [$base_name] Preserving existing status: $status_file (${old_status:-unreadable/empty})." >&2
        echo "[WARNING] [$base_name] THIS invocation failed; prior status was not changed." >&2
        return 0
    fi

    if ! (set -C; printf 'FAILED\n' > "$status_file"); then
        echo "[ERROR] [$base_name] Could not create FAILED status: $status_file" >&2
        return 1
    fi
)

job_failed() {
    local base_name="$1"
    local reason="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAIL  $base_name — $reason" >&2
    record_failed_status "$base_name" || true
    return 1
}

run_one() (
    local phylip_file="$1"
    local base_name job_tmp='' seq_copy tree_copy tree_input_name
    local f fname dest
    local -a sources=() destinations=()
    local -A seen_destinations=()

    base_name=$(basename "$phylip_file" .phylip)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $base_name"

    if ! job_tmp=$(mktemp -d "${TMPDIR:-/tmp}/nhphyml_${base_name}_XXXXXX"); then
        job_failed "$base_name" 'could not create temporary working directory'
        return 1
    fi
    trap 'rm -rf -- "$job_tmp"' EXIT

    seq_copy="$job_tmp/${base_name}.phylip"
    tree_input_name=$(basename "$TREE_FILE")
    tree_copy="$job_tmp/$tree_input_name"

    if ! cp -- "$phylip_file" "$seq_copy" || ! cp -- "$TREE_FILE" "$tree_copy"; then
        job_failed "$base_name" 'could not copy input alignment or tree'
        return 1
    fi

    # Scientific settings and nhPhyml invocation are IDENTICAL to v2.
    if ! (
        cd "$job_tmp"
        "$NHPHYML_BIN" \
            -sequences="$seq_copy" \
            -tree="$tree_copy" \
            -format=i \
            -positions=3 \
            -tstv=e \
            -rates=1 \
            -alpha=e \
            -topology=k \
            -outseqs=n \
            -eqfreq=lim \
            -numeqfreq=5 \
            -treefile="${base_name}.treefile"
    ); then
        job_failed "$base_name" 'nhPhyml returned a nonzero exit status'
        return 1
    fi

    # Build and validate the COMPLETE filename mapping before publishing any
    # result. dotglob includes unexpected hidden output files too.
    shopt -s nullglob dotglob
    for f in "$job_tmp"/*; do
        fname=${f##*/}
        case "$fname" in
            "${base_name}.phylip"|"$tree_input_name")
                continue  # temporary INPUT copies, not analysis results
                ;;
        esac

        if [[ ! -f "$f" ]]; then
            echo "[ERROR] [$base_name] Unexpected non-regular output cannot be transferred: $f" >&2
            job_failed "$base_name" 'untransferable generated output'
            return 1
        fi

        case "$fname" in
            "${base_name}.treefile")
                dest="${base_name}.treefile"
                ;;
            "${base_name}.phylip_nhPhyml"*)
                dest="${base_name}_nhPhyml${fname#"${base_name}.phylip_nhPhyml"}"
                ;;
            *)
                dest="${base_name}__${fname}"
                echo "[WARNING] [$base_name] Unexpected output: '$fname' -> '$dest'" >&2
                ;;
        esac

        if [[ -n "${seen_destinations[$dest]+present}" ]]; then
            echo "[ERROR] [$base_name] Two job outputs map to the same filename '$dest': '${seen_destinations[$dest]}' and '$fname'." >&2
            job_failed "$base_name" 'duplicate destination names within job'
            return 1
        fi
        seen_destinations["$dest"]="$fname"
        sources+=("$f")
        destinations+=("$dest")
    done
    shopt -u nullglob dotglob

    # Publish in a separate subshell, so its output-directory lock is released
    # before job_failed() tries to write FAILED.status. The staging directory
    # is always removed on return (including a conflict/error).
    publish_outputs() (
        local stage_dir='' lock_fd i path published_dest
        local -a published=()
        trap 'if [[ -n "$stage_dir" ]]; then rm -rf -- "$stage_dir"; fi' EXIT

        if ! exec {lock_fd}<"$OUTPUT_DIR" || ! flock -x "$lock_fd"; then
            echo "[ERROR] [$base_name] Cannot obtain output-directory lock." >&2
            return 1
        fi

        # Include status in the upfront collision check; no results are moved
        # if ANY destination exists (including a dangling symlink).
        for dest in "${destinations[@]}" "${base_name}.status"; do
            path="$OUTPUT_DIR/$dest"
            if [[ -e "$path" || -L "$path" ]]; then
                echo "[ERROR] [$base_name] Output filename conflict: $path" >&2
                return 1
            fi
        done

        if ! stage_dir=$(mktemp -d "$OUTPUT_DIR/.nhphyml-stage.XXXXXXXXXX"); then
            echo "[ERROR] [$base_name] Cannot create staging directory in $OUTPUT_DIR" >&2
            return 1
        fi

        for i in "${!sources[@]}"; do
            if ! cp -p -- "${sources[i]}" "$stage_dir/${destinations[i]}"; then
                echo "[ERROR] [$base_name] Failed staging: ${sources[i]}" >&2
                return 1
            fi
        done
        if ! printf 'SUCCESS\n' > "$stage_dir/${base_name}.status"; then
            echo "[ERROR] [$base_name] Failed staging SUCCESS status." >&2
            return 1
        fi

        # A hard link fails atomically when a destination already exists; it
        # never clobbers. Links are made on one FS after stage copies succeed.
        for dest in "${destinations[@]}" "${base_name}.status"; do
            path="$OUTPUT_DIR/$dest"
            if ln -- "$stage_dir/$dest" "$path"; then
                published+=("$path")
            else
                echo "[ERROR] [$base_name] Could not publish '$path' (possible concurrent conflict)." >&2
                for published_dest in "${published[@]}"; do
                    if ! rm -f -- "$published_dest"; then
                        echo "[ERROR] [$base_name] Could not roll back: $published_dest" >&2
                    fi
                done
                return 1
            fi
        done
    )

    if ! publish_outputs; then
        job_failed "$base_name" 'output transfer failed (conflict or filesystem error)'
        return 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE  $base_name"
    return 0
)

export INPUT_DIR TREE_FILE OUTPUT_DIR NHPHYML_BIN RUN_MARKERS_DIR
export -f record_failed_status job_failed run_one

# Every worker writes a unique per-invocation marker, completely independent of
# preexisting .status files. Returning 1 (not 255) lets xargs schedule other jobs.
xargs_rc=0
if find "$INPUT_DIR" -maxdepth 1 -type f -name '*.phylip' -print0 |
    xargs -0 -n 1 -P "$THREADS" bash -c '
        if run_one "$1"; then
            result=SUCCESS
        else
            result=FAILED
        fi
        if ! mktemp "$RUN_MARKERS_DIR/${result}.XXXXXXXXXX" >/dev/null; then
            echo "[ERROR] Could not record current-invocation result for: $1" >&2
            exit 1
        fi
        if [[ "$result" == SUCCESS ]]; then exit 0; else exit 1; fi
    ' _; then
    :
else
    xargs_rc=$?
fi

SUCCESS_COUNT=$(find "$RUN_MARKERS_DIR" -maxdepth 1 -type f -name 'SUCCESS.*' -printf '.' | wc -c)
FAILED_COUNT=$(find "$RUN_MARKERS_DIR" -maxdepth 1 -type f -name 'FAILED.*' -printf '.' | wc -c)
SUCCESS_COUNT=$(echo "$SUCCESS_COUNT" | tr -d '[:space:]')
FAILED_COUNT=$(echo "$FAILED_COUNT" | tr -d '[:space:]')
NOT_RUN_COUNT=$(( PHYLIP_COUNT - SUCCESS_COUNT - FAILED_COUNT ))

echo "-----------------------------------"
echo "Completed nhPhyml runs (this invocation)."
echo "Successful: $SUCCESS_COUNT"
echo "Failed    : $FAILED_COUNT"
echo "Not recorded / not run: $NOT_RUN_COUNT"
echo "Outputs   : $OUTPUT_DIR"

if [ "$FAILED_COUNT" -gt 0 ] || [ "$NOT_RUN_COUNT" -ne 0 ] || [ "$xargs_rc" -ne 0 ]; then
    exit 1
fi
