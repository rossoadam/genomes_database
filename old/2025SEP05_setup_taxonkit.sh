#!/usr/bin/env bash
# Setup NCBI taxonomy dump for taxonkit and run a quick test.
# Usage: bash setup_taxonkit.sh

set -euo pipefail

# --- sanity checks ---
if ! command -v curl >/dev/null 2>&1; then
  echo "Error: curl not found. Please install curl and re-run." >&2
  exit 1
fi

if ! command -v tar >/dev/null 2>&1; then
  echo "Error: tar not found. Please install tar and re-run." >&2
  exit 1
fi

if ! command -v taxonkit >/dev/null 2>&1; then
  echo "Error: taxonkit not found in PATH. Activate your env (e.g., 'conda activate <env>') and re-run." >&2
  exit 1
fi

# --- prepare dirs ---
TAXONKIT_HOME="${HOME}/.taxonkit"
mkdir -p "${TAXONKIT_HOME}"

# Use a dedicated temp dir so we can clean it up safely
TMPDIR="$(mktemp -d)"
TAR_GZ="${TMPDIR}/taxdump.tar.gz"

cleanup() {
  rm -rf "${TMPDIR}"
}
trap cleanup EXIT

# --- download and extract ---
echo "Downloading NCBI taxonomy dump..."
curl -L --fail "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" -o "${TAR_GZ}"

echo "Extracting selected files..."
tar -xzf "${TAR_GZ}" -C "${TMPDIR}" names.dmp nodes.dmp delnodes.dmp merged.dmp

# --- install into ~/.taxonkit ---
echo "Copying taxonomy files into ${TAXONKIT_HOME} ..."
cp "${TMPDIR}/names.dmp" "${TMPDIR}/nodes.dmp" "${TMPDIR}/delnodes.dmp" "${TMPDIR}/merged.dmp" "${TAXONKIT_HOME}/"

# --- quick test ---
echo "Running a quick taxonkit test (taxid 190870)..."
echo 190870 \
  | taxonkit lineage \
  | taxonkit reformat -f "{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F

echo "Done. Taxonkit taxonomy data installed in ${TAXONKIT_HOME}."

