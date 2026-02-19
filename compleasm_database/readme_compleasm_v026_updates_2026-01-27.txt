COMPLEASM v0.2.6 ENVIRONMENT — UPDATE LOG (since prior README)
Last updated: 2026-01-27
Host: macOS (osx-arm64 / Apple Silicon)
Conda env: /Users/rossoaa/miniconda3/envs/compleasm_v026


SUMMARY OF WHAT CHANGED
1) Added missing external dependencies into the conda env
   - Installed HMMER (hmmer 3.4; provides `hmmsearch`)
   - Installed miniprot (miniprot 0.18)
   - Result: compleasm no longer depends on user-level fallback binaries for these tools.

2) Verified tool resolution inside the active environment
   - which compleasm  => /Users/rossoaa/miniconda3/envs/compleasm_v026/bin/compleasm
   - which miniprot   => /Users/rossoaa/miniconda3/envs/compleasm_v026/bin/miniprot
   - which hmmsearch  => /Users/rossoaa/miniconda3/envs/compleasm_v026/bin/hmmsearch

3) Confirmed odb10 lineage directory is present and complete
   - BUSCO lineage root (passed via -L):
       /Users/rossoaa/projects/mass_predicts_dna_dynamics/pipeline/data/busco_lineages/mb_downloads
   - odb10 lineage folder exists:
       sauropsida_odb10/
   - Required files confirmed present (examples):
       lengths_cutoff, scores_cutoff, hmms/, refseq_db.faa.gz, dataset.cfg

4) Validated the stable workflow for odb10 on a real genome
   - `miniprot` produces a non-empty GFF (miniprot_output.gff ~1.1 GB for Podarcis muralis)
   - `compleasm analyze` succeeds on that GFF using sauropsida_odb10

   Example successful analyze output (Podarcis muralis, sauropsida_odb10):
     S:94.08%, 7037
     D:3.09%, 231
     F:0.49%, 37
     I:0.00%, 0
     M:2.34%, 175
     N:7480

5) Notes on workflow behavior observed during debugging
   - When `compleasm run` was executed earlier without a valid lineage file path, it errored on missing `lengths_cutoff`.
     This was resolved by ensuring `-L` points to the directory that CONTAINS `sauropsida_odb10/`.
   - A robust and debuggable approach is to separate the pipeline into:
       (A) miniprot alignment → GFF
       (B) compleasm analyze  → BUSCO-style summary
     This matches how we validated odb10 functionality on Podarcis muralis.


CURRENT ENVIRONMENT STATE (FROM CONDA LIST)
Key packages relevant to this workflow:
  - compleasm  0.2.6        (pypi)
  - miniprot   0.18         (bioconda)   [now inside env]
  - hmmer      3.4          (bioconda)   [now inside env, provides hmmsearch]
  - python     3.10.19
  - pandas     2.3.3
  - numpy      2.2.5
  - matplotlib 3.10.8


WHAT DID NOT CHANGE
- Compleasm version remains 0.2.6 (legacy / stable target).
- Lineage storage layout remains under:
    pipeline/data/busco_lineages/mb_downloads/
- Script-level metadata parsing logic for selecting genomes from genomes_metadata.csv is unchanged
  (only dependency resolution and the validated run/analyze usage were updated).


QUICK CHECK COMMANDS
Activate and verify tools:
  conda activate compleasm_v026
  which compleasm
  which miniprot
  which hmmsearch

Verify lineage contents:
  ls -lh /Users/rossoaa/projects/mass_predicts_dna_dynamics/pipeline/data/busco_lineages/mb_downloads/sauropsida_odb10

Smoke test (two-step, recommended):
  # 1) generate GFF using miniprot (example path)
  /Users/rossoaa/miniconda3/envs/compleasm_v026/bin/miniprot --trans -u -I --outs=0.85 -t 7 --gff \
    /path/to/genome.fna \
    /path/to/mb_downloads/sauropsida_odb10/refseq_db.faa.gz \
    > /path/to/output/miniprot_output.gff

  # 2) analyze the GFF with compleasm
  compleasm analyze -g /path/to/output/miniprot_output.gff -o /tmp/compleasm_analyze_test -l sauropsida_odb10 -L /path/to/mb_downloads -t 7


END
