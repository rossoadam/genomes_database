2026OCT6

conda activate liftoff_2
(liftoff_2) rossoaa@BIOL-135085 sql % python 00a_natural_history_standardized_and_name_audit_v6.py \
/Users/rossoaa/projects/genomes \
        --manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
        --genome-size /Users/rossoaa/projects/genomes/records/natural_history/reptile_summary/species_summary.csv \
        --traits-oskyrko /Users/rossoaa/projects/genomes/records/natural_history/oskyrko_et_al_2024_repttraits.csv \
        --mass-title /Users/rossoaa/projects/genomes/records/natural_history/title_at_al_data.csv \
        --mass-ji /Users/rossoaa/projects/genomes/records/natural_history/ji_et_al.tsv \
        --thermal-bennett /Users/rossoaa/projects/genomes/records/natural_history/bennett_et_al_2019_globaltherm_upload_10_11_17.csv
[OK] Loaded 85 focal species from /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv (manifest)
[INFO] Existing /Users/rossoaa/projects/genomes/records/natural_history_standardized/natural_history_standardized.tsv has no species_pk column (legacy output); canonical IDs will be initialized by 00a
[INFO] 00a generated canonical species_pk for 85 unresolved focal species (1..85) in deterministic species_normalized order
[OK] Canonical species_pk resolved for all 85 focal species
[INFO] Authority: 00a owns species_pk; natural_history_standardized.tsv is the canonical downstream handoff
[OK] species_name_audit columns match the required SQL order
[OK] species_name_audit_pk is non-null, numeric, unique, and sequential from 1..25,496
[OK] All 425 non-null species_pk values are numeric and integer-valued
[OK] species_pk status rules pass: matched and missing_trait are populated; unmatched_raw_trait is NULL
[OK] Each of the 85 focal normalized species maps to exactly one canonical species_pk
[OK] species_name_audit species_pk values agree with the authoritative focal-species mapping
[OK] Wrote /Users/rossoaa/projects/genomes/records/natural_history_standardized/natural_history_standardized.tsv (85 focal species; genome_size=17, mass_preferred=84, ct_max=19, ct_min=11)
[OK] Wrote /Users/rossoaa/projects/genomes/records/sql_tsvs/species_name_audit.tsv (rows=25,496, matched=203, unmatched_raw_trait=25071, missing_trait=222)
[INFO] species_name_audit.tsv SQL NULL representation: \\N
[INFO] 00a owns canonical species_pk assignment; natural_history_standardized.tsv is the authoritative handoff to 00b
[DONE] Natural-history standardization and species-name audit complete.
[INFO] Canonical join key: species_normalized
[INFO] Preferred mass priority: oskyrko > title > ji

python 00b_build_starter_sql_tsvs_v17.py \
	--genomes-dir /Users/rossoaa/projects/genomes \
	--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
	--natural-history /Users/rossoaa/projects/genomes/records/natural_history_standardized/natural_history_standardized.tsv \
    --species-name-audit /Users/rossoaa/projects/genomes/records/sql_tsvs/species_name_audit.tsv \
    --output-dir /Users/rossoaa/projects/genomes/records/sql_tsvs/ \
    --all-accessions

[OK] Wrote SQL TSVs:
  /Users/rossoaa/projects/genomes/records/sql_tsvs/natural_history.tsv	85 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/species_name_audit.tsv	25496 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv	85 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv	114705 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/sequence_type_audit.tsv	85 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/analysis_run.tsv	1 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/window_set.tsv	340 data rows
  /Users/rossoaa/projects/genomes/records/sql_tsvs/window_chunk_manifest.tsv	1360 data rows

[WARN] Sequence type audit rows needing review:
  GCA_003113815.1	status=missing_ncbi_stats	chromosomes=0/\N	components=16536/16536
  GCA_014337955.2	status=review	chromosomes=23/23	components=3514/3513
  GCA_020142125.1	status=review	chromosomes=18/17	components=5291/5291
  GCA_030015325.2	status=review	chromosomes=16/16	components=37/36
  GCA_030412105.1	status=review	chromosomes=19/19	components=30/29
  GCA_038048745.1	status=review	chromosomes=19/19	components=273/272
  GCA_039797435.1	status=review	chromosomes=19/19	components=27/26
  GCA_039880765.1	status=review	chromosomes=18/18	components=179/178
  GCA_041380405.1	status=review	chromosomes=16/16	components=21/20
  GCA_041722995.2	status=review	chromosomes=36/16	components=172/171
  GCA_046524025.1	status=review	chromosomes=19/19	components=222/221
  GCA_051529865.1	status=review	chromosomes=19/19	components=1262/1261
  GCA_053572275.1	status=review	chromosomes=23/23	components=3148/3147
  GCA_055504975.1	status=review	chromosomes=16/16	components=17/16
  GCA_055824665.1	status=review	chromosomes=18/18	components=111/110
  GCA_947247035.1	status=review	chromosomes=19/19	components=384/383
  GCA_947686815.1	status=review	chromosomes=20/20	components=2149/2148
  GCA_964106915.2	status=review	chromosomes=20/20	components=211/210
  GCA_964188175.1	status=review	chromosomes=20/20	components=294/293
  GCA_964188305.1	status=review	chromosomes=20/20	components=507/506
  GCA_964252035.1	status=review	chromosomes=20/20	components=326/325
  GCA_964265115.1	status=review	chromosomes=23/23	components=61/60
  GCA_964270895.1	status=review	chromosomes=20/20	components=349/348
  GCA_964273705.1	status=review	chromosomes=18/18	components=487/486
  GCA_964340495.1	status=review	chromosomes=19/19	components=111/110
  GCA_964659585.1	status=review	chromosomes=18/18	components=210/209
  GCA_965112155.1	status=review	chromosomes=20/20	components=290/289
  GCA_965194845.1	status=review	chromosomes=18/18	components=535/534
  GCA_965280105.1	status=review	chromosomes=15/15	components=44/43
  GCA_965636955.1	status=review	chromosomes=21/21	components=78/77
  GCF_009819535.1	status=review	chromosomes=20/20	components=29/28
  GCF_028583425.1	status=review	chromosomes=19/19	components=76/75
  GCF_032191835.1	status=review	chromosomes=21/21	components=661/660
  GCF_051106095.1	status=review	chromosomes=18/18	components=90/89
  GCF_951804945.1	status=review	chromosomes=19/19	components=51/50
  GCF_964188315.1	status=review	chromosomes=20/20	components=185/184
  GWHFILM00000000.1	status=missing_ncbi_stats	chromosomes=40/\N	components=40/\N

python 00c_build_tsvs_from_gff.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
    --metadata /Users/rossoaa/projects/genomes/records/genomes_metadata.csv

WARNING: GCA_003113815.1: path_to_gff is empty; accession skipped.
WARNING: GCA_003400415.2: path_to_gff is empty; accession skipped.
WARNING: GCA_014337955.2: path_to_gff is empty; accession skipped.
WARNING: GCA_016801065.1: path_to_gff is empty; accession skipped.
WARNING: GCA_019473425.1: path_to_gff is empty; accession skipped.
WARNING: GCA_023653725.1: path_to_gff is empty; accession skipped.
WARNING: GCA_030015325.2: path_to_gff is empty; accession skipped.
WARNING: GCA_033807585.1: path_to_gff is empty; accession skipped.
WARNING: GCA_037367245.1: path_to_gff is empty; accession skipped.
WARNING: GCA_038048745.1: path_to_gff is empty; accession skipped.
WARNING: GCA_039707465.1: path_to_gff is empty; accession skipped.
WARNING: GCA_039880765.1: path_to_gff is empty; accession skipped.
WARNING: GCA_040285375.1: path_to_gff is empty; accession skipped.
WARNING: GCA_041380405.1: path_to_gff is empty; accession skipped.
WARNING: GCA_041722995.2: path_to_gff is empty; accession skipped.
WARNING: GCA_042257475.1: path_to_gff is empty; accession skipped.
WARNING: GCA_043643385.1: path_to_gff is empty; accession skipped.
WARNING: GCA_046524025.1: path_to_gff is empty; accession skipped.
WARNING: GCA_047301725.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050042745.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050231175.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050613815.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051473905.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051529865.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051940855.1: path_to_gff is empty; accession skipped.
WARNING: GCA_052054735.1: path_to_gff is empty; accession skipped.
WARNING: GCA_053572275.1: path_to_gff is empty; accession skipped.
WARNING: GCA_053574215.1: path_to_gff is empty; accession skipped.
WARNING: GCA_054791115.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055504975.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055773955.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055824665.1: path_to_gff is empty; accession skipped.
WARNING: GCA_947247035.1: path_to_gff is empty; accession skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56703', 'gene-PODLI_1NCB018006'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56697', 'gene-PODLI_1NCB001829-2'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56712', 'gene-PODLI_1NCB019844'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56712', 'gene-PODLI_1NCB021970'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56721', 'gene-PODLI_1NCB023742'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56721', 'gene-PODLI_1NCB023742'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56727', 'gene-PODLI_1NCB025569'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56742', 'gene-PODLI_1NCB033971'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56742', 'gene-PODLI_1NCB034247'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56742', 'gene-PODLI_1NCB035688'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56746', 'gene-PODLI_1NCB037447'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56746', 'gene-PODLI_1NCB039080'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56749', 'gene-PODLI_1NCB040292'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56750', 'gene-PODLI_1NCB003740'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56752', 'gene-PODLI_1NCB008200'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56753', 'gene-PODLI_1NCB009169'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56753', 'gene-PODLI_1NCB009364'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56760', 'gene-PODLI_1NCB042169'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56760', 'gene-PODLI_1NCB043676'); duplicate gene row skipped.
WARNING: GCA_947686815.1: duplicate (sequence_pk, gene_id) ('56755', 'gene-PODLI_1NCB013526'); duplicate gene row skipped.
WARNING: GCA_964106635.2: path_to_gff is empty; accession skipped.
WARNING: GCA_964106915.2: path_to_gff is empty; accession skipped.
WARNING: GCA_964188175.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964188305.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964234825.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964252035.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964265115.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964270895.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964273705.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964340495.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964659585.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965112155.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965113305.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965153305.2: path_to_gff is empty; accession skipped.
WARNING: GCA_965194845.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965280105.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965636955.1: path_to_gff is empty; accession skipped.
Wrote 1,881,213 gene rows to: /Users/rossoaa/projects/genomes/records/sql_tsvs/genes.tsv
Wrote 1,194,613 IGS rows to: /Users/rossoaa/projects/genomes/records/sql_tsvs/igs.tsv
Wrote 85 gene summary rows to: /Users/rossoaa/projects/genomes/records/sql_tsvs/gene_summary.tsv
Wrote 85 IGS summary rows to: /Users/rossoaa/projects/genomes/records/sql_tsvs/igs_summary.tsv

mysql -u root -p gc3_dynamics < /Users/rossoaa/Library/CloudStorage/OneDrive-UTArlington/Documents/UTA/Research-001/Research/Genome_Structure_Isochores/gc3_dynamics_9_reduced.sql

conda deactivate
conda activate mysql
python 01_load_starter_sql_tsvs_v5.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password * \
--dry-run

[INFO] Running fast preflight checks...
[OK] Python executable: /Users/rossoaa/miniconda3/envs/mysql/bin/python
[OK] Python environment: conda:mysql (/Users/rossoaa/miniconda3/envs/mysql)
[OK] Python dependency: pymysql 2.2.8
[INFO] Ignoring non-schema helper file(s):
  /Users/rossoaa/projects/genomes/records/sql_tsvs/sequence_type_audit.tsv
  /Users/rossoaa/projects/genomes/records/sql_tsvs/window_chunk_manifest.tsv
[OK] Required table inputs found and readable
[OK] TSV headers match expected schema
[OK] MySQL connection successful: root@localhost:3306/gc3_dynamics
[OK] Required MySQL tables found
[OK] MySQL table columns match loader requirements
[OK] LOCAL INFILE enabled for client/server
[OK] LOAD DATA SQL parameterization check passed
[OK] Fast preflight passed
[INFO] Running full TSV relationship/duplicate validation...
[OK] Full TSV data validation passed.
  natural_history	85 rows	natural_history.tsv
  genomes	85 rows	genomes.tsv
  sequences	114705 rows	sequences.tsv
  genes	1881213 rows	genes.tsv
  igs	1194613 rows	igs.tsv
  gene_summary	85 rows	gene_summary.tsv
  igs_summary	85 rows	igs_summary.tsv
  species_name_audit	25496 rows	species_name_audit.tsv
  analysis_run	1 rows	analysis_run.tsv
  window_set	340 rows	window_set.tsv
  genomic_windows	40885963 rows	340 file(s)
  gc_window_stats	40885963 rows	340 file(s)
  sequence_summary	458684 rows	340 file(s)
  genome_summary	340 rows	340 file(s)

[INFO] Load order:
  natural_history (1 file(s))
  genomes (1 file(s))
  sequences (1 file(s))
  genes (1 file(s))
  igs (1 file(s))
  gene_summary (1 file(s))
  igs_summary (1 file(s))
  species_name_audit (1 file(s))
  analysis_run (1 file(s))
  window_set (1 file(s))
  genomic_windows (340 file(s))
  gc_window_stats (340 file(s))
  sequence_summary (340 file(s))
  genome_summary (340 file(s))

[OK] No SQL changes made.

conda deactivate
conda activate bio
python 03_raw_ortholog_validity_batch.py \
/Users/rossoaa/projects/genomes \
/Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv

Resolved 85 Compleasm metadata rows.
Missing accession roots in Compleasm metadata: 0
Processed GCA_003113815.1: 5122/5929 passed
Processed GCA_003400415.2: 4880/5300 passed
Processed GCA_009733165.1: 4542/5597 passed
Processed GCA_014337955.2: 5658/5902 passed
Processed GCA_016801065.1: 5562/5808 passed
Processed GCA_019473425.1: 4562/5649 passed
Processed GCA_020142125.1: 5233/5616 passed
Processed GCA_021292165.1: 4771/5616 passed
Processed GCA_023653725.1: 5281/5663 passed
Processed GCA_024294585.1: 5441/5751 passed
Processed GCA_030015325.2: 5433/5610 passed
Processed GCA_030412105.1: 5747/5954 passed
Processed GCA_030440675.1: 5688/5907 passed
Processed GCA_030867105.1: 5625/5866 passed
Processed GCA_033807585.1: 5433/5718 passed
Processed GCA_037367245.1: 4919/5662 passed
Processed GCA_038048745.1: 5098/5349 passed
Processed GCA_039707465.1: 5492/5732 passed
Processed GCA_039797435.1: 5480/5756 passed
Processed GCA_039880765.1: 5343/5635 passed
Processed GCA_040285375.1: 5732/5937 passed
Processed GCA_041380405.1: 5761/5950 passed
Processed GCA_041722995.2: 5688/5888 passed
Processed GCA_042257475.1: 5705/5923 passed
Processed GCA_043643385.1: 3362/5666 passed
Processed GCA_046524025.1: 5502/5737 passed
Processed GCA_047301725.1: 5519/5768 passed
Processed GCA_050042745.1: 5580/5862 passed
Processed GCA_050231175.1: 5530/5762 passed
Processed GCA_050613815.1: 5529/5756 passed
Processed GCA_051312515.2: 5743/5953 passed
Processed GCA_051473905.1: 5755/5938 passed
Processed GCA_051529865.1: 5711/5948 passed
Processed GCA_051903775.2: 5315/5536 passed
Processed GCA_051940855.1: 5393/5556 passed
Processed GCA_052054735.1: 4390/5193 passed
Processed GCA_053572275.1: 5584/5849 passed
Processed GCA_053574215.1: 5685/5920 passed
Processed GCA_054791115.1: 4964/5124 passed
Processed GCA_055504975.1: 5703/5911 passed
Processed GCA_055773955.1: 5727/5927 passed
Processed GCA_055824665.1: 5190/5470 passed
Processed GCA_947247035.1: 5426/5702 passed
Processed GCA_947686815.1: 5709/5934 passed
Processed GCA_964106635.2: 5049/5446 passed
Processed GCA_964106915.2: 5244/5666 passed
Processed GCA_964188175.1: 5275/5672 passed
Processed GCA_964188305.1: 5661/5855 passed
Processed GCA_964234825.1: 5030/5452 passed
Processed GCA_964252035.1: 5617/5800 passed
Processed GCA_964265115.1: 5534/5769 passed
Processed GCA_964270895.1: 5296/5673 passed
Processed GCA_964273705.1: 5457/5696 passed
Processed GCA_964340495.1: 5513/5740 passed
Processed GCA_964659585.1: 5521/5740 passed
Processed GCA_965112155.1: 5286/5661 passed
Processed GCA_965113305.1: 5084/5464 passed
Processed GCA_965153305.2: 5050/5467 passed
Processed GCA_965194845.1: 5479/5749 passed
Processed GCA_965280105.1: 5713/5919 passed
Processed GCA_965636955.1: 5755/5953 passed
Processed GCF_009769535.1: 4755/5541 passed
Processed GCF_009819535.1: 5398/5889 passed
Processed GCF_019175285.1: 5122/5759 passed
Processed GCF_021028975.2: 5043/5622 passed
Processed GCF_023053635.1: 5765/5970 passed
Processed GCF_027172205.1: 5307/5697 passed
Processed GCF_027244095.1: 5766/5948 passed
Processed GCF_028583425.1: 5717/5926 passed
Processed GCF_028640845.1: 5356/5618 passed
Processed GCF_029931775.1: 5742/5946 passed
Processed GCF_030035675.1: 5659/5912 passed
Processed GCF_031021105.1: 5158/5532 passed
Processed GCF_032191835.1: 5486/5807 passed
Processed GCF_035046505.1: 5526/5729 passed
Processed GCF_035149785.1: 5625/5858 passed
Processed GCF_035594765.1: 5652/5850 passed
Processed GCF_037176765.1: 5535/5731 passed
Processed GCF_049243985.1: 5676/5917 passed
Processed GCF_051106095.1: 5739/5932 passed
Processed GCF_951804945.1: 5766/5949 passed
Processed GCF_963506605.1: 5630/5854 passed
Processed GCF_964188315.1: 5755/5936 passed
Processed GCF_964194415.1: 5441/5694 passed
INFO GWHFILM00000000.1: Duplicated/Duplicate Compleasm orthologs are allowed by the triploid accession exception; raw CDS checks still apply. Use 04_build_compleasm_feature_tsvs_v5.py to include passing duplicates in summaries/features.
Processed GWHFILM00000000.1: 5457/5763 passed

Finished raw ortholog CDS QC.
Accessions requested: 85
Compleasm rows resolved: 85
Missing accession roots in metadata: 0
Rows skipped because cds_fasta was missing: 0
Rows skipped because FASTA was empty: 0
Total sequences: 488332
Passed raw CDS QC: 459153
Failed raw CDS QC: 29179
Failed length % 3 check: 22451
Failed internal stop check: 27968
Failed invalid base check: 0
Sequences with terminal stop codon: 346231
Detailed QC table written to: /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/body-size-dna-dynamics_manifest_raw_ortholog_validity.tsv
Summary QC table written to: /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/body-size-dna-dynamics_manifest_raw_ortholog_validity_summary.tsv

conda deactivate
conda activate liftoff_2
python 04_build_compleasm_feature_tsvs_v7.py \
--genomes /Users/rossoaa/projects/genomes \
--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
--orthodb /Users/rossoaa/projects/genomes/records/compleasm/mb_downloads/sauropsida_odb12.2025-07-01.tar.gz \
--sequences-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv \
--genomes-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv \
--ortholog-validation /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/body-size-dna-dynamics_manifest_raw_ortholog_validity.tsv


Loaded 85 accessions from manifest
Resolved 85 Compleasm metadata rows
Loaded ortholog-validation entries for 85 accession/accession-root keys
Missing accession roots in Compleasm metadata: 0
Processed GCA_964252035.1: 5617 valid Single orthologs in CDS summaries; 5617 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44707 genomic introns; 16851 genomic flank rows
Processed GCA_964188305.1: 5661 valid Single orthologs in CDS summaries; 5661 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45311 genomic introns; 16983 genomic flank rows
Processed GCF_009819535.1: 5398 valid Single orthologs in CDS summaries; 5398 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42744 genomic introns; 16194 genomic flank rows
Processed GCA_947686815.1: 5709 valid Single orthologs in CDS summaries; 5709 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45633 genomic introns; 17127 genomic flank rows
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_964270895.1: 5296 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5296 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_964106915.2: 5244 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5244 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_965112155.1: 5286 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5286 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028744054.1|gene:gene-ENPEP' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_964234825.1: 5030 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5030 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_964106635.2: 5049 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5049 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028744054.1|gene:gene-ENPEP' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCF_027172205.1: 5307 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5307 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_964188175.1: 5275 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5275 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_965153305.2: 5050 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5050 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
INFO GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv; retaining the ortholog and QC-passing CDS GC values with blank genomic fields
Processed GCA_965113305.1: 5084 valid Single orthologs in CDS summaries; 0 linked to genomic sequences; 5084 unlinked ortholog rows retain CDS GC; 0 genomic introns; 0 genomic flank rows
Processed GCF_963506605.1: 5630 valid Single orthologs in CDS summaries; 5630 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45053 genomic introns; 16890 genomic flank rows
Processed GCA_041722995.2: 5688 valid Single orthologs in CDS summaries; 5688 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45836 genomic introns; 17064 genomic flank rows
Processed GCF_028640845.1: 5356 valid Single orthologs in CDS summaries; 5356 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42842 genomic introns; 16068 genomic flank rows
Processed GCA_051312515.2: 5743 valid Single orthologs in CDS summaries; 5743 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46016 genomic introns; 17229 genomic flank rows
Processed GCF_035594765.1: 5652 valid Single orthologs in CDS summaries; 5652 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45108 genomic introns; 16956 genomic flank rows
Processed GCF_037176765.1: 5535 valid Single orthologs in CDS summaries; 5535 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44101 genomic introns; 16605 genomic flank rows
Processed GCA_042257475.1: 5705 valid Single orthologs in CDS summaries; 5705 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45493 genomic introns; 17115 genomic flank rows
Processed GCA_053572275.1: 5584 valid Single orthologs in CDS summaries; 5584 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43905 genomic introns; 16752 genomic flank rows
Processed GCA_014337955.2: 5658 valid Single orthologs in CDS summaries; 5658 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44891 genomic introns; 16974 genomic flank rows
Processed GCA_053574215.1: 5685 valid Single orthologs in CDS summaries; 5685 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45247 genomic introns; 17055 genomic flank rows
Processed GCA_023653725.1: 5281 valid Single orthologs in CDS summaries; 5281 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42000 genomic introns; 15843 genomic flank rows
Processed GCF_035149785.1: 5625 valid Single orthologs in CDS summaries; 5625 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45012 genomic introns; 16875 genomic flank rows
Processed GCA_965280105.1: 5713 valid Single orthologs in CDS summaries; 5713 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45792 genomic introns; 17139 genomic flank rows
Processed GCA_046524025.1: 5502 valid Single orthologs in CDS summaries; 5502 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44236 genomic introns; 16506 genomic flank rows
Processed GCA_043643385.1: 3362 valid Single orthologs in CDS summaries; 3362 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 22112 genomic introns; 10086 genomic flank rows
Processed GCA_039797435.1: 5480 valid Single orthologs in CDS summaries; 5480 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43643 genomic introns; 16440 genomic flank rows
Processed GCA_003400415.2: 4880 valid Single orthologs in CDS summaries; 4880 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 37042 genomic introns; 14640 genomic flank rows
Processed GCA_030015325.2: 5433 valid Single orthologs in CDS summaries; 5433 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43524 genomic introns; 16299 genomic flank rows
Processed GCA_030412105.1: 5747 valid Single orthologs in CDS summaries; 5747 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46014 genomic introns; 17241 genomic flank rows
Processed GCA_054791115.1: 4964 valid Single orthologs in CDS summaries; 4964 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 39674 genomic introns; 14892 genomic flank rows
Processed GCA_052054735.1: 4390 valid Single orthologs in CDS summaries; 4390 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 33149 genomic introns; 13170 genomic flank rows
Processed GCA_050231175.1: 5530 valid Single orthologs in CDS summaries; 5530 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44337 genomic introns; 16590 genomic flank rows
Processed GCF_023053635.1: 5765 valid Single orthologs in CDS summaries; 5765 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46280 genomic introns; 17295 genomic flank rows
Processed GCF_031021105.1: 5158 valid Single orthologs in CDS summaries; 5158 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 40929 genomic introns; 15474 genomic flank rows
Processed GCF_028583425.1: 5717 valid Single orthologs in CDS summaries; 5717 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45676 genomic introns; 17151 genomic flank rows
Processed GCF_029931775.1: 5742 valid Single orthologs in CDS summaries; 5742 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46045 genomic introns; 17226 genomic flank rows
Processed GCA_030440675.1: 5688 valid Single orthologs in CDS summaries; 5688 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45410 genomic introns; 17064 genomic flank rows
Processed GCF_027244095.1: 5766 valid Single orthologs in CDS summaries; 5766 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46319 genomic introns; 17298 genomic flank rows
Processed GCA_964340495.1: 5513 valid Single orthologs in CDS summaries; 5513 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44284 genomic introns; 16539 genomic flank rows
Processed GCF_032191835.1: 5486 valid Single orthologs in CDS summaries; 5486 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43411 genomic introns; 16458 genomic flank rows
Processed GCA_019473425.1: 4562 valid Single orthologs in CDS summaries; 4562 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 33932 genomic introns; 13686 genomic flank rows
Processed GCA_033807585.1: 5433 valid Single orthologs in CDS summaries; 5433 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43494 genomic introns; 16299 genomic flank rows
Processed GCA_051529865.1: 5711 valid Single orthologs in CDS summaries; 5711 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45587 genomic introns; 17133 genomic flank rows
Processed GCA_040285375.1: 5732 valid Single orthologs in CDS summaries; 5732 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45866 genomic introns; 17196 genomic flank rows
Processed GCA_965194845.1: 5479 valid Single orthologs in CDS summaries; 5479 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43932 genomic introns; 16437 genomic flank rows
Processed GCA_964265115.1: 5534 valid Single orthologs in CDS summaries; 5534 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44379 genomic introns; 16602 genomic flank rows
Processed GCA_965636955.1: 5755 valid Single orthologs in CDS summaries; 5755 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46037 genomic introns; 17265 genomic flank rows
Processed GCA_009733165.1: 4542 valid Single orthologs in CDS summaries; 4542 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 33623 genomic introns; 13626 genomic flank rows
Processed GCA_964273705.1: 5457 valid Single orthologs in CDS summaries; 5457 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43549 genomic introns; 16371 genomic flank rows
Processed GCA_964659585.1: 5521 valid Single orthologs in CDS summaries; 5521 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44495 genomic introns; 16563 genomic flank rows
Processed GCF_049243985.1: 5676 valid Single orthologs in CDS summaries; 5676 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45289 genomic introns; 17028 genomic flank rows
Processed GCA_037367245.1: 4919 valid Single orthologs in CDS summaries; 4919 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 37230 genomic introns; 14757 genomic flank rows
Processed GCA_020142125.1: 5233 valid Single orthologs in CDS summaries; 5233 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 40525 genomic introns; 15699 genomic flank rows
Processed GCF_951804945.1: 5766 valid Single orthologs in CDS summaries; 5766 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46271 genomic introns; 17298 genomic flank rows
Processed GCF_964188315.1: 5755 valid Single orthologs in CDS summaries; 5755 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46049 genomic introns; 17265 genomic flank rows
Processed GCF_051106095.1: 5739 valid Single orthologs in CDS summaries; 5739 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45957 genomic introns; 17217 genomic flank rows
Processed GCA_050613815.1: 5529 valid Single orthologs in CDS summaries; 5529 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44654 genomic introns; 16587 genomic flank rows
Processed GCA_039707465.1: 5492 valid Single orthologs in CDS summaries; 5492 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44149 genomic introns; 16476 genomic flank rows
Processed GCF_030035675.1: 5659 valid Single orthologs in CDS summaries; 5659 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45221 genomic introns; 16977 genomic flank rows
Processed GCA_016801065.1: 5562 valid Single orthologs in CDS summaries; 5562 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44164 genomic introns; 16686 genomic flank rows
Processed GCF_019175285.1: 5122 valid Single orthologs in CDS summaries; 5122 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 38752 genomic introns; 15366 genomic flank rows
Processed GCA_021292165.1: 4771 valid Single orthologs in CDS summaries; 4771 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 36312 genomic introns; 14313 genomic flank rows
Processed GCA_039880765.1: 5343 valid Single orthologs in CDS summaries; 5343 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42752 genomic introns; 16029 genomic flank rows
Processed GCF_021028975.2: 5043 valid Single orthologs in CDS summaries; 5043 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 38108 genomic introns; 15129 genomic flank rows
Processed GCA_041380405.1: 5761 valid Single orthologs in CDS summaries; 5761 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 46047 genomic introns; 17283 genomic flank rows
Processed GCA_051473905.1: 5755 valid Single orthologs in CDS summaries; 5755 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45930 genomic introns; 17265 genomic flank rows
Processed GCF_009769535.1: 4755 valid Single orthologs in CDS summaries; 4755 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 36827 genomic introns; 14265 genomic flank rows
Processed GCF_035046505.1: 5526 valid Single orthologs in CDS summaries; 5526 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44213 genomic introns; 16578 genomic flank rows
Processed GCA_050042745.1: 5580 valid Single orthologs in CDS summaries; 5580 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44090 genomic introns; 16740 genomic flank rows
Processed GCA_047301725.1: 5519 valid Single orthologs in CDS summaries; 5519 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44384 genomic introns; 16557 genomic flank rows
Processed GCA_024294585.1: 5441 valid Single orthologs in CDS summaries; 5441 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43472 genomic introns; 16323 genomic flank rows
Processed GCA_051940855.1: 5393 valid Single orthologs in CDS summaries; 5393 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43607 genomic introns; 16179 genomic flank rows
Processed GCA_038048745.1: 5098 valid Single orthologs in CDS summaries; 5098 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 40905 genomic introns; 15294 genomic flank rows
Processed GCA_055504975.1: 5703 valid Single orthologs in CDS summaries; 5703 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45788 genomic introns; 17109 genomic flank rows
Processed GCA_947247035.1: 5426 valid Single orthologs in CDS summaries; 5426 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43194 genomic introns; 16278 genomic flank rows
Processed GCA_030867105.1: 5625 valid Single orthologs in CDS summaries; 5625 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 44783 genomic introns; 16875 genomic flank rows
Processed GCA_051903775.2: 5315 valid Single orthologs in CDS summaries; 5315 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42281 genomic introns; 15945 genomic flank rows
Processed GCA_055773955.1: 5727 valid Single orthologs in CDS summaries; 5727 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 45823 genomic introns; 17181 genomic flank rows
Processed GCA_055824665.1: 5190 valid Single orthologs in CDS summaries; 5190 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 42558 genomic introns; 15570 genomic flank rows
Processed GCF_964194415.1: 5441 valid Single orthologs in CDS summaries; 5441 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43484 genomic introns; 16323 genomic flank rows
Processed GCA_003113815.1: 5122 valid Single orthologs in CDS summaries; 5122 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 38449 genomic introns; 15366 genomic flank rows
Processed GWHFILM00000000.1: 5457 valid Single or Duplicated/Duplicate orthologs in CDS summaries; 5457 linked to genomic sequences; 0 unlinked ortholog rows retain CDS GC; 43665 genomic introns; 16371 genomic flank rows

Finished building Compleasm feature TSVs.
Output directory: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features
sauropsida_odb12 rows: 6118
orthologs rows: 498949
ortholog_source_audit rows: 498949
Ortholog rows retaining CDS GC without sequence_pk: 46621
ortholog_source_audit.tsv is a helper file; do not load it as a v9 SQL table.
ortholog_summary rows: 85
intron_compleasm rows: 3271603
intron_compleasm_summary rows: 85
flanks_compleasm rows: 1237596
flank_sets_compleasm rows: 3
flank_compleasm_summary rows: 255

conda deactivate
conda activate mysql
python 05_to_sql_compleasm_features_v5.py \
  --input-dir /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features \
  --starter-tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
  --tsv-check-only \
  --check-nhphyml-mappings

python 05_to_sql_compleasm_features_v5.py \
  --input-dir /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features \
  --starter-tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
  --tsv-check-only \
  --check-nhphyml-mappings

[INFO] Running TSV-only validation; MySQL is not used
[OK] Python executable: /Users/rossoaa/miniconda3/envs/mysql/bin/python (3.14.4)
[INFO] Ignoring flank_summary.tsv as a compatibility alias; using flank_compleasm_summary.tsv.
[OK] Required Compleasm TSV inputs found and readable
[OK] TSV headers match gc3_dynamics_9_reduced (with explicit legacy alias mapping, if used)
[INFO] sauropsida_odb12: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/sauropsida_odb12.tsv
[INFO] orthologs: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/orthologs.tsv
[INFO] ortholog_summary: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/ortholog_summary.tsv
[INFO] intron_compleasm: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/intron_compleasm.tsv
[INFO] intron_compleasm_summary: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/intron_compleasm_summary.tsv
[INFO] flanks_compleasm: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/flanks_compleasm.tsv
[INFO] flank_sets_compleasm: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/flank_sets_compleasm.tsv
[INFO] flank_compleasm_summary: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features/flank_compleasm_summary.tsv
[INFO] genomes: /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv
[INFO] sequences: /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv
[INFO] natural_history: /Users/rossoaa/projects/genomes/records/sql_tsvs/natural_history.tsv
[INFO] species_name_audit: /Users/rossoaa/projects/genomes/records/sql_tsvs/species_name_audit.tsv
[INFO] Parsing/normalizing TSV rows and building disk-backed validation indexes...
[OK] sauropsida_odb12: 6,118 incoming rows; accepted=6,118, rejected during parsing=0
[OK] flank_sets_compleasm: 3 incoming rows; accepted=3, rejected during parsing=0
[INFO] orthologs: parsed 250,000 rows
[OK] orthologs: 498,949 incoming rows; accepted=498,949, rejected during parsing=0
[OK] ortholog_summary: 85 incoming rows; accepted=85, rejected during parsing=0
[INFO] intron_compleasm: parsed 250,000 rows
[INFO] intron_compleasm: parsed 500,000 rows
[INFO] intron_compleasm: parsed 750,000 rows
[INFO] intron_compleasm: parsed 1,000,000 rows
[INFO] intron_compleasm: parsed 1,250,000 rows
[INFO] intron_compleasm: parsed 1,500,000 rows
[INFO] intron_compleasm: parsed 1,750,000 rows
[INFO] intron_compleasm: parsed 2,000,000 rows
[INFO] intron_compleasm: parsed 2,250,000 rows
[INFO] intron_compleasm: parsed 2,500,000 rows
[INFO] intron_compleasm: parsed 2,750,000 rows
[INFO] intron_compleasm: parsed 3,000,000 rows
[INFO] intron_compleasm: parsed 3,250,000 rows
[OK] intron_compleasm: 3,271,603 incoming rows; accepted=3,271,603, rejected during parsing=0
[OK] intron_compleasm_summary: 85 incoming rows; accepted=85, rejected during parsing=0
[INFO] flanks_compleasm: parsed 250,000 rows
[INFO] flanks_compleasm: parsed 500,000 rows
[INFO] flanks_compleasm: parsed 750,000 rows
[INFO] flanks_compleasm: parsed 1,000,000 rows
[OK] flanks_compleasm: 1,237,596 incoming rows; accepted=1,237,596, rejected during parsing=0
[OK] flank_compleasm_summary: 255 incoming rows; accepted=255, rejected during parsing=0
[OK] genomes: 85 incoming rows; accepted=85, rejected during parsing=0
[OK] sequences: 114,705 incoming rows; accepted=114,705, rejected during parsing=0
[OK] natural_history: 85 incoming rows; accepted=85, rejected during parsing=0
[OK] species_name_audit: 25,496 incoming rows; accepted=25,496, rejected during parsing=0
[INFO] Indexed sauropsida_odb12 for offline validation
[INFO] Indexed orthologs for offline validation
[INFO] Indexed ortholog_summary for offline validation
[INFO] Indexed intron_compleasm for offline validation
[INFO] Indexed intron_compleasm_summary for offline validation
[INFO] Indexed flanks_compleasm for offline validation
[INFO] Indexed flank_sets_compleasm for offline validation
[INFO] Indexed flank_compleasm_summary for offline validation
[INFO] Indexed genomes for offline validation
[INFO] Indexed sequences for offline validation
[INFO] Indexed natural_history for offline validation
[INFO] Indexed species_name_audit for offline validation
[OK] sauropsida_odb12:odb12_id: duplicate rows=0
[OK] orthologs:ortholog_pk: duplicate rows=0
[OK] orthologs:sequence_pk,odb12_id: duplicate rows=0
[OK] ortholog_summary:ortholog_summary_pk: duplicate rows=0
[OK] ortholog_summary:genome_pk: duplicate rows=0
[OK] intron_compleasm:intron_pk: duplicate rows=0
[OK] intron_compleasm:ortholog_pk,intron_id: duplicate rows=0
[OK] intron_compleasm_summary:intron_summary_pk: duplicate rows=0
[OK] intron_compleasm_summary:genome_pk: duplicate rows=0
[OK] flanks_compleasm:flank_pk: duplicate rows=0
[OK] flanks_compleasm:ortholog_pk,flank_set_pk: duplicate rows=0
[OK] flank_sets_compleasm:flank_set_pk: duplicate rows=0
[OK] flank_sets_compleasm:upstream_bp,downstream_bp: duplicate rows=0
[OK] flank_compleasm_summary:flank_summary_pk: duplicate rows=0
[OK] flank_compleasm_summary:genome_pk,flank_set_pk: duplicate rows=0
[OK] genomes:genome_pk: duplicate rows=0
[OK] genomes:accession_id: duplicate rows=0
[OK] sequences:sequence_pk: duplicate rows=0
[OK] sequences:genome_pk,sequence_id: duplicate rows=0
[OK] natural_history:species_pk: duplicate rows=0
[OK] natural_history:species_normalized: duplicate rows=0
[OK] species_name_audit:species_name_audit_pk: duplicate rows=0
[OK] orthologs.sequence_pk -> sequences.sequence_pk: rows=498,949, non-NULL=447,256, distinct=3,588, matched=447,256, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=51,693
[OK] orthologs.genome_pk -> genomes.genome_pk: rows=498,949, non-NULL=498,949, distinct=85, matched=498,949, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] orthologs.odb12_id -> sauropsida_odb12.odb12_id: rows=498,949, non-NULL=498,949, distinct=6,092, matched=498,949, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] intron_compleasm.ortholog_pk -> orthologs.ortholog_pk: rows=3,271,603, non-NULL=3,271,603, distinct=383,814, matched=3,271,603, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] flanks_compleasm.ortholog_pk -> orthologs.ortholog_pk: rows=1,237,596, non-NULL=1,237,596, distinct=412,532, matched=1,237,596, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] flanks_compleasm.flank_set_pk -> flank_sets_compleasm.flank_set_pk: rows=1,237,596, non-NULL=1,237,596, distinct=3, matched=1,237,596, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] flank_compleasm_summary.flank_set_pk -> flank_sets_compleasm.flank_set_pk: rows=255, non-NULL=255, distinct=3, matched=255, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] ortholog_summary.genome_pk -> genomes.genome_pk: rows=85, non-NULL=85, distinct=85, matched=85, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] ortholog_summary.species_pk -> natural_history.species_pk: rows=85, non-NULL=85, distinct=85, matched=85, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] intron_compleasm_summary.genome_pk -> genomes.genome_pk: rows=85, non-NULL=85, distinct=85, matched=85, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] intron_compleasm_summary.species_pk -> natural_history.species_pk: rows=85, non-NULL=85, distinct=85, matched=85, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] flank_compleasm_summary.genome_pk -> genomes.genome_pk: rows=255, non-NULL=255, distinct=85, matched=255, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] flank_compleasm_summary.species_pk -> natural_history.species_pk: rows=255, non-NULL=255, distinct=85, matched=255, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] sequences.genome_pk -> genomes.genome_pk: rows=114,705, non-NULL=114,705, distinct=85, matched=114,705, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] genomes.species_pk -> natural_history.species_pk: rows=85, non-NULL=85, distinct=85, matched=85, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=0
[OK] species_name_audit.species_pk -> natural_history.species_pk: rows=25,496, non-NULL=425, distinct=85, matched=425, unavailable=0, references to rejected parents=0, IDs absent from source=0, absence inconclusive=0 NULL skipped=25,071
[OK] orthologs: genome_pk differs from the referenced sequence's genome_pk; comparable=447,256, mismatches=0, uncomparable=51,693
[OK] intron_compleasm: parent_id differs from referenced ortholog odb12_id; comparable=3,271,603, mismatches=0, uncomparable=0
[OK] ortholog_summary: species_pk differs from the referenced genome's species_pk; comparable=85, mismatches=0, uncomparable=0
[OK] intron_compleasm_summary: species_pk differs from the referenced genome's species_pk; comparable=85, mismatches=0, uncomparable=0
[OK] flank_compleasm_summary: species_pk differs from the referenced genome's species_pk; comparable=255, mismatches=0, uncomparable=0
[INFO] Checking nhPhyML reference mappings (not run inputs or result files)...
[INFO] nhPhyML species mapping route: audit.species_pk -> genomes.species_pk
[INFO] Audit rows usable for the selected route=425, NULL name/ID excluded from that route=25,071
[INFO] nhPhyML references: resolved ortholog rows=498,949, focal genomes=85, resolved species names=85, audit names without a genome=0
[WARN] Offline text comparisons are exact; live SQL collations, string lengths, existing SQL records and upload safety must still be checked in MySQL mode.
[OK] TSV relationship validation passed. No source TSV or SQL data changed.
[OK] Reference mappings passed the requested nhPhyML readiness checks (reference scope only).
  
conda deactivate
conda activate bio
python 06_nhpyhml_to_tsvs_v4.py \
--sql-tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--output-dir

  -h, --help            show this help message and exit
  --sql-tsv-dir SQL_TSV_DIR
                        Directory containing orthologs.tsv, species_name_audit.tsv, sequences.tsv, and genomes.tsv. Default: /genomes/records/sql_tsvs
  --ortholog-alignments-root ORTHOLOG_ALIGNMENTS_ROOT
                        Root used for gc3_per_gene.tsv and detailed.tsv run mappings. Default follows the paths supplied in the request.
  --summary-alignments-root SUMMARY_ALIGNMENTS_ROOT
                        Root used for master_output.tsv run mappings. Default follows the paths supplied in the request.
  --output-dir OUTPUT_DIR


*Note to self, it may make sense to make a final validation script that takes the schema as an input*

=============================================================================================================================================================================================
=============================================================================================================================================================================================
=============================================================================================================================================================================================
=============================================================================================================================================================================================
=============================================================================================================================================================================================
=============================================================================================================================================================================================

2026SEP01 - 2026SEP04
python 00a_natural_history_standardized_and_name_audit.py \
/Users/rossoaa/projects/genomes \
	--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
	--genome-size /Users/rossoaa/projects/genomes/records/natural_history/reptile_summary/species_summary.csv \
	--traits-oskyrko /Users/rossoaa/projects/genomes/records/natural_history/oskyrko_et_al_2024_repttraits.csv \
	--mass-title /Users/rossoaa/projects/genomes/records/natural_history/title_at_al_data.csv \
	--mass-ji /Users/rossoaa/projects/genomes/records/natural_history/ji_et_al.tsv \
	--thermal-bennett /Users/rossoaa/projects/genomes/records/natural_history/bennett_et_al_2019_globaltherm_upload_10_11_17.csv

python 00b_build_starter_sql_tsvs_v15.py \
	--genomes-dir /Users/rossoaa/projects/genomes \
	--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
	--natural-history /Users/rossoaa/projects/genomes/records/natural_history_standardized/natural_history_standardized.tsv \
    --species-name-audit /Users/rossoaa/projects/genomes/records/sql_tsvs/species_name_audit.tsv \
    --output-dir /Users/rossoaa/projects/genomes/records/sql_tsvs/ \
    --all-accessions

python 00c_build_tsvs_from_gff.py \
	../../../genomes \
	--debug

here are some of the outputs to the terminal:
WARNING: GCA_003113815.1: path_to_gff is empty; accession skipped.
WARNING: GCA_003400415.2: path_to_gff is empty; accession skipped.
WARNING: GCA_014337955.2: path_to_gff is empty; accession skipped.
WARNING: GCA_016801065.1: path_to_gff is empty; accession skipped.
WARNING: GCA_019473425.1: path_to_gff is empty; accession skipped.
WARNING: GCA_021292165.1: path_to_gff is empty; accession skipped.
WARNING: GCA_023653725.1: path_to_gff is empty; accession skipped.
WARNING: GCA_030015325.2: path_to_gff is empty; accession skipped.
WARNING: GCA_030412105.1: path_to_gff is empty; accession skipped.
WARNING: GCA_030440675.1: path_to_gff is empty; accession skipped.
WARNING: GCA_033807585.1: path_to_gff is empty; accession skipped.
WARNING: GCA_037367245.1: path_to_gff is empty; accession skipped.
WARNING: GCA_038048745.1: path_to_gff is empty; accession skipped.
WARNING: GCA_039707465.1: path_to_gff is empty; accession skipped.
WARNING: GCA_039880765.1: path_to_gff is empty; accession skipped.
WARNING: GCA_040285375.1: path_to_gff is empty; accession skipped.
WARNING: GCA_041380405.1: path_to_gff is empty; accession skipped.
WARNING: GCA_041722995.2: path_to_gff is empty; accession skipped.
WARNING: GCA_042257475.1: path_to_gff is empty; accession skipped.
WARNING: GCA_043643385.1: path_to_gff is empty; accession skipped.
WARNING: GCA_046524025.1: path_to_gff is empty; accession skipped.
WARNING: GCA_047301725.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050042745.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050231175.1: path_to_gff is empty; accession skipped.
WARNING: GCA_050613815.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051312515.2: path_to_gff is empty; accession skipped.
WARNING: GCA_051473905.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051529865.1: path_to_gff is empty; accession skipped.
WARNING: GCA_051940855.1: path_to_gff is empty; accession skipped.
WARNING: GCA_052054735.1: path_to_gff is empty; accession skipped.
WARNING: GCA_053572275.1: path_to_gff is empty; accession skipped.
WARNING: GCA_053574215.1: path_to_gff is empty; accession skipped.
WARNING: GCA_054791115.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055504975.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055773955.1: path_to_gff is empty; accession skipped.
WARNING: GCA_055824665.1: path_to_gff is empty; accession skipped.
WARNING: GCA_947247035.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964106635.2: path_to_gff is empty; accession skipped.
WARNING: GCA_964106915.2: path_to_gff is empty; accession skipped.
WARNING: GCA_964188175.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964188305.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964234825.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964252035.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964265115.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964270895.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964273705.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964340495.1: path_to_gff is empty; accession skipped.
WARNING: GCA_964659585.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965112155.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965113305.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965153305.2: path_to_gff is empty; accession skipped.
WARNING: GCA_965194845.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965280105.1: path_to_gff is empty; accession skipped.
WARNING: GCA_965636955.1: path_to_gff is empty; accession skipped.

WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 77348 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 143920 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 187742 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 230134 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 264540 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 280136 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 301296 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 319696 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 328594 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 338212 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 357668 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 365380 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 378160 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 388230 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 466018 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 533122 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 576312 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 619276 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 657802 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 673160 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 695386 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 713540 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 725508 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 735604 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 755880 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 761112 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 774534 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 785630 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 862900 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 930386 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 972576 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1015430 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1051130 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1071062 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1089304 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1097984 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1107560 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1127298 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1132328 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff: line 1142912 has fewer than 9 GFF fields; skipping.

After inspecting looks like these lines are blank:

GFF="/Users/rossoaa/projects/genomes/GWHFILM00000000.1/GWHFILM00000000.1_genomic.gff"

awk -F '\t' '
!/^#/ && NF < 9 {
    if ($0 == "")
        print NR, "BLANK"
    else
        print NR, "NON-BLANK:", $0
}' "$GFF"
77348 BLANK
143920 BLANK
187742 BLANK
230134 BLANK
264540 BLANK
280136 BLANK
301296 BLANK
319696 BLANK
328594 BLANK
338212 BLANK
357668 BLANK
365380 BLANK
378160 BLANK
388230 BLANK
466018 BLANK
533122 BLANK
576312 BLANK
619276 BLANK
657802 BLANK
673160 BLANK
695386 BLANK
713540 BLANK
725508 BLANK
735604 BLANK
755880 BLANK
761112 BLANK
774534 BLANK
785630 BLANK
862900 BLANK
930386 BLANK
972576 BLANK
1015430 BLANK
1051130 BLANK
1071062 BLANK
1089304 BLANK
1097984 BLANK
1107560 BLANK
1127298 BLANK
1132328 BLANK
1142912 BLANK

Here is my validation test for the igs:

(base) rossoaa@BIOL-135085 synthetic_genome_v1 % ls -o
total 48
-rw-r--r--  1 rossoaa  1162 Sep  4 11:58 readme.txt
-rw-r--r--  1 rossoaa   199 Sep  4 11:58 test_genomes_metadata.csv
-rw-r--r--  1 rossoaa    26 Sep  4 11:58 test_manifest.csv
-rw-r--r--  2 rossoaa   108 Sep  4 11:58 test_sequences.tsv
-rw-r--r--  1 rossoaa   112 Sep  4 11:40 TEST000000001.1_genomic.fna
-rw-r--r--  1 rossoaa  1253 Sep  4 11:42 TEST000000001.1_genomic.gff
(base) rossoaa@BIOL-135085 synthetic_genome_v1 % python 00c_build_tsvs_from_gff.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_manifest.csv \
    --metadata /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_genomes_metadata.csv \
    --sequences /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_sequences.tsv \
    --output /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_igs.tsv \
    --debug
python: can't open file '/Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/00c_build_tsvs_from_gff.py': [Errno 2] No such file or directory
(base) rossoaa@BIOL-135085 synthetic_genome_v1 % head test.igs.tsv
head: test.igs.tsv: No such file or directory
(base) rossoaa@BIOL-135085 synthetic_genome_v1 % head test_igs.tsv
igs_pk	sequence_pk	igs_seq_start	igs_seq_end	length	gc
1	9999991	10	19	10	1.0
2	9999991	38	47	10	0.0
3	9999991	66	75	10	0.4

these values match what we expect

I edited the 00c script so that the debug function is removed since it seems to work fine:

(liftoff_2) rossoaa@BIOL-135085 sql % python 00c_build_tsvs_from_gff.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_manifest.csv \
    --metadata /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_genomes_metadata.csv \
    --sequences /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_sequences.tsv \
    --output /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_igs.tsv
WARNING: /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/TEST000000001.1_genomic.gff: line 7 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/TEST000000001.1_genomic.gff: line 12 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/TEST000000001.1_genomic.gff: line 17 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/TEST000000001.1_genomic.gff: line 22 has fewer than 9 GFF fields; skipping.
WARNING: /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/TEST000000001.1_genomic.gff: line 27 has fewer than 9 GFF fields; skipping.

Next I need to update the 01_load_starter_sql_tsvs_v3.py script so that it uploads the igs.tsv table as well as the others.

But before doing that I need to confirm which gffs are missing from my local genome repository
in /genomes/ the following accessions are missing gffs:
	GCA_021292165.1

(ncbi_datasets) rossoaa@BIOL-135085 genomes_database % python genomes_file_recovery_v1.py ../../genomes GCA_021292165.1 --source ncbi 

I am going to try the refactored genomes_file_recovery_v1.py using only the GCA_021292165.1
- that looked like it worked.

Now I am doing:
(ncbi_datasets) rossoaa@BIOL-135085 genomes_database % python genomes_file_recovery_v1.py ../../genomes squamata --source ncbi

I ran out of space so I now need to edit genomes.py to compress genomes by default.
1) currently running the temporary script - to compress all but one project manifest
2) I was considering running multiple genomes.py multiple times but,
instead I asked chat to make the following revisions so that it would be a one time run:
	--reconcile-manifest PATH
	Direct accession-by-accession NCBI/GWH repository checking
	Automatic decompression/activation of manifest genomes
	Missing FASTA and GFF recovery without redownloading valid files
	Manifest assemblies left uncompressed
	Metadata synchronization across CSV, JSONL, species index, and source_metadata.json
	Preservation of prior storage/provenance history
	reconcile_manifest provenance events
	Terminal progress and final summary
	Independent failure handling so one bad accession does not stop the manifest
	No duplicate catalog entries
	Existing NCBI/GWH/storage/rebuild functionality preserved

run:
python genomes_file_recovery_v3.py /Users/rossoaa/projects/genomes \
  --reconcile-manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv

done.

Manifest reconciliation complete.
  Accessions checked: 85
  FASTAs downloaded: 5
  FASTAs decompressed: 0
  GFFs downloaded: 2
  GFFs decompressed: 0
  Repository GFF unavailable/not recovered: 50
  Failed accessions: 1
    - GWHFILM00000000.1

The failure on GWH was expected - I originally downloaded the genome, gff, obtained metadata from the webpag and mannually entered the metadata

Now that I have all the possible gffs for these genomes - I need to rerun 00c_build_tsvs_from_gff.py

python 00c_build_tsvs_from_gff.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_manifest.csv \
    --metadata /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_genomes_metadata.csv \
    --sequences /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_sequences.tsv \
    --output /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_igs.tsv \
    --output-genes /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_genes.tsv \
    --output-gene-summary /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_gene_summary.tsv \
    --output-igs-summary /Users/rossoaa/projects/genomes/records/test_data/synthetic_genome_v1/test_igs_summary.tsv

I did not check the gc values but the files themselves look fine. Now running for I think 35 genomes with gffs in manifest

python 00c_build_tsvs_from_gff.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
    --metadata /Users/rossoaa/projects/genomes/records/genomes_metadata.csv

/Users/rossoaa/projects/genomes/records/phylogenetic_analysis/promoted_products/tree_sets/t1.0_e1_o5_with_i_braminus/3_rooted_nhphyml_input.nwk
/Users/rossoaa/projects/genomes/records/compleasm/alignments/t1.0_e1_o35_with_i_braminus/04_concatenated/concatenated.phy
/Users/rossoaa/projects/genomes/records/phylogenetic_analysis/promoted_products/tree_sets/t1.0_e1_o5_with_i_braminus/1_unrooted.treefile
/Users/rossoaa/projects/genomes/records/phylogenetic_analysis/promoted_products/tree_sets/t1.0_e1_o5_with_i_braminus/2_rooted_pgls_input.nwk
/Users/rossoaa/projects/genomes/records/phylogenetic_analysis/promoted_products/tree_sets/t1.0_e1_o5_with_i_braminus/3_rooted_nhphyml_input.nwk
/Users/rossoaa/projects/genomes/records/phylogenetic_analysis/promoted_products/tree_sets/t1.0_e1_o5_with_i_braminus/4_gc3_divergence.nwk

(base) rossoaa@BIOL-135085 sql_tsvs % awk -F '\t' '$4 > 0' igs_summary.tsv | wc -l           
      36

looks like 35 gffs have made it to the igs_summmary.tsv
next I think I need to create the following two tables by hand
	tree_file x
	nhphyml_runs x
I'll need to take a closer look at the following to see how they can be automated
	nhphyml_ortholog_results write prompt in notes x
	nhphyml_species_summary write prompt in notes x

okay, I've written the prompts for the two above tables. Now I need to run the 01 loader:

conda activate mysql
python 01_load_starter_sql_tsvs_v5.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password **** \
--dry-run

Please write a prompt for 01_load_starter_sql_tsvs_v6.py that includes the quick 

python 01_load_starter_sql_tsvs_v6.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password \
--preflight-only

looks good

python 01_load_starter_sql_tsvs_v6.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password \
--skip-data-validation

ran into error

python 01_load_starter_sql_tsvs_v7.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password \
--preflight-only

since validation worked.. I think I can skip uploading for now...

python 01_load_starter_sql_tsvs_v7.py \
--tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
--mysql-db gc3_dynamics \
--mysql-user root \
--mysql-password \
--skip-data-validation

conda activate bio
python 03_raw_ortholog_validity_batch.py \
/Users/rossoaa/projects/genomes \
/Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv

conda activate liftoff_2
python 04_build_compleasm_feature_tsvs_v7.py \
--genomes /Users/rossoaa/projects/genomes \
--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
--orthodb /Users/rossoaa/projects/genomes/records/compleasm/mb_downloads/sauropsida_odb12.2025-07-01.tar.gz \
--sequences-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv \
--genomes-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv \
--ortholog-validation /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/body-size-dna-dynamics_manifest_raw_ortholog_validity.tsv \

Finished building Compleasm feature TSVs.
Output directory: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features
sauropsida_odb12 rows: 6118
orthologs rows: 498949
ortholog_source_audit rows: 498949
Ortholog rows retaining CDS GC without sequence_pk: 46621
ortholog_source_audit.tsv is a helper file; do not load it as a v9 SQL table.
ortholog_summary rows: 85
intron_compleasm rows: 3271603
intron_compleasm_summary rows: 85
flanks_compleasm rows: 1237596
flank_sets_compleasm rows: 3
flank_compleasm_summary rows: 255

python 06_nhphyml_to_tsvs_v2.py \
  --sql-tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs
  


2026AUG31-----------------------------------

python 00a_make_natural_history_standardized.py \
/Users/rossoaa/projects/genomes \
	--manifest /Users/rossoaa/projects/genomes/records/project_manifests/body-size-dna-dynamics_manifest.csv \
	--genome-size /Users/rossoaa/projects/genomes/records/natural_history/reptile_summary/species_summary.csv \
	--traits-oskyrko /Users/rossoaa/projects/genomes/records/natural_history/oskyrko_et_al_2024_repttraits.csv \
	--mass-title /Users/rossoaa/projects/genomes/records/natural_history/title_at_al_data.csv \
	--mass-ji /Users/rossoaa/projects/genomes/records/natural_history/ji_et_al.tsv \
	--thermal-bennett /Users/rossoaa/projects/genomes/records/natural_history/bennett_et_al_2019_globaltherm_upload_10_11_17.csv \



everything below this is prior to 2026AUG31

  python 00_build_starter_sql_tsvs.py \
    --genomes-dir /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv \
    --natural-history /Users/rossoaa/projects/genomes/records/natural_history/natural_history.tsv \
    --species-name-audit /Users/rossoaa/projects/genomes/records/natural_history/species_name_audit.tsv

  python 01_load_starter_sql_tsvs.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_test \
    --mysql-user root
    
I did this for test data, I don't think I've uploaded all the windows tables.

python 00_build_starter_sql_tsvs_v14.py \
    --genomes-dir /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv \
    --natural-history /Users/rossoaa/projects/genomes/records/natural_history/natural_history.tsv \
    --species-name-audit /Users/rossoaa/projects/genomes/records/natural_history/species_name_audit.tsv
    
 python 03_raw_ortholog_validity_translate_dna.py /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_cds_compleasm.fasta

Finished raw ortholog CDS QC.
Total sequences: 5929
Passed raw CDS QC: 5122
Failed raw CDS QC: 807
Failed length % 3 check: 249
Failed internal stop check: 797
Failed invalid base check: 0
Sequences with terminal stop codon: 4319
QC table written to: /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_cds_compleasm_raw_ortholog_validity.tsv


(bio) rossoaa@ sql % python 03_raw_ortholog_validity_batch.py /Users/rossoaa/projects/genomes /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv
Loaded 84 accessions from manifest: /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv
Resolved 84 Compleasm metadata rows.
Missing accession roots in Compleasm metadata: 0
Processed GCA_003113815.1: 5122/5929 passed
Processed GCA_003400415.2: 4880/5300 passed
Processed GCA_009733165.1: 4542/5597 passed
Processed GCA_014337955.2: 5658/5902 passed
Processed GCA_016801065.1: 5562/5808 passed
Processed GCA_019473425.1: 4562/5649 passed
Processed GCA_020142125.1: 5233/5616 passed
Processed GCA_021292165.1: 4771/5616 passed
Processed GCA_023653725.1: 5281/5663 passed
Processed GCA_024294585.1: 5441/5751 passed
Processed GCA_030015325.2: 5433/5610 passed
Processed GCA_030412105.1: 5747/5954 passed
Processed GCA_030440675.1: 5688/5907 passed
Processed GCA_030867105.1: 5625/5866 passed
Processed GCA_033807585.1: 5433/5718 passed
Processed GCA_037367245.1: 4919/5662 passed
Processed GCA_038048745.1: 5098/5349 passed
Processed GCA_039707465.1: 5492/5732 passed
Processed GCA_039797435.1: 5480/5756 passed
Processed GCA_039880765.1: 5343/5635 passed
Processed GCA_040285375.1: 5732/5937 passed
Processed GCA_041380405.1: 5761/5950 passed
Processed GCA_041722995.2: 5688/5888 passed
Processed GCA_042257475.1: 5705/5923 passed
Processed GCA_043643385.1: 3362/5666 passed
Processed GCA_046524025.1: 5502/5737 passed
Processed GCA_047301725.1: 5519/5768 passed
Processed GCA_050042745.1: 5580/5862 passed
Processed GCA_050231175.1: 5530/5762 passed
Processed GCA_050613815.1: 5529/5756 passed
Processed GCA_051312515.2: 5743/5953 passed
Processed GCA_051473905.1: 5755/5938 passed
Processed GCA_051529865.1: 5711/5948 passed
Processed GCA_051903775.2: 5315/5536 passed
Processed GCA_051940855.1: 5393/5556 passed
Processed GCA_052054735.1: 4390/5193 passed
Processed GCA_053572275.1: 5584/5849 passed
Processed GCA_053574215.1: 5685/5920 passed
Processed GCA_054791115.1: 4964/5124 passed
Processed GCA_055504975.1: 5703/5911 passed
Processed GCA_055773955.1: 5727/5927 passed
Processed GCA_055824665.1: 5190/5470 passed
Processed GCA_947247035.1: 5426/5702 passed
Processed GCA_947686815.1: 5709/5934 passed
Processed GCA_964106635.2: 5049/5446 passed
Processed GCA_964106915.2: 5244/5666 passed
Processed GCA_964188175.1: 5275/5672 passed
Processed GCA_964188305.1: 5661/5855 passed
Processed GCA_964234825.1: 5030/5452 passed
Processed GCA_964252035.1: 5617/5800 passed
Processed GCA_964265115.1: 5534/5769 passed
Processed GCA_964270895.1: 5296/5673 passed
Processed GCA_964273705.1: 5457/5696 passed
Processed GCA_964340495.1: 5513/5740 passed
Processed GCA_964659585.1: 5521/5740 passed
Processed GCA_965112155.1: 5286/5661 passed
Processed GCA_965113305.1: 5084/5464 passed
Processed GCA_965153305.2: 5050/5467 passed
Processed GCA_965194845.1: 5479/5749 passed
Processed GCA_965280105.1: 5713/5919 passed
Processed GCA_965636955.1: 5755/5953 passed
Processed GCF_009769535.1: 4755/5541 passed
Processed GCF_009819535.1: 5398/5889 passed
Processed GCF_019175285.1: 5122/5759 passed
Processed GCF_021028975.2: 5043/5622 passed
Processed GCF_023053635.1: 5765/5970 passed
Processed GCF_027172205.1: 5307/5697 passed
Processed GCF_027244095.1: 5766/5948 passed
Processed GCF_028583425.1: 5717/5926 passed
Processed GCF_028640845.1: 5356/5618 passed
Processed GCF_029931775.1: 5742/5946 passed
Processed GCF_030035675.1: 5659/5912 passed
Processed GCF_031021105.1: 5158/5532 passed
Processed GCF_032191835.1: 5486/5807 passed
Processed GCF_035046505.1: 5526/5729 passed
Processed GCF_035149785.1: 5625/5858 passed
Processed GCF_035594765.1: 5652/5850 passed
Processed GCF_037176765.1: 5535/5731 passed
Processed GCF_049243985.1: 5676/5917 passed
Processed GCF_051106095.1: 5739/5932 passed
Processed GCF_951804945.1: 5766/5949 passed
Processed GCF_963506605.1: 5630/5854 passed
Processed GCF_964188315.1: 5755/5936 passed
Processed GCF_964194415.1: 5441/5694 passed

Finished raw ortholog CDS QC.
Accessions requested: 84
Compleasm rows resolved: 84
Missing accession roots in metadata: 0
Rows skipped because cds_fasta was missing: 0
Rows skipped because FASTA was empty: 0
Total sequences: 482569
Passed raw CDS QC: 453696
Failed raw CDS QC: 28873
Failed length % 3 check: 22233
Failed internal stop check: 27676
Failed invalid base check: 0
Sequences with terminal stop codon: 342203
Detailed QC table written to: /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/mass_predicts_dna_dynamics_with_s_punctatus_manifest_raw_ortholog_validity.tsv
Summary QC table written to: /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/mass_predicts_dna_dynamics_with_s_punctatus_manifest_raw_ortholog_validity_summary.tsv

(liftoff_2) rossoaa@ sql % python extract_introns_for_one_ortholog.py --genome /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna --full_table /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/full_table.tsv --busco_id 414176at8457 --species sphenodon_punctatus
Status: Genome FASTA indexed successfully: /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna

Finished extracting introns.
BUSCO ID: 414176at8457
Species: sphenodon_punctatus
Introns extracted: 78
TSV written to: /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns.tsv
FASTA written to: /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns.fasta
(liftoff_2) rossoaa@ sql % python validate_intron_slicing.py --genome /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna --introns_tsv /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns.tsv
Status: Genome FASTA indexed successfully: /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna

Finished validating intron slicing.
Introns validated: 78
Length math OK: 78/78
Sequence length match: 78/78
Sequence match: 78/78
Validation TSV written to: /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns_slicing_validation.tsv

(liftoff_2) rossoaa@ sql % python validate_intron_slicing.py \ 
        --genome /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna \
        --introns_tsv /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns.tsv \
        --boundary_window 10
Status: Genome FASTA indexed successfully: /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna

Finished validating intron slicing.
Introns validated: 78
Length math OK: 78/78
Sequence length match: 78/78
Sequence match: 78/78
Validation TSV written to: /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_414176at8457_introns_slicing_validation.tsv

(liftoff_2) rossoaa@ sql % python 04_build_compleasm_feature_tsvs_v4.py --genomes /Users/rossoaa/projects/genomes --manifest /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv --orthodb /Users/rossoaa/projects/genomes/records/compleasm/mb_downloads/sauropsida_odb12.2025-07-01.tar.gz --sequences-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv --genomes-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv --ortholog-validation /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/mass_predicts_dna_dynamics_with_s_punctatus_manifest_raw_ortholog_validity.tsv
 --test
Loaded 84 accessions from manifest
Resolved 84 Compleasm metadata rows
Loaded ortholog-validation entries for 168 accession/accession-root keys
Missing accession roots in Compleasm metadata: 0
Processed GCA_003113815.1: 5122 valid Single orthologs for summaries/features
Processed GCA_003400415.2: 4880 valid Single orthologs for summaries/features
Processed GCA_009733165.1: 4542 valid Single orthologs for summaries/features
Processed GCA_014337955.2: 5658 valid Single orthologs for summaries/features
Processed GCA_016801065.1: 5562 valid Single orthologs for summaries/features
Processed GCA_019473425.1: 4562 valid Single orthologs for summaries/features
Processed GCA_020142125.1: 5233 valid Single orthologs for summaries/features
Processed GCA_021292165.1: 4771 valid Single orthologs for summaries/features
Processed GCA_023653725.1: 5281 valid Single orthologs for summaries/features
Processed GCA_024294585.1: 5441 valid Single orthologs for summaries/features
Processed GCA_030015325.2: 5433 valid Single orthologs for summaries/features
Processed GCA_030412105.1: 5747 valid Single orthologs for summaries/features
Processed GCA_030440675.1: 5688 valid Single orthologs for summaries/features
Processed GCA_030867105.1: 5625 valid Single orthologs for summaries/features
Processed GCA_033807585.1: 5433 valid Single orthologs for summaries/features
Processed GCA_037367245.1: 4919 valid Single orthologs for summaries/features
Processed GCA_038048745.1: 5098 valid Single orthologs for summaries/features
Processed GCA_039707465.1: 5492 valid Single orthologs for summaries/features
Processed GCA_039797435.1: 5480 valid Single orthologs for summaries/features
Processed GCA_039880765.1: 5343 valid Single orthologs for summaries/features
Processed GCA_040285375.1: 5732 valid Single orthologs for summaries/features
Processed GCA_041380405.1: 5761 valid Single orthologs for summaries/features
Processed GCA_041722995.2: 5688 valid Single orthologs for summaries/features
Processed GCA_042257475.1: 5705 valid Single orthologs for summaries/features
Processed GCA_043643385.1: 3362 valid Single orthologs for summaries/features
Processed GCA_046524025.1: 5502 valid Single orthologs for summaries/features
Processed GCA_047301725.1: 5519 valid Single orthologs for summaries/features
Processed GCA_050042745.1: 5580 valid Single orthologs for summaries/features
Processed GCA_050231175.1: 5530 valid Single orthologs for summaries/features
Processed GCA_050613815.1: 5529 valid Single orthologs for summaries/features
Processed GCA_051312515.2: 5743 valid Single orthologs for summaries/features
Processed GCA_051473905.1: 5755 valid Single orthologs for summaries/features
Processed GCA_051529865.1: 5711 valid Single orthologs for summaries/features
WARNING: GCA_051903775.2 10114at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10114at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10114at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10114at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10180at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10180at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10180at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10180at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10546at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 10546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11023at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11023at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11023at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11023at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11378at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11378at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11378at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 114629at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 114629at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 114629at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 114629at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11487at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11501at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11501at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11501at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11501at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11684at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11684at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11684at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 11684at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12065at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12065at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12065at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12065at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1210at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1210at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1210at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1210at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12147at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12147at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12147at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12147at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12505at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12505at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12505at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 12505at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 126036at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 126036at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 126036at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 126036at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13242at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13242at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13242at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13242at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13315at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13315at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13315at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 13315at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 144024at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 144024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 144024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 144024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1442at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1442at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1442at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 1442at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15225at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15406at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15406at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15406at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15406at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15598at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15652at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15652at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15652at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15652at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15870at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 15870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16248at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16248at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16248at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16248at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16253at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16253at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16253at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16253at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16418at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 164369at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 164369at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 164369at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 164369at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16744at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16744at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 16744at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17049at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17097at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17097at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17097at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17097at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17370at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17370at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17370at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17370at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17670at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17838at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17838at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17838at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 17838at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18208at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18208at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18208at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18208at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18363at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18404at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18404at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18404at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18404at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 18490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 185888at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 185888at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 185888at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 185888at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19111at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19111at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19111at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19111at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19346at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19346at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19346at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19346at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 193488at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 193488at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 193488at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 193488at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19566at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19566at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19566at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19566at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19670at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 19670at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 20190at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 20190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 20190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 20190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 211910at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 211910at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 211910at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 211910at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21334at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21334at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21334at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21334at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21572at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21572at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21572at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 21572at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2163at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2163at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2163at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2163at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22515at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22515at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22515at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 226181at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 226181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 226181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 226181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22722at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 22722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23099at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23137at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23137at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23137at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23137at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 235058at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 235058at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 235058at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 235058at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23694at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 237189at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 237189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 237189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 237189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23846at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 23846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 241722at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 241722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 241722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 241722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24221at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24221at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24221at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24221at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24689at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24689at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24689at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24689at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 246900at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 246900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 246900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 246900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2470at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24751at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24751at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24751at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24751at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24939at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24939at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24939at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 24939at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25084at8457: intron extraction failed: "Sequence 'JAIFHB020000003.1' not found in FASTA"
WARNING: GCA_051903775.2 25084at8457: flank extraction failed: "Sequence 'JAIFHB020000003.1' not found in FASTA"
WARNING: GCA_051903775.2 25084at8457: flank extraction failed: "Sequence 'JAIFHB020000003.1' not found in FASTA"
WARNING: GCA_051903775.2 25084at8457: flank extraction failed: "Sequence 'JAIFHB020000003.1' not found in FASTA"
WARNING: GCA_051903775.2 251337at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 251337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 251337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 251337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2519at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2519at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2519at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2519at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25201at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25201at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25201at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25201at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 252865at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 252865at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 252865at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 252865at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25365at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25365at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25365at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 25365at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2537at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2537at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2537at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2537at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 254199at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 254199at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 254199at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 254199at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255420at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255420at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255420at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255420at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255828at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255828at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255828at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 255828at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 257067at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 257067at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 257067at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 259467at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 259467at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 259467at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 259467at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2608at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26181at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26181at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 261876at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 261876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 261876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 261876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26250at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26250at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26250at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26250at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26319at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 263942at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 263942at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 263942at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 263942at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26397at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26397at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26397at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 26397at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 265413at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 265413at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 265413at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 265413at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2662at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2662at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2662at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 2662at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 267240at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 267240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 267240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 267240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 268021at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 268021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 268021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 268021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27183at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27183at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27183at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27183at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 272003at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 272003at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 272003at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 272003at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27317at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27317at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27317at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27317at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27362at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 27362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28024at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28024at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 281729at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 281729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 281729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 281729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28517at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28517at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28517at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28517at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 285706at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 285706at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 285706at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 285706at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28874at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28874at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28874at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28874at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28913at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28913at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28913at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 28913at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 292396at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 292396at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 292396at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 292396at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 295986at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 295986at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 295986at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 295986at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 296230at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 296230at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 296230at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 296230at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 298354at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 298354at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 298354at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 298354at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29899at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29899at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29899at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29899at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29953at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29953at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29953at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29953at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29954at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29954at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 29954at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306410at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306410at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306410at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306410at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306928at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306928at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306928at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 306928at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 30759at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 30759at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 30759at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 30759at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 307920at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 307920at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 307920at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 307920at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310279at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310279at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310279at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310279at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310961at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310961at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310961at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 310961at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31101at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31101at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31101at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31101at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31308at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31308at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31308at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31308at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313469at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313469at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313469at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313469at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313850at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313850at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313850at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 313850at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 316216at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 316216at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 316216at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 316216at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31715at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31715at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31715at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 31715at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 317462at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 317462at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 317462at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322815at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322815at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322815at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322815at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3256at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3256at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3256at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3256at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32707at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32707at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32707at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 32707at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327612at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327709at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327709at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327709at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327709at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327839at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327839at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327839at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 327839at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 328608at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 328608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 328608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 328608at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330372at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330372at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330372at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330372at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330448at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330448at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330448at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 330448at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 331795at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 331795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 331795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 331795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33455at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33455at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33455at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33455at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 335345at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 335345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 335345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 335345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33576at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 33576at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 338656at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 338656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 338656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 338656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 341436at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 341436at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 341436at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 34154at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 34154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 34154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 34154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 344977at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 344977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 344977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 344977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 346851at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 346851at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 346851at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 346851at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 347357at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 347357at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 347357at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 347357at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 349075at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 349075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 349075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 349075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35055at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35322at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35322at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3532at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353791at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353791at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353791at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353791at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353956at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353956at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353956at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 353956at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 354270at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 354270at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 354270at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 354270at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 355556at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 355556at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 355556at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 355556at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35673at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35673at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 35673at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 359800at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 359800at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 359800at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 359800at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 360213at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 360213at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 360213at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 360213at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36139at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36139at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36139at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36139at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362729at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362729at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362837at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 362837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363049at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363560at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363560at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363560at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 363560at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36433at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36433at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36433at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36489at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36489at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36489at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36489at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3654at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 3654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365543at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365543at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365543at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365543at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36578at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 36578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365921at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365921at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365921at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 365921at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 366989at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 366989at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 366989at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 366989at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 369178at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 369178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 369178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 369178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370398at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370671at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370671at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370671at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370671at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370724at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370724at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370724at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 370724at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371247at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371699at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371699at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371699at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 371699at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372154at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372154at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372189at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372755at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372755at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372755at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 372755at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37380at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37380at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37380at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37385at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37385at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37385at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37385at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37559at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37559at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37559at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37559at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 375980at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 375980at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 375980at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 375980at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37676at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37676at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37676at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 37676at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378333at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378896at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378896at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378896at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 378896at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 379708at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 379708at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 379708at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 379708at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 380996at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 380996at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 380996at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 380996at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 381898at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 381898at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 381898at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 381898at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 383047at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 383047at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 383047at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 383047at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 385236at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 385236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 385236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 385236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386182at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386182at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386182at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386182at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386368at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386368at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386368at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 386368at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387345at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387615at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 387615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389021at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389021at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389482at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389482at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389482at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389482at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389562at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389562at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389562at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 389562at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 390654at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 390654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 390654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 390654at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392345at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392345at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392686at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 392686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39268at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39268at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39268at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39268at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394419at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394419at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394419at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394419at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394637at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394637at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394637at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394637at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394897at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394897at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394897at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 394897at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 395880at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 395880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 395880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 395880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39598at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 39598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 396130at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 396130at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 396130at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 396130at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397203at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397203at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397203at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397203at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397496at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397496at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397496at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 397496at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 399134at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 399134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 399134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 399134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401351at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401351at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401351at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401351at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40142at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40142at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40142at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40142at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401532at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 401532at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 402568at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 402568at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 402568at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 402568at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40398at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40398at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40434at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40434at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40434at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404622at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404622at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404622at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404622at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404971at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404971at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404971at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 404971at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40520at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40520at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40520at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 40520at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 409527at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 409527at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 409527at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41055at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41055at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41217at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41217at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41217at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41217at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41257at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41257at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41257at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41301at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41301at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41301at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41301at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414071at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414071at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414071at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414071at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414422at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414422at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414422at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414422at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414679at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414679at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414679at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 414679at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415141at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415141at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415141at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415141at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415569at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415569at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415569at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415569at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41565at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41565at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41565at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 41565at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415784at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415784at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415784at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 415784at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416195at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416195at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416195at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416195at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416209at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416209at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416209at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416209at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416904at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416904at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416904at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 416904at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417006at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417006at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417006at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417006at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417533at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417533at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417533at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417533at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417664at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417664at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417664at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417664at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417757at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417757at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417757at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417757at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417809at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417809at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417809at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 417809at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418487at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418700at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418700at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418700at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418700at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418876at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 418876at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419134at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419134at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419166at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419166at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419166at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419166at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419222at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419222at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419222at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419222at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419529at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419529at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419529at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419529at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419633at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419633at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419633at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419633at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419889at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419889at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419889at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 419889at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420008at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420008at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420008at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420008at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420375at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420738at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420738at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420738at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 420738at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421445at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421445at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421445at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421445at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421585at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421585at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421585at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421585at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421962at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421962at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421962at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 421962at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422786at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422786at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422786at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422786at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422959at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422959at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422959at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 422959at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423327at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423327at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423327at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423327at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423796at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423796at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423796at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423796at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423837at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423895at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423900at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 423900at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424099at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424099at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424236at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424236at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424271at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424271at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424271at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424271at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424335at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424335at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424335at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 424335at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425178at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425178at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425190at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425795at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425795at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425858at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425858at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425858at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425858at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425863at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425863at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425863at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425863at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425945at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425945at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425945at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 425945at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426066at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426066at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426066at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426066at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426240at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426240at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426299at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426299at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426299at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426299at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426333at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426333at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426374at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426374at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426374at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426374at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426539at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426539at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426539at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426539at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426722at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426722at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426912at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426912at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426912at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 426912at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427112at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427330at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427330at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427330at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 427330at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428022at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428022at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428022at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428022at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428027at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428285at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428285at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428285at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428285at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428538at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428866at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428866at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428866at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 428866at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429225at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429225at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429319at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429337at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429337at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429615at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429615at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429737at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429737at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429737at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429737at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429857at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429857at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429857at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 429857at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430063at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430292at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430292at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430292at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430292at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430509at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430509at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430509at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430509at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43063at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43063at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430698at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 430698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431081at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431081at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431081at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431081at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431319at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431319at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431425at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431425at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431425at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431425at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431470at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431470at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431487at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431487at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431610at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431610at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431610at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431610at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431686at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 431686at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432206at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432206at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432206at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432206at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432480at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432480at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432480at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432480at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432511at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432511at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432511at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432511at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432546at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432546at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432598at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 432598at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433040at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433040at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433040at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433040at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433112at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433112at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433391at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433391at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433391at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433391at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433578at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433578at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433770at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433770at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433770at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 433770at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434309at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434309at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434309at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434309at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434386at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434386at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434386at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434386at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434497at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434497at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434497at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434497at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434837at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 434837at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 435530at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 435530at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 435530at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43604at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43604at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43604at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43604at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436075at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436075at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436355at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436355at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436355at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436355at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436363at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436363at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43641at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43641at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43641at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43641at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436735at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436735at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436735at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436735at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 436880at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43698at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43698at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43764at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43764at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43764at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43764at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437936at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437936at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437936at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437936at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437987at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437987at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437987at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437997at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437997at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437997at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 437997at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438016at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438016at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438016at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438016at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438049at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438570at8457: intron extraction failed: "Sequence 'JAIFHB020000090.1' not found in FASTA"
WARNING: GCA_051903775.2 438570at8457: flank extraction failed: "Sequence 'JAIFHB020000090.1' not found in FASTA"
WARNING: GCA_051903775.2 438570at8457: flank extraction failed: "Sequence 'JAIFHB020000090.1' not found in FASTA"
WARNING: GCA_051903775.2 438570at8457: flank extraction failed: "Sequence 'JAIFHB020000090.1' not found in FASTA"
WARNING: GCA_051903775.2 438618at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438643at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438643at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438643at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438643at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438657at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438657at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438657at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 438657at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43911at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43911at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43911at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 43911at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439224at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439224at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439224at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439224at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439692at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439692at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439692at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 439692at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44200at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44200at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44200at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44200at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44548at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44548at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44548at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 44548at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45362at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45362at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45731at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45731at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45731at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45731at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45825at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45825at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45825at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45825at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45840at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45840at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45840at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 45840at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46153at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46153at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46153at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46651at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46651at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46651at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 46651at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47352at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47352at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47352at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47352at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47418at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47418at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47636at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47636at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47636at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47636at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47688at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47688at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47688at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 47688at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48189at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48189at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48227at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48227at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48227at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 48227at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 4975at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 4975at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 4975at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 4975at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49782at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49782at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49782at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49782at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49947at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49947at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49947at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 49947at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5001at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5001at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5001at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5001at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5172at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5172at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5172at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5172at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52190at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52190at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52215at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52215at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52215at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52215at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52296at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52296at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52296at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52296at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5238at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5238at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5238at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5238at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5247at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5247at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52535at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52535at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52535at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52535at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52597at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52597at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52597at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52597at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52905at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52905at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52905at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 52905at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5318at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5318at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5318at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5318at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5324at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5324at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5324at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5324at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53490at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53490at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53522at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53522at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53522at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53522at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53575at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53575at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53575at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53575at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53656at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53656at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53694at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53694at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53846at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 53846at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54375at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54375at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54648at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54648at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54648at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54648at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5478at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5478at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5478at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5478at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54946at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54946at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54946at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54946at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54977at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 54977at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55175at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55175at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55175at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55175at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5554at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5554at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55710at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55710at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55710at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 55710at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5600at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5600at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5600at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56027at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56027at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5602at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5602at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56135at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56135at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56135at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 56135at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57521at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57521at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57521at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57521at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57667at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57667at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57701at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57701at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57701at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 57701at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5789at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5789at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5789at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 5789at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58843at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58843at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58843at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58843at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58849at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58849at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58849at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 58849at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59553at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59553at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59553at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59553at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59593at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59593at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59593at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 59593at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 60138at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 60138at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 60138at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 60138at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 612at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 612at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 61325at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 61325at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 61325at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 61325at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 6260at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 6260at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 6260at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 6260at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 75895at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 75895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 75895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 75895at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 7734at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 7734at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 7734at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 7734at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 79228at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 79228at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 79228at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 79228at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8620at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8620at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8620at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8620at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8803at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8803at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8803at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 8803at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9010at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9010at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9010at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9010at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 91618at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 91618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 91618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 91618at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 94870at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 94870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 94870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 94870at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9538at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 9538at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 95887at8457: intron extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 95887at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 95887at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
WARNING: GCA_051903775.2 95887at8457: flank extraction failed: "Sequence 'CM122935.2' not found in FASTA"
Processed GCA_051903775.2: 5315 valid Single orthologs for summaries/features
Processed GCA_051940855.1: 5393 valid Single orthologs for summaries/features
Processed GCA_052054735.1: 4390 valid Single orthologs for summaries/features
Processed GCA_053572275.1: 5584 valid Single orthologs for summaries/features
Processed GCA_053574215.1: 5685 valid Single orthologs for summaries/features
Processed GCA_054791115.1: 4964 valid Single orthologs for summaries/features
Processed GCA_055504975.1: 5703 valid Single orthologs for summaries/features
Processed GCA_055773955.1: 5727 valid Single orthologs for summaries/features
Processed GCA_055824665.1: 5190 valid Single orthologs for summaries/features
Processed GCA_947247035.1: 5426 valid Single orthologs for summaries/features
Processed GCA_947686815.1: 5709 valid Single orthologs for summaries/features
Processed GCA_964106635.2: 0 valid Single orthologs for summaries/features
Processed GCA_964106915.2: 0 valid Single orthologs for summaries/features
Processed GCA_964188175.1: 0 valid Single orthologs for summaries/features
Processed GCA_964188305.1: 5661 valid Single orthologs for summaries/features
Processed GCA_964234825.1: 0 valid Single orthologs for summaries/features
Processed GCA_964252035.1: 5617 valid Single orthologs for summaries/features
Processed GCA_964265115.1: 5534 valid Single orthologs for summaries/features
Processed GCA_964270895.1: 0 valid Single orthologs for summaries/features
Processed GCA_964273705.1: 5457 valid Single orthologs for summaries/features
Processed GCA_964340495.1: 5513 valid Single orthologs for summaries/features
Processed GCA_964659585.1: 5521 valid Single orthologs for summaries/features
Processed GCA_965112155.1: 0 valid Single orthologs for summaries/features
Processed GCA_965113305.1: 0 valid Single orthologs for summaries/features
Processed GCA_965153305.2: 0 valid Single orthologs for summaries/features
Processed GCA_965194845.1: 5479 valid Single orthologs for summaries/features
Processed GCA_965280105.1: 5713 valid Single orthologs for summaries/features
Processed GCA_965636955.1: 5755 valid Single orthologs for summaries/features
Processed GCF_009769535.1: 4755 valid Single orthologs for summaries/features
Processed GCF_009819535.1: 5398 valid Single orthologs for summaries/features
Processed GCF_019175285.1: 5122 valid Single orthologs for summaries/features
Processed GCF_021028975.2: 5043 valid Single orthologs for summaries/features
Processed GCF_023053635.1: 5765 valid Single orthologs for summaries/features
Processed GCF_027172205.1: 0 valid Single orthologs for summaries/features
Processed GCF_027244095.1: 5766 valid Single orthologs for summaries/features
Processed GCF_028583425.1: 5717 valid Single orthologs for summaries/features
Processed GCF_028640845.1: 5356 valid Single orthologs for summaries/features
Processed GCF_029931775.1: 5742 valid Single orthologs for summaries/features
Processed GCF_030035675.1: 5659 valid Single orthologs for summaries/features
Processed GCF_031021105.1: 5158 valid Single orthologs for summaries/features
Processed GCF_032191835.1: 5486 valid Single orthologs for summaries/features
Processed GCF_035046505.1: 5526 valid Single orthologs for summaries/features
Processed GCF_035149785.1: 5625 valid Single orthologs for summaries/features
Processed GCF_035594765.1: 5652 valid Single orthologs for summaries/features
Processed GCF_037176765.1: 5535 valid Single orthologs for summaries/features
Processed GCF_049243985.1: 5676 valid Single orthologs for summaries/features
Processed GCF_051106095.1: 5739 valid Single orthologs for summaries/features
Processed GCF_951804945.1: 5766 valid Single orthologs for summaries/features
Processed GCF_963506605.1: 5630 valid Single orthologs for summaries/features
Processed GCF_964188315.1: 5755 valid Single orthologs for summaries/features
Processed GCF_964194415.1: 5441 valid Single orthologs for summaries/features

Finished building Compleasm feature TSVs.
Output directory: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features
sauropsida_odb12 rows: 6118
orthologs rows: 493173
ortholog_summary rows: 75
intron_compleasm rows: 3224709
intron_compleasm_summary rows: 75
flanks_compleasm rows: 1220013
flank_sets_compleasm rows: 3
flank_compleasm_summary rows: 225

rossoaa@ sql % python 04_build_compleasm_feature_tsvs_v4.py --genomes /Users/rossoaa/projects/genomes --manifest /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv --orthodb /Users/rossoaa/projects/genomes/records/compleasm/mb_downloads/sauropsida_odb12.2025-07-01.tar.gz --sequences-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/sequences.tsv --genomes-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/genomes.tsv --ortholog-validation /Users/rossoaa/projects/genomes/records/compleasm/records/raw_ortholog_validity/mass_predicts_dna_dynamics_with_s_punctatus_manifest_raw_ortholog_validity.tsv
Loaded 84 accessions from manifest
Resolved 84 Compleasm metadata rows
Loaded ortholog-validation entries for 84 accession/accession-root keys
Missing accession roots in Compleasm metadata: 0
Processed GCA_964252035.1: 5617 valid Single orthologs for summaries/features
Processed GCA_964188305.1: 5661 valid Single orthologs for summaries/features
Processed GCF_009819535.1: 5398 valid Single orthologs for summaries/features
Processed GCA_947686815.1: 5709 valid Single orthologs for summaries/features
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_964270895.1: genome_pk=52 sequence_id='Podarcis_filfolensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv
Processed GCA_964270895.1: 5296 valid Single orthologs for summaries/features
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_964106915.2: genome_pk=46 sequence_id='Podarcis_gaigeae|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv
Processed GCA_964106915.2: 5244 valid Single orthologs for summaries/features
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv
DEBUG GCA_965112155.1: genome_pk=56 sequence_id='Podarcis_liolepis|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv
Processed GCA_965112155.1: 5286 valid Single orthologs for summaries/features
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028726806.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028744054.1|gene:gene-ENPEP' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_964234825.1: genome_pk=49 sequence_id='Podarcis_melisellensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv
Processed GCA_964234825.1: 5030 valid Single orthologs for summaries/features
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv
DEBUG GCA_964106635.2: genome_pk=45 sequence_id='Podarcis_pityusensis|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv
Processed GCA_964106635.2: 5049 valid Single orthologs for summaries/features
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028744054.1|gene:gene-ENPEP' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028723607.1|gene:gene-SHPRH' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCF_027172205.1: genome_pk=67 sequence_id='Podarcis_raffonei|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
Processed GCF_027172205.1: 5307 valid Single orthologs for summaries/features
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028724289.1|gene:gene-RASGRP3' not found in sequences.tsv
DEBUG GCA_964188175.1: genome_pk=47 sequence_id='Podarcis_siculus|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv
Processed GCA_964188175.1: 5275 valid Single orthologs for summaries/features
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv
DEBUG GCA_965153305.2: genome_pk=58 sequence_id='Podarcis_tiliguerta|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv
Processed GCA_965153305.2: 5050 valid Single orthologs for summaries/features
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028734048.1|gene:gene-LOC114599228' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028726805.1|gene:gene-RIOX2' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028718810.1|gene:gene-GSTZ1' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028711322.1|gene:gene-DGKE' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028744562.1|gene:gene-TMEM144' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028731351.1|gene:gene-C5H3orf33' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028737233.1|gene:gene-CLUL1' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028747750.1|gene:gene-APPL2' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028724288.1|gene:gene-RASGRP3' not found in sequences.tsv
DEBUG GCA_965113305.1: genome_pk=57 sequence_id='Podarcis_vaucheri|rna-XM_028729768.1|gene:gene-RBP4' not found in sequences.tsv
Processed GCA_965113305.1: 5084 valid Single orthologs for summaries/features
Processed GCF_963506605.1: 5630 valid Single orthologs for summaries/features
Processed GCA_041722995.2: 5688 valid Single orthologs for summaries/features
Processed GCF_028640845.1: 5356 valid Single orthologs for summaries/features
Processed GCA_051312515.2: 5743 valid Single orthologs for summaries/features
Processed GCF_035594765.1: 5652 valid Single orthologs for summaries/features
Processed GCF_037176765.1: 5535 valid Single orthologs for summaries/features
Processed GCA_042257475.1: 5705 valid Single orthologs for summaries/features
Processed GCA_053572275.1: 5584 valid Single orthologs for summaries/features
Processed GCA_014337955.2: 5658 valid Single orthologs for summaries/features
Processed GCA_053574215.1: 5685 valid Single orthologs for summaries/features
Processed GCA_023653725.1: 5281 valid Single orthologs for summaries/features
Processed GCF_035149785.1: 5625 valid Single orthologs for summaries/features
Processed GCA_965280105.1: 5713 valid Single orthologs for summaries/features
Processed GCA_046524025.1: 5502 valid Single orthologs for summaries/features
Processed GCA_043643385.1: 3362 valid Single orthologs for summaries/features
Processed GCA_039797435.1: 5480 valid Single orthologs for summaries/features
Processed GCA_003400415.2: 4880 valid Single orthologs for summaries/features
Processed GCA_030015325.2: 5433 valid Single orthologs for summaries/features
Processed GCA_030412105.1: 5747 valid Single orthologs for summaries/features
Processed GCA_054791115.1: 4964 valid Single orthologs for summaries/features
Processed GCA_052054735.1: 4390 valid Single orthologs for summaries/features
Processed GCA_050231175.1: 5530 valid Single orthologs for summaries/features
Processed GCF_023053635.1: 5765 valid Single orthologs for summaries/features
Processed GCF_031021105.1: 5158 valid Single orthologs for summaries/features
Processed GCF_028583425.1: 5717 valid Single orthologs for summaries/features
Processed GCF_029931775.1: 5742 valid Single orthologs for summaries/features
Processed GCA_030440675.1: 5688 valid Single orthologs for summaries/features
Processed GCF_027244095.1: 5766 valid Single orthologs for summaries/features
Processed GCA_964340495.1: 5513 valid Single orthologs for summaries/features
Processed GCF_032191835.1: 5486 valid Single orthologs for summaries/features
Processed GCA_019473425.1: 4562 valid Single orthologs for summaries/features
Processed GCA_033807585.1: 5433 valid Single orthologs for summaries/features
Processed GCA_051529865.1: 5711 valid Single orthologs for summaries/features
Processed GCA_040285375.1: 5732 valid Single orthologs for summaries/features
Processed GCA_965194845.1: 5479 valid Single orthologs for summaries/features
Processed GCA_964265115.1: 5534 valid Single orthologs for summaries/features
Processed GCA_965636955.1: 5755 valid Single orthologs for summaries/features
Processed GCA_009733165.1: 4542 valid Single orthologs for summaries/features
Processed GCA_964273705.1: 5457 valid Single orthologs for summaries/features
Processed GCA_964659585.1: 5521 valid Single orthologs for summaries/features
Processed GCF_049243985.1: 5676 valid Single orthologs for summaries/features
Processed GCA_037367245.1: 4919 valid Single orthologs for summaries/features
Processed GCA_020142125.1: 5233 valid Single orthologs for summaries/features
Processed GCF_951804945.1: 5766 valid Single orthologs for summaries/features
Processed GCF_964188315.1: 5755 valid Single orthologs for summaries/features
Processed GCF_051106095.1: 5739 valid Single orthologs for summaries/features
Processed GCA_050613815.1: 5529 valid Single orthologs for summaries/features
Processed GCA_039707465.1: 5492 valid Single orthologs for summaries/features
Processed GCF_030035675.1: 5659 valid Single orthologs for summaries/features
Processed GCA_016801065.1: 5562 valid Single orthologs for summaries/features
Processed GCF_019175285.1: 5122 valid Single orthologs for summaries/features
Processed GCA_021292165.1: 4771 valid Single orthologs for summaries/features
Processed GCA_039880765.1: 5343 valid Single orthologs for summaries/features
Processed GCF_021028975.2: 5043 valid Single orthologs for summaries/features
Processed GCA_041380405.1: 5761 valid Single orthologs for summaries/features
Processed GCA_051473905.1: 5755 valid Single orthologs for summaries/features
Processed GCF_009769535.1: 4755 valid Single orthologs for summaries/features
Processed GCF_035046505.1: 5526 valid Single orthologs for summaries/features
Processed GCA_050042745.1: 5580 valid Single orthologs for summaries/features
Processed GCA_047301725.1: 5519 valid Single orthologs for summaries/features
Processed GCA_024294585.1: 5441 valid Single orthologs for summaries/features
Processed GCA_051940855.1: 5393 valid Single orthologs for summaries/features
Processed GCA_038048745.1: 5098 valid Single orthologs for summaries/features
Processed GCA_055504975.1: 5703 valid Single orthologs for summaries/features
Processed GCA_947247035.1: 5426 valid Single orthologs for summaries/features
Processed GCA_030867105.1: 5625 valid Single orthologs for summaries/features
Processed GCA_051903775.2: 5315 valid Single orthologs for summaries/features
Processed GCA_055773955.1: 5727 valid Single orthologs for summaries/features
Processed GCA_055824665.1: 5190 valid Single orthologs for summaries/features
Processed GCF_964194415.1: 5441 valid Single orthologs for summaries/features
Processed GCA_003113815.1: 5122 valid Single orthologs for summaries/features

Finished building Compleasm feature TSVs.
Output directory: /Users/rossoaa/projects/genomes/records/sql_tsvs/compleasm_features
sauropsida_odb12 rows: 6118
orthologs rows: 493173
ortholog_summary rows: 84
intron_compleasm rows: 3227938
intron_compleasm_summary rows: 84
flanks_compleasm rows: 1221225
flank_sets_compleasm rows: 3
flank_compleasm_summary rows: 252
zsh: command not found: --test



 [TABLE TOTAL] genome_summary: 336 row(s)
[INFO] Running post-load orphan checks...
[INFO] Database row counts after load:
  natural_history	84 rows
  genomes	84 rows
  sequences	114665 rows
  species_name_audit	39639 rows
  analysis_run	1 rows
  window_set	336 rows
  genomic_windows	39473164 rows
  gc_window_stats	39473164 rows
  sequence_summary	458524 rows
  genome_summary	336 rows

[OK] Load committed successfully.
