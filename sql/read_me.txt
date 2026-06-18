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
