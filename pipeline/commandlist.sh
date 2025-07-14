#!/bin/bash
subfamily = THE1C

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/01progressive_annotation_mending.py \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline \
-s "$subfamily"

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/02e_value_calculation.py \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-m /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/ \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline \
-x /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline \
-s "$subfamily"

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/03age_assignment.py \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-e /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-m /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/ \
-x /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-s "$subfamily"

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/04ltr_patch.py \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-a /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/ \
-m ltr-int \
-s "$subfamily"


python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/02e_value_calculation.py \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-m /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/ \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/internal_id \
-x /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/e_value_simple \
-S /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/LTRs_with_internal.txt

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/02e_value_calculation.py \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-m /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/ \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/ \
-x /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/e_value_simple \
-s final_LTR_sample_3000

python3 /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/03age_assignment.py \
-d /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/ \
-e /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/e_value_simple \
-r /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv \
-m /rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/ \
-x /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt \
-o /home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/teatime_simple \
-s final_LTR_sample_3000