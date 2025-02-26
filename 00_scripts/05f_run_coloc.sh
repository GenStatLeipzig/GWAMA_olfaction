#!/bin/bash

# olfaction gwas regions
file=coloc_region_grch37.txt
base=/storage03/ekoenig/olfaction

i=$(( $SLURM_ARRAY_TASK_ID ))

line=$(head -n $i $file | tail -n 1)
line=($line)
chr=${line[1]}
region=${line[0]}

# qtl prefixes
brain_eQTL=${base}/Yang/eQTL/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr${chr}
blood_eQTL=${base}/Yang/eQTL/cage_eqtl_data/CAGE.sparse
brain_mQTL=${base}/Yang/mQTL/Brain-mMeta/Brain-mMeta
blood_mQTL=${base}/Yang/mQTL/EUR/EUR_chr${chr}
brain_sQTL=${base}/Yang/sQTL/BrainMeta_cis_sqtl_summary/BrainMeta_cis_sQTL_chr${chr}

# N data from https://yanglab.westlake.edu.cn/software/smr/
qtls=( $brain_eQTL $blood_eQTL $brain_mQTL $blood_mQTL $brain_sQTL )
N_samples=( 2865 2765 1160 3701 2865 )
labels=( brain_eQTL blood_eQTL brain_mQTL blood_mQTL brain_sQTL )

for i in {0..4}
do
    qtl=${qtls[$i]}
    N=${N_samples[$i]}
    label=${labels[$i]}
    echo "Running coloc for region $region and qtl $label with samples size $N on chr$chr"
        Rscript 05_coloc.R $region ${qtl}_region${region}.txt.gz $chr $N \
            results/coloc_region${region}_${label}.txt \
            credible_set/cs_region${region}_${label} \
            $label
done

# just do once and add header line
# chr	probe_id	probe_bp	gene
cat ${base}/Yang/*/*/*.epi | cut -f 1,4,5,6 | sort -n | uniq > epi_mapping.txt
