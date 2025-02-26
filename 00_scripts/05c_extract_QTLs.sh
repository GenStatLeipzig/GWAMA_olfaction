#!/bin/bash

base=/storage03/ekoenig/olfaction
smr=/home/ekoenig/Software/smr-1.3.1-linux-x86_64/smr
qtl_file=qtl_positions.txt
rm $qtl_file
touch $qtl_file

# added the GRCh37 positions for the GRCh38 regions manually
file=coloc_region_grch37.txt

while read line
do
    line=($line)

    chr=${line[1]}
    start=${line[4]}
    end=${line[5]}
    region=${line[0]}

    echo "extracting data for locus $region and region $chr:${start}-${end} (GRCh37)"

    # extract sum stats of region in text format for the different QTL inputs
    brain_eQTL=${base}/Yang/eQTL/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr${chr}
    blood_eQTL=${base}/Yang/eQTL/cage_eqtl_data/CAGE.sparse
    brain_mQTL=${base}/Yang/mQTL/Brain-mMeta/Brain-mMeta
    blood_mQTL=${base}/Yang/mQTL/EUR/EUR_chr${chr}
    brain_sQTL=${base}/Yang/sQTL/BrainMeta_cis_sqtl_summary/BrainMeta_cis_sQTL_chr${chr}

    qtls=( $brain_eQTL $blood_eQTL $brain_mQTL $blood_mQTL $brain_sQTL )

    for qtl in "${qtls[@]}"
      do
          echo "getting data for QTL $qtl"
          # write SNP list for this region
          Rscript 03_make_snplist.R ${qtl}.esi $chr $start $end ${qtl}_region${region}_snplist.txt

          # use smr tool to extract regions --query <p> min p-value to include. put 1 to get all snps
          $smr --query 1 --beqtl-summary $qtl --extract-snp ${qtl}_region${region}_snplist.txt \
            --out ${qtl}_region${region}

          # QTLs are in GRCh37, get unique variants to lift to GRCh38
          cut -f 2-3 ${qtl}_region${region}.txt | sort -n | uniq >> $qtl_file

          # compress for faster reading in R
          bgzip ${qtl}_region${region}.txt
      done

done <"$file"

sort -n $qtl_file | uniq > tmp.txt
mv tmp.txt $qtl_file

echo "done"

