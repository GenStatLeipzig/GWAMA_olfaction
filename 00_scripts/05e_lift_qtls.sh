#!/bin/bash

# 1. create fake vcf
n=$(wc -l < qtl_positions.txt)
printf '%0.s.\n' $(seq 1 $n) > id.txt
printf '%0.s100\n' $(seq 1 $n) > qual.txt
printf '%0.sA\n' $(seq 1 $n) > ref.txt
printf '%0.sT\n' $(seq 1 $n) > alt.txt

paste qtl_positions.txt id.txt ref.txt alt.txt qual.txt id.txt id.txt > tmp.vcf
cat header.vcf tmp.vcf > qtl.grch37.vcf
bcftools sort -Oz -o lift/qtl.grch37.vcf.gz qtl.grch37.vcf
tabix lift/qtl.grch37.vcf.gz

rm qtl.grch37.vcf
rm id.txt
rm qual.txt
rm ref.txt
rm alt.txt


# bcftools lifover does not work since alleles are not ref/alt, but EA/AA so bcftools complains
# therefore, use vcf-liftover
for chr in 4 5 6 11 13 14
do
    bash /home/ekoenig/Software/vcf-liftover.sh \
    /home/ekoenig/NGSpipeline/data/GRCh38/GRCh37_to_GRCh38.chain \
    lift/qtl.grch37.vcf.gz \
    lift/qtl.chr${chr}.grch38.vcf.gz \
    $chr nostop fast tmp/

    #3. create map of lifting
    gunzip -c lift/qtl.chr${chr}.grch38.vcf.gz | grep -v "^#" | cut -f 1-3 | uniq > lift/qtl.mapping.chr${chr}.vcfliftover.txt
done


#colnames: chr pos38 pos37 ref alt

