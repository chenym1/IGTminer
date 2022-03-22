#!/usr/bin/env bash

#     IGTminer - 01.reannotation.sh
#     Copyright (C) Yongming Chen
#     Contact: chen_yongming@126.com
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

set -e

if [ -z "$1" ]; then
        echo "The arguments is empty!"
        exit 1
fi

gettime() {
        echo -e `date '+%Y-%m-%d %H:%M:%S ... '`
}

Usage () {
        echo ""
        echo "Tool:         IGTminer reannotation"
        echo "Description:  Reannotation of nucleus organelle genes (NOGs)"
        echo ""
        echo "Usage:        IGTminer reannotation [options]"
        echo ""
        echo "Options:"
        echo "  -h          Show this message and exit"
        echo "  -g <file>   Gff3"
        echo "  -f <file>   Nucleus genome fasta"
        echo "  -c <file>   Chromosome name file"
        echo "  -r <file>   Custom protein sequence file"
        echo "  -m <str>    If parameter '-r' is empty, 
               please select plastid (pt) or mitochondrion (mt) [default pt]"
        echo "  -p <file>   Organelle genome fasta"
        echo "  -d <int>    Maximum distance for clustering NOGs [default 140000]"
        echo "  -n <str>    Prefix name of output file [default NOG]"
        echo "  -t <int>    Number of threads [default 12]"
        echo ""
        echo "Author: Chen,Yongming; chen_yongming@126.com"
        echo ""
        exit 1
}

prot=""
name="NOG"
num_threads=12
protraw="pt"
distance=140000

while getopts "hg:f:c:r:m:p:d:n:t:" opt
do
    case $opt in
        h)
                Usage
                exit 1
                ;;
        g)
                gff3=$OPTARG
                ;;
        f)
                fa=$OPTARG
                ;;
        c)
                chrlist=$OPTARG
                ;;
        r)
                prot=$OPTARG
                ;;
        m)
                protraw=$OPTARG
                ;;
        p)
                origin_genomefasta=$OPTARG
                ;;
        d)
                distance=$OPTARG
                ;;
        n)
                name=$OPTARG
                ;;
        t)
                num_threads=$OPTARG
                ;;
        ?)
                echo "Unknow argument!"
                exit 1
                ;;
        esac
done

#
logo () {
	echo ""
	echo "   ==============================="
	echo "  ||                             ||"
	echo "  ||          IGTminer           ||"
	echo "  ||       Version: v1.0.0       ||"
	echo "  ||                             ||"
	echo "   ==============================="
	echo ""
}
logo

gettime() {
	echo -e `date '+%Y-%m-%d %H:%M:%S ... '`
}

dec=`echo $(dirname $(readlink -f "$0"))`

if [ "${prot}"x = x ];then
  if [ "${protraw}"x = "pt"x ];then
    prot=${dec}"/pt_prot.fa"
  elif [ "${protraw}"x = "mt"x ];then
    prot=${dec}"/mt_prot.fa"
  else
    echo "Unknow -m argument"
    exit 1
  fi
fi

echo `gettime`"Build input database..."
# get seq
mkdir ${name}_out
${dec}/../include/gffread ${gff3} -g ${fa} -x ${name}_out/${name}.cds.fa
python ${dec}/get_longest_gene.py -i ${name}_out/${name}.cds.fa -f ${gff3} -n ${name}_out/${name}.longest.cds.fa
makeblastdb -in ${name}_out/${name}.longest.cds.fa -parse_seqids -hash_index -dbtype nucl -out ${name}_out/cds_db
${dec}/../include/gffread ${gff3} -g ${fa} -y ${name}_out/${name}.pep.fa
python ${dec}/get_longest_gene.py -i ${name}_out/${name}.pep.fa -f ${gff3} -n ${name}_out/${name}.longest.pep.fa
makeblastdb -in ${name}_out/${name}.longest.pep.fa -parse_seqids -hash_index -dbtype prot -out ${name}_out/prot_db
${dec}/../include/gff2bed --input=gff < ${gff3} | gawk -vOFS="\t" '{if($8=="gene")print $1,$2,$3,$4,$5,$6}' | sed 's/gene://g' > ${name}_out/${name}.gene.bed
makeblastdb -in ${fa} -dbtype nucl -title ${name} -parse_seqids -hash_index -out ${name}_out/fa_db

echo `gettime`"Reannotate nucleus organelle genes..."

# remove raw annotation
cat ${origin_genomefasta} | seqkit sliding -s 2000 -W 2000 > ${name}_out/${name}_query.fa
blastn -query ${name}_out/${name}_query.fa -db ${name}_out/cds_db -evalue 1e-5 -num_threads ${num_threads} -outfmt 6 -out ${name}_out/${name}_removepep.blast
cat ${name}_out/${name}_removepep.blast | cut -f2 | sort | uniq > ${name}_out/${name}.match.gene
python ${dec}/filter_bed.py ${name}_out/${name}.match.gene ${name}_out/${name}.gene.bed > ${name}_out/${name}.clean.bed

# split beds
for chr in `cat ${chrlist} | sed 's/,/\n/g'`;do
  python ${dec}/split_bed_by_chromosome_group.py -l ${name}_out/${name}.clean.bed -m ${chr} > ${name}_out/${name}_${chr}.bed
done

# reannotation of organelle genes
${dec}/../include/exonerate -t $origin_genomefasta -q $prot --model protein2genome --showtargetgff --maxintron 2000 > ${name}_out/${name}.exonerate.out

python ${dec}/exonerate_gff3.origin.py ${name}_out/${name}.exonerate.out > ${name}_out/${name}.reann_pt.gff3

${dec}/../include/gffread -y ${name}_out/${name}.reann_pt.pep.fa -g $origin_genomefasta  ${name}_out/${name}.reann_pt.gff3

blastp -query $prot -subject ${name}_out/${name}.reann_pt.pep.fa -evalue 1e-5 \
  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend sstart send evalue bitscore' > ${name}_out/${name}.blast.out
cat ${name}_out/${name}.blast.out | grep -v '#' | gawk '{if($5>=80)print}'| cut -f2 | sort | uniq > ${name}_out/${name}.high_coverage.gene

python ${dec}/filter_unique.py ${name}_out/${name}.reann_pt.pep.fa ${name}_out/${name}.reann_pt.gff3 ${name}_out/${name}.high_coverage.gene > ${name}_out/${name}.origin.gff3

${dec}/../include/gffread -y ${name}_out/${name}.origin.pep.fa -g $origin_genomefasta  ${name}_out/${name}.origin.gff3
sed -i 's/transcript://g' ${name}_out/${name}.origin.pep.fa
${dec}/stat_length.py ${name}_out/${name}.origin.pep.fa > ${name}_out/${name}.origin.pep.fa.fai
${dec}/../include/gffread -x ${name}_out/${name}.origin.cds.fa -g $origin_genomefasta  ${name}_out/${name}.origin.gff3
sed -i 's/transcript://g' ${name}_out/${name}.origin.cds.fa
${dec}/stat_length.py ${name}_out/${name}.origin.cds.fa > ${name}_out/${name}.origin.cds.fa.fai

# reannotation of nucleus genes
cds=${name}_out/${name}.origin.cds.fa
prot=${name}_out/${name}.origin.pep.fa

blastn -query $cds -db ${name}_out/fa_db -evalue 1e-5 -num_threads $num_threads -outfmt 6 -out ${name}_out/${name}.cds.blast

python ${dec}/filter_blast_hits.py ${fa}.fai ${name}_out/${name}.cds.blast | bedtools sort | bedtools merge | \
  gawk -vOFS="\t" '{print $1,$2,$3,$1"__"$2+1"__"$3,".","+"}' > ${name}_out/${name}.cds_blast.bed

bedtools getfasta -name -fi $fa -bed ${name}_out/${name}.cds_blast.bed -fo ${name}_out/${name}.cds_blast.fa

${dec}/../include/exonerate -t ${name}_out/${name}.cds_blast.fa -q $prot --model protein2genome --showtargetgff --maxintron 2000 --percent 75 > ${name}_out/${name}.exonerate.out

python ${dec}/exonerate_gff3.py ${name}_out/${name}.exonerate.out ${name} > ${name}_out/${name}.reann_pt.gff3

${dec}/../include/gffread -y ${name}_out/${name}.reann_pt.raw.pep.fa -g ${name}_out/${name}.cds_blast.fa ${name}_out/${name}.reann_pt.gff3

makeblastdb -in ${name}_out/${name}.reann_pt.raw.pep.fa -parse_seqids -hash_index -dbtype prot -out ${name}_out/${name}_blast_db/${name}
blastp -query $prot -db ${name}_out/${name}_blast_db/${name} -num_threads $num_threads -evalue 1e-5 -max_target_seqs 10000 \
  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend sstart send evalue bitscore' > ${name}_out/${name}.blast.out
cat ${name}_out/${name}.blast.out | grep -v '#' | gawk '{if($5>=80&&$3>=70)print}'| cut -f2 | sort | uniq > ${name}_out/${name}.high_coverage.gene

python ${dec}/filter_coding_gene.py \
  ${name}_out/${name}.reann_pt.raw.pep.fa \
  ${name}_out/${name}.reann_pt.gff3 \
  ${name}_out/${name}.high_coverage.gene \
  ${prot}.fai > ${name}_out/${name}.reann_pt.clean.gff3

python ${dec}/gff3_to_ref.py ${name}_out/${name}.reann_pt.clean.gff3 > ${name}_out/${name}.nucleus.gff3

${dec}/../include/gff2bed --input=gff < ${name}_out/${name}.nucleus.gff3 | gawk -vOFS="\t" '{if($8=="gene")print $1,$2,$3,$4,$5,$6}' > ${name}_out/${name}.nucleus.bed

for chr in `cat ${chrlist} | sed 's/,/\n/g'`;do
  python ${dec}/split_bed_by_chromosome_group.py -l ${name}_out/${name}.nucleus.bed -m ${chr} > ${name}_out/${name}_${chr}.nucleus.bed
  cat ${name}_out/${name}_${chr}.nucleus.bed ${name}_out/${name}_${chr}.bed | bedtools sort > ${name}_out/${name}_${chr}.total.bed
  Rscript ${dec}/cluster.R ${name}_${chr} ${name}_out/${name}_${chr}.total.bed ${name}_out/${name}_${chr}.nucleus.bed $distance ${name}_out/${name}_${chr}.cluster
done

echo `gettime`"Write output files..."

# move to file
mkdir ${name}_reannotation
mv ${name}_out/${name}_*.cluster ./${name}_reannotation/
mv ${name}_out/${name}.nucleus.gff3 ./${name}_reannotation/${name}.NOGs.gff3

mkdir ${name}_toCollinearity
for chr in `cat ${chrlist} | sed 's/,/\n/g'`;do
  mv ${name}_out/${name}_${chr}.bed ./${name}_toCollinearity/
done
mv ${name}_out/*.longest.pep.fa ./${name}_toCollinearity/
mv ${name}_out/prot_db* ./${name}_toCollinearity/
rm -rf ${name}_out/
