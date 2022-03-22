#!/usr/bin/env bash

#     IGTminer - 02.collinearity.sh
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
        echo "Tool:         IGTminer collinearityNetwork"
        echo "Description:  Pairwise collinear block search and construction of NOG map"
        echo ""
        echo "Usage:        IGTminer collinearityNetwork [options]"
        echo ""
        echo "Options:"
        echo "  -h          Show this message and exit"
        echo "  -c <file>   Chromosome name path"
        echo "  -n <file>   Prefix name of output file [default IGTminerOutput]" 
        echo "  -e          E-value of BLASTP [default 1e-5]"
        echo "  -d <int>    Extent of flanking regions to search using MCscan [default 8]"
        echo "  -s <int>    Minimum number of anchors in a cluster using MCscan [default 10]"
        echo "  -t          Number of threads for blast [default 12]"
        echo ""
        echo "Author: Chen,Yongming; chen_yongming@126.com"
        echo ""
        exit 1
}

output="IGTminerOutput"
evalue=1e-5
distance=8
cluster=10
num_threads=12

while getopts "hc:n:e:d:s:t:" opt
do
    case $opt in
        h)
                Usage
                exit 1
                ;;
        c)
                chrlist=$OPTARG
                ;;
        n)
                output=$OPTARG
                ;;
        e)
                evalue=$OPTARG
                ;;
        d)
                distance=$OPTARG
                ;;
        s)
                cluster=$OPTARG
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

# header 
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
#

dec=`echo $(dirname $(readlink -f "$0"))`

if [[ ! -d $output ]];then
  mkdir $output
fi

echo `gettime`"Merge NOGC into single site..."
# get clean nucleus bed
for spec in `cat $chrlist`;do
  name=`echo $spec | gawk -vFS="_" '{print $1}'`
	(cat ${name}_toCollinearity/${spec}.bed; python ${dec}/make_bed.py ${name}_reannotation/${spec}.cluster) | bedtools sort > ${output}/${spec}.bed
done

# blastp
echo `gettime`"All-vs-all blast..."
for name1 in `cat $chrlist | cut -f1`;do
  for name2 in `cat $chrlist | cut -f1`;do
    if [[ $name1 != $name2 ]];then
      spec1=`echo $name1 | gawk -vFS="_" '{print $1}'`
      spec2=`echo $name2 | gawk -vFS="_" '{print $1}'`
      if [[ ! -f ${output}/${spec1}_${spec2}.blast ]];then
        blastp -query ${spec1}_toCollinearity/${spec1}.longest.pep.fa -db ${spec2}_toCollinearity/prot_db -evalue $evalue -num_threads $num_threads -outfmt 6 -out ${output}/${spec1}_${spec2}.blast
      fi
    fi
  done
done

# split blast file
echo `gettime`"Split blast output..."
for spec1 in `cat $chrlist`;do
  for spec2 in `cat $chrlist`;do
    if [[ $spec1 != $spec2 ]];then
      name1=`echo $spec1 | gawk -vFS="_" '{print $1}'`
      chr1=`echo $spec1 | gawk -vFS="_" '{print $2}'`
      name2=`echo $spec2 | gawk -vFS="_" '{print $1}'`
      chr2=`echo $spec2 | gawk -vFS="_" '{print $2}'`
      python ${dec}/split_hits_by_chromosome_group.py -i ${output}/${name1}_${name2}.blast \
        -l ${name1}_toCollinearity/${spec1}.bed \
        -f ${name2}_toCollinearity/${spec2}.bed \
        -m ${chr1},${chr2} > ${output}/${spec1}.${spec2}.last
    fi
  done
done

# MCscan
echo `gettime`"Collinearity block search..."

echo -en > ${output}_total_links.txt

for i in `cat $chrlist | cut -f1`;do
  for j in `cat $chrlist | cut -f1`;do
    if [[ $i != $j ]];then
      spec1=${i}
      spec2=${j}
      name1=`echo $spec1 | gawk -vFS="_" '{print $1}'`
      name2=`echo $spec2 | gawk -vFS="_" '{print $1}'`
      python ${dec}/make_blast_hits.py ${name1}_reannotation/${spec1}.cluster ${name2}_reannotation/${spec2}.cluster >> ${output}/${spec1}.${spec2}.last
      cd ${output}
      python -m jcvi.compara.catalog ortholog --dist=${distance} --min_size=${cluster} --no_strip_names ${spec1} ${spec2} --cscore=.99
      python -m jcvi.compara.synteny mcscan ${spec1}.bed ${spec1}.${spec2}.lifted.anchors --iter=1 -o ${spec1}.${spec2}.i1.blocks
      rm ${spec1}.${spec2}.last
      cat ${spec1}.${spec2}.i1.blocks | grep ${spec1} | grep ${spec2} >> ../${output}_total_links.txt
      cd ..
    fi
  done
done

# MCL
echo `gettime`"Construct collinearity network..."

${dec}/../include/mcl ${output}_total_links.txt --abc -I 2 -o ${output}_total.mcl 2>/dev/null
python ${dec}/parser_mcl.py ${chrlist} ${output}_total.mcl > ${output}_NOG_raw_matrix.txt

#
echo `gettime`"Write output file..."
Rscript ${dec}/parser_matrix.R ${output}_NOG_raw_matrix.txt ${chrlist} ${output}

rm ${output}_NOG_raw_matrix.txt ${output}_total_links.txt ${output}_total.mcl
rm -rf $output
