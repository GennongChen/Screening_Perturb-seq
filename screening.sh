#!/bin/bash
#install: https://www.youtube.com/watch?v=0FLWSp1tdTI
#Comparison between samples:https://www.youtube.com/watch?v=iIyP2w44ric
#count/test/pathway/plot/mle

#install mageck 0.5.9.2 
#cd /home/chengennong/biosoft && gunzip liulab-mageck-0.5.9.2.tar.gz && tar -xvf liulab-mageck-0.5.9.2.tar && mv liulab-mageck-ef7c39474ed0 mageck
#cd mageck
#python setup.py install --prefix=/home/chengennong/biosoft/mageck
#cp -rf /home/chengennong/biosoft/mageck/lib/python3.6/site-packages/mageck /home/chengennong/miniconda3/lib/python3.6/site-packages
#export PATH="/home/chengennong/biosoft/mageck/bin/:$PATH"

#test
#mageck --help
#cd /home/chengennong/biosoft/mageck/demo/demo1
#sh run.sh
#cd /home/chengennong/biosoft/mageck/demo/demo2
#sh runmageck.sh

#metatable
##screening main
cut -f4,17 /home/chengennong/ncbi/public/data/SRP319477/filereport_read_run_SRP319477_tsv.txt | \
awk 'NR>1' | \
awk 'BEGIN {FS="_"; print "Cun\tStategy\tLibrary\tSet\tDonor\tMarker\tExpression"} {if ($3=="Plasmid") {$4="Control";$5="Control"}; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' \
> /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/screening_sample_info.txt


#library
##screening main
for table_i in `ls /home/chengennong/code-manual/vscode/ylp/Science2022/sgRNA_library/*txt`
do 
    library_i=${table_i}.library
    awk 'BEGIN {FS="_"} {printf $1 "\t" $2 "\n"}' $table_i | \
    awk '{gene=$1; if (!a[$1]++ == 0) {num++; print $1 "_" $2 "\t" $2 "\t" $3}  else {num=1; print  $1 "_" $2 "\t" $2 "\t" $3}}' > ${library_i}
done  


#1 Reads were aligned to the appropriate reference library using MAGeCK version 0.5.9.2 (45) using the –trim-5 22, 23, 24, 25, 26, 28, 29, 30 argument to remove the staggered 5′ adapter.



library_dir=/home/chengennong/code-manual/vscode/ylp/Science2022/sgRNA_library
screening_meta_table=/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/screening_sample_info.txt
out_dir=/home/chengennong/project/ylp/Science2022/screening/mageck_count
mkdir -p $out_dir
for Library_i in `awk 'NR>1' $screening_meta_table |cut -f3|sort -u`
do
    grep "$Library_i" $screening_meta_table
    library=${library_dir}/`grep "$Library_i" $screening_meta_table|cut -f2,3|sed 's/\t/./'|sort -u`.txt.library
    out_prefix=${out_dir}/`grep "$Library_i" $screening_meta_table|cut -f2,3|sed 's/\t/./'|sort -u`
    sample_label=`grep "$Library_i" $screening_meta_table|cut -f3-6|sed 's/\t/_/g'|xargs echo|sed 's/ /,/g'`
    fastq=`grep "$Library_i" $screening_meta_table|cut -f1|awk '{print  "/home/chengennong/ncbi/public/data/SRP319477/" $1 ".fastq.gz"}'|xargs echo`
    mageck count \
    -l $library \
    -n $out_prefix \
    --sample-label $sample_label \
    --fastq $fastq \
    --trim-5 22,23,24,25,26,28,29,30
done


#2 Next, raw read counts across both library sets were normalized to the total read count in each sample, and each of the matching samples across two sets were merged to generate a single normalized read count table.
/home/chengennong/anaconda3/envs/R/bin/Rscript /home/chengennong/code-manual/vscode/ylp/Science2022/script_dir/nor_screen_counts.R


#3 Normalized read counts in high versus low bins were compared using mageck test with -–norm-method none, -–paired, and –-control-sgrna options, pairing samples by donor and using nontargeting sgRNAs as controls, respectively. 
#mageck test -k demo.count.txt -t L1 -c CTRL -n demo
#mageck test -k demo.count.txt -t L_r1,L_r2 -c C_r1,C_r2 -n demo

for marker in IL2 IFNG
do 
    #CRISPRa 2-marker 1-notarget-library 4-sample-each_time 2-result
    mageck test \
        -k /home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRa.CalabreseSet.normalized_count.txt \
        -t CalabreseSet_Donor1_${marker}_high,CalabreseSet_Donor2_${marker}_high \
        -c CalabreseSet_Donor1_${marker}_low,CalabreseSet_Donor2_${marker}_low \
        --norm-method none \
        --paired \
        --control-sgrna /home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRa.CalabreseSet.no_target.txt \
        -n /home/chengennong/project/ylp/Science2022/screening/mageck_test/CalabreseSet_${marker}
    #CRISPRi 2-marker 1-notarget-library 4-sample-each_time 2-result
    mageck test \
        -k /home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRi.DolcettoSet.normalized_count.txt \
        -t DolcettoSet_Donor1_${marker}_high,DolcettoSet_Donor2_${marker}_high \
        -c DolcettoSet_Donor1_${marker}_low,DolcettoSet_Donor2_${marker}_low \
        --norm-method none \
        --paired \
        --control-sgrna /home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRi.DolcettoSet.no_target.txt \
        -n /home/chengennong/project/ylp/Science2022/screening/mageck_test/DolcettoSet_${marker}
done

#4 Gene hits were classified as having a median absolute log2-fold change >0.5 and a false discovery rate (FDR) <0.05. For supplemental CD4+ screens (fig. S9), reads were aligned to the full Calabrese A and B library in a single reference file.
#statistic & draw result
for i in `ls /home/chengennong/code-manual/vscode/ylp/Science2022/script_dir/fig_script_dir/* |grep -v "#" |grep fig[12]`
do /home/chengennong/anaconda3/envs/R/bin/Rscript $i
done

#5 For the supplemental CD4+ IFN-g screen, which was sorted and sequenced as two technical replicates, normalized counts were averaged across technical replicates before analysis with mageck test.
#mageck count -l library.txt -n demo --sample-label L1,CTRL  --fastq test1.fastq test2.fastq
#stategy library set donor marker expression