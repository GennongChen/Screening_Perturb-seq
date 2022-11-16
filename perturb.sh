#merge sample fastq & index file
mkdir -p /home/chengennong/project/ylp/Science2022/perturb/raw_fq
for raw_id in `cat /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|cut -f1|grep -v Run`
do 
    new_id=`grep $raw_id /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|cut -f2`
    #merge R1
    echo ${new_id}_S1_L001_R1_001.fastq.gz && ls /home/chengennong/ncbi/public/data/SRP350148|grep $raw_id|grep _1
    ls /home/chengennong/ncbi/public/data/SRP350148|grep $raw_id|grep _1| \
        xargs -I {} cat /home/chengennong/ncbi/public/data/SRP350148/{} >> \
        /home/chengennong/project/ylp/Science2022/perturb/raw_fq/${new_id}_S1_L001_R1_001.fastq.gz
    #merge R2
    echo ${new_id}_S1_L001_R2_001.fastq.gz && ls /home/chengennong/ncbi/public/data/SRP350148|grep $raw_id|grep _2
    ls /home/chengennong/ncbi/public/data/SRP350148|grep $raw_id|grep _2| \
        xargs -I {} cat /home/chengennong/ncbi/public/data/SRP350148/{} >> \
        /home/chengennong/project/ylp/Science2022/perturb/raw_fq/${new_id}_S1_L001_R2_001.fastq.gz
done

#count
##create sgRNA reference
/home/chengennong/anaconda3/envs/R/bin/Rscript /home/chengennong/code-manual/vscode/ylp/Science2022/script_dir/mk_sgrna_feature_ref_csv.R

##create library.csv & scripts
fastqs_dir=/home/chengennong/project/ylp/Science2022/perturb/raw_fq
mRNA_ref_dir=/home/chengennong/project/single_cell_fyz/refdata-gex-GRCh38-2020-A
counts_dir=/home/chengennong/project/ylp/Science2022/perturb/counts
feature_ref=/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/feature_ref.csv
mkdir -p $counts_dir /home/chengennong/project/ylp/Science2022/perturb/scripts
for new_id in `cat /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|grep mRNA|cut -f2`
do 
    sample=`echo $new_id|cut -f2,3 -d_`
    guide_id=`grep $sample /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|grep guide_$sample|cut -f2`
    lib_strategy1=`grep $new_id /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|cut -f3`
    lib_strategy2=`grep $guide_id /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|cut -f3`
    sample_library_csv=/home/chengennong/project/ylp/Science2022/perturb/scripts/${sample}.csv
    echo "fastqs,sample,library_type" > $sample_library_csv
    echo ${fastqs_dir},${new_id},${lib_strategy1} >> $sample_library_csv
    echo ${fastqs_dir},${guide_id},${lib_strategy2} >> $sample_library_csv
    cat << EOF > /home/chengennong/project/ylp/Science2022/perturb/scripts/${sample}.sh
    cd $counts_dir
    /home/chengennong/biosoft/single_cell/cellranger-7.0.0/cellranger count --id=${sample} \
    --localcores=18 \
    --libraries=$sample_library_csv \
    --transcriptome=$mRNA_ref_dir \
    --feature-ref=$feature_ref
EOF
done


##qsub mRNA
cd /home/chengennong/project/ylp/Science2022/perturb/scripts/
for i in `ls *sh`
do qsub -l nodes=1:ppn=18 $i
done


#aggr
##create library.csv & scripts
counts_dir=/home/chengennong/project/ylp/Science2022/perturb/counts
for condition in nostim stim 
do
    out_dir=/home/chengennong/project/ylp/Science2022/perturb/aggr/$condition
    mkdir -p $out_dir
    aggr_library_csv=$out_dir/aggr_library.csv
    echo "sample_id,molecule_h5" > $aggr_library_csv
    for new_id in `cat /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|grep mRNA|grep _${condition}_|cut -f2`
    do 
        sample_id=`echo $new_id|cut -f2,3 -d_`
        echo ${sample_id},${counts_dir}/${sample_id}/outs/molecule_info.h5
    done  >> $aggr_library_csv
    cat << EOF > /home/chengennong/project/ylp/Science2022/perturb/scripts/${condition}_aggr.sh
    cd /home/chengennong/project/ylp/Science2022/perturb/aggr
    /home/chengennong/biosoft/single_cell/cellranger-7.0.0/cellranger aggr \
        --localcores 24 \
        --id=$condition \
        --csv=$aggr_library_csv \
        --normalize=mapped
EOF
done
##qsub
for i in `ls *_aggr.sh`
do qsub -l nodes=1:ppn=18 $i
done

#donor
##souporcell/
#conda activate souporcell
#export PATH="/home/chengennong/anaconda3/envs/souporcell/bin/:$PATH" 
ref_fa=/home/chengennong/project/single_cell_fyz/refdata-gex-GRCh38-2020-A/fasta/genome.fa
for new_id in `cat /home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_sample_info.txt|grep mRNA|cut -f2`
do 
    sample_id=`echo $new_id|cut -f2,3 -d_`
    echo ${sample_id}
    bam=/home/chengennong/project/ylp/Science2022/perturb/counts/$sample_id/outs/possorted_genome_bam.bam
    bc_csv=/home/chengennong/project/ylp/Science2022/perturb/counts/$sample_id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
    output_dir_name=/home/chengennong/project/ylp/Science2022/perturb/donor/souporcell/$sample_id
    cat << EOF > /home/chengennong/project/ylp/Science2022/perturb/scripts/${sample_id}_souporcell.sh
    export PATH="/home/chengennong/anaconda3/envs/souporcell/bin/:\$PATH" 
    cd /home/chengennong/project/ylp/Science2022/perturb/souporcell
    /home/chengennong/biosoft/souporcell/souporcell_pipeline.py \
        -i $bam \
        -b $bc_csv \
        -f $ref_fa \
        -t 18 \
        -o $output_dir_name \
        -k 2
EOF
done
##qsub souporcell
cd /home/chengennong/project/ylp/Science2022/perturb/scripts
for i in `ls *_souporcell.sh`
do qsub -l nodes=1:ppn=18 $i
done

##apigenome vcf-match-sample-ids
perl /home/chengennong/biosoft/apigenome/scripts/vcf-match-sample-ids \
    -vcf1 STR \
    -vcf2 STR

/home/chengennong/biosoft/single_cell/cellranger-7.0.0/cellranger multi -h