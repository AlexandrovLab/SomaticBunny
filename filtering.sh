#!/bin/bash
# Post-EVC for consensus calling and filtering 

############################
## ALL YOU NEED TO CHANGE ##
############################

WORK_DIR='/YOUR/IDEAL/OUTPUT/FOLDER'
BAM_PATH='/YOUR/RESULT/DIR/RECALIBRATE'
HOME_PATH='/YOUR/RESULT/DIR'
tlod=10
sevs=13

#####################
## Copy over files ##
#####################
# Starting the party ...
mkdir $WORK_DIR
cd $WORK_DIR

# Copy over snv files
echo Copying files over...

mkdir snvs_filtered
mkdir indels_filtered
mkdir special_mutations


for sample in $(cat ../sample.txt);
do

gunzip $HOME_PATH/STRELKA/${sample}.somatic_indels.vcf.gz
gunzip $HOME_PATH/STRELKA/${sample}.somatic_snvs.vcf.gz
gunzip $HOME_PATH/SAGE/${sample}.sage.vcf.gz

done



for sample in $(cat ../sample.txt);
do
cp $HOME_PATH/Mutect2/${sample}_mutect2_filtered.vcf ./${sample}_mutect.vcf

cp $HOME_PATH/MuSE2/${sample}.vcf ./${sample}_muse_snv.vcf

cp $HOME_PATH/STRELKA/${sample}.somatic_indels.vcf ./${sample}_strelka_indel.vcf
cp $HOME_PATH/STRELKA/${sample}.somatic_snvs.vcf ./${sample}_strelka_snv.vcf

cp $HOME_PATH/SAGE/${sample}.sage.vcf ./${sample}_sage.vcf

done

#######################
## Fix Mutect2 files ##
#######################
# Remove multiallelic mutations and seperate biallelic mutations
# Seperate snv and indels
echo Fixing Mutect2 files...

for f in *_mutect.vcf;
do
normal_name=$(grep normal_sample $f | cut -d= -f2)
tumor_name=$(grep tumor_sample $f | cut -d= -f2)
swapped_file=$(basename $f .vcf).swapped.vcf
multiallelic_file=$(basename $f .vcf).multiMutation.txt
biallelic_file=$(basename $f .vcf).biallelic.vcf
indel_file=$(basename $f .vcf)_indel.vcf
snv_file=$(basename $f .vcf)_snv.vcf
separated_snvs_file=$(basename $f .vcf).snv.separated.vcf

#change this if additional steps are added
final_vcf=${separated_snvs_file}
final_file=$(basename $f .vcf)_snv.vcf


## Swap columns
#col10 is the first column
if [ "$(grep "CHROM" $f | cut -f10)" == "${tumor_name}" ]
then
        echo Wrong order. Swapping columns...
        awk 'BEGIN{OFS="\t";}; { t = $10; $10 = $11; $11 = t; print; }' $f > ${swapped_file}
else
        echo Correct order. Copying to sample_swapped.vcf...
        #make a copy to match naming
        cp $f ${swapped_file}
fi


## Multiallele processing
#if col5 has commas, it is multiallelic - we ignore multiallelic
grep -v "#" ${swapped_file} | awk 'BEGIN{OFS="\t";}; $5 ~ /,/ { print }' > special_mutations/${multiallelic_file}
grep -v "#" ${swapped_file} | awk 'BEGIN{OFS="\t";}; ! ($5 ~ /,/) { print }' > ${biallelic_file}


#indels
awk 'BEGIN{OFS="\t";}; length($4) != length($5) { print }' ${biallelic_file} > ${indel_file}
#snvs
awk 'BEGIN{OFS="\t";}; length($4) == length($5) { print }' ${biallelic_file} > ${snv_file}

awk 'BEGIN{OFS="\t";};
        {if (length($4) > 1) {
                totalLen = length($4);
                ref=$4;
                alt=$5;
                pos=$2;
                for (i = 1; i <= totalLen; i++) {
                        new_ref=substr(ref, i, 1);
                        new_alt=substr(alt, i, 1);
                        $4 = new_ref
                        $5 = new_alt
                        $2 = pos + i - 1
                        print
                }
        }
        else {print}
}' ${snv_file} > ${separated_snvs_file}

## Generate final file
grep "#" ${swapped_file} > ${final_file}
cat ${final_vcf} >> ${final_file}
rm $f
done

rm *.swapped.vcf *.biallelic.vcf *.snv.separated.vcf


####################
## Fix MuSE files ##
####################
# Seperate biallelic mutations to specical mutations and swap tumor normal score columns
# No additional filtering for WGS based on MuSE2 recommandation from GitHub
echo Fixing MuSE files...

for f in *_muse_snv.vcf;
do fa=`basename $f _muse_snv.vcf`;
cat $f|awk 'length($5)>1'> special_mutations/${fa}_muse_snvMulti.txt;
cat $f|awk '{OFS="\t";print $1,$2,".",$4,substr($5,1,1),$6,$7,$8,$9,$11,$10}'>a;
mv a $f;
done

#######################
## Fix SAGE files ##
#######################
# Seperate biallelic mutations
# Seperate snv and indels
echo Fixing SAGE files...

for f in *_sage.vcf;
do
biallelic_file=$(basename $f .vcf).biallelic.vcf
indel_file=$(basename $f .vcf)_indel.vcf
snv_file=$(basename $f .vcf)_snv.vcf
separated_snvs_file=$(basename $f .vcf).snv.separated.vcf

cp $f ${biallelic_file}

#change this if additional steps are added
final_vcf=${separated_snvs_file}
final_file=$(basename $f .vcf)_snv.vcf

#indels
awk 'BEGIN{OFS="\t";}; length($4) != length($5) { print }' ${biallelic_file} > ${indel_file}
#snvs
awk 'BEGIN{OFS="\t";}; length($4) == length($5) { print }' ${biallelic_file} > ${snv_file}

awk 'BEGIN{OFS="\t";};
        {if (length($4) > 1) {
                totalLen = length($4);
                ref=$4;
                alt=$5;
                pos=$2;
                for (i = 1; i <= totalLen; i++) {
                        new_ref=substr(ref, i, 1);
                        new_alt=substr(alt, i, 1);
                        $4 = new_ref
                        $5 = new_alt
                        $2 = pos + i - 1
                        print
                }
        }
        else {print}
}' ${snv_file} > ${separated_snvs_file}

## Generate final file

grep "#" ${biallelic_file} > ${final_file}
cat ${final_vcf} >> ${final_file}
rm $f
done

rm *.biallelic.vcf *.snv.separated.vcf

############################
## Collect filtered files ##
############################
#Select only the PASS-ed mutations
for f in *indel.vcf; do cat $f|grep PASS|grep -v "#">indels_filtered/$f;done
for f in *snv.vcf; do cat $f|grep PASS|grep -v "#">snvs_filtered/$f;done

echo Merging INDELs..
cd indels_filtered
mkdir 2outof3 3outof3

# Select mutations called by at least 2/3 callers and filter out chromosomes rather than 1-22 and XY
for f in *_mutect_indel.vcf;do fa=`basename $f _mutect_indel.vcf`;cat ${fa}_*|cut -f1-5|sort|uniq -c|awk '$1>1'|awk '{OFS="\t";print $2,$3,$4,$5,$6}'| grep -v "_" | grep -v chrM>2outof3/${fa}_2outof3.vcf;done

for f in *_mutect_indel.vcf;do fa=`basename $f _mutect_indel.vcf`;cat ${fa}_*|cut -f1-5|sort|uniq -c|awk '$1>2'|awk '{OFS="\t";print $2,$3,$4,$5,$6}'| grep -v "_" | grep -v chrM>3outof3/${fa}_3outof3.vcf;done

mkdir mutect_indels sage_indels strelka_indels
mv *mutect_indel.vcf mutect_indels
mv *strelka_indel.vcf strelka_indels
mv *sage_indel.vcf sage_indels


echo Merging SNVs...
cd ../snvs_filtered
mkdir 2outof4 3outof4 4outof4

# Select mutations called by at least 2/3/4 callers and filter out chromosomes rather than 1-22 and XY
for f in *_mutect_snv.vcf;do fa=`basename $f _mutect_snv.vcf`;cat ${fa}_*|cut -f1-5|sort|uniq -c|awk '$1>1'|awk '{OFS="\t";print $2,$3,$4,$5,$6}'| awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t'| grep -v "_" | grep -v chrM>2outof4/${fa}_2outof4.vcf;done

for f in *_mutect_snv.vcf;do fa=`basename $f _mutect_snv.vcf`;cat ${fa}_*|cut -f1-5|sort|uniq -c|awk '$1>2'|awk '{OFS="\t";print $2,$3,$4,$5,$6}' | awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t'| grep -v "_" | grep -v chrM>3outof4/${fa}_3outof4.vcf;done

for f in *_mutect_snv.vcf;do fa=`basename $f _mutect_snv.vcf`;cat ${fa}_*|cut -f1-5|sort|uniq -c|awk '$1>3'|awk '{OFS="\t";print $2,$3,$4,$5,$6}'| awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t'| grep -v "_" | grep -v chrM>4outof4/${fa}_4outof4.vcf;done

mkdir mutect_snvs muse_snvs strelka_snvs sage_snvs
mv *muse_snv.vcf muse_snvs
mv *mutect_snv.vcf mutect_snvs
mv *strelka_snv.vcf strelka_snvs
mv *sage_snv.vcf sage_snvs


################
## Annotation ##
################
echo starting annotation in 2outof4 snv folder ...

cd $WORK_DIR/snvs_filtered/2outof4
mkdir tmp

## Run Bseq form DKFZ
echo running bseq....
# source /tscc/nfs/home/tiy002/anaconda3/etc/profile.d/conda.sh
source ~/.bashrc
conda activate dkfz

for f in *_2outof4.vcf;do
fa=`basename $f _2outof4.vcf`
cat /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/bseq_header $f > ${f}.tmp
/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/DKFZBiasFilter/scripts/biasFilter.py \
${f}.tmp \
$BAM_PATH/${fa}_tumor_recal.bam \
/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/GRCh38_ref/GRCh38.d1.vd1.fa \
${fa}_2outof4_bseq.vcf \
--tempFolder=$WORK_DIR/snvs_filtered/2outof4/tmp \
--mapq=1 \
--baseq=1;
grep -v "#" ${fa}_2outof4_bseq.vcf >a; mv a ${fa}_2outof4_bseq.vcf
rm ${f}.tmp
done

conda deactivate
rm -r tmp


## Annotate callers
## snvs
echo Annotating callers in snvs...

for f in *_2outof4_bseq.vcf;do
fa=`basename $f _2outof4_bseq.vcf`;
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="mt"}else{$6=$6",mt"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../mutect_snvs/${fa}_mutect_snv.vcf $f>${fa}_match1.vcf
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="sa"}else{$6=$6",sa"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../sage_snvs/${fa}_sage_snv.vcf ${fa}_match1.vcf>${fa}_match2.vcf
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="st"}else{$6=$6",st"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../strelka_snvs/${fa}_strelka_snv.vcf ${fa}_match2.vcf>${fa}_match3.vcf
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="ms"}else{$6=$6",ms"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../muse_snvs/${fa}_muse_snv.vcf ${fa}_match3.vcf>${fa}_match4.vcf
cat ${fa}_match4.vcf|awk '{OFS="\t";a=substr($8,3,length($8));$9=a"|"$10"|"$11"|"$12"|"$13;print $1,$2,$3,$4,$5,$6,$7,".",$9}'>${fa}_bseq_annotated.vcf
rm *match*
done


## Annotate additional filters
## SNVs
for f in *_bseq_annotated.vcf;
do
fa=`basename $f _bseq_annotated.vcf`;

## annotate lowAF, lowTLOD (mutect), lowSomaticEVS (strelka) filters

cat $f | \
awk -F"\t" -v tlod="$tlod" '{
    OFS=FS
    split($9,caller,"|")
    # Search through all fields for TLOD
    for(j=1; j<=length(caller); j++) {
        split(caller[j],t,";")
        for(i in t) {
            if(t[i] ~ /^TLOD=/) {
                split(t[i],ta,"=")
                ta[2]<tlod?($7=="PASS"?$7="lowTLOD":$7=$7";lowTLOD"):$7=$7
                break
            }
        }
    }
    print
}' | \
awk -F"\t" -v sevs="$sevs" '{
    OFS=FS
    split($9,caller,"|")
    # Search through all fields for SomaticEVS
    for(j=1; j<=length(caller); j++) {
        split(caller[j],t,";")
        for(i in t) {
            if(t[i] ~ /^SomaticEVS=/) {
                split(t[i],ta,"=")
                ta[2]<sevs?($7=="PASS"?$7="lowSomaticEVS":$7=$7";lowSomaticEVS"):$7=$7
                break
            }
        }
    }
    print
}' > ${f}.tmp1;

#save mutations called by at least 2 variant callers
awk -F"\t" '{OFS=FS; if($7=="lowTLOD"&&length($6)>5) {$7="PASS";gsub("mt,", "", $6);print} else if($7=="lowSomaticEVS"&&length($6)>5) {$7="PASS";gsub("st,", "", $6);gsub(",st", "", $6);print} else if($7=="lowTLOD;lowSomaticEVS"&&length($6)==11) {$7="PASS";$6="sa,ms";print} else {print}}' ${f}.tmp1 > ${f}.tmp2


## filter mutations annotated by mutect as "panel_of_normal"
#grep mutations annotated as "panel_of_normal" in mutect variant calling
grep -v "#" $WORK_DIR/snvs_filtered/mutect_snvs/${fa}_mutect_snv.vcf | grep panel_of_normals > ${fa}_mutect_snv_PON.vcf
pon=$(cat ${fa}_mutect_snv_PON.vcf | wc -l)

if [ $pon -eq 0 ]
then
        cp ${f}.tmp2 ${f}.tmp3
else
        awk -F"\t" 'NR==FNR{a[$1$2$4$5];next} NR>FNR{if($1$2$4$5 in a){$7=="PASS"?$7="panel_of_normals":$7=$7";panel_of_normals";print $0}else{print $0}}' OFS='\t' ${fa}_mutect_snv_PON.vcf ${f}.tmp2 > ${f}.tmp3
fi
mv ${f}.tmp3 ${fa}_snv_final_annotated.vcf
rm *tmp*
rm *mutect_snv_PON.vcf

# Cleaning up and removing low-quality samples
{
    # Print headers first
    cat ${fa}_snv_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep "^#"
    
    # Then print sorted data
    cat ${fa}_snv_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep -v "^#" | \
    sort -k1,1V -k2,2n
} > ${fa}_PASSed.vcf

done



## indels
echo Annotating callers in indels...
cd $WORK_DIR/indels_filtered/2outof3

for f in *_2outof3.vcf
do
awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t' $f > a; mv a $f;
done

for f in *_2outof3.vcf;do
fa=`basename $f _2outof3.vcf`;
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="mt"}else{$6=$6",mt"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../mutect_indels/${fa}_mutect_indel.vcf $f>${fa}_match1.vcf
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="sa"}else{$6=$6",sa"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../sage_indels/${fa}_sage_indel.vcf ${fa}_match1.vcf>${fa}_match2.vcf
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="st"}else{$6=$6",st"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' ../strelka_indels/${fa}_strelka_indel.vcf ${fa}_match2.vcf>${fa}_match3.vcf
cat ${fa}_match3.vcf|awk '{OFS="\t";$9=".|"$10"|"$11"|"$12;print $1,$2,$3,$4,$5,$6,$7,$8,$9}'>${fa}_indel_annotated.vcf
rm *match*
done

for f in *_indel_annotated.vcf;
do
fa=`basename $f _indel_annotated.vcf`;

## annotate lowAF, lowTLOD, lowSomaticEVS filters
cat $f | \
awk -F"\t" -v tlod="$tlod" '{
    OFS=FS
    split($9,caller,"|")
    # Search through all fields for TLOD
    for(j=1; j<=length(caller); j++) {
        split(caller[j],t,";")
        for(i in t) {
            if(t[i] ~ /^TLOD=/) {
                split(t[i],ta,"=")
                ta[2]<tlod?($7=="PASS"?$7="lowTLOD":$7=$7";lowTLOD"):$7=$7
                break
            }
        }
    }
    print
}' | \
awk -F"\t" -v sevs="$sevs" '{
    OFS=FS
    split($9,caller,"|")
    # Search through all fields for SomaticEVS
    for(j=1; j<=length(caller); j++) {
        split(caller[j],t,";")
        for(i in t) {
            if(t[i] ~ /^SomaticEVS=/) {
                split(t[i],ta,"=")
                ta[2]<sevs?($7=="PASS"?$7="lowSomaticEVS":$7=$7";lowSomaticEVS"):$7=$7
                break
            }
        }
    }
    print
}' > ${f}.tmp1;

#save mutations called by at least 2 variant callers
cat ${f}.tmp1|awk -F"\t" '{OFS=FS; if($7=="lowTLOD"&&length($6)>5) {$7="PASS";gsub("mt,", "", $6);print} else if($7=="lowSomaticEVS"&&length($6)>5) {$7="PASS";gsub("st,", "", $6);gsub(",st", "", $6);print} else if($7=="lowTLOD;lowSomaticEVS"&&length($6)==11) {$7="PASS";$6="sa,ms";print} else {print}}' > ${f}.tmp2

## filter mutations annotated by mutect as "panel_of_normal"
#grep mutations annotated as "panel_of_normal" in mutect variant calling
grep -v "#" ../../${fa}_mutect_indel.vcf | grep panel_of_normals > ${fa}_mutect_indel_PON.vcf
pon=$(cat ${fa}_mutect_indel_PON.vcf | wc -l)

if [ $pon -eq 0 ]
then
        cp ${f}.tmp2 ${f}.tmp3
else
        awk -F"\t" 'NR==FNR{a[$1$2$4$5];next} NR>FNR{if($1$2$4$5 in a){$7=="PASS"?$7="panel_of_normals":$7=$7";panel_of_normals";print $0}else{print $0}}' OFS='\t' ${fa}_mutect_indel_PON.vcf ${f}.tmp2 > ${f}.tmp3
fi
mv ${f}.tmp3 ${fa}_indel_final_annotated.vcf
rm *tmp*
rm *mutect_indel_PON.vcf

# Cleaning up and removing low-quality samples
{
    # Print headers first
    cat ${fa}_indel_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep "^#"
    
    # Then print sorted data
    cat ${fa}_indel_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep -v "^#" | \
    sort -k1,1V -k2,2n
} > ${fa}_PASSed.vcf

done
