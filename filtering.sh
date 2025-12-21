#!/bin/bash
# Post-EVC for consensus calling and filtering 

######################
## ACCEPT ARGUMENTS ##
######################

SAMPLE_NAME=$1
MUTECT2_VCF=$2
MUSE_VCF=$3
STRELKA_SNV=$4
STRELKA_INDEL=$5
SAGE_VCF=$6
TUMOR_BAM=$7
DATABASE_PATH=$8
REF_FASTA=$9
tlod=10
sevs=13

#####################
## Copy over files ##
#####################
# Starting the party ...
WORK_DIR=$PWD

mkdir -p snvs_filtered
mkdir -p indels_filtered
mkdir -p special_mutations

# Decompress if needed and copy files
if [[ $STRELKA_INDEL == *.gz ]]; then
    gunzip -c $STRELKA_INDEL > ${SAMPLE_NAME}_strelka_indel.vcf
else
    cp $STRELKA_INDEL ${SAMPLE_NAME}_strelka_indel.vcf
fi

if [[ $STRELKA_SNV == *.gz ]]; then
    gunzip -c $STRELKA_SNV > ${SAMPLE_NAME}_strelka_snv.vcf
else
    cp $STRELKA_SNV ${SAMPLE_NAME}_strelka_snv.vcf
fi

if [[ $SAGE_VCF == *.gz ]]; then
    gunzip -c $SAGE_VCF > ${SAMPLE_NAME}_sage.vcf
else
    cp $SAGE_VCF ${SAMPLE_NAME}_sage.vcf
fi

cp $MUTECT2_VCF ${SAMPLE_NAME}_mutect.vcf
cp $MUSE_VCF ${SAMPLE_NAME}_muse_snv.vcf

#######################
## Fix Mutect2 files ##
#######################
# Remove multiallelic mutations and seperate biallelic mutations
# Seperate snv and indels
echo Fixing Mutect2 files...

f="${SAMPLE_NAME}_mutect.vcf"
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
final_file=${SAMPLE_NAME}_mutect_snv.vcf


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

rm ${swapped_file} ${biallelic_file} ${separated_snvs_file}

####################
## Fix MuSE files ##
####################
# Seperate biallelic mutations to specical mutations and swap tumor normal score columns
# No additional filtering for WGS based on MuSE2 recommandation from GitHub
echo Fixing MuSE files...

f="${SAMPLE_NAME}_muse_snv.vcf"
cat $f|awk 'length($5)>1'> special_mutations/${SAMPLE_NAME}_muse_snvMulti.txt
cat $f|awk '{OFS="\t";print $1,$2,".",$4,substr($5,1,1),$6,$7,$8,$9,$11,$10}'> temp_muse
mv temp_muse $f

#######################
## Fix SAGE files ##
#######################
# Seperate biallelic mutations
# Seperate snv and indels
echo Fixing SAGE files...

f="${SAMPLE_NAME}_sage.vcf"
biallelic_file=${SAMPLE_NAME}_sage.biallelic.vcf
indel_file=${SAMPLE_NAME}_sage_indel.vcf
snv_file=${SAMPLE_NAME}_sage_snv.vcf
separated_snvs_file=${SAMPLE_NAME}_sage.snv.separated.vcf

cp $f ${biallelic_file}

#change this if additional steps are added
final_vcf=${separated_snvs_file}
final_file=${SAMPLE_NAME}_sage_snv.vcf

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

rm ${biallelic_file} ${separated_snvs_file}

############################
## Collect filtered files ##
############################
#Select only the PASS-ed mutations
# Process indels
for f in *_indel.vcf; do 
    cat $f | grep PASS | grep -v "#" > indels_filtered/$f
done

# Process SNVs
for f in *_snv.vcf; do 
    cat $f | grep PASS | grep -v "#" > snvs_filtered/$f
done

############################
## Merge and Annotate SNVs ##
############################
echo "Merging and annotating SNVs..."
cd snvs_filtered
mkdir -p 2outof4 mutect_snvs muse_snvs strelka_snvs sage_snvs

# Merge SNVs - 2 out of 4 callers
cat ${SAMPLE_NAME}_* | cut -f1-5 | sort | uniq -c | awk '$1>1' | \
    awk '{OFS="\t";print $2,$3,$4,$5,$6}' | \
    awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t' | \
    grep -v "_" | grep -v chrM > 2outof4/${SAMPLE_NAME}_2outof4.vcf

# Move individual caller files
mv *mutect_snv.vcf mutect_snvs/
mv *muse_snv.vcf muse_snvs/
mv *strelka_snv.vcf strelka_snvs/
mv *sage_snv.vcf sage_snvs/

cd 2outof4
mkdir -p tmp

echo "Running bseq..."
# Add bseq header and run bias filter
cat ${DATABASE_PATH}/EVC_nextflow/bseq_header ${SAMPLE_NAME}_2outof4.vcf > ${SAMPLE_NAME}_2outof4.vcf.tmp
${DATABASE_PATH}/EVC_nextflow/DKFZBiasFilter/scripts/biasFilter.py \
    ${SAMPLE_NAME}_2outof4.vcf.tmp \
    $TUMOR_BAM \
    ${REF_FASTA} \
    ${SAMPLE_NAME}_2outof4_bseq.vcf \
    --tempFolder=$WORK_DIR/snvs_filtered/2outof4/tmp \
    --mapq=1 \
    --baseq=1

grep -v "#" ${SAMPLE_NAME}_2outof4_bseq.vcf > temp_bseq
mv temp_bseq ${SAMPLE_NAME}_2outof4_bseq.vcf
rm ${SAMPLE_NAME}_2outof4.vcf.tmp
rm -r tmp

echo "Annotating callers in SNVs..."
f="${SAMPLE_NAME}_2outof4_bseq.vcf"

# Annotate with each caller
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="mt"}else{$6=$6",mt"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../mutect_snvs/${SAMPLE_NAME}_mutect_snv.vcf $f > ${SAMPLE_NAME}_match1.vcf

awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="sa"}else{$6=$6",sa"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../sage_snvs/${SAMPLE_NAME}_sage_snv.vcf ${SAMPLE_NAME}_match1.vcf > ${SAMPLE_NAME}_match2.vcf

awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="st"}else{$6=$6",st"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../strelka_snvs/${SAMPLE_NAME}_strelka_snv.vcf ${SAMPLE_NAME}_match2.vcf > ${SAMPLE_NAME}_match3.vcf

awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="ms"}else{$6=$6",ms"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../muse_snvs/${SAMPLE_NAME}_muse_snv.vcf ${SAMPLE_NAME}_match3.vcf > ${SAMPLE_NAME}_match4.vcf

cat ${SAMPLE_NAME}_match4.vcf | awk '{OFS="\t";a=substr($8,3,length($8));$9=a"|"$10"|"$11"|"$12"|"$13;print $1,$2,$3,$4,$5,$6,$7,".",$9}' > ${SAMPLE_NAME}_bseq_annotated.vcf

rm *match*

# Apply filters
f="${SAMPLE_NAME}_bseq_annotated.vcf"

cat $f | \
awk -F"\t" -v tlod="$tlod" '{
    OFS=FS
    split($9,caller,"|")
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
}' > ${f}.tmp1

awk -F"\t" '{OFS=FS; if($7=="lowTLOD"&&length($6)>5) {$7="PASS";gsub("mt,", "", $6);print} else if($7=="lowSomaticEVS"&&length($6)>5) {$7="PASS";gsub("st,", "", $6);gsub(",st", "", $6);print} else if($7=="lowTLOD;lowSomaticEVS"&&length($6)==11) {$7="PASS";$6="sa,ms";print} else {print}}' ${f}.tmp1 > ${f}.tmp2

# Filter panel of normals
grep -v "#" $WORK_DIR/${SAMPLE_NAME}_mutect_snv.vcf | grep panel_of_normals > ${SAMPLE_NAME}_mutect_snv_PON.vcf || touch ${SAMPLE_NAME}_mutect_snv_PON.vcf
pon=$(cat ${SAMPLE_NAME}_mutect_snv_PON.vcf | wc -l)

if [ $pon -eq 0 ]; then
    cp ${f}.tmp2 ${f}.tmp3
else
    awk -F"\t" 'NR==FNR{a[$1$2$4$5];next} NR>FNR{if($1$2$4$5 in a){$7=="PASS"?$7="panel_of_normals":$7=$7";panel_of_normals";print $0}else{print $0}}' OFS='\t' ${SAMPLE_NAME}_mutect_snv_PON.vcf ${f}.tmp2 > ${f}.tmp3
fi

mv ${f}.tmp3 ${SAMPLE_NAME}_snv_final_annotated.vcf
rm -f *tmp* *mutect_snv_PON.vcf

# Generate final PASS file
{
    cat ${SAMPLE_NAME}_snv_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep "^#"
    
    cat ${SAMPLE_NAME}_snv_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep -v "^#" | \
    sort -k1,1V -k2,2n
} > ${SAMPLE_NAME}_PASSed.vcf


###############################
## Merge and Annotate INDELs ##
###############################
echo "Merging and annotating INDELs..."
cd $WORK_DIR/indels_filtered
mkdir -p 2outof3 mutect_indels sage_indels strelka_indels

# Merge INDELs - 2 out of 3 callers
cat ${SAMPLE_NAME}_* | cut -f1-5 | sort | uniq -c | awk '$1>1' | \
    awk '{OFS="\t";print $2,$3,$4,$5,$6}' | \
    grep -v "_" | grep -v chrM > 2outof3/${SAMPLE_NAME}_2outof3.vcf

# Add VCF columns
awk '{$6=".";$7="PASS";$8=".";$9=".";print}' OFS='\t' 2outof3/${SAMPLE_NAME}_2outof3.vcf > temp_indel
mv temp_indel 2outof3/${SAMPLE_NAME}_2outof3.vcf

# Move individual caller files
mv *mutect_indel.vcf mutect_indels/
mv *strelka_indel.vcf strelka_indels/
mv *sage_indel.vcf sage_indels/

cd 2outof3

f="${SAMPLE_NAME}_2outof3.vcf"

# Annotate with each caller
awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="mt"}else{$6=$6",mt"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../mutect_indels/${SAMPLE_NAME}_mutect_indel.vcf $f > ${SAMPLE_NAME}_match1.vcf

awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="sa"}else{$6=$6",sa"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../sage_indels/${SAMPLE_NAME}_sage_indel.vcf ${SAMPLE_NAME}_match1.vcf > ${SAMPLE_NAME}_match2.vcf

awk -F "\t" 'NR==FNR{a[$1$2$3$4$5]=$8;next} NR>FNR{if($1$2$3$4$5 in a){OFS="\t";if($6=="."){$6="st"}else{$6=$6",st"};print $0"\t"a[$1$2$3$4$5]}else{OFS="\t";print $0"\t."}}' \
    ../strelka_indels/${SAMPLE_NAME}_strelka_indel.vcf ${SAMPLE_NAME}_match2.vcf > ${SAMPLE_NAME}_match3.vcf

cat ${SAMPLE_NAME}_match3.vcf | awk '{OFS="\t";$9=".|"$10"|"$11"|"$12;print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ${SAMPLE_NAME}_indel_annotated.vcf

rm *match*

# Apply filters
f="${SAMPLE_NAME}_indel_annotated.vcf"

cat $f | \
awk -F"\t" -v tlod="$tlod" '{
    OFS=FS
    split($9,caller,"|")
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
}' > ${f}.tmp1

cat ${f}.tmp1 | awk -F"\t" '{OFS=FS; if($7=="lowTLOD"&&length($6)>5) {$7="PASS";gsub("mt,", "", $6);print} else if($7=="lowSomaticEVS"&&length($6)>5) {$7="PASS";gsub("st,", "", $6);gsub(",st", "", $6);print} else if($7=="lowTLOD;lowSomaticEVS"&&length($6)==11) {$7="PASS";$6="sa,ms";print} else {print}}' > ${f}.tmp2

# Filter panel of normals
grep -v "#" $WORK_DIR/${SAMPLE_NAME}_mutect_indel.vcf | grep panel_of_normals > ${SAMPLE_NAME}_mutect_indel_PON.vcf || touch ${SAMPLE_NAME}_mutect_indel_PON.vcf
pon=$(cat ${SAMPLE_NAME}_mutect_indel_PON.vcf | wc -l)

if [ $pon -eq 0 ]; then
    cp ${f}.tmp2 ${f}.tmp3
else
    awk -F"\t" 'NR==FNR{a[$1$2$4$5];next} NR>FNR{if($1$2$4$5 in a){$7=="PASS"?$7="panel_of_normals":$7=$7";panel_of_normals";print $0}else{print $0}}' OFS='\t' ${SAMPLE_NAME}_mutect_indel_PON.vcf ${f}.tmp2 > ${f}.tmp3
fi

mv ${f}.tmp3 ${SAMPLE_NAME}_indel_final_annotated.vcf
rm -f *tmp* *mutect_indel_PON.vcf

# Generate final PASS file
{
    cat ${SAMPLE_NAME}_indel_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep "^#"
    
    cat ${SAMPLE_NAME}_indel_final_annotated.vcf | \
    awk '$7 == "PASS" {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
    awk -v OFS="\t" '$1=$1' | \
    grep -v "^#" | \
    sort -k1,1V -k2,2n
} > ${SAMPLE_NAME}_PASSed.vcf

echo "Processing complete for ${SAMPLE_NAME}"
