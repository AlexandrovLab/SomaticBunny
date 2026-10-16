set -eu

suffix_1="mutect2"
suffix_2="strelka"
suffix_3="muse" 

echo "INFO/DP DP_${suffix_1}" > rename_mutect.txt
echo "INFO/DP DP_${suffix_2}" > rename_strelka.txt
echo "INFO/SOMATIC SOMATIC_${suffix_2}" >> rename_strelka.txt
echo "INFO/SOMATIC SOMATIC_${suffix_3}" > rename_muse.txt
bcftools annotate --rename-annots rename_mutect.txt $1 -o output_mutect.vcf
bcftools annotate --rename-annots rename_strelka.txt $2 -o output_strelka.vcf
bcftools annotate --rename-annots rename_muse.txt $3 -o output_muse.vcf

