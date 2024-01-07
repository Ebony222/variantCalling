# directories
reads="/home/ebony/Vaccitech/variants/reads"
data="/home/ebony/Vaccitech/variants/data"
aligned_reads="/home/ebony/Vaccitech/variants/aligned_reads"
results="/home/ebony/Vaccitech/variants/results"
ref="/home/ebony/Vaccitech/variants/refGenome/refGENOME.fasta"
fastqc="/home/ebony/variantCalling_GATK/FastQC/fastqc"



echo "Checking QC using fastQC"
# ----------------
# STEP 1: using fastQC to check the quality of the reads
# -----------------
${fastqc} ${reads}/illumina_R1.fq.gz -o ${reads}/
${fastqc} ${reads}/illumina_R2.fq.gz -o ${reads}/
echo "fastQC finished running!"



echo "adapter trimming with trimmomatic"
# ----------------
# STEP 2: Removal of hypothetical adapter sequence with trimmomatic
# ----------------
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 ${reads}/illumina_R1.fq.gz ${reads}/illumina_R2.fq.gz ${data}/illumina_trimmed_R1.paired.fq.gz ${data}/illumina_trimmed_R1.unpaired.R1.fq.gz ${data}/illumina_trimmed_R2.paired.fq.gz ${data}/illumin_trimmed_R2.unpaired.fz.gz TRAILING:10 -phred33
echo "Trimmomatic finished running!"




echo "Checking QC of trimmed reads"
# ------------------
# STEP 3: using FastQC to check the quality of trimmmed reads
# -------------------
{fastqc} ${data}/illumina_trimmed_R1.paired.fq.gz -o ${data}/
{fastqc} ${data}/illumina_trimmed_R2.paired.fq.gz -o ${data}/
echo "fastQC finished running!"




echo "creating ref index"
# -------------------------
# STEP 4a: using BWA for alignment
# ------------------------
# Create the ref index
bwa index ${ref}




echo "alignment of reads with bwa"
# -------------------------
# STEP 4b: Alignment of reads
# ------------------------
bwa mem -t 4 -R "@RG\tID:illumina\tPL:ILLUMINA\\tSM:illumina" ${ref} ${reads}/illumina_R1.fq.gz ${reads}/illumina_R2.fq.gz > ${aligned_reads}/illumina.paired.sam
echo "alignment completed!"




echo "load freebayes"
# ----------------------
# STEP 5: using freebayes to identify variants
# ----------------------
 module load freebayes



echo "conversion of SAM to BAM file"
# ----------------------
# STEP 5a: convert sam to bam file
# ----------------------
samtools view -b -o ${aligned_reads}/illumina.bam ${aligned_reads}/illumina.paired.sam



echo "sorting reads position in BAM file"
# ----------------------
# STEP 5b: sorting reads position
# ----------------------
samtools sort -o ${aligned_reads}/illumina_sorted.bam ${aligned_reads}/illumina.bam




echo "creating BAM index file"
# ----------------------
# STEP 5c: sorting reads position
# ----------------------
bamtools index -in ${aligned_reads}/illumina_sorted.bam



echo "identifying variants with freebayes"
# ----------------------
# STEP 5d: identify variants
# ----------------------
freebayes -f  ${ref}  -b ${aligned_reads}/illumina_sorted.bam > ${results}/illumina.vcf
echo "variant calling completed!"




echo "compress vcf file"
# ----------------------
# STEP 5e: compress vcf with bgzip
# ----------------------
gzip ${results}/illumina.vcf 




echo "create tabix of vcf file"
# download tabix
# sudo apt-get install tabix
tabix -p vcf ${results}/illumina.vcf.gz



echo "download rtg tools"
# download rtg tools
# wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip
# unzip rtg-tools-3.12.1-linux-x64.zip



echo "statistics of VCF" 
rtg-tools-3.12.1/rtg vcfstats ${results}/illumina.vcf.gz

