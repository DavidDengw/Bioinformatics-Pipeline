import os
import subprocess

# Reference Varaiables
reference = './myPythonPipeline/data/hg19.fa'
read = ['0', './myPythonPipeline/data/NA12878-20k_R1.fastq.gz','./myPythonPipeline/data/NA12878-20k_R2.fastq.gz']
sample_name = 'NA12878-20k'
output_dir = './myPythonPipeline/Output/'

#add GATK to path
os.environ["PATH"] += os.pathsep + "/usr/local/bin/gatk-4.6.0.0"


def run_command(command):
    """Run a shell command and handle errors."""
    try: 
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        print(f"Output: {result.stdout} Done")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        raise Exception(f"Command failed: {command}")


def fastqc(input_fastq, output_dir):
    """Construct shell for generating QC statistics for FASTQ"""
    command = f"fastqc {input_fastq} -o {output_dir}"
    print(f"Running FastQC for {input_fastq}")
    run_command(command)


def bwa_align(reference, read1, read2, output_sam, rg_id, rg_pl, rg_sm):
    """Construct shell command for bwa alignment"""
    read_group = f"@RG\\tID:{rg_id}\\tPL:{rg_pl}\\tSM:{rg_sm}"
    command = f"bwa mem -R '{read_group}' {reference} {read1} {read2} > {output_sam}"
    print("Running BWA alignment")
    run_command(command)


def sam_bam_sort(input_sam, output_bam):
    """Construct shell command for sorting sam_file while converting to bam_file"""
    command = f" samtools sort {input_sam} > {output_bam}"
    print(f"Sorting and Converting {input_sam} to {output_bam}")
    run_command(command)


def bam_qc(input_bam, output_flagstats):
    """Construct shell command for generating QC stats for bam_file"""
    command = f"samtools flagstat {input_bam} > {output_flagstats}"
    print(f"Running BAM QC for {input_bam}")
    run_command(command)


def variant_call(reference, input_bam, output_vcf):
    """Constructing shell command for Generating Pile Ups and Variant Calling"""
    command = f"gatk HaplotypeCaller -R {reference} -I {input_bam} -O {output_vcf}"
    print(f"Calling variants for {input_bam} using GATK HaplotypeCaller")
    run_command(command)


def filter_false_positives(input_vcf, output_vcf):
    """Constructing shell command for applying false-positive filter(quality score less than 50)"""
    command = f"bcftools filter -s Q50 -e 'QUAL<50' {input_vcf} -o {output_vcf}"
    print(f"Filtering false positives in {input_vcf}")
    run_command(command)



def main():

    if not os.path.exists(output_dir): # Create output directory
        os.makedirs(output_dir)  

    #Variables of file paths
    sam_file = os.path.join(output_dir, f"{sample_name}.sam")
    sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
    vcf_file = os.path.join(output_dir, f"{sample_name}.vcf")
    filtered_vcf_file = os.path.join(output_dir, f"{sample_name}.filtered.vcf")        
    qc_file = os.path.join(output_dir, f"{sample_name}.qc.txt")
    fastqc_dir = os.path.join(output_dir, "fastqc")
        

    #FastQC on the reads
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)
    fastqc(read[1], fastqc_dir)
    fastqc(read[2], fastqc_dir)

    #Read group information
    rg_id = "1"
    rg_pl = "illumina"
    rg_sm = sample_name

    #Align readss to reference sequence
    bwa_align(reference, read[1], read[2], sam_file, rg_id, rg_pl, rg_sm)

    #Sorts sam_file while converting to sorted_bam_file
    sam_bam_sort(sam_file, sorted_bam_file)

    #BAM Quality Control Statistics using Flagstat
    bam_qc(sorted_bam_file, qc_file)

    #Variant Calling using GATK
    variant_call(reference, sorted_bam_file, vcf_file)

    #Applying False-Positive Filter for quality score less than 50
    filter_false_positives(vcf_file, filtered_vcf_file)

    #Print results
    print(f"Results for {sample_name}:")
    with open(filtered_vcf_file, 'r') as vcf:
        print(vcf.read())
    with open(qc_file, 'r') as qc:
        print(qc.read())


main()