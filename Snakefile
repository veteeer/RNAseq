import os
import glob
import re
import logging

logger = logging.getLogger(__name__)

wildcard_constraints:
    sample = r"[^/.]+",
    lane   = r"L\d+",
    rgroup = r"R[12]"

configfile: "config/config.yaml"

INPUT_DIR= config['input']['input_folder']
OUT_DIR = config['output']
REF_DIR = config['reference']
LOG_DIR   = os.path.join(OUT_DIR, "logs")
QC_DIR    = os.path.join(OUT_DIR, "qc")


def filename_to_sample(name):
    '''
    Extracts sample name from fastq filename 
    :param: name: str, fastq filename
    :return: sample_name: str, sample name 
    '''
    basename = os.path.basename(name)

    for ext in config['input']['fastq_extensions']:
        if ext in basename: 
            basename = basename.split(ext)[0]
            break

    if '_L' in basename:
        # Find position of _L followed by 3 digits and _R1 or _R2
        match = re.search(r'_L\d{3}_R[12]', basename)
        if match: 
            return basename[:match.start()]

    if '_R1' in basename:
        return basename.split('_R1')[0]
    elif '_R2' in basename:
        return basename.split('_R2')[0]

    return basename 

def get_sample_names(input_folder): 
    '''
    Detects all sample names in input folder 
    :param: input_folder: str, given input folder
    :return: final_samples: list of str, samples
    '''
    files = []
    for ext in config['input']['fastq_extensions']:
        files += glob.glob(os.path.join(input_folder, f'*{ext}'))

    sample_names = [filename_to_sample(name) for name in files]  

    # find sample names that have at least 2 fastq 
    sample_counts = {}
    for s in sample_names: 
        if s in sample_counts: 
            sample_counts[s] += 1
        else: 
            sample_counts[s] = 1

    unpaired = [s for s in sample_counts if sample_counts[s] == 1]
    if len(unpaired) > 0:
        logger.warning(f"Unpaired data will not be processed: {sorted(unpaired)}")
    
    final_samples = [s for s in sample_counts if sample_counts[s] > 1]
    return final_samples


def get_fastq_input(wildcards, rgroup):
    '''
    Gets R1 or R2 fastq file for the sample 
    :param: wildcards: Snakemake wildcards
    :param: input_folder: str, folder of input data
    :param: rgroup: str, R1 or R2
    :return: str, fastq file
    '''
    sample = wildcards.sample
    all_files = []
    for ext in config['input']['fastq_extensions']:
        pattern = os.path.join(INPUT_DIR, f'{sample}*{rgroup}*{ext}')
        all_files += glob.glob(pattern)
    return sorted(all_files)[0]

def get_star_input(wildcards, rgroup):
    if config['to_trim']:
        return os.path.join(OUT_DIR, f'{wildcards.sample}_{rgroup}.trimmed.fastq.gz')
    else:
        return get_fastq_input(wildcards, rgroup)

samples = get_sample_names(INPUT_DIR)
logger.info(f"Detected {len(samples)} paired samples: {samples}")

rule all:
    input:
        qc = expand(os.path.join(OUT_DIR, "qc","{sample}_{rgroup}_fastqc.{ext}"), sample=samples, rgroup=['R1', 'R2'], ext=['html', 'zip']),
        quant = expand(os.path.join(OUT_DIR, "{sample}/quant.sf"), sample=samples),
        distr = expand(os.path.join(OUT_DIR, "qc", "{sample}.read_distribution.txt"), sample=samples),
        cov = expand(os.path.join(OUT_DIR, "qc", "{sample}.geneBodyCoverage.txt"), sample=samples)

rule fastqc:
    input:
        r1 = lambda wildcards: get_fastq_input(wildcards, 'R1'),
        r2 = lambda wildcards: get_fastq_input(wildcards, 'R2')
    output:
        html_r1 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_R1_fastqc.html"),
        html_r2 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_R2_fastqc.html"),
        zip_r1 = temp(os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_R1_fastqc.zip")),
        zip_r2 = temp(os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_R2_fastqc.zip"))  
    params:
        outdir=os.path.join(OUT_DIR, "qc", "fastqc")
    threads: 
        config["threads"]["fastqc"]
    log:
        os.path.join(LOG_DIR, "fastqc/{sample}.log")
    benchmark: 
        "benchmarks/{sample}.fastqc.tsv"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        """

rule trim:
    input:
        r1 = lambda wildcards: get_fastq_input(wildcards, 'R1'),
        r2 = lambda wildcards: get_fastq_input(wildcards, 'R2')
    output:
        trimmed_r1 = temp(os.path.join(OUT_DIR, "{sample}_R1.trimmed.fastq.gz")),
        trimmed_r2 = temp(os.path.join(OUT_DIR, "{sample}_R2.trimmed.fastq.gz")),
        html = os.path.join(OUT_DIR, "qc", "fastp", "{sample}_fastp.html"),
        json = os.path.join(OUT_DIR, "qc", "fastp", "{sample}_fastp.json")
    threads: 
        config["threads"]["fastp"]
    log:
        os.path.join(LOG_DIR, "trim/{sample}.log")
    benchmark:  
        "benchmarks/{sample}.fastp.tsv"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.trimmed_r1} -O {output.trimmed_r2} \
        -h {output.html} -j {output.json} -x -5 -3 \
        --thread {threads} > {log} 2>&1
        """

rule star_index:
    input:
        genome = os.path.join(REF_DIR, "genome.fa"),
        annot = os.path.join(REF_DIR, "genes.gtf")
    output:
        index = directory(os.path.join(REF_DIR, 'STAR_index'))
    params:
        readLength = config['star']['read_length']
    threads:
        config['threads']['star_index']
    resources:
        mem_mb = int(config['mem_mb']['star'])
    shell:
        """
        mkdir -p {REF_DIR}/STAR_index

        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {REF_DIR} \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.annot} \
        --sjdbOverhang {params.readLength}
        """

rule star:
    input:
        r1 = lambda wildcards: get_star_input(wildcards, 'R1'),
        r2 = lambda wildcards: get_star_input(wildcards, 'R2'),
        index = os.path.join(REF_DIR, 'STAR_index')
    output:
        bam = temp(os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Aligned.out.bam")),
        transcriptome_bam = temp(os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Aligned.toTranscriptome.out.bam")),
        log = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Log.out"),
        sj = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.SJ.out.tab")
    params:
        outdir = os.path.join(OUT_DIR, "star_align/{sample}")
        prefix = os.path.join(OUT_DIR, "star_align/{sample}/{sample}.")
    threads: 
        config['threads']['star']
    resources:
        mem_mb = int(config['mem_mb']['star'])
    log:
        os.path.join(LOG_DIR, 'star/{sample}.log')
    benchmark:  "benchmarks/{sample}.star.tsv"
    shell:
        """
        mkdir -p {params.outdir}

        STAR --runMode alignReads \
        --quantMode TranscriptomeSAM \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix  {params.prefix} \
        --outSAMtype BAM Unsorted > {log} 2>&1 
        """

rule samtools:
    input:
        bam = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Aligned.out.bam")
    output:
        sorted_bam = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.sorted.markdup.bam"),
        sorted_bai = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.sorted.markdup.bam.bai")
    params:
        prefix = os.path.join(OUT_DIR, "star_align", "{sample}/{sample}")
    threads:
        config['threads']['samtools']
    benchmark: "benchmarks/{sample}.samtools.tsv"
    shell:
        """
        samtools collate -@ {threads} -Ou {input.bam} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {params.prefix}.sorted.markdup.bam
        samtools index -@ {threads} {params.prefix}.sorted.markdup.bam
        """

rule salmon:
    input:
        bam = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Aligned.toTranscriptome.out.bam"),
        transcriptome = os.path.join(REF_DIR, "transcriptome.fa")
    output:
        quant = os.path.join(OUT_DIR, "salmon_counts", "{sample}/quant.sf")
    params:
        outdir = os.path.join(OUT_DIR, "salmon_counts", "{sample}")
    threads:
        config['threads']['salmon']        
    log:
        os.path.join(LOG_DIR, "salmon/{sample}.log")
    benchmark:
        "benchmarks/{sample}.salmon.tsv"
    shell:
        """
        salmon quant -a {input.bam} -t {input.transcriptome} -o {params.outdir} \
        --seqBias --gcBias --posBias \
        -p {threads} --quiet 2> {log}
        """

# QC 
rule rustqc:
    input:
        bam = os.path.join(OUT_DIR,  "star_align", "{sample}", "{sample}.sorted.markdup.bam"),
        gtf = os.path.join(REF_DIR, "genes.gtf")
    output:
        rseqc = directory(os.path.join(OUT_DIR, "qc", "rustqc/{sample}"))
    params:
        prefix = os.path.join(OUR_DIR, "qc/rustqc/{sample}")
        yaml = config['rustqc']['params']
    threads:
        config['threads']['rustqc']
    log:
        os.path.join(LOG_DIR, "rustqc/{sample.log}")
    benchmark: "benchmarks/{sample}.rustqc.tsv"
    shell:
        """
        rustqc {input.bam} --gtf {input.gtf} -p -o {params.prefix} --sample-name {wildcards.sample} \
        -t {threads} -q --flat-output -c {params.yaml} 2> {log}
        """  

rule multiqc:
    input:
        fastqc = expand(os.path.join(QC_DIR, "fastqc","{sample}_{rgroup}_fastqc.zip"),sample=samples, rgroup=['R1', 'R2']),
        fastp = expand(os.path.join(QC_DIR, "fastp","{sample}_fastp.json"),nsample=samples) if config['to_trim'] else [],
        rustqc = os.path.join(QC_DIR, "rustqc"),
        star = expand(os.path.join(OUT_DIR, "star_align", "{sample}/{sample}.Log.final.out"), sample=samples),
        salmon = expand(os.path.join(OUT_DIR, "{sample}/quant.sf"), sample=samples)
    output:
        report = os.path.join(QC_DIR, "multiqc_report.html")
    params:
        indirs = f"{OUT_DIR}",
        outdir = QC_DIR
    log:
        os.path.join(LOG_DIR, "multiqc.log")
    shell:
        """
        multiqc {params.indirs} -o {params.outdir} -n multiqc_report.html -f > {log} 2>&1
        """
