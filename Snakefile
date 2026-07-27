import os
import logging

logger = logging.getLogger(__name__)

configfile: "config/config.yaml"

include: "utils.smk"

rule all:
    input:
        bam = expand(os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.sorted.markdup.bam"), sample=samples),
        quant = expand(os.path.join(OUT_DIR, "salmon_counts", "{sample}", "quant.sf"), sample=samples),
        multiqc_report = os.path.join(OUT_DIR, "qc/multiqc_report.html")

rule fastqc:
    input:
        r1 = lambda wildcards: UNITS[wildcards.sample][wildcards.lane]["R1"],
        r2 = lambda wildcards: UNITS[wildcards.sample][wildcards.lane]["R2"]
    output:
        html_r1 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_{lane}_R1_fastqc.html"),
        html_r2 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_{lane}_R2_fastqc.html"),
        zip_r1 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_{lane}_R1_fastqc.zip"),
        zip_r2 = os.path.join(OUT_DIR, "qc", "fastqc", "{sample}_{lane}_R2_fastqc.zip")  
    params:
        outdir=os.path.join(OUT_DIR, "qc", "fastqc")
    threads: 
        config["threads"]["fastqc"]
    log:
        os.path.join(LOG_DIR, "fastqc/{sample}_{lane}.log")
    benchmark: 
        "benchmarks/{sample}_{lane}.fastqc.tsv"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        """

rule trim:
    input:
        r1 = lambda wildcards: UNITS[wildcards.sample][wildcards.lane]["R1"],
        r2 = lambda wildcards: UNITS[wildcards.sample][wildcards.lane]["R2"]
    output:
        trimmed_r1 = temp(os.path.join(OUT_DIR, "{sample}_{lane}_R1.trimmed.fastq.gz")),
        trimmed_r2 = temp(os.path.join(OUT_DIR, "{sample}_{lane}_R2.trimmed.fastq.gz")),
        html = os.path.join(OUT_DIR, "qc", "fastp", "{sample}_{lane}_fastp.html"),
        json = os.path.join(OUT_DIR, "qc", "fastp", "{sample}_{lane}_fastp.json")
    threads: 
        config["threads"]["fastp"]
    log:
        os.path.join(LOG_DIR, "trim/{sample}_{lane}.log")
    benchmark:  
        "benchmarks/{sample}_{lane}.fastp.tsv"
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
        mkdir -p {output.index}

        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.index} \
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
        log = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.Log.final.out"),
        sj = os.path.join(OUT_DIR, "star_align", "{sample}", "{sample}.SJ.out.tab")
    params:
        outdir = os.path.join(OUT_DIR, "star_align/{sample}"),
        prefix = os.path.join(OUT_DIR, "star_align/{sample}/{sample}."),
        r1 = lambda wildcards, input: ",".join(input.r1),
        r2 = lambda wildcards, input: ",".join(input.r2)
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
        --readFilesIn {params.r1} {params.r2} \
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
        rustqc = directory(os.path.join(OUT_DIR, "qc", "rustqc/{sample}"))
    params:
        prefix = os.path.join(OUT_DIR, "qc/rustqc/{sample}"),
        yaml = config['rustqc']['params']
    threads:
        config['threads']['rustqc']
    log:
        os.path.join(LOG_DIR, "rustqc/{sample}.log")
    benchmark: "benchmarks/{sample}.rustqc.tsv"
    shell:
        """
        rustqc rna {input.bam} --gtf {input.gtf} -p -o {params.prefix} --sample-name {wildcards.sample} \
        -t {threads} -q --flat-output -c {params.yaml} > {log} 2>&1
        """  

rule multiqc:
    input: 
        files = get_multiqc_input
    output:
        report = os.path.join(OUT_DIR, "qc/multiqc_report.html")
    params:
        indir = OUT_DIR,
        outdir = os.path.join(OUT_DIR, "qc")
    log:
        os.path.join(LOG_DIR, "multiqc.log")
    shell:
        """
        multiqc {params.indir} -o {params.outdir} -n multiqc_report.html -f > {log} 2>&1
        """
