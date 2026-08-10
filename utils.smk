import os
import glob
import re

INPUT_DIR= config['input']['input_folder']
OUT_DIR = config['output']
REF_DIR = config['reference']
LOG_DIR   = os.path.join(OUT_DIR, "logs")

wildcard_constraints:
    sample = r"[^/.]+",
    lane   = r"L\d+",
    rgroup = r"R[12]"

def strip_fastq_ext(path):
    """
    Returns basename without fastq extension
    """
    b = os.path.basename(path)
    for ext in config['input']['fastq_extensions']:
        if b.endswith(ext):
            return b.split(ext)[0]

def filename_to_sample(name):
    '''
    Extracts sample name, lane, read group from fastq filename 
    :param: name: str, fastq filename
    :return: sample_name: str, lane: str, read group: str
    '''
    basename = strip_fastq_ext(name)    
    # sample with multi lanes
    m = re.search(r'_(L\d{3})_(R[12])', basename)
    if m:
        return basename[:m.start()], m.group(1), m.group(2)
    
    # without lanes, set default lane to L000
    m = re.search(r'_(R[12])', basename)
    if m:
        return basename[:m.start()], 'L000', m.group(1)


def get_units(input_folder):
    '''
    Find paired fastq files and group them by sample and lane.
    :return: dict: {sample: {lane: {"R1": path, "R2": path}}}
    '''
    files = []
    for ext in config['input']['fastq_extensions']:
        files += glob.glob(os.path.join(input_folder, f'*{ext}'))

    units = {}
    for f in sorted(files):
        sample, lane, read = filename_to_sample(f)

        if sample not in units:
            units[sample] = {}
        if lane not in units[sample]:
            units[sample][lane] = {}
        units[sample][lane][read] = f

    clean = {}
    for sample, lanes in units.items():
        for lane, reads in lanes.items():
            if 'R1' in reads and 'R2' in reads:
                clean.setdefault(sample, {})[lane] = reads
            else:
                missing = {'R1', 'R2'} - set(reads)
                logger.warning(f"Incomplete pair for {sample} {lane} (missing {missing}); lane skipped")
    return clean


UNITS   = get_units(INPUT_DIR)
samples = sorted(UNITS.keys())

logger.info(f"Detected {len(samples)} samples")
for s in samples:
    logger.info(f"  {s}: lanes {sorted(UNITS[s].keys())}")


def get_star_input(wildcards, read):
    """
    Return all R1 or R2 fastq files for sample.
    If trimming enabled in configfile, returns trimmed reads. Otherwise orginal fastq files.
    """
    lanes = sorted(UNITS[wildcards.sample].keys())
    if config['to_trim']:
        return expand(os.path.join(OUT_DIR, "{sample}_{lane}_{read}.trimmed.fastq.gz"), sample=wildcards.sample, lane=lanes, read=read)
    return [UNITS[wildcards.sample][lane][read] for lane in lanes]

def all_fastqc():
    out = []
    for s in samples:
        for lane in UNITS[s]:
            for r in ['R1', 'R2']:
                out.append(os.path.join(OUT_DIR, "qc/fastqc", f"{s}", f"{s}_{lane}_{r}_fastqc.zip"))
    return out

def all_fastp():
    out = []
    for s in samples:
        for lane in UNITS[s]:
            out.append(os.path.join(OUT_DIR, "qc/fastp", f"{s}", f"{s}_{lane}_fastp.json"))
    return out

def get_multiqc_input(wildcards):
    """
    Return all QC and analysis outputs for MultiQC.
    """
    inputs = []
    inputs += all_fastqc()
    if config['to_trim']:
        inputs += all_fastp()
    inputs += expand(os.path.join(OUT_DIR, "star_output", "{sample}", "{sample}.Log.final.out"), sample=samples)
    inputs += expand(os.path.join(OUT_DIR, "salmon_counts", "{sample}", "quant.sf"), sample=samples)
    inputs += expand(os.path.join(OUT_DIR, "qc/rustqc", "{sample}"), sample=samples)
    return inputs