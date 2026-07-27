Required reference files:
- `genome.fa`
- `genes.gtf`
- `genome.fa.fai`
- `transcriptome.fa`
- `STAR_index/`

Download the genome FASTA and GTF annotation from the documented source.
Generate transcriptome FASTA with: `gffread genes.gtf -g genome.fa -w transcriptome.fa`.
The STAR index is generated automatically when the pipeline runs.

