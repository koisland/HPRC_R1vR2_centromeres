Does the following:
* Align R1 (query) assemblies to R2 (target). We need alignment coordinates to be in R2.
* Format CDR bed
* Liftover to R1 coordinates
* Determine breaks in CDR based on following criteria:
    *  If only one mapping must be close to start or end of contig and is incomplete
    *  If only multiple mapping.
* Count and plot. Also Fisher's exact test per chrom. Also do FDR control.

```bash
which bedtools impg
snakemake -np
```
