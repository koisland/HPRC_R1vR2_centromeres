from os.path import join

OUTDIR = "results"
LOGDIR = "logs"
BMKDIR = "benchmarks"

H1_ALN_CFG = {
    "ref": {
        "HG00673_H1_release_2": "fasta/R2/HG00673_chr16_HG00673#1#JAHBBZ020000008.1:15687200-22504109.fa.gz"
    },
    "sm": {
        "HG00673_H1_release_1": [
            "fasta/R1/HG00673_chr16_HG00673#1#JAHBBZ010000141.1:16117034-17742360.fa.gz",
            "fasta/R1/HG00673_rc-chr16_HG00673#1#JAHBBZ010000125.1:1-1560945.fa.gz",
            "fasta/R1/HG00673_rc-chr16_HG00673#1#JAHBBZ010000081.1:1-3527707.fa.gz",
        ]
    },
    "temp_dir": join(OUTDIR, "temp"),
    "output_dir": OUTDIR,
    "logs_dir": LOGDIR,
    "benchmarks_dir": BMKDIR,
    "aln_threads": 4,
    "aln_mem": "4GB",
    "mm2_opts": "-x asm5 -K 8G",
}
H2_ALN_CFG = {
    "ref": {
        "HG00673_H2_release_2": "fasta/R2/HG00673_chr16_HG00673#2#JAHBBY020000044.1:15855519-24453678.fa.gz"
    },
    "sm": {
        "HG00673_H2_release_1": [
            "fasta/R1/HG00673_rc-chr16_HG00673#2#JAHBBY010000128.1:15857114-17666371.fa.gz",
            "fasta/R1/HG00673_chr16_HG00673#2#JAHBBY010000121.1:1-2516994.fa.gz",
            "fasta/R1/HG00673_rc-chr16_HG00673#2#JAHBBY010000088.1:1-4360130.fa.gz",
        ]
    },
    "temp_dir": join(OUTDIR, "temp"),
    "output_dir": OUTDIR,
    "logs_dir": LOGDIR,
    "benchmarks_dir": BMKDIR,
    "aln_threads": 4,
    "aln_mem": "4GB",
    "mm2_opts": "-x asm5 -K 8G",
}


# Align assemblies to reference.
module align_asm_to_ref_h1:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="minimal",
        )
    config:
        H1_ALN_CFG


module align_asm_to_ref_h2:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="minimal",
        )
    config:
        H2_ALN_CFG


use rule * from align_asm_to_ref_h1 as h1_*


use rule * from align_asm_to_ref_h2 as h2_*


rule all:
    input:
        rules.h1_all.input,
        rules.h2_all.input,
    default_target: True
