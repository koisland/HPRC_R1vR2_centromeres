from os.path import join


OUTDIR = "results"
LOGDIR = "logs"
BMKDIR = "benchmarks"

HG00673_H1_ALN_CFG = {
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
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}
HG00673_H2_ALN_CFG = {
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
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}
HG01071_H1_ALN_CFG = {
    "ref": {
        "HG01071_H1_release_2": "/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/CenMAP/results/8-humas_annot/seq/HG00621_R2_chr6_HG00621#1#CM087851.1:57498547-64828711.fa"
    },
    "sm": {"HG01071_H1_release_1": []},
    "temp_dir": join(OUTDIR, "temp"),
    "output_dir": OUTDIR,
    "logs_dir": LOGDIR,
    "benchmarks_dir": BMKDIR,
    "aln_threads": 4,
    "aln_mem": "4GB",
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}
HG01071_H2_ALN_CFG = {
    "ref": {
        "HG01071_H2_release_2": "/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/CenMAP/results/8-humas_annot/seq/HG00621_R2_chr6_HG00621#2#CM087865.1:57362978-66193602.fa"
    },
    "sm": {"HG01071_H2_release_1": []},
    "temp_dir": join(OUTDIR, "temp"),
    "output_dir": OUTDIR,
    "logs_dir": LOGDIR,
    "benchmarks_dir": BMKDIR,
    "aln_threads": 4,
    "aln_mem": "4GB",
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}


# Align assemblies to reference.
module HG00673_H1_align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        HG00673_H1_ALN_CFG


module HG00673_H2_align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        HG00673_H2_ALN_CFG


module HG01071_H1_align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        HG01071_H1_ALN_CFG


module HG01071_H2_align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        HG01071_H2_ALN_CFG


use rule * from HG00673_H1_align_asm_to_ref as HG00673_H1_*


use rule * from HG00673_H2_align_asm_to_ref as HG00673_H2_*


use rule * from HG01071_H1_align_asm_to_ref as HG01071_H1_*


use rule * from HG01071_H2_align_asm_to_ref as HG01071_H2_*


rule all:
    input:
        rules.HG00673_H1_all.input,
        rules.HG00673_H2_all.input,
        rules.HG01071_H1_all.input,
        rules.HG01071_H2_all.input,
    default_target: True
