wcs = glob_wildcards("cenplot/{ctg}.yaml")


rule replot:
    input:
        "cenplot/{ctg}.yaml",
    output:
        directory("plots/{ctg}"),
    shell:
        """
        cenplot draw -t {input} -d {output}
        """


rule all:
    input:
        expand(rules.replot.output, ctg=wcs.ctg),
    default_target: True
