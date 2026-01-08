from os.path import join
from collections import defaultdict


OUTPUT_DIR = config.get("output_dir", ".")
CENMAP_DIR = config.get("cenmap_dir", "/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/results")
# {sample: {dtype: { R1: [], R2: []}}}
DATA = defaultdict(lambda: defaultdict(dict))
# {sample: {release: []}}
NF_DATA = defaultdict(lambda: defaultdict(list))
RELEASES = ["R1", "R2"]

with open(config["manifest"]) as fh:
    # Skip header.
    next(fh)
    for line in fh:
        sample, asm_h1, asm_h2, bam, release = line.strip().split(",")
        DATA[sample]["asm"][release] = [asm_h1, asm_h2]
        DATA[sample]["bam"][release] = [bam]

with open(config["nucflag_r1_manifest"]) as fh:
    next(fh)
    for line in fh:
        sample_id, dtype, urls, uri = line.strip().split("\t")
        if dtype != "bedfile":
            continue
        NF_DATA[sample_id]["R1"].append(uri)


with open(config["nucflag_r2_manifest"]) as fh:
    next(fh)
    for line in fh:
        sample_id, haplotype, _, uri = line.strip().split(",")
        NF_DATA[sample_id]["R2"].append(uri)


wildcard_constraints:
    sm="|".join(DATA.keys()),
    release="|".join(RELEASES),


rule download_assemblies:
    output:
        asm=join(OUTPUT_DIR, "data", "{sm}_{release}.fa.gz"),
        asm_fai=join(OUTPUT_DIR, "data", "{sm}_{release}.fa.gz.fai"),
    params:
        uris=lambda wc: DATA[wc.sm]["asm"][wc.release],
    shell:
        """
        for uri in {params.uris}; do
            aws s3 --no-sign-request cp "${{uri}}" - | zcat | bgzip >> {output.asm}
        done
        samtools faidx {output.asm}
        """

rule download_nucflag:
    output:
        bed=join(OUTPUT_DIR, "data", "{sm}_{release}_nucflag.bed"),
    params:
        uris=lambda wc: NF_DATA[wc.sm][wc.release],
    shell:
        """
        for uri in {params.uris}; do
            aws s3 --no-sign-request cp "${{uri}}" - >> {output.bed}
        done
        """

rule filter_hor_array_length:
    input:
        lengths=join(CENMAP_DIR, "{release}/results/{release}/final/bed/all_AS-HOR_lengths.bed"),
        fai=join(CENMAP_DIR, "{release}/results/{release}/2-concat_asm/{sm}-asm-renamed-reort.fa.fai"),
        nucflag=rules.download_nucflag.output,
    output:
        final_bed=join(OUTPUT_DIR, "data", "{release}_centromeres", "{sm}_AS-HOR_lengths.bed"),
    params:
        script="scripts/nucflag/filter_arrays.py"
    shell:
        """
        # Add the contig length to the array length bed.
        # If reverse complemented, convert the coordinates back to the original coordinates.
        # Intersect with bedtools taking any centromere's HOR arrays without any overlap with NucFlag non-HET calls.
        # Format back to original lengths bed.
        python {params.script} <(join \
            <(grep {wildcards.sm} {input.lengths} | awk -v OFS="\\t" '{{
                og=$1;
                split($1, ctg_name, ":");
                is_rc=(og ~ "rc-");
                $1=ctg_name[1]
                print $0, is_rc, og
            }}' | sort -k1,1) \
            <(sort -k1,1 {input.fai} | cut -f1,2) | \
        awk -v OFS="\\t" '{{
            is_rc=$(NF-2);
            og=$(NF-1);
            len=$(NF);
            if (is_rc) {{
                st=len-$3;
                end=len-$2
            }} else {{
                st=$2;
                end=$3;
            }};
            split($1, ctg_name, "_");
            print ctg_name[3], st, end, og, $2, $3, $4, $5, $6
        }}' | \
        bedtools intersect -a - -b <(grep -v HET {input.nucflag}) -loj) > {output}
        """

rule merge_hor_array_lengths:
    input:
        expand(rules.filter_hor_array_length.output, sm=DATA.keys(), release="{release}")
    output:
        final_bed=join(OUTPUT_DIR, "data", "{release}_centromeres", "all_AS-HOR_lengths.bed"),
    shell:
        """
        cat {input} > {output}
        """

rule download:
    input:
        expand(rules.download_assemblies.output, sm=DATA.keys(), release=RELEASES),
        expand(rules.download_nucflag.output, sm=DATA.keys(), release=RELEASES),
    default_target:
        True

rule filter_cens:
    input:
        expand(rules.filter_hor_array_length.output, sm=DATA.keys(), release=RELEASES),
        expand(rules.merge_hor_array_lengths.output, release=RELEASES),
