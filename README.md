# HPRC R2 Centromere Analysis
Using `CenMAP` commit `ca012f7c80642446c25836f8144ea8106a42fbd1`.

## Getting Started
Clone and setup env.
```bash
git clone https://github.com/koisland/HPRC_R1vR2_centromeres.git
cd HPRC_R1vR2_centromeres
conda env create -f env.yaml --name hprc_r1vr2_cens
conda activate hprc_r1vr2_cens
```

Clone `CenMAP`. 
```bash
git clone https://github.com/logsdon-lab/CenMAP \
--revision ca012f7c80642446c25836f8144ea8106a42fbd1 \
--depth 1 \
--recursive
```

Download assemblies and NucFlag data.
```bash
which aws
snakemake -np download \
-s scripts/nucflag/data.smk \
--config manifest=config/hprc_all_manifest.csv nucflag_r1_manifest=config/nucflag_r1.tsv nucflag_r2_manifest=config/nucflag_r2.csv
```

## Run CenMAP
```bash
pushd CenMAP
# Generate config.
./cenmap --generate-config > r1_config.yaml
# Fill samples and concat_asm.input_dir
# Comment out cdr_finder and nucflag
# Then run.
./cenmap -c r1_config.yaml -j 30 --workflow-profile workflow/profiles/lpc/
# Repeat for R2

# Then return to main dir.
popd
```

> [!NOTE]
> We are using the repo's script not the bioconda script.

## Filter centromeres NucFlag results
```bash
which bedtools python
python -c "import polars"

snakemake -np filter_cens \
-s scripts/nucflag/data.smk \
--config manifest=config/hprc_all_manifest.csv nucflag_r1_manifest=config/nucflag_r1.tsv nucflag_r2_manifest=config/nucflag_r2.csv
```

## Count complete centromeres
```bash
python scripts/complete_counts/count_complete_cens.py \
-a /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/results/R1/results/R1/final/bed/all_AS-HOR_lengths.bed \
-b /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/results/R2/results/R2/final/bed/all_AS-HOR_lengths.bed \
-o figures/complete_centromeres.png
```

Manual changes:
* R1: Remove the following:
    * All chr3. Only one array complete.
    * Sample HG03492. No R2 equivalent.
* R2: Remove the following incomplete scaffolds:
    * HG00621_rc-chr2_HG00621#2#JAHBCC020000017.1:2-689335
    * HG00438_chrX_HG00438#2#JAHBCA020000020.1:2-703696

### And accurate centromeres
```bash
python scripts/complete_counts/count_complete_cens.py \
-a /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R1_centromeres/all_AS-HOR_lengths.bed \
-b /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R2_centromeres/all_AS-HOR_lengths.bed \
-o figures/complete_and_accurate_centromeres.png
```

## HOR array length for complete centromeres
```bash
python scripts/live_arrays/plot_hor_length.py \
-a /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/results/R1/results/R1/final/bed/all_AS-HOR_lengths.bed \
-b /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/results/R2/results/R2/final/bed/all_AS-HOR_lengths.bed \
-o figures/hor_array_length_complete_centromeres.png
```

### And accurate centromeres
```bash
python scripts/live_arrays/plot_hor_length.py \
-a /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R1_centromeres/all_AS-HOR_lengths.bed \
-b /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R2_centromeres/all_AS-HOR_lengths.bed \
-o figures/hor_array_length_complete_and_accurate_centromeres.png
```

Manual changes:
* Same as [Count complete centromeres](#count-complete-centromeres) but without filter:
    * All chr3. Only one array complete.
