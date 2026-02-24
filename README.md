# HPRC R2 Centromere Analysis
Using `CenMAP` commit `80d553f9e88172d4fc77ce889c3637041dbc5677`.

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
--revision 80d553f9e88172d4fc77ce889c3637041dbc5677 \
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
-o figures/hor_array_length_complete_centromeres.png \
--mode total
```

### And accurate centromeres
```bash
python scripts/live_arrays/plot_hor_length.py \
-a /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R1_centromeres/all_AS-HOR_lengths.bed \
-b /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/data/R2_centromeres/all_AS-HOR_lengths.bed \
-o figures/hor_array_length_complete_and_accurate_centromeres.png \
--mode total
```

Manual changes:
* Same as [Count complete centromeres](#count-complete-centromeres) but without filter:
    * All chr3. Only one array complete.

## Figure 2: Release 1 and CDR breaks
Show example where R2's improved contiguity allows characterization of the centromere dip region.

```bash
pushd scripts/figure_2
```

### CenMAP
Run R1 and R2 with special branch of `CenMAP` allowing omitting entropy bed:
* https://github.com/logsdon-lab/CenMAP/tree/feature/omit-entropy-filter

Clone repo.
```bash
git clone git@github.com:logsdon-lab/CenMAP.git --recursive --revision 45a0fd326cbfbdb8669ff28658594f05f4fea8ee
pushd CenMAP
```

Do the following:
1. Link data for sample HG00673 and HG01071. Both assembly and ONT data should be downloaded before.
2. Run CenMAP for R1 and R2 on chr16 for HG00673. Also chr6 HG01071 H1.
```bash
# Done via ../cenmap.sh
bash ../link_data.sh HG00673
bash ../link_data.sh HG01071
snakemake -p --configfile ../config_chr16.yaml --workflow-profile workflow/profiles/lpc_all -j 50
snakemake -p --configfile ../config_chr6.yaml --workflow-profile workflow/profiles/lpc_all -j 50
popd
```

### SafFire and CenPlot
```bash
snakemake -p --workflow-profile workflow/profiles/lpc_all -j 4
popd
```

### CDR breakpoints by chromosome
Determine number of centromeres breaking in CDR.
```bash
pushd scripts/cdr_breakpoints
```

Does the following:
* Align R1 (query) assemblies to R2 (target). We need alignment coordinates to be in R2.
* Format CDR bed
* Liftover to R1 coordinates
* Determine breaks in CDR based on following criteria:
    *  If only one mapping must be close to start or end of contig and is incomplete
    *  If only multiple mapping.
* Count and plot. Also Fisher's exact test per chrom. Also do FDR control.

```bash
snakemake -p --workflow-profile workflow/profiles/lpc_all -j 4
popd
```

### Final
Merged all previous outputs in `InkScape`. Legend and additional formatting by Glennis.
