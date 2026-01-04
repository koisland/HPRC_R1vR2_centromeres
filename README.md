# HPRC R2 Centromere Analysis
Using `CenMAP` commit `ca012f7c80642446c25836f8144ea8106a42fbd1`.

## Getting Started
```bash
git clone https://github.com/logsdon-lab/CenMAP \
--revision ca012f7c80642446c25836f8144ea8106a42fbd1 \
--depth 1 \
--recursive
cd CenMAP
```

## Get NucFlag results
```bash
```
> WIP

## Count complete centromeres
```bash
python exp/scripts/complete_counts/count_complete_cens.py \
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
> WIP

## HOR array length
> WIP

Manual changes:
* Recalculate chr19 array length due to chimeric HOR breaking array. Reduce threshold on required number of monomers in HOR.
