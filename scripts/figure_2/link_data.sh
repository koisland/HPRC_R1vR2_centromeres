#!/bin/bash

set -euo pipefail

sm="${1}"

indir_asm_r1="/project/logsdon_shared/data/HPRC/assemblies/hifiasm_year1/${sm}"
indir_asm_r2="/project/logsdon_shared/data/HPRC/assemblies/hifiasm/${sm}"
indir_ont="/project/logsdon_shared/data/HPRC/ont/${sm}"

outdir_asm="/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/CenMAP/data/assemblies"
outdir_ont="/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/CenMAP/data/ont"
mkdir -p "${outdir_ont}" "${outdir_asm}"

ln -s "${indir_asm_r1}" "${outdir_asm}/${sm}_R1"
ln -s "${indir_asm_r2}" "${outdir_asm}/${sm}_R2"
ln -s "${indir_ont}" "${outdir_ont}/${sm}_R1"
ln -s "${indir_ont}" "${outdir_ont}/${sm}_R2"
