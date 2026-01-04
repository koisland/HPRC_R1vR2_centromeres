#!/bin/bash

set -euo pipefail

R1_ASM_DIR="/project/logsdon_shared/data/HPRC/assemblies/hifiasm_year1/"
R2_ASM_DIR="/project/logsdon_shared/data/HPRC/assemblies/hifiasm/"
INPUT_ASM_DIR="/project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/data/asm"
INPUT_R1_ASM_DIR="${INPUT_ASM_DIR}/R1"
INPUT_R2_ASM_DIR="${INPUT_ASM_DIR}/R2"
mkdir -p "${INPUT_R1_ASM_DIR}" "${INPUT_R2_ASM_DIR}"

# No HG03492
while read -r line; do
    sample=$(basename $line .done)
    echo "On ${sample}"
    
    ln -s "${R1_ASM_DIR}/${sample}" "${INPUT_R1_ASM_DIR}/${sample}"
    ln -s "${R2_ASM_DIR}/${sample}" "${INPUT_R2_ASM_DIR}/${sample}"

done < <(find /project/logsdon_shared/data/HPRC/assemblies/hifiasm_year1/ -name "*.done")
