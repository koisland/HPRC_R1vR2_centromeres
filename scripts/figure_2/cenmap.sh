#!/bin/bash

set -euo pipefail

git clone git@github.com:logsdon-lab/CenMAP.git \
    --recursive \
    --revision 45a0fd326cbfbdb8669ff28658594f05f4fea8ee

bash /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/link_data.sh HG01071
bash /project/logsdon_shared/projects/HPRC/CenMAP-R1vR2/exp/scripts/figure_2/link_data.sh HG00673

pushd CenMAP
snakemake -np \
    --configfile ../config_chr6.yaml \
    -j 30 \
    --workflow-profile workflow/profiles/lpc_all
snakemake -np \
    --configfile ../config_chr16.yaml \
    -j 30 \
    --workflow-profile workflow/profiles/lpc_all
popd
