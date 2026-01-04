#!/bin/bash

set -euo pipefail

./cenmap --generate-config > r2_config.yaml

find data/asm/R2/ -mindepth 1 -type l -exec basename {} \; > R2_samples.txt

./cenmap -c r2_config.yaml -j 30 --workflow-profile workflow/profiles/lpc_all/
