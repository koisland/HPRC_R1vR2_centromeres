#!/bin/bash

set -euo pipefail

./cenmap --generate-config > r1_config.yaml

find data/asm/R1/ -mindepth 1 -type l -exec basename {} \; > R1_samples.txt

# Update config with samples.
# Remove all options and run minimally.

./cenmap -c r1_config.yaml -j 30 --workflow-profile workflow/profiles/lpc_all/
