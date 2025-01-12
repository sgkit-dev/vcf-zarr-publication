#!/bin/bash

if [ ! -d SAIGE ]; then
    git clone https://github.com/Will-Tyler/SAIGE.git
fi

export PATH="/opt/miniconda3/bin:$PATH"

if ! conda list -n RSAIGE; then
    CONDA_SUBDIR=osx-64 conda env create -f env.yml
fi

conda activate RSAIGE

R CMD build SAIGE
R CMD INSTALL SAIGE
