#!/bin/bash

export PATH="/opt/miniconda3/bin:$PATH"

rm -rf SAIGE
rm SAIGE_1.3.6.tar.gz
rm chr21_10_4.*
rm chr21_10_5.*
conda env remove -n RSAIGE -y
