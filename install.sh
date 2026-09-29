#!/bin/bash
# Install genome-decoded's software into the repo: ./env (+ ./env-x86 on Apple Silicon) and HIBAG.
# Needs conda/mamba (Miniforge recommended). Reference data is separate: ./download_refs.sh
# Usage: ./install.sh            (re-running is safe; existing envs are updated)
set -euo pipefail
cd "$(dirname "$0")"
MAMBA=$(command -v mamba || command -v conda || true)
[ -n "$MAMBA" ] || { echo "Need mamba or conda (https://github.com/conda-forge/miniforge)"; exit 1; }
OS=$(uname -s) ARCH=$(uname -m)
echo "== platform: $OS $ARCH"

if [ -d env ]; then "$MAMBA" env update -p ./env -f environment.yml; else "$MAMBA" env create -p ./env -f environment.yml; fi

if [ "$OS" = Linux ]; then
    "$MAMBA" install -y -p ./env -c conda-forge -c bioconda plink2=2.0.0a.6.9 admixture=1.3.0 bioconductor-hibag
    ln -sfn env env-x86          # scripts call env-x86/bin/plink2; on Linux that is the same env
elif [ "$OS" = Darwin ] && [ "$ARCH" = arm64 ]; then
    if ! arch -x86_64 /usr/bin/true 2>/dev/null; then
        echo "Rosetta 2 is required for plink2/admixture: softwareupdate --install-rosetta"; exit 1
    fi
    if [ -d env-x86 ]; then CONDA_SUBDIR=osx-64 "$MAMBA" env update -p ./env-x86 -f environment-x86.yml
    else CONDA_SUBDIR=osx-64 "$MAMBA" env create -p ./env-x86 -f environment-x86.yml; fi
    echo "subdir: osx-64" > env-x86/.condarc
    # HIBAG (Bioconductor): into ./r-lib for the system R if there is one, else into a conda R inside ./env.
    if command -v Rscript >/dev/null; then
        mkdir -p r-lib
        R_LIBS_USER=$PWD/r-lib Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib="r-lib", repos="https://cloud.r-project.org"); .libPaths(c("r-lib", .libPaths())); BiocManager::install("HIBAG", lib="r-lib", ask=FALSE, update=FALSE)'
    else
        "$MAMBA" install -y -p ./env -c conda-forge r-base r-biocmanager compilers
        env/bin/Rscript -e 'BiocManager::install("HIBAG", ask=FALSE, update=FALSE)'
    fi
else
    "$MAMBA" install -y -p ./env -c conda-forge -c bioconda plink2=2.0.0a.6.9 admixture=1.3.0
    ln -sfn env env-x86
    echo "NOTE: install HIBAG for your R (BiocManager::install('HIBAG')) into ./r-lib"
fi

echo "== checks"
env/bin/plink --version | head -1
env-x86/bin/plink2 --version | head -1
env/bin/bcftools --version | head -1
env/bin/snakemake --version
env/bin/python -c "import pandas, scipy, pysam; print('python libs ok')"
if [ -x env/bin/Rscript ] && env/bin/Rscript -e 'library(HIBAG)' >/dev/null 2>&1; then echo "HIBAG ok (conda R)"
elif R_LIBS_USER=$PWD/r-lib Rscript -e 'library(HIBAG)' >/dev/null 2>&1; then echo "HIBAG ok (r-lib)"
else echo "WARNING: HIBAG not found; the HLA step will fail"; fi
echo "== done. Next: ./download_refs.sh"
