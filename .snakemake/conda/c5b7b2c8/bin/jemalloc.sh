#!/bin/sh

prefix=/Users/biplab/Documents/rna_seq_workflow/.snakemake/conda/c5b7b2c8
exec_prefix=/Users/biplab/Documents/rna_seq_workflow/.snakemake/conda/c5b7b2c8
libdir=${exec_prefix}/lib

DYLD_INSERT_LIBRARIES=${libdir}/libjemalloc.2.dylib
export DYLD_INSERT_LIBRARIES
exec "$@"
