#
# This spec file is read by gfortran when linking.
# It is used to specify the libraries we need to link in, in the right
# order.
#

%rename lib liborig
*lib: -lquadmath -lm %(libgcc) %(liborig) -rpath /Users/biplab/Documents/rna_seq_workflow/.snakemake/conda/c8a34079/lib
