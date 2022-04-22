#!/usr/bin/env python

import os
import pandas as pd
import logging
from Bio import SeqIO

##################################################
# LOGGING
##################################################
# logger
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=logging.DEBUG,
    format='[%(asctime)s] %(name)s %(levelname)s: %(message)s'
)
logger = logging.getLogger(__file__)

# Reading the data
data=pd.read_csv(snakemake.input.TSV, sep="\t", header=0, usecolss="contig_id").values

# Saving each complete contig as a fasta file
with open(snakemake.input.FNA, "r") as ifile:
    for record in SeqIO.parse(ifile, "fasta"):
        if record.id in data:
            with open(os.path.join(snakemake.params.bins, "%.fna" % record.id), "w") as ofile:
                SeqIO.write(record, ofile, "fasta")

