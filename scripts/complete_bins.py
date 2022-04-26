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
    filename=str(snakemake.log[0]),
    filemode="w",
    level=logging.DEBUG,
    format='[%(asctime)s] %(name)s %(levelname)s: %(message)s'
)
logger = logging.getLogger(__file__)

# Reading the data
col=['contig_id']
data=pd.read_csv(snakemake.input.TSV, sep="\t", header=0, usecols=col).values

# Saving each complete contig as a fasta file
with open(snakemake.input.FNA[0], "r") as ifile:
    for record in SeqIO.parse(ifile, "fasta"):
        if record.id in data:
            with open(os.path.join(snakemake.params.bins, "%.fna" % record.id), "w") as ofile:
                SeqIO.write(record, ofile, "fasta")