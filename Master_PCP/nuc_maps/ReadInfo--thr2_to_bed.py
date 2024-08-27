import subprocess
import pysam
import sys
import pandas as pd
import gzip
from itertools import islice, combinations
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


# ---> INPUT file and parameters
path= ''
read_info= '../read_info.filt.thr2.txt'


#----> OUTPUT files


#df1 = pd.read_table(path + read_info, sep=' ', names=['bc','chr','start','length','rID', 'midpoint', 'diff'])
df_reads = pd.read_table(read_info, names=['chr', 'read_start', 'read_stop', 'length', 'read_midpoint', 'read_strand', 'bc', 'read_rID', 'read_delta'])

# Determine the tagged extremity as the 5' of R1
df_reads['read_tag_extremity'] = np.where(df_reads['read_strand'] == '+',
                               df_reads['read_start'],
                               df_reads['read_stop'])

df_reads['value'] = np.where(
    df_reads['length'] < 260,
    df_reads['read_midpoint'],
    np.where(
        df_reads['read_strand'] == '+',
        df_reads['read_tag_extremity'] + 83,
        df_reads['read_tag_extremity'] - 83
    )
)

df_reads['value+1'] = df_reads['value'] +1

df_bed = df_reads[['chr','value', 'value+1']]

df_bed.to_csv('read_info.mixed.bed', header=None, index=None, sep='\t')