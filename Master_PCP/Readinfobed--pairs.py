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
read_info= 'read_info.filt.thr2.txt'


#----> OUTPUT files

plt_size_distribution = read_info + '_size-ditribution.pdf'

plt_fragment_per_CCF_total = path + read_info + '_fragments_per_CCF_total.pdf'
plt_fragment_per_CCF_nosingle = path + read_info + '_fragments_per_CCF_nosingle.pdf'

tojuice = path + read_info +'.tojuice'
juicer = path + read_info +'.juicer'


current_time = datetime.now().time()
print("\n",
    'program starting at: ',current_time, "\n",
    "\n",
    "---------------------------------------------------------------")


#----------------------------------------------------------------------
#----------------------------------------------------------------------

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


current_time = datetime.now().time()
print("\n",
    'Generating the size distribution graph at: ', current_time, "\n"
    "---------------------------------------------------------------")
#make a graph about the size distribution and save it as PDF
df_reads['length'].hist(bins=500, range=[1,1000])
plt.savefig(plt_size_distribution)
plt.close()

# create a new df with only the columns we want
df_hic = df_reads[['bc','chr','value','length']]

#df_hic = df_hic.dropna(subset=['midpoint'])

# Get it back to the fromat BC  chr-pos-length
df_hic['concatenated'] = df_hic.apply(lambda row: '-'.join([str(row['chr']), str(row['value']), str(row['length'])]), axis=1)
df_hic2 = df_hic[['bc','concatenated']]

current_time = datetime.now().time()
print("\n",
    'Grouping by barcode ', current_time, "\n"
    "---------------------------------------------------------------")

# Groupby back on the barcode
df_before_ccf = df_hic2.groupby('bc').agg(lambda x: list(x)).apply(list).reset_index()

df_before_ccf['number of fragments'] = df_before_ccf['concatenated'].str.len()

df_after_ccf = df_before_ccf[df_before_ccf["number of fragments"] > 1]


current_time = datetime.now().time()
print("\n",
    'Generating the fragment per ccf graph', current_time, "\n"
    "---------------------------------------------------------------")

df_before_ccf['number of fragments'].hist(bins=50, range=[1,50])
plt.savefig(plt_fragment_per_CCF_total)
plt.close()

df_after_ccf['number of fragments'].hist(bins=50, range=[1,50])
plt.savefig(plt_fragment_per_CCF_nosingle)
plt.close()


juice = df_after_ccf
juice.rename(columns={'concatenated': 'locus'})

juice.to_csv(tojuice, header=None, index=None, sep='\t')

juice = pd.read_table(tojuice, names=['barcode', 'locus', 'number of ccfs'])

# convert the column into a list
juice['locus'] = juice.locus.apply(lambda x: x[1:-1].split(','))


# convert dataframe columns into lists
arr = juice['locus'].values.tolist()



# convert roman numerals to arabic numbers
chromosome_conversion = {"'chr1": "chrI",
                         " 'chr1": "chrI", 
                        "'chr2": "chrII",
                        " 'chr2": "chrII",
                        "'chr3": "chrIII",
                        " 'chr3": "chrIII",
                        "'chr4": "chrIV",
                        " 'chr4": "chrIV",
                        "'chr5": "chrV",
                        " 'chr5": "chrV",
                        "'chr6": "chrVI",
                        " 'chr6": "chrVI",
                        "'chr7": "chrVII",
                        " 'chr7": "chrVII",
                        "'chr8": "chrVIII",
                        " 'chr8": "chrVIII",
                        "'chr9": "chrIX",
                        " 'chr9": "chrIX",
                        "'chr10": "chrX",
                        " 'chr10": "chrX",
                        "'chr11": "chrXI",
                        " 'chr11": "chrXI",
                        "'chr12": "chrXII",
                        " 'chr12": "chrXII",
                        "'chr13": "chrXIII",
                        " 'chr13": "chrXIII",
                        "'chr14": "chrXIV",
                        " 'chr14": "chrXIV",
                        "'chr15": "chrXV",
                        " 'chr15": "chrXV",
                        "'chr16": "chrXVI",
                        " 'chr16": "chrXVI"}

# format data to create hic file using juicer pre
# short format with score
#
# 0 chrI 1020 1000 16 chrII 2030 2000 0.333
#

current_time = datetime.now().time()
print("\n",
    'Generating the juicer file', current_time, "\n"
    "---------------------------------------------------------------")

## open both output juicer files at the same time (this method will also close them automatically)
with open(juicer, 'w') as juicer:
    ## loop through list
    for i in range(len(arr)):
        ## for each barcode, collect combinations of 2
        comb = combinations(arr[i], 2)
        ## iterate over each individual combination (pair)
        for j in comb:
            ## set chr1 equal to the first element of the locus value split by the dash
            chr1 = j[0].split('-')[0]
            ## if the chr1 value is present in any of the keys in the chromosome_conversion dictionary, replace it with its value 
            if chr1 in chromosome_conversion:
                chr1 = chromosome_conversion[chr1]
            ## set pos1 equal to the second element of the locus value split by the dash
            pos1 = j[0].split('-')[1]
            ## set strand1 equal to the third element of the locus value split by the dash
            str1 = j[0].split('-')[2]
            ## replace the lagging quotation mark with nothing
            str1 = str1.replace("'", '')
            ## set fragment1
            frag1 = "0"
            ## set chr2 equal to the first element of the locus value split by the dash
            chr2 = j[1].split('-')[0]
            ## if the chr2 value is present in any of the keys in the chromosome_conversion dictionary, replace it with its value
            if chr2 in chromosome_conversion:
                chr2 = chromosome_conversion[chr2]
            ## set pos2 equal to the second element of the locus value split by the dash
            pos2 = j[1].split('-')[1]
            ## set strand2 equal to the third element of the locus value split by the dash
            str2 = j[1].split('-')[2]
            ## replace the lagging quotation mark with nothing
            str2 = str2.replace("'", '')
            ## set fragment2
            frag2 = "1"

            ## write the output to the juicer file formatted by separating each value by a space and indenting after each row
            juicer.write(str(str1) + ' ' + chr1 + ' ' + str(pos1) + ' ' + str(frag1) + ' ' + str(str2) + ' ' + chr2 + ' ' + str(pos2) + ' ' + str(frag2) +'\n')


current_time = datetime.now().time()
print("\n",
    'Program done at: ', current_time, "\n"
    "---------------------------------------------------------------""\n",
    "\n")