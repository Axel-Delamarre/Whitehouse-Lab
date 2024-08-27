import pysam
import pandas as pd
import numpy as np


#### INPUT FILES #####
seed_read_info = 'seed_read_info.bed'

read_info = 'read_info.r.s.bed'


#### OUTPUT FILES #####
seedreach = 'seedreach_distribution_15.txt'

seed_reach_bed = 'seed_reach_15.ns.bed'

####### SETTINGS ######

Upper_limit_fragment_number = 20

print("\n",
    'Upper_limit_fragment_number is set at ',Upper_limit_fragment_number, "\n",
    "\n",
    "---------------------------------------------------------------")

#######################

df_preseed = pd.read_table(seed_read_info, names=['seed_chr', 'seed_start', 'seed_stop', 'seed_length', 'seed_midpoint', 'seed_strand', 'BC', 'seed_rID', 'seed_delta'])

# calculate the real extremity at which ligated the seed as the 3' end of R1
df_preseed['seed_extremity'] = np.where(df_preseed['seed_strand'] == '+',
                               df_preseed['seed_stop'],
                               df_preseed['seed_start'])

# select for the one that have an acceptable size
df_select_preseed = df_preseed[(df_preseed['seed_length'] > 15) & (df_preseed['seed_length'] < 1000)]

# make a new DF replacing the seed pos by the real seed pos
columns_desired = ['seed_chr', 'seed_start', 'seed_stop', 'seed_length', 'seed_midpoint', 'seed_strand', 'BC', 'seed_rID', 'seed_extremity']
df_seed = df_select_preseed[columns_desired]

df_reads = pd.read_table(read_info, names=['read_chr', 'read_start', 'read_stop', 'read_length', 'read_midpoint', 'read_strand', 'BC', 'read_rID', 'read_delta'])

# Determine the tagged extremity as the 5' of R1
df_reads['read_tag_extremity'] = np.where(df_reads['read_strand'] == '+',
                               df_reads['read_start'],
                               df_reads['read_stop'])

df_reads['value'] = np.where(
    df_reads['read_length'] < 260,
    df_reads['read_midpoint'],
    np.where(
        df_reads['read_strand'] == '+',
        df_reads['read_tag_extremity'] + 83,
        df_reads['read_tag_extremity'] - 83
    )
)

# group the reads on the BC
df_reads_groupby = df_reads.groupby('BC').agg(lambda x: list(x)).apply(list).reset_index()

# measure the amount of fragments per BC
df_reads_groupby['number of fragments'] = df_reads_groupby['read_rID'].str.len()

# select the complexes that have more than 1 but less than 20 reads 
df_multi_reads_groupby = df_reads_groupby[(df_reads_groupby["number of fragments"] > 1) & (df_reads_groupby["number of fragments"] < Upper_limit_fragment_number)]

# select the complexes that have a BC present in the seeds DF
df_seeded = df_multi_reads_groupby[df_multi_reads_groupby['BC'].isin(df_seed['BC'])]
# to know how many are not seeded
df_not_seeded = df_multi_reads_groupby[~df_multi_reads_groupby['BC'].isin(df_seed['BC'])]

# select the seeds whose BC is present in the seeded
df_seedin =  df_seed[df_seed['BC'].isin(df_seeded['BC'])]
# to Know how many seeds are not in the df_seedin
df_noseedin =  df_seed[~df_seed['BC'].isin(df_seeded['BC'])]

with open('report.txt', 'w') as file:
    # Write an introductory sentence
    file.write("Numbers on the >1 CCFs:\n\n")
    
    # Write the number of rows for each DataFrame
    file.write(f"number_of_rows_df_multi_reads_groupby: {len(df_multi_reads_groupby):,} rows, Number of CCFs\n")
    file.write(f"number_of_rows_df_seed: {len(df_seed):,} rows, Number of Seeds\n")
    file.write(f"number_of_rows_df_seeded: {len(df_seeded):,} rows, Number of CCFs that have a seed\n")
    file.write(f"number_of_rows_df_not_seeded: {len(df_not_seeded):,} rows, Number of Orphans CCF\n")
    file.write(f"number_of_rows_df_seedin: {len(df_seedin):,} rows, Number of of seed that have a CCF\n")
    file.write(f"number_of_rows_df_noseedin: {len(df_noseedin):,} rows, Number of Childless seed \n")



df_merge = pd.merge(df_seedin, df_seeded, on='BC', how='outer')

df_exp = df_merge.explode(['read_chr', 'read_start', 'read_stop', 'read_length', 'read_midpoint', 'read_strand', 'read_rID', 'read_tag_extremity', 'value'])

# Keep only the lines where the seed and read chr arethe same
df_cis = df_exp[df_exp['seed_chr'] == df_exp['read_chr']]

#remove the lines where the seed is the read
df_ns = df_cis[df_cis['seed_rID'] != df_cis['read_rID']]

#calculate the seed-reach
df_ns['distance']=df_ns['value']-df_ns['seed_extremity']


####################### Create a CSV file that can be used to make the seed-reach plot in Prism #########################

data = df_ns['distance']

hist, bins = np.histogram(data, bins=6000, range=[-3000,3000])

df = pd.DataFrame({'bin_edges': bins[:-1], 'bin_counts': hist})

df.to_csv(seedreach, sep='\t', index=False)

####################### Create a bed with 1 line per interaction

#df_ns['seed_end'] = np.where(df_ns['seed_strand'] == '+',
 #                              df_preseed['seed_real_pos'] - df_preseed['seed_length'],
 #                              df_preseed['seed_real_pos'] + df_preseed['seed_length'])


df_bed = df_ns[['seed_chr', 'seed_start', 'seed_stop', 'seed_length', 'seed_midpoint', 'seed_strand', 'BC', 'seed_rID', 'seed_extremity', 'read_chr', 'read_start', 'read_stop', 'read_length', 'read_midpoint', 'read_strand', 'read_rID', 'read_tag_extremity', 'distance']]

#df_bed['distance'] = df_bed['distance'].abs()

df_bed.to_csv(seed_reach_bed, sep='\t', index=False, header=False)






