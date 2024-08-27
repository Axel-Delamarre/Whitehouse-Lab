import pandas as pd
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from mpl_interactions import ioff, panhandler, zoom_factory
import re
import numpy as np
import mplcursors
from matplotlib.widgets import MultiCursor
import mpld3
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns
import subprocess
import sys
import glob
import os
from pathlib import Path
#import proplot as pplt

########### SET CURRENT DIRECTORY #############

#plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['font.family'] = 'Arial' # Change to 'Helvetica', 'Times New Roman', etc. 
#plt.rcParams['font.size'] = 8 # Set a readable font size

colorp = "Blues"

read_info = 'read_info.chrIV-100k-114k.bed'

setgridsize = 200
setvmin = 0
setvmax = 50
setalpha = 1
ylimit = 450
yminlim = 0
xlimmin = 100000
xlimmax = 113600
figwidth = 5
figheight = 2


## get current directory
currentPath = os.path.dirname(__file__)

## change directory to subfolders just in case
#sub = os.chdir(currentPath)

#print(sub)

###########################################################

########### INPUT FILES #############

## store bed_file in variable
for bed_file in Path(currentPath).rglob(read_info):
    bed_file = bed_file.name

bed_file = currentPath + '/' + bed_file

## store .gff file with origin of replication coordinates in variable
#for gff_new_oor2 in Path(currentPath).rglob('OriDB_all_origins_saccer1_to_W303_redo.confirmed.gff'):
#    gff_new_oor2 = gff_new_oor2.name

#gff_new_oor2 = currentPath + '/' + gff_new_oor2

## store new genes .csv file with gene standard names and coordinates in variable
for new_genes in Path(currentPath).rglob('verified_ORF_simple.bed'):
    new_genes = new_genes.name

new_genes = currentPath + '/' + new_genes


###### Extra Bed ####
#for new_BS in Path(currentPath).rglob('Gal4BS_GAL2.r.bed'):
#    new_BS = new_BS.name

#new_BS = currentPath + '/' + new_BS


##### Extra Bed ####
#for new_TSS in Path(currentPath).rglob('W3_Nuc_Plus1_Jiang-Pugh2009.r.bed'):
#    new_TSS = new_TSS.name

#new_TSS = currentPath + '/' + new_TSS

###########################################################

########### INPUT VARIABLES #############

# input variables
#chromosome = input("Enter the chromosome: ")
chromosome = "chrIV"
bed_column = "strand"
track1_feature_type = "gene"

###########################################################

########### COLOR AND CHR DICTIONARIES #############

# All features have to be written and assigned a color.
color_lookup = {
    '+': "black",
    '-': "black"

}

# All features have to be written and assigned a color.
#color_lookup_ARS = {
#    '+': "orange",
#    '-': "orange_r"

#}

## dictionary for changing chr values to .bed format chr names
#chroms = {'chr1.chr1.final': "chr1", 
#          'chr2.chr2.final': "chr2", 
#          'chr3.chr3.final': "chr3",
#         'chr4.chr4.final': "chr4",
#         'chr5.chr5.final': "chr5",
#         'chr6.chr6.final': "chr6",
#         'chr7.chr7.final': "chr7",
#         'chr8.chr8.final': "chr8",
#         'chr9.chr9.final': "chr9",
#         'chr10.chr10.final': "chr10",
#         'chr11.chr11.final': "chr11",
#         'chr12.chr12.final': "chr12",
#         'chr13.chr13.final': "chr13",
#         'chr14.chr14.final': "chr14",
#         'chr15.chr15.final': "chr15",
#         'chr16.chr16.final': "chr16"}

# using dictionary to convert specific columns for new_genes .bed file
convert_dict = {'start': int,
                'stop': int,
                'standard_name': str,
                'strand1': str,
                'length': float}

###########################################################

########### PROCESS GENOMIMC .BED FILE (MAIN PLOT) #############

## import the .bed file
genome = pd.read_table(bed_file, header=None, comment='#')

## select specific columns
genome = genome[[0, 1, 2]]

## rename those columns
genome = genome.rename(columns={0: "chr", 1: "start", 2: "stop"})

## only select the records based on the chr we are interested in (set at the top)
genome = genome[genome["chr"] == chromosome]

## create length column
genome["length"] = genome["stop"] - genome["start"]

## set x to the largest value in the stop column (for x axis)
x = genome["stop"].max()

## set y to the largest value in the length column (for y axis)
y = ylimit

## create middle column
genome['middle'] = genome['stop'] - (genome['length'] / 2)

## create list of all the values
lst = genome.values.tolist()

###########################################################

########### PROCESS NEW .CSV FILE WITH CORRECT GENE NAMES AND COORDINATES #############

# read file as pandas df
new_genes3 = pd.read_table(new_genes, header=None, comment='#')

## rename the columns
new_genes3 = new_genes3.rename(columns={0: "chr_x", 1: "start", 2: "stop", 3: "standard_name",
                                        4: "length", 5: "strand1"})

new_genes3["length"] = new_genes3["stop"] - new_genes3['start']

## drop the first row because it is just the column names
new_genes3 = new_genes3.drop(new_genes3.index[0])

## only select the records based on the chr we are interested in (set at the top)
new_genes3 = new_genes3[new_genes3["chr_x"] == chromosome]

## drop the unnecessary columns
new_genes3 = new_genes3.drop(columns = ['chr_x'])

## rearrange the columns for the list
new_genes3 = new_genes3[['start', 'stop', 'standard_name', 'strand1', 'length']]

## use convert_dict to change types of each column
new_genes3 = new_genes3.astype(convert_dict)

## create pandas series that only contains positive strand (red)
blue_bar = new_genes3[new_genes3["strand1"] == '+']

## create pandas series that only contains negative strand (blue)
darkblue_bar = new_genes3[new_genes3["strand1"] == '-']

## create a list of the positive strand genes
pos_genes = blue_bar.values.tolist()

## create a list of the negative strand genes
neg_genes = darkblue_bar.values.tolist()

# print(new_genes)


########### PROCESS NEW .CSV FILE WITH CORRECT GENE NAMES AND COORDINATES #############

# read file as pandas df
#new_BS3 = pd.read_table(new_BS, header=None, comment='#')

## rename the columns
#new_BS3 = new_BS3.rename(columns={0: "chr_x", 1: "start", 2: "stop", 3: "standard_name",
#                                        4: "length", 5: "strand1"})

#new_BS3["length"] = new_BS3["stop"] - new_BS3['start']

## drop the first row because it is just the column names
#new_BS3 = new_BS3.drop(new_BS3.index[0])

## only select the records based on the chr we are interested in (set at the top)
#new_BS33 = new_BS3[new_BS3["chr_x"] == chromosome]

## drop the unnecessary columns
#new_BS3 = new_BS3.drop(columns = ['chr_x'])

## rearrange the columns for the list
#new_BS3 = new_BS3[['start', 'stop', 'standard_name', 'strand1', 'length']]

## use convert_dict to change types of each column
#new_BS3 = new_BS3.astype(convert_dict)

## create pandas series that only contains positive strand (red)
#blue_barBS = new_BS3[new_BS3["strand1"] == '+']

## create pandas series that only contains negative strand (blue)
#darkblue_barBS = new_BS3[new_BS3["strand1"] == '-']

## create a list of the positive strand genes
#pos_BS = blue_barBS.values.tolist()

## create a list of the negative strand genes
#neg_BS = darkblue_barBS.values.tolist()

# print(new_genes)

########### PROCESS NEW .CSV FILE WITH CORRECT GENE NAMES AND COORDINATES #############

# read file as pandas df
#new_TSS3 = pd.read_table(new_TSS, header=None, comment='#')

## rename the columns
#new_TSS3 = new_TSS3.rename(columns={0: "chr_x", 1: "start", 2: "stop", 3: "standard_name",
#                                        4: "length", 5: "strand1"})

#new_TSS3["length"] = new_TSS3["stop"] - new_TSS3['start']

## drop the first row because it is just the column names
#new_TSS3 = new_TSS3.drop(new_TSS3.index[0])

## only select the records based on the chr we are interested in (set at the top)
#new_TSS3 = new_TSS3[new_TSS3["chr_x"] == chromosome]

## drop the unnecessary columns
#new_TSS3 = new_TSS3.drop(columns = ['chr_x'])

## rearrange the columns for the list
#new_TSS3 = new_TSS3[['start', 'stop', 'standard_name', 'strand1', 'length']]

## use convert_dict to change types of each column
#new_TSS3 = new_TSS3.astype(convert_dict)

## create pandas series that only contains positive strand (red)
#blue_barTSS = new_TSS3[new_TSS3["strand1"] == '+']

## create pandas series that only contains negative strand (blue)
#darkblue_barTSS = new_TSS3[new_TSS3["strand1"] == '-']

## create a list of the positive strand genes
#pos_TSS = blue_barTSS.values.tolist()

## create a list of the negative strand genes
#neg_TSS = darkblue_barTSS.values.tolist()


###########################################################

########### ORIGIN OF REPLICATION #############

## read file as pandas df
#gff_new_oor2 = pd.read_table(gff_new_oor2, header=None, comment='#')

## extract certain columns and label them
#gff_new_oor2 = gff_new_oor2[[0, 3, 4, 6, 8, 9]]
#gff_new_oor2 = gff_new_oor2.rename(columns={0: "chr", 3: "start", 4: "stop", 6: "strand",
 #                                        8: "description", 9:"name"})

# assign the positive and negative strands with blue and red respectfully  
#gff_new_oor2["colors"] = gff_new_oor2["strand"].apply(lambda x: color_lookup_ARS[x])

## for every row in the chr column, convert it to the key in chroms dictionary
#gff_new_oor2["chr"] = gff_new_oor2["chr"].apply(lambda x: chroms[x])

## only select the records based on the chr we are interested in (set at the top)
#gff_new_oor2 = gff_new_oor2[gff_new_oor2["chr"] == chromosome]

## create standard name column with pattern searching description column extracting something like gene=GENE323
#gff_new_oor2['gene_name'] = gff_new_oor2['description'].str.extract('(gene\S[A-Z]+[0-9]+)', expand=True)

## remove the gene= from each row
#gff_new_oor2['gene_name'] = gff_new_oor2['gene_name'].str.extract('([A-Z]+[0-9]+)', expand=True)

## drop columns for plotting
#gff_new_oor2 = gff_new_oor2.drop(columns = ['chr', 'strand', 'description', 'name'])

## create length column
#gff_new_oor2['length'] = gff_new_oor2['stop'] - gff_new_oor2['start']

## rearrange the columns for the list
#gff_new_oor2 = gff_new_oor2[['start', 'stop', 'length', 'colors', 'gene_name']]

## create a list of the nascent genes for plotting
#gff_new_oor2_list = gff_new_oor2.values.tolist()

###########################################################

########### INITITATE PLOTS #############

## create figure and axis with 5 subplots that share the x axis but not the y axis
fig, axs = plt.subplots(3, gridspec_kw={'height_ratios': [30, 1, 1]}, sharex=True, sharey=False)

## use the seaborn-whitegrid style
#plt.style.use('seaborn-whitegrid')

## show cursor when hovering over figure on the plot
#multi = MultiCursor(fig.canvas, (axs[0], axs[1], axs[2]), color='yellow', lw=1,
#                    horizOn=False, vertOn=True)


## features function for plotting (may not need)
def features(df):
    df['width'] = df['stop'] - df['start']
    for chrom, group in df.groupby('chr'):
        yrange = (0, 1)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(xranges, yrange, facecolors=group["colors"])


###########################################################

########### SUBPLOT 0 - HISTOGRAM SHOWING COUNTS OF MIDDLE (x) AND LENGTH (y) #############

## histogram of middle point with custom bins
#sns.histplot(data=genome, x="middle", kde=True, color="lightblue", ax=axs[0], bins=100)

## plot average line
#axs[0].axvline(genome['middle'].mean(), color='k', linestyle='dashed', linewidth=1)

## to have plots extend to the boarder of the figure
#axs[0].axis('tight')

## set title for the first subplot
#axs[0].set_title(chromosome)

###########################################################

########### SUBPLOT 1 - HEXBIN SHOWING COUNTS OF MIDDLE (x) AND LENGTH (y) #############

## hexbin that is like a heatmap but hexagon shaped for each bin
hb = axs[0].hexbin(genome['middle'], genome['length'], gridsize=setgridsize, cmap=colorp, 
vmin = setvmin, vmax = setvmax, edgecolors = "white", linewidths = 0, mincnt = 0.5, alpha = setalpha)
#hb_s = axs[1].hexbin(genome['stop'], genome['length'], gridsize=300, cmap="ocean_r", 
#vmin = 0, vmax = 10, edgecolors = "white", linewidths = 0.1, mincnt = 0.5, alpha = 0.8)


## create colorbar for the hexbin plot only
cb = fig.colorbar(hb, ax=axs[0])


## create label for the colorbar
cb.set_label('Counts')

## to have plots extend to the boarder of the figure
axs[0].axis('tight')

## set the y axis to go as far as the y variable set above
axs[0].set_ylim(yminlim, y)

## set the y axis label
axs[0].set_ylabel("Length (bp)")

# Format the x-axis labels as thousands
import matplotlib.ticker as ticker
formatter = ticker.FuncFormatter(lambda x, pos: f'{int(x/1000)}')
axs[0].xaxis.set_major_formatter(formatter)


###########################################################

########### SUBPLOT 2 - RECTANGLES SHOWING SHADED IN AREA OF START AND END OF EACH ORIGIN OF REP  #############

## for each element in the list, create a rectangle
#for s in range(len(gff_new_oor2_list)):
#    box_size = 1
#    axs[2].add_patch(patches.Rectangle((gff_new_oor2_list[s][0], 0), gff_new_oor2_list[s][2], 1,
#                                       color='orange', lw=1, label = gff_new_oor2_list[s][4]))

## have plot extend to the boarder of the figure
#axs[2].axis('tight')

## have axis start at 0 and extend to where last chromosome ends on genome
#axs[2].set_xlim(0, x)

## set y axis ticks to nothing
#axs[2].get_yaxis().set_ticks([])

###########################################################

#for g in range(len(pos_TSS)):
#    box_size = 1
#    axs[2].add_patch(patches.Rectangle((pos_TSS[g][0], 0), pos_TSS[g][4], 1,
#                                       color='forestgreen', lw=1, label=pos_TSS[g][2]))

## have plot extend to the boarder of the figure
#axs[2].axis('tight')

## have axis start at 0 and extend to where last chr ends on genome
#axs[2].set_xlim(xlimmin, x)

## set ticks to nothing
#axs[2].get_yaxis().set_ticks([])


########### SUBPLOT 3 - RECTANGLES SHOWING SHADED IN AREA OF START AND END OF EACH GENE ON + STRAND  #############

for g in range(len(pos_genes)):
    box_size = 1
    axs[1].add_patch(patches.Rectangle((pos_genes[g][0], 0), pos_genes[g][4], 1,
                                       color='darkgray', lw=1, label=pos_genes[g][2]))

## have plot extend to the boarder of the figure
#axs[1].axis('tight')

## have axis start at 0 and extend to where last chr ends on genome
axs[1].set_xlim(xlimmin, xlimmax)

## set ticks to nothing
axs[1].get_yaxis().set_ticks([])

###########################################################

########### SUBPLOT 4 - RECTANGLES SHOWING SHADED IN AREA OF START AND END OF EACH GENE ON - STRAND  #############

for n in range(len(neg_genes)):
    box_size = 1
    axs[2].add_patch(patches.Rectangle((neg_genes[n][0], 0), neg_genes[n][4], 1,
                                       color='darkgray', lw=1, label=neg_genes[n][2]))

## have plot extend to the boarder of the figure
axs[2].axis('tight')

##have axis start at 0 and extend to where last chr ends on genome
axs[2].set_xlim(xlimmin, xlimmax)

## label x axis
axs[2].set_xlabel("Chr IV (Kbp)")

## set ticks to nothing
axs[2].get_yaxis().set_ticks([])


###########################################################

#for n in range(len(neg_TSS)):
#    box_size = 1
#    axs[5].add_patch(patches.Rectangle((neg_TSS[n][0], 0), neg_TSS[n][4], 1,
 #                                      color='darkgreen', lw=1, label=neg_TSS[n][2]))

## have plot extend to the boarder of the figure
#axs[5].axis('tight')

##have axis start at 0 and extend to where last chr ends on genome
#axs[5].set_xlim(xlimmin, x)

## label x axis
#axs[4].set_xlabel("Position (bp)")

## set ticks to nothing
#axs[5].get_yaxis().set_ticks([])

###########################################################

########### SUBPLOT 3 - RECTANGLES SHOWING SHADED IN AREA OF START AND END OF EACH GENE ON + STRAND  #############

#for g in range(len(pos_BS)):
#    box_size = 1
#    axs[6].add_patch(patches.Rectangle((pos_BS[g][0], 0), pos_BS[g][4], 1,
#                                       color='crimson', lw=1, label=pos_BS[g][2]))

## have plot extend to the boarder of the figure
#axs[6].axis('tight')

## label x axis
#axs[6].set_xlabel("Position (bp)")

## have axis start at 0 and extend to where last chr ends on genome
#axs[6].set_xlim(xlimmin, x)

## set ticks to nothing
#axs[6].get_yaxis().set_ticks([])

###########################################################

########### SUBPLOT 4 - RECTANGLES SHOWING SHADED IN AREA OF START AND END OF EACH GENE ON - STRAND  #############

#for n in range(len(neg_BS)):
 #   box_size = 1
  #  axs[6].add_patch(patches.Rectangle((neg_BS[n][0], 0), neg_BS[n][4], 1,
   #                                    color='darkgreen', lw=1, label=neg_BS[n][2]))

## have plot extend to the boarder of the figure
#axs[6].axis('tight')

##have axis start at 0 and extend to where last chr ends on genome
#axs[6].set_xlim(0, x)

## label x axis
#axs[6].set_xlabel("Position (bp)")

## set ticks to nothing
#axs[6].get_yaxis().set_ticks([])


###########################################################

########### CREATE HOVER-OVER INTERACTIVENESS  #############

import threading

def update_annotation(sel):
    """ Update the annotation for the currently hovered item. """
    label = sel.artist.get_label()
    sel.annotation.set_text(label)
    sel.annotation.set(position=(1, 10), anncoords="offset points")
    sel.annotation.get_bbox_patch().set(fc="grey", alpha=.5)
    sel.annotation.arrow_patch.set(arrowstyle="simple", fc="grey", alpha=.5)
    sel.annotation.set_visible(True)

def hide_annotation(sel):
    """ Hide the annotation after a delay. """
    def delayed_hide():
        # Delay in seconds before hiding the annotation
        delay = 2
        threading.Timer(delay, lambda: sel.annotation.set_visible(False)).start()

    delayed_hide()

# Create mplcursor objects for interactivity
cursor = mplcursors.cursor(axs[1].patches, hover=True)
cursor1 = mplcursors.cursor(axs[2].patches, hover=True)
#cursor2 = mplcursors.cursor(axs[3].patches, hover=True)
#cursor3 = mplcursors.cursor(axs[4].patches, hover=True)

# Connect the update and hide functions to the appropriate events
cursor.connect("add", update_annotation)
cursor.connect("remove", hide_annotation)
cursor1.connect("add", update_annotation)
cursor1.connect("remove", hide_annotation)
#cursor2.connect("add", update_annotation)
#cursor2.connect("remove", hide_annotation)
#cursor3.connect("add", update_annotation)
#cursor3.connect("remove", hide_annotation)


########### CREATE AND PLOT THE LEGEND  #############

## create legend elements
#legend_elements = [Line2D([0], [0], color='forestgreen', lw=2, label='+'),
#                   Line2D([0], [0], color='darkgreen', lw=2, label='-')]
#                   Line2D([0], [0], color='orange', lw=2, label='origin of replication')]

## plot the legend
#plt.legend(handles=legend_elements, bbox_to_anchor=(1, 33), frameon=False, ncol=3)

###########################################################

########### ADJUST SPACING BETWEEN SUBPLOTS AND FIGURE SIZE FOR BEST RESOLUTION #############

## ensure the subplots are connected with no whitespace in between them
plt.subplots_adjust(wspace=0, hspace=0)

plt.rcParams.update({'font.size': 6})

## set the width of the entire figure (all subplots)
fig.set_figwidth(figwidth)

## set the height of the entire figure (all subplots)
fig.set_figheight(figheight)

## align all axis after implementing colorbar (colorbar will overlap other supblots by defualt - here's how to fix)

## store the position of subplot [0]
pos = axs[0].get_position()

## store the position of subplot [1]
pos1 = axs[0].get_position()
pos2 = axs[1].get_position()
pos3 = axs[2].get_position()
#pos4 = axs[4].get_position()
#pos5 = axs[5].get_position()
#pos6 = axs[6].get_position()
#axs[0].set_position([pos1.x0, pos.y0, pos1.width, pos.height])
axs[1].set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
axs[2].set_position([pos1.x0, pos3.y0, pos1.width, pos3.height])
#axs[4].set_position([pos1.x0, pos4.y0, pos1.width, pos4.height])
#axs[5].set_position([pos1.x0, pos5.y0, pos1.width, pos5.height])
#axs[6].set_position([pos1.x0, pos6.y0, pos1.width, pos6.height])

## save the plot
plt.savefig('heax_local_'+read_info+'.png', dpi=600, bbox_inches="tight")

## show the plot
#plt.show()
