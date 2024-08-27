
# extrarequirement:

bedgraphtobigwig http://hgdownload.soe.ucsc.edu/admin/exe/



# Ensure that the "read_info.filt.thr2.txt" is present in the parent folder

# Run the following command:

source source.txt

# The '80' file adds and substract 40bp to the midpoint, producing a 'smoothed' nucleosome profile

# The .py script is easily editable to include extra information in produced bedfile. This is useful if one want to select reads on size, strand or midpoint position.