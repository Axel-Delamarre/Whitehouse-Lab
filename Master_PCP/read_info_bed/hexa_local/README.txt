
# Run an awk command to select the region of interest and plot it on an hexaplot
# $1 is chromosome
# $2 is start
# $3 is stop

# There is a need to modify the python script (hexaplot_local.py) accordingly as well !!!

# Read_info.bed --> read_info.subset of interest.bed 

# example:

awk '$1=="chrIV" && $2>100000 && $3<115000' ../read_info.r.s.bed > read_info.chrIV-100k-115k.bed

# then run the hexaplot script

python hexaplot_local.py
