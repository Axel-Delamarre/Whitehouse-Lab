import pysam
import pandas as pd

path= ''
bam= '../deduplicated.bam'
out_ccf= 'read_info.bed'

# open the BAM file
bamfile = pysam.AlignmentFile(path + bam, 'rb')

# initialize an empty list to store the read ID and sequence
ccf = path + out_ccf

# loop through each read in the BAM file
with open(ccf, "w") as new:
    for read in bamfile.fetch():
        # Check if the read is R1
        if read.is_read1:

        # extract the read ID and sequence
            read_id = '@'+ read.query_name
            pair_length = abs(read.template_length)

            if read.is_forward: # assign 16
                    strand = "+"
                    start = 1 + read.reference_start
                    end = start + pair_length
                    midpoint = start + pair_length//2
                    
            else:
                    strand = "-"
                    end = read.reference_end - 1
                    start = end - pair_length
                    midpoint = end - pair_length//2

            chromosome = read.reference_name

            #start_position = 1 + read.reference_start

            read_id_parts = read_id.split('_')

        # append the read ID parts and sequence to the read_info list
            new.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome, start, end, pair_length, midpoint, strand, read_id_parts[1], read_id_parts[0]))

# close the BAM file
bamfile.close()
