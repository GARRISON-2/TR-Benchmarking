import os
import itertools 
import Levenshtein as lv
from helpers.readers import *

# set directory variables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH


def mainloop(bed_file, file1, file2):
    print("\n\n\n")
    # open all input and output files
    with open(os.path.join(DATA_DIR, bed_file) , 'r') as bf, \
    open(os.path.join(DATA_DIR, file1) , 'r') as f1, \
    open(os.path.join(DATA_DIR, file2) , 'r') as f2:
        # File Prep
        file_stack = [bed_file, f1, f2]         # CURRENTLY REPRESENTS INPUT FILE STACK

        # create list of reader objects to keep track of individual file info
        reader_list = [BEDReader(bf)]
        for i, file in enumerate(file_stack[1:]):
                reader_list.append(VCFReader(file))

        bed = reader_list[0]

        # move past header data in vcf files
        for i, rdr in enumerate(reader_list[1:]):
            rdr.skipMetaData()

        overlap_count = 0
        cur_pos = bed.pos_info
        bed.read() # move bed file ahead by one line
        ahead_pos = bed.pos_info

        # loop until the bed file has reached its end
        while not bed.end_state:
          
            # cycle through the bed and vcf files 
            for i, reader in enumerate(reader_list):   
                end_of_cycle = ( i == len(reader_list) -1 ) 

            
                if type(reader) != BEDReader:              
                   
                    # if position is ahead of current bed pos or with ahead pos
                    if ((reader.pos_info["start"] > cur_pos["end"]) & \
                            (reader.pos_info["chrom"] == cur_pos["chrom"])) | \
                                reader.end_state:

                        reader.pause = True # pause current vcf from moving to the next line
                        
                    else:
                        reader.pause = False

                        # check against 1 ahead bed position
                        if (reader.pos_info["end"] > ahead_pos["start"]) & (reader.pos_info["chrom"] == ahead_pos["chrom"]) & (reader.pos_info["chrom"] == cur_pos["chrom"]):
                            print(f"{reader.name} POS: {list(reader.pos_info.values())}, BEDPOS1: {list(cur_pos.values())}, BEDPOS2: {list(ahead_pos.values())}")
                            overlap_count += 1

                          
                if end_of_cycle:
                    cur_pos = bed.pos_info
                    bed.read() 
                    ahead_pos = bed.pos_info

                # read in and set the new current line
                if not type(reader) == BEDReader:
                    reader.read()


        print(f"Test Done, {overlap_count} overlaps")




mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.strdust.vcf", "HG001.strkit.vcf")
