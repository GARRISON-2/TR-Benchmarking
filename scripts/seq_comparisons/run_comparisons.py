import os
import itertools 
import Levenshtein as lv
from readers import *

# set directory variables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH


def compareString(line1, line2):
    # compare the bed info and the refs (for now) 
    #print(f"{seq1.cur_line[3]},    {seq2.cur_line[3]}")
    return_dict = {}

    dist = lv.distance(line1.ref, line2.ref)

    # get the positional difference and calculate the lv distance based only on positions in the actual reference
    start_diff = line2.pos_info["start"] - line1.pos_info["start"]
    end_diff = len(line1.ref) - line2.pos_info["end"] - line1.pos_info["end"]

    dist2 = lv.distance(line1.ref[start_diff: end_diff], line2.ref[start_diff: end_diff])
    
    return_dict = {
        "dist": dist,
        "dist2": dist2
    }

    return return_dict 


def checkIdx(idx):
    try:
        return int(idx)
    except ValueError:
        return 0



def mainloop(bed_file, file1, file2):
    # print("Enter catalog path:")
    # print("Enter 1st VCF path:")
    # print("Enter 2nd VCF path")

 
    # open all input and output files
    with open(os.path.join(DATA_DIR, bed_file) , 'r') as bf, \
    open(os.path.join(DATA_DIR, file1) , 'r') as f1, \
    open(os.path.join(DATA_DIR, file2) , 'r') as f2, \
    open(os.path.join(LOCAL_DATA, "reference-comp.tsv"), "w") as bof, \
    open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w") as vof:
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

        # Write metadata to output file


        # set column headers for output tsv files
        bof.write("#FILE\tBEDPOS\tVCFPOS\tDIST\n")
        vof.write(f"#BEDPOS\tTOTDIST\tfile1\tfile2\n")


        # make list to keep track of which files have reached their ends
        f_states = [obj.end_state for obj in reader_list]


        # Main Operations loop
        # loop until all files have reached their final line
        #while not all(f_states):

        # loop until the bed file has reached its end
        while not bed.end_state:
            bof_out_str = ""
            vof_out_str = ""
            gt_list = []
            compare_vcf = True # boolean value to disable/enable running the vcf comparisons as needed

            # cycle through the bed and vcf files and perform operations on the current line
            for i, reader in enumerate(reader_list):   
                end_of_cycle = ( i == len(reader_list) -1 ) 

                # bed alignment checks
                #if type(reader) != BEDReader:




                # run vcf specific operations 
                if type(reader) != BEDReader:                   
                    # print(f"File[{i}]{reader.pos_info[1]}:{reader_list[0].pos_info[1]} - {reader.pos_info[1] > reader_list[0].pos_info[1]}") # DEBUGGING - REMOVE FOR LAUNCH
                    
                    # if current file position is ahead of bed position range, then pause operations
                    if ((reader.pos_info["start"] > bed.pos_info["end"]) & (reader.pos_info["chrom"] == bed.pos_info["chrom"])) | \
                            reader.end_state:    

                        reader.pause = True # pause current vcf from moving to the next line
                        compare_vcf = False # pause vcf comparisons for the rest of the cycle

                    else:
                        reader.pause = False

                        # BED Ref Comparison Operations
                        # compare current vcf ref with bed ref
                        bed_dist = lv.distance(bed.ref, reader.ref)


                        # VCF Comparison Operations
                        # if at least 2 of the vcf files have not hit their end, 
                        # then proceed with VCF comparisons
                        if (not sum(f_states[1:]) > 2) & compare_vcf: # FIX FOR READABILITY BEFORE LAUNCH

                            # check indices to make sure they are ints (eg. accounting for 1|., etc.)
                            gt_idx = (checkIdx(reader.cur_line[-2][0]), checkIdx(reader.cur_line[-2][2]))
                            alleles = [reader.ref, *reader.alt]

                            # build genotype based on gt_idx indices
                            gt = [alleles[gt_idx[0]], alleles[gt_idx[1]]]

                            # add to list of genotypes for comparison between VCFs
                            gt_list.append(gt)                     

                            if end_of_cycle:

                                # get distances for both permutations of comparing alleles
                                vert_dist = lv.distance(gt_list[0][0], gt_list[1][0]) + lv.distance(gt_list[0][1], gt_list[1][1])
                                cross_dist = lv.distance(gt_list[0][1], gt_list[1][0]) + lv.distance(gt_list[0][0], gt_list[1][1])

                                gt_dist = vert_dist

                                # if the allele order is mismatched between the genotypes
                                if vert_dist > cross_dist:
                                    # swap the order of the alleles in the second gt to ensure 
                                    # they can be compared properly
                                    #gt_list[1][0], gt_list[1][1] = gt_list[1][1], gt_list[1][0]

                                    # use the cross distance since it is less,
                                    # which means it is more likely to be the correct comparison ordering
                                    gt_dist = cross_dist


                                vof_out_str = f"{list(bed.pos_info.values())}\t{gt_dist}\t{gt_list[0]}\t{gt_list[1]}\n"
                         

                        bof_out_str += (f"{reader.name}\t{list(bed.pos_info.values())}\t{list(reader.pos_info.values())}\t{bed_dist}\n")
                
        
                # read in and set the new current line
                if not type(reader) == BEDReader:
                    reader.read()

            # wait to move the bed file's current line until
            # all other files have been gone through
            bed.read()

            bof.write(bof_out_str)
            vof.write(vof_out_str)


        print("\n---PROGRAM COMPLETE---\n")




mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.strdust.vcf", "HG001.strkit.vcf")
