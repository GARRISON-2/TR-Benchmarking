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
    open(os.path.join(LOCAL_DATA, "reference-comp.tsv"), "w") as rof, \
    open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w") as vof:
        #---File Prep---
        # CURRENTLY REPRESENTS FILE STACK
        file_stack = [bed_file, f1, f2] 

        # create list of reader objects that 
        # help to keep the individual file information together
        file_list = [BEDReader(bf)]
        for i, file in enumerate(file_stack[1:]):
                file_list.append(VCFReader(file))


        # move past header information in vcf files
        for i, fi_info in enumerate(file_list[1:]):
            fi_info.skipMetaData()

        # Write metadata to output file

        # set column headers for output tsv files
        rof.write("#FILE\tBEDPOS\tVCFPOS\tDIST\n")
        vof.write(f"#BEDPOS\tTOTDIST\tfile1\tfile2\n")

        # set specific bed file variable for readability
        bed = file_list[0]



        # make list to easily keep track of which files have reached their ends
        f_states = [obj.end_state for obj in file_list]

        #---Main Operations loop---
        # loop until all files have reached their final line
        #while not all(f_states):
        # loop until the bed file has reached its end
        while not bed.end_state:
            rof_out_str = ""
            vof_out_str = ""
            gt_list = []
            vcf_comp = True

            # cycle through the bed and vcf files and perform operations on the current line
            for i, fi_info in enumerate(file_list):   
                last_cycle = ( i == len(file_list) -1 )
                

                # position alignment check
                # if i == 0:
                #     out_str += (f"BED Ref ({fi_info.pos_info}): {cur_line[3]}\n")


                # run vcf specific operations 
                if type(fi_info) != BEDReader:                   
                    # print(f"File[{i}]{fi_info.pos_info[1]}:{file_list[0].pos_info[1]} - {fi_info.pos_info[1] > file_list[0].pos_info[1]}") # DEBUGGING - REMOVE FOR LAUNCH
                    
                    # if current file position is ahead of bed position range, then pause operations
                    if ((fi_info.pos_info["start"] > bed.pos_info["end"]) & \
                            (fi_info.pos_info["chrom"] == bed.pos_info["chrom"]) | \
                                fi_info.end_state):    

                        fi_info.pause = True
                        vcf_comp = False

                    else:
                        fi_info.pause = False

                        #---BED Ref Comparison Operations---
                        # compare current vcf ref with bed ref
                        #diff_dict = compareString(bed, fi_info)
                        bed_dist = lv.distance(bed.ref, fi_info.ref)

                        #---VCF Comparison Operations---
                        # at least 2 of the vcf files have not hit their end, 
                        # then proceed with VCF-only operations
                        bool = (not sum(f_states[1:]) > 2) & vcf_comp # FIX BEFORE LAUNCH
                        if bool:

                            # check indices to make sure they are ints (eg. accounting for 1|., etc.)
                            gt_idx = (checkIdx(fi_info.cur_line[-2][0]), checkIdx(fi_info.cur_line[-2][2]))
                            alleles = [fi_info.ref, *fi_info.alt]

                            # build genotype based on gt_idx indices
                            gt = [alleles[gt_idx[0]], alleles[gt_idx[1]]]

                            # add to list of genotypes for comparison between VCFs
                            gt_list.append(gt)                     

                            # if it is the last cycle
                            if last_cycle:
                                vert_order = lv.distance(gt_list[0][0], gt_list[1][0]) + lv.distance(gt_list[0][1], gt_list[1][1])

                                cross_order = lv.distance(gt_list[0][1], gt_list[1][0]) + lv.distance(gt_list[0][0], gt_list[1][1])

                                gt_dist = vert_order

                                # if the allele order is mismatched between the genotypes
                                if vert_order > cross_order:
                                    # swap the order of the alleles in the second gt to ensure 
                                    # they can be compared properly
                                    #gt_list[1][0], gt_list[1][1] = gt_list[1][1], gt_list[1][0]

                                    # use the cross distance since it is less,
                                    # which means it is more likely to be the correct comparison ordering
                                    gt_dist = cross_order

                                
                                vof_out_str = f"{list(bed.pos_info.values())}\t{gt_dist}\t{gt_list[0]}\t{gt_list[1]}\n"

                            

                        rof_out_str += (f"{fi_info.name}\t{list(bed.pos_info.values())}\t{list(fi_info.pos_info.values())}\t{bed_dist}\n")
                
                    # wait to move the bed file to the next line until 
                    # all other files have been gone through
                    if last_cycle:
                        bed.nextLine()

                # move to next line in file
                # reader object automatically formats the line and sets it to cur_line parameter
                if not type(fi_info) == BEDReader:
                    fi_info.nextLine()

            # output results from single cycle through the files
            rof.write(rof_out_str)
            vof.write(vof_out_str)


        print("\n---PROGRAM COMPLETE---\n")



mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.strdust.vcf", "HG001.strkit.vcf")
