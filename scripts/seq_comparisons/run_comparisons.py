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
        # CURRENTLY REPRESENTS FILE STACK
        file_stack = [bed_file, f1, f2] 

        # Write metadata to output file


        # set column headers for output tsv files
        rof.write("#FILE\tBEDPOS\tVCFPOS\tDIST\n")
        vof.write(f"#BEDPOS\tTOTDIST\tfile1\tfile2\n")

        # make list of open file objects for easier management
        file_list = [BEDReader(bf)]
        for i, file in enumerate(file_stack[1:]):
                file_list.append(VCFReader(file))

        bed = file_list[0]

        # move past header information
        for i, fi_info in enumerate(file_list[1:]):
            fi_info.skipMetaData()


        f_states = [obj.end_state for obj in file_list]

        # loop until all files have reached their final line
        #while not all(f_states):
        while not bed.end_state:
        #for j in range(10):
            rof_out_str = ""
            vof_out_str = ""
            gt_list = []
    
            for i, fi_info in enumerate(file_list):   
                last_cycle = ( i == len(file_list) -1 ) 

                # position alignment check
                # if i == 0:
                #     out_str += (f"BED Ref ({fi_info.pos_info}): {cur_line[3]}\n")


                # if not the bed file
                if type(fi_info) != BEDReader:                   
                    # print(f"File[{i}]{fi_info.pos_info[1]}:{file_list[0].pos_info[1]} - {fi_info.pos_info[1] > file_list[0].pos_info[1]}") # DEBUGGING - REMOVE FOR LAUNCH
                    
                    # run vcf specific operations

                    # if current file position is ahead of bed position range, then pause operations
                    if ((fi_info.pos_info["start"] > bed.pos_info["end"]) & \
                            (fi_info.pos_info["chrom"] == bed.pos_info["chrom"]) ):    

                        fi_info.pause = True

                    else:
                        fi_info.pause = False

                        # compare current vcf ref with bed ref
                        #diff_dict = compareString(bed, fi_info)
                        bed_dist = lv.distance(bed.ref, fi_info.ref)

                        # at least 2 of the vcf files have not hit their end, 
                        # then proceed with VCF-only operations
                        if (sum(f_states[1:]) > 2):
                            # check indices to make sure they are ints (eg. accounting for 1|., etc.)
                            gt_idx = (checkIdx(fi_info.cur_line[-2][0]), checkIdx(fi_info.cur_line[-2][2]))
                            alleles = [fi_info.ref, *fi_info.alt]
                            # build genotype based on gt_idx indices
                            gt = [alleles[gt_idx[0]], alleles[gt_idx[1]]]

                            # add to list of genotypes for comparison between VCFs
                            gt_list.append(gt)                     

                            # if it is the last cycle and 
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
                

                    if last_cycle:
                        bed.nextLine()

                # move to next line in file
                # reader object automatically formats the line and sets it to cur_line parameter
                if not type(fi_info) == BEDReader:
                    fi_info.nextLine()

            # output results from single cycle through the files
            rof.write(rof_out_str)
            vof.write(vof_out_str)


        print("PROGRAM COMPLETE")



mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.strdust.vcf", "HG001.strkit.vcf")
