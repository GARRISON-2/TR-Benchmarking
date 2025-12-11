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



def compareGt(gt1, gt2):
    
    # get distances for both permutations of comparing alleles
    vert_dist = lv.distance(gt1[0], gt2[0]) + lv.distance(gt1[1], gt2[1])
    cross_dist = lv.distance(gt1[1], gt2[0]) + lv.distance(gt1[0], gt2[1])

    gt_dist = vert_dist

    # if the allele order is mismatched between the genotypes
    if vert_dist > cross_dist:
        # swap the order of the alleles in the second gt to ensure 
        # they can be compared properly
        #gt_list[1][0], gt_list[1][1] = gt_list[1][1], gt_list[1][0]

        # use the cross distance since it is less,
        # which means it is more likely to be the correct comparison ordering
        gt_dist = cross_dist

    return gt_dist


def stateCheck(rl):
    end_states = [obj.end_state for obj in rl]
    pse_states = [obj.pause for obj in rl]

    return ((len(end_states) - sum(end_states)) >= 2) & ((len(pse_states) - sum(pse_states)) >= 2)




def mainloop(bed_file, file1, file2):
    # print("Enter catalog path:")
    # print("Enter 1st VCF path:")
    # print("Enter 2nd VCF path")

 
    # open all input and output files
    with open(os.path.join(DATA_DIR, bed_file) , 'r') as bf, \
    open(os.path.join(DATA_DIR, file1) , 'r') as f1, \
    open(os.path.join(DATA_DIR, file2) , 'r') as f2, \
    open(os.path.join(LOCAL_DATA, "bed-comp.tsv"), "w") as bof, \
    open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w") as vof:
        # File Prep
        file_stack = [f1, f2, bf]         # CURRENTLY REPRESENTS INPUT FILE STACK

        # create list of reader objects to keep track of individual file info
        vcf_rdrs = []
        offsets = [[-1, 0], [-1, 0]]
        for i in range(len(file_stack)):
            if file_stack[i].name.endswith("vcf"):
                vcf_rdrs.append(VCFReader(file_stack[i], 
                                        start_offset=offsets[i][0], 
                                        end_offset=offsets[i][1]))
            elif file_stack[i].name.endswith("bed"):
                bed = BEDReader(file_stack[i])
        

        # move past header data in vcf files
        for i, rdr in enumerate(vcf_rdrs):
            if type(rdr) != BEDReader:
                rdr.skipMetaData()

        # Write metadata to output file


        # set column headers for output tsv files
        bof.write("#FILE\tBEDPOS\tVCFPOS\tDIST\n")
        vof.write(f"#CHROM\tSTART\tEND\tPOSDIST\tGTDIST\tfile1\tfile2\n")


        # Main Operations loop
        # loop until the bed file has reached its end
        while not bed.end_state:
        #while not all(end_states): or loop till all files are done
            bof_out_str = ""
            vof_out_str = ""
            gt_list = []
            compare_vcf = True # boolean value to disable/enable running the vcf comparisons as needed


            # bed alignment checks
            #if type(reader) != BEDReader:



            # cycle through the vcf files and perform operations on the current line
            for i, reader in enumerate(vcf_rdrs): 
                last_file = ( i == len(vcf_rdrs) -1 )
           
                # if current file position is ahead of bed position range, then pause operations
                if ((reader.pos_info["start"] > bed.pos_info["end"]) & (reader.pos_info["chrom"] == bed.pos_info["chrom"])) | \
                        reader.end_state:    

                    reader.pause = True # pause current vcf from moving to the next line

                else:
                    reader.pause = False

                    # BED Ref Comparison Operations
                    # compare current vcf ref with bed ref
                    start_diff = abs(reader.pos_info["start"] - bed.pos_info["start"])
                    end_diff = abs(reader.pos_info["end"] - bed.pos_info["end"])

                    bed_dist = start_diff + end_diff

                    end_num = sum([obj.end_state for obj in vcf_rdrs]) # number of files that have ended
                    pse_num = sum([obj.pause for obj in vcf_rdrs]) # number of files that are paused

                    bof_out_str += (f"{reader.name}\t{list(bed.pos_info.values())}\t{list(reader.pos_info.values())}\t{bed_dist}\n")


                    # Build Current Genotype
                    reader.buildGt()


                if stateCheck(vcf_rdrs) & last_file:

                    # compare 2 genotypes
                    gt_dist = compareGt(vcf_rdrs[0].genotype, vcf_rdrs[1].genotype)

                    vcf_start_diff = abs(vcf_rdrs[0].pos_info["start"] - vcf_rdrs[1].pos_info["start"])
                    vcf_end_diff = abs(vcf_rdrs[0].pos_info["end"] - vcf_rdrs[1].pos_info["end"])

                    pos_dist = vcf_start_diff + vcf_end_diff 

                    vof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}\t{pos_dist}\t{gt_dist}\t{vcf_rdrs[0].genotype}\t{vcf_rdrs[1].genotype}\n"    
            

            # read current lines for all files
            [rdr.read() for rdr in vcf_rdrs]
            bed.read()

            bof.write(bof_out_str)
            vof.write(vof_out_str)


    print("\n---PROGRAM COMPLETE---\n")




mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.PAW79146.haplotagged.URfix.atarva.vcf", "HG001.PAW79146.haplotagged.URfix.strkit.vcf")
