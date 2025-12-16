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



def mainloop(bed_file, vcf_list):

    vcf_rdrs = []
    comps = 0
    offsets = [[-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0]]

    # File Prep
    # create list of reader objects to keep track of individual file info
    for i in range(len(vcf_list)):
                # create VCFReader object for file. VCFReader automatically opens the file
                vcf_rdrs.append(VCFReader(os.path.join(DATA_DIR, vcf_list[i]), 
                                    start_offset=offsets[i][0], 
                                    end_offset=offsets[i][1]))
                
                # move past header data in vcf files
                vcf_rdrs[i].skipMetaData()

    bed = BEDReader(os.path.join(DATA_DIR, bed_file))
    

    # print("Comparisons:\n1 (default): BED-VCF & VCF-VCF\n 2 : BED-VCF\n 3 : VCF-VCF")
    # #while comps > 3 and comps < 1:
    # comps = int(input(": "))
    comps = 1

    if comps == 1 or comps == 2:
        bof = open(os.path.join(LOCAL_DATA, "bed-comp.tsv"), "w")
        bof.write("#FILE\tBEDPOS\tVCFPOS\tDIST\n")

    if comps == 1 or comps == 3:
        vof = open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w")
        vof.write(f"#CHROM\tSTART\tEND\tPOSDIST\tGTDIST\tfile1\tfile2\tfile3\n")


    # Write metadata to output file



    # Main Operations loop
    # loop until the bed file has reached its end
    while not bed.end_state:
    #while not all(end_states): or loop till all files are done
        bof_out_str = ""
        vof_out_str = ""
        gt_dist_list = []
        ps_dist_list = []


        # cycle through the vcf files and perform operations on the current line
        for i, reader in enumerate(vcf_rdrs): 
            last_file = ( i == len(vcf_rdrs) -1 )

        
            # if current file position is ahead of bed position range, then pause operations
            if ((reader.pos_info["start"] > bed.pos_info["end"]) and (reader.pos_info["chrom"] == bed.pos_info["chrom"])) or \
                    reader.end_state:    

                reader.pause = True # pause current vcf from moving to the next line

            else:
                reader.pause = False

                # check to see if the vcf has fallen behind the bed
                while ((reader.pos_info["end"] < bed.pos_info["start"]) and (reader.pos_info["chrom"] == bed.pos_info["chrom"])):
                    print(f"{reader.name} skipping {list(reader.pos_info.values())}")
                    reader.read()

            
            if comps == 1 or comps == 2:

                if reader.pause:
                    bed_dist = None
                else:
                    # BED Ref Comparison Operations
                    # compare current vcf ref with bed ref
                    start_diff = abs(reader.pos_info["start"] - bed.pos_info["start"])
                    end_diff = abs(reader.pos_info["end"] - bed.pos_info["end"])

                    bed_dist = start_diff + end_diff

                bof_out_str = (f"{reader.name}\t{list(bed.pos_info.values())}\t{list(reader.pos_info.values())}\t{bed_dist}\n")
                bof.write(bof_out_str)

            if stateCheck(vcf_rdrs) and (comps == 1 or comps == 3): # and last_file
                gtd_sub = []
                psd_sub = []
                for other_reader in vcf_rdrs[i+1:]:
                    # compare 2 genotypes
                    gtd = compareGt(reader.genotype, other_reader.genotype)
                    gtd_sub.append(gtd)

                    # calculate difference in positions between vcf files
                    vcf_start_diff = abs(reader.pos_info["start"] - other_reader.pos_info["start"])
                    vcf_end_diff = abs(reader.pos_info["end"] - other_reader.pos_info["end"])
                    psd_sub.append(vcf_start_diff + vcf_end_diff)

                gt_dist_list.append(gtd_sub)  
                ps_dist_list.append(psd_sub)

                
       # vof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}\t{psd_sub}\t{gt_dist_list}\t{reader.genotype}\t{other_reader.genotype}\n"                  
        vof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}\t{psd_sub}\t{gt_dist_list}\n"                  
        vof.write(vof_out_str)

        # read current lines for all files
        [rdr.read() for rdr in vcf_rdrs] 
        [rdr.buildGt() for rdr in vcf_rdrs]         
        bed.read()

           


    # close all files
    bed.close_file()
    for rdr in vcf_rdrs:
        rdr.close_file()
    if comps == 1 or comps == 2:
        bof.close() 
    if comps == 1 or comps == 3:
        vof.close() 


    print("\n\n---PROGRAM COMPLETE---\n")



mainloop("test-isolated-vc-catalog.atarva.bed",
         ["HG001.PAW79146.haplotagged.URfix.atarva.vcf", 
         "HG001.strdust.vcf",
         "HG001.PAW79146.haplotagged.URfix.strkit.vcf",
         #"HG001.PAW79146.haplotagged.URfix.longTR.vcf",
         #"HG001.PAW79146.haplotagged.URfix.vcf",
         #"HG001.PAW79146.haplotagged.URfix.vamos.vcf",
         #"medaka_to_ref.TR.vcf"
         ])
