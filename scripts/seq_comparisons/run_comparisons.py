import os
import sys
from contextlib import ExitStack
from readers import *
from utils import *

# set directory variables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH




def mainloop(bed_file, vcf_list):
    vcf_rdrs = []
    comps = 0

    with ExitStack() as stack: 

        # File Prep
        bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))

        # create list of vcf reader objects
        for i, vcf_info in enumerate(vcf_list):
            try:
                # create VCFReader object for file.
                vcf_rdrs.append(SC_VCFReader(os.path.join(DATA_DIR, vcf_info[0]), 
                                start_offset=vcf_info[1][0], 
                                end_offset=vcf_info[1][1]))


                # add vcf to exit stack which puts it under control of the with statement
                stack.enter_context(vcf_rdrs[i])

            except FileIOError as e:
                sys.exit(f"\nERROR\nExiting program due to file error: {e}")

            # move past header data
            vcf_rdrs[i].skipMetaData()


            # TEMP FIX
            # HG001.PAW79146.haplotagged.URfix.vcf & HG001.PAW79146.haplotagged.URfix.vamos.vcf are currently both position only
            if i == 3 or i == 4:
                vcf_rdrs[i].pos_only = True


            # build first line's genotype and save it      
            vcf_rdrs[i].safeRead()
            vcf_rdrs[i].buildGt()


        # open file and put it into the exit stack
        bof = stack.enter_context(open(os.path.join(LOCAL_DATA, "bed-comp.tsv"), "w")) 
        vof = stack.enter_context(open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w"))
        

        # Configure headers for output files
        bof.write("#CHROM\tSTART\tEND")
        vof.write(f"#CHROM\tSTART\tEND")
        temp1 = ""
        temp2 = ""
        for i, rdr_1 in enumerate(vcf_rdrs):
            for j, rdr_j in enumerate(vcf_rdrs[i+1:], i+1):
                temp1 += f"\tPSDIST_START_{i}-{j}\tPSDIST_END_{i}-{j}"
                temp2 += f"\tGTDIST_ALL1_{i}-{j}\tGTDIST_ALL2_{i}-{j}"
            bof.write(f"\tBDDIST_START_{i}\tBDDIST_END_{i}")
        bof.write("\n")
        vof.write(temp1 + temp2 + "\n")


        # Write metadata to output file
        


        # Main Operations loop       
        while not bed.end_state: # loop until the bed file has reached its end

            bof_out_str = f"{bed.pos_info['chrom']}\t{bed.pos_info['start']}\t{bed.pos_info['end']}"
            vof_out_str = f"{bed.pos_info['chrom']}\t{bed.pos_info['start']}\t{bed.pos_info['end']}"
            gtd_sub_str = ""
            psd_sub_str = ""


            # cycle through all vcf files and ensure they are synced to the bed before performing comparisons
            for i, reader in enumerate(vcf_rdrs): 

                # run alignment logic to ensure VCF and bed file are synced to the same positions
                reader.syncToBed(bed.pos_info)              


            # Run comparisons on each VCF
            for i, reader in enumerate(vcf_rdrs):

                # VCF-BED Comparisons
                if reader.pause or reader.end_state:
                    bof_out_str += "\tNA\tNA"
                else:
                    start_diff = reader.pos_info["start"] - bed.pos_info["start"]
                    end_diff = reader.pos_info["end"] - bed.pos_info["end"]
                    
                    if start_diff > 500 or end_diff > 500:
                        print(f"WARNING: Large Positional difference at {bed.pos_info}")

                    bof_out_str += f"\t{start_diff}\t{end_diff}"


                # VCF-VCF Comparisons
                for other_reader in vcf_rdrs[i+1:]:
                    # if both readers are not paused or ended
                    if stateCheck(reader) and stateCheck(other_reader):
                        # calculate levenshtein distance of genotypes
                        gt_diff, a1_diff, a2_diff = compareGt(reader.genotype, other_reader.genotype)
                        gtd_sub_str += f"\t{a1_diff}\t{a2_diff}"
                    
                        # calculate difference in positions between vcf files
                        vcf_start_diff = reader.pos_info["start"] - other_reader.pos_info["start"]
                        vcf_end_diff = reader.pos_info["end"] - other_reader.pos_info["end"]
                        psd_sub_str += f"\t{vcf_start_diff}\t{vcf_end_diff}"
                    else:
                        gtd_sub_str += "\tNA\tNA"
                        psd_sub_str += "\tNA\tNA"


            # write data to output files
            bof.write(bof_out_str + "\n")   
            vof.write(vof_out_str + psd_sub_str + gtd_sub_str + "\n")


            # read lines for all files
            for rdr in vcf_rdrs:
                rdr.safeRead()

                if not rdr.end_state:
                    rdr.buildGt()   

            bed.read() 

            

    # End of Program checks
    for rdr in vcf_rdrs:
        # if any vcfs are still not at their end, then they are likely out of order
        if not rdr.end_state:
            print(f"\nWARNING: BED file finished before {rdr.path}.\nPossible chromosome ordering error.")

        # if any lines were skipped in the file, print a warning
        if rdr.skip_num > 0:
            print(f"\nWARNING: {rdr.skip_num} lines skipped in {rdr.path}")


    print("\n\n---PROGRAM COMPLETE---\n")



mainloop("test-isolated-vc-catalog.atarva.bed.gz",
        [
        ["HG001.PAW79146.haplotagged.URfix.atarva.vcf", [0,0]], # the file name, and the offset amount (eg. [-1, 0] for 1 based inclusive)
        ["HG001.strdust.vcf.gz", [0,0]],
        ["HG001.PAW79146.haplotagged.URfix.strkit.vcf", [0,0]],
       # ["HG001.PAW79146.haplotagged.URfix.longTR.vcf.gz", [0,0]],
        ["HG001.PAW79146.haplotagged.URfix.vcf", [0,0]],
        ["HG001.PAW79146.haplotagged.URfix.vamos.vcf", [0,0]],
       # ["medaka_to_ref.TR.vcf", [0,0]]
       # ["test_cases.vcf", [0,0]]
        ])