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
    offsets = [[-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0]]
    offsets = [[0,0]] * 7

    with ExitStack() as stack: 

        # File Prep
        bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))

        # create list of vcf reader objects
        for i, vcf in enumerate(vcf_list):
          
            try:
                # create VCFReader object for file.
                vcf_rdrs.append(SC_VCFReader(os.path.join(DATA_DIR, vcf), 
                                start_offset=offsets[i][0], 
                                end_offset=offsets[i][1]))


                # add vcf to exit stack which puts it under control of the with statement
                stack.enter_context(vcf_rdrs[i])

            except FileIOError as e:
                sys.exit(f"\nERROR\nExiting program due to file error: {e}")

            # move past header data
            vcf_rdrs[i].skipMetaData()

            # build first line's genotype and save it
            try: 
                vcf_rdrs[i].safeRead()
                vcf_rdrs[i].buildGt()
            except (FileReadError, VCFFormatError) as e:
                sys.exit(f"\nERROR\n{e}\nFrom file: {vcf_rdrs[i].path}")


        
        # print("Comparisons:\n1 (default): BED-VCF & VCF-VCF\n 2 : BED-VCF\n 3 : VCF-VCF")
        # #while comps > 3 and comps < 1:
        # comps = int(input(": "))
        comps = 1

        if comps == 1 or comps == 2:
            bof = stack.enter_context(open(os.path.join(LOCAL_DATA, "bed-comp.tsv"), "w")) # open file and put it into the exit stack
            bof.write("#CHROM\tSTART\tEND")

        if comps == 1 or comps == 3:
            vof = stack.enter_context(open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w"))
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

            bof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}"
            vof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}"
            gtd_sub_str = ""
            psd_sub_str = ""


            # cycle through the vcf files and perform operations on the current line
            for i, reader in enumerate(vcf_rdrs): 
                skip_count = 0

                # check VCF chromosome ordering
                reader.checkOrder()

                # Run Alingment Checks
                chrom_match = (reader.pos_info["chrom"] == bed.pos_info["chrom"])

                # if vcf position is behind the bed, or the vcf chrom is behind, loop until the vcf catches up         
                while (((reader.pos_info["end"] < bed.pos_info["start"]) and chrom_match) or (reader.pos_info["chrom"] < bed.pos_info["chrom"])) \
                    and not reader.end_state:

                    #if skip_count == 0:
                    #    print(f"\nWARNING: Skipping lines starting at {list(reader.pos_info.values())} from {reader.path}")

                    # move the file line forward until it is no longer behind, or the end of the file is reached
                    try: 
                        vcf_rdrs[i].safeRead()
                        vcf_rdrs[i].checkOrder()
                    except (FileReadError, VCFFormatError) as e:
                        sys.exit(f"\nERROR\n{e}\nFrom file: {vcf_rdrs[i].path}")

                    reader.skip_num += 1


                # if current vcf position is ahead of bed position range, or if the vcf chrom is ahead, then pause operations
                if ((reader.pos_info["start"] > bed.pos_info["end"]) and chrom_match) or \
                    (reader.pos_info["chrom"] > bed.pos_info["chrom"]):

                    reader.pause = True # pause current vcf from being able to move to the next line or run comparisons

                # if the current vcf is not ahead or the chromosomes dont match, continue moving forward
                else: 
                    reader.pause = False

                # VCF-BED Comparisons
                if comps == 1 or comps == 2:

                    if reader.pause or reader.end_state:
                        bof_out_str += "\tNA\tNA"
                    else:
                        # BED Ref Comparison Operations
                        # compare current vcf ref with bed ref
                        start_diff = reader.pos_info["start"] - bed.pos_info["start"]
                        end_diff = reader.pos_info["end"] - bed.pos_info["end"]
                        
                        if start_diff > 500 or end_diff > 500:
                            print(f"WARNING: Large Positional difference at {bed.pos_info}")

                        bof_out_str += f"\t{start_diff}\t{end_diff}"
                    

            # VCF-VCF Comparisons
            for i, reader in enumerate(vcf_rdrs):
                if (comps == 1 or comps == 3): 

                    for other_reader in vcf_rdrs[i+1:]:
                        # if both readers are not paused or ended
                        if stateCheck(reader) and stateCheck(other_reader):
                            # calculate lv distance of genotypes
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
            try:
                for rdr in vcf_rdrs:
                    rdr.safeRead()

                    if not rdr.end_state:
                        rdr.buildGt()   

                bed.read() 
            except (FileReadError, VCFFormatError, BEDFormatError) as e:
                sys.exit(f"\nERROR\n{e}\nFrom file: {vcf_rdrs[i].path}")
            

            
    for rdr in vcf_rdrs:
        # if any vcfs are still not at their end, then they are likely out of order
        if not rdr.end_state:
            print(f"\nWARNING: BED file finished before {reader.path}.\nPossible chromosome ordering error.")


        if rdr.skip_num > 0:
            print(f"\nWARNING: {rdr.skip_num} lines skipped in {rdr.path}")


    print("\n\n---PROGRAM COMPLETE---\n")



mainloop("test-isolated-vc-catalog.atarva.bed",
        [
        "HG001.PAW79146.haplotagged.URfix.atarva.vcf", 
        "HG001.strdust.vcf.gz",
        "HG001.PAW79146.haplotagged.URfix.strkit.vcf",
        "HG001.PAW79146.haplotagged.URfix.longTR.vcf.gz",
        "HG001.PAW79146.haplotagged.URfix.vcf",
        "HG001.PAW79146.haplotagged.URfix.vamos.vcf",
        "test_cases.vcf"
        "medaka_to_ref.TR.vcf"
        ])