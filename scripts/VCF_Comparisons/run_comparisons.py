import os
import sys
from contextlib import ExitStack
from helpers.readers import *
from helpers.comp_readers import *
from helpers.utils import *

# set directory variables for easy file i/o while testing
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'local_data')
LOCAL = os.path.join(PROJ_ROOT, 'scripts\\VCF_Comparisons')
#LOCAL = os.path.join(PROJ_ROOT, 'local_data')



def mainloop():
    bed_file = "test-isolated-vc-catalog.atarva.bed.gz"
    
    vcf_list = [ # the file name, and the offset amount (eg. [-1, 0] for 1 based inclusive), and bool for whether is is position only
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.atarva.vcf"), [0, 0], False], 
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.strdust.sorted.vcf"), [0,0], False],
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.longTR.sorted.vcf"), [0,0], False],
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.straglr.vcf"), [0,0], True],
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.vamos.vcf"), [0,0], True],
        [os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.strkit.vcf"), [0,0], False],
        [os.path.join(DATA_DIR, "medaka_to_ref.TR.sorted.vcf"), [0,0], False]
       # ["test_cases.vcf", [0,0], False]
        ]



    
    with ExitStack() as stack: 
        vcf_rdrs = []
        # output strings
        header_start = f"#CHROM\tSTART\tEND"
        bof_header = ""
        ps_header = ""
        gt_header = ""

        # create bed reader and enter the file into the stack, putting it under control of the with statement
        bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))
        bed.skipMetaData()

        # create list of vcf reader objects
        for i, vcf_info in enumerate(vcf_list):
            vcf_rdrs.append(setupVCFReader(vcf=vcf_info[0], 
                                           offset=vcf_info[1], 
                                           pos_only=vcf_info[2],
                                           stk=stack))



        # open output files and put it into the exit stack
        bof = stack.enter_context(open(os.path.join(DATA_DIR, "bed-comp.tsv"), "w")) 
        vof = stack.enter_context(open(os.path.join(DATA_DIR, "vcf-comp.tsv"), "w"))
        

        # Write metadata to output file
        for i, rdr in enumerate(vcf_rdrs):
            file_name = getFileName(vcf_rdrs[i].path)

            bof.write(f"##FILE_{i}=<Name={file_name}>\n")
            vof.write(f"##FILE_{i}=<Name={file_name}>\n")



        # Configure headers for output files
        for i in range(len(vcf_rdrs)):
            for j in range(i+1, len(vcf_rdrs)):
                ps_header += f"\tPSDIST_START_{i}-{j}\tPSDIST_END_{i}-{j}"
                gt_header += f"\tLVDIST_ALL1_{i}-{j}\tLVDIST_ALL2_{i}-{j}"
            bof_header += f"\tBDDIST_START_{i}\tBDDIST_END_{i}"
        bof.write(header_start + bof_header + "\n")
        vof.write(header_start + ps_header + gt_header + "\n")


        # Main Operations loop       
        while not bed.end_state: # loop until the bed file has reached its end
            bof_out_str = f"{bed.chrom}\t{bed.pos}\t{bed.end_pos}"
            vof_out_str = f"{bed.chrom}\t{bed.pos}\t{bed.end_pos}"
            gtd_sub_str = ""
            psd_sub_str = ""


            # cycle through all vcf files and ensure they are synced to the bed
            [reader.syncToBed(bed) for reader in vcf_rdrs]
                                      

            # Run comparisons on each VCF
            for i, reader in enumerate(vcf_rdrs):

                # VCF-BED Comparisons
                if reader.pause or reader.end_state: # if the vcf skipped the current line or has ended
                    bof_out_str += "\tNA\tNA"
                else:
                    # POSDIST: compre vcf ref position with bed
                    start_diff = reader.pos - bed.pos
                    end_diff = reader.end_pos - bed.end_pos
                    
                    if start_diff > 500 or end_diff > 500:
                        print(f"WARNING: Large Positional difference at {bed.pos} from {reader.path}")

                    bof_out_str += f"\t{start_diff}\t{end_diff}"


                # VCF-VCF Comparisons
                for other_reader in vcf_rdrs[i+1:]:
                    # if both readers are not paused or ended
                    if stateCheck(reader) and stateCheck(other_reader):
                        # LVDIST: calculate levenshtein distance of genotypes
                        gt_diff, a1_diff, a2_diff = compareGt(reader.genotype, other_reader.genotype)
                        gtd_sub_str += f"\t{a1_diff}\t{a2_diff}"
                    
                        # POSDIST: calculate difference in positions between vcf files
                        vcf_start_diff = reader.pos - other_reader.pos
                        vcf_end_diff = reader.end_pos - other_reader.end_pos
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



mainloop()