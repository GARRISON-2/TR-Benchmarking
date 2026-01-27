import os
import time
from contextlib import ExitStack
from helpers.readers import BEDReader
from helpers.utils import *
from helpers.constants import *



def mainloop():
    # PROGRAM SETTINGS/VARIABLES
    '''    vcf_list = []
    with os.scandir(DATA_DIR) as dir_entries:
        for entry in dir_entries:
            if entry.is_file():
                if entry.path.endswith(".vcf") or entry.path.endswith(".vcf.gz"):
                    if entry.path.lower().contains("straglr"):
                        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.straglr.sorted.vcf"), [0,0], SPECIAL_CASE.STRAGLR],
                    elif entry.path.lower().contains("vamos"):
                        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.vamos.sorted.vcf"), [0,0], SPECIAL_CASE.VAMOS],
                    else:
                        vcf_list.append([entry.path, [0,0], None])
    '''
    
    # Input File Paths
    bed_path = os.path.join(DATA_DIR, "BED_files\\benchmark-catalog-v2.vamos.bed")   
    vcf_list = [
        # the file name, the offset amount (eg. [-1, 0] for 1 based inclusive), and whether it needs special parsing parameters    
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.atarva.sorted.vcf"), [-1, 0], None], 
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.strdust.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.longTR.sorted.vcf"), [-1,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.straglr.sorted.vcf"), [0,0], SPECIAL_CASE.STRAGLR],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.vamos.sorted.vcf"), [0,0], SPECIAL_CASE.VAMOS],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.strkit.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.medaka_to_ref.TR.sorted.vcf"), [-1,0], None]
        ]
    # Output File Names
    bed_comp_file = 'bed-comp.tsv'
    position_comp_file = 'pos-comp.tsv'
    levenshtein_comp_file = 'lev-comp.tsv'
    length_comp_file = 'len-comp.tsv'
    # Program Options
    trim_alleles = False


    with ExitStack() as stack: 
        vcf_rdrs = []

        # create bed reader and enter the file into the stack
        bed = stack.enter_context(BEDReader(bed_path))
        bed.read()
        bed.skipMetaData()

        # setup list of vcf reader objects, and put them in file stack
        for i, vcf_info in enumerate(vcf_list):
            vcf_rdrs.append(setupVCFReader(vcf=vcf_info[0], 
                                           offset=vcf_info[1], 
                                           sc=vcf_info[2],
                                           stk=stack))
            
            # build first line's genotype and save it      
            vcf_rdrs[i].VCFParse()
            vcf_rdrs[i].buildGtData()


        # open output files and put them into the exit stack
        bdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, bed_comp_file), "w")) 
        pdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, position_comp_file), "w"))
        lvdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, levenshtein_comp_file), "w"))
        ldof = stack.enter_context(open(os.path.join(OUTPUT_DIR, length_comp_file), "w"))
        
        # write metadata to output files
        bdof_meta, pdof_meta, lvdof_meta, ldof_meta = setupMetadata(vcf_rdrs, header_only=False)
        bdof.write(bdof_meta)
        pdof.write(pdof_meta)
        lvdof.write(lvdof_meta)
        ldof.write(ldof_meta)
        

        str_time = time.perf_counter()

        # Main Operations loop       
        while bed.cur_line: # loop until BED file reaches end
            bed_pos_str = f"{bed.chrom}\t{bed.pos}\t{bed.end_pos}"
            bdof_out_str = bed_pos_str
            pdof_out_str = bed_pos_str
            lvdof_out_str = bed_pos_str
            ldof_out_str = bed_pos_str


            if bed.prev_line is None or bed.chrom != bed.prev_line[0]:
                print(f"Comparing {bed.chrom}")


            # cycle through all vcf files and ensure they are synced to the bed
            [reader.syncToBed(bed) for reader in vcf_rdrs]
                                      
            # Run comparisons on each VCF
            for i, reader in enumerate(vcf_rdrs):

                # VCF-BED Comparisons
                if reader.pause or reader.end_state: # if the vcf skipped the current line or has ended
                    bdof_out_str += "\tNA\tNA"
                else:
                    # BDDIST: compre vcf ref position with bed
                    start_diff = bed.pos - reader.pos
                    end_diff = bed.end_pos - reader.end_pos
                    
                    # add start_diff and end_diff to object
                    if trim_alleles:
                        for allele in reader.gt_data:
                            # trim only if the current tool has widened the position
                            allele.start_trim = start_diff if start_diff > 0 else 0
                            allele.end_trim = -end_diff if end_diff > 0 else 0

                    bdof_out_str += f"\t{start_diff}\t{end_diff}"


                # VCF-VCF Comparisons
                for other_reader in vcf_rdrs[i + 1:]:
                    # if both readers are not paused or ended
                    if stateCheck(reader) and stateCheck(other_reader):
                        # LVDIST: calculate levenshtein distance of alleles between vcf files
                        gt_lvdiff, a1_lvdiff, a2_lvdiff, order = compareGt(reader.gt_data, 
                                                                    other_reader.gt_data,
                                                                    comp_method=COMP_METHOD.LEVENSHTEIN,
                                                                    trim=trim_alleles)
                        lvdof_out_str += f"\t{a1_lvdiff}\t{a2_lvdiff}"
                    
                        # LENDIST: calculate difference in allele lengths between vcf files
                        gt_ldiff, a1_ldiff, a2_ldiff, order = compareGt(reader.gt_data, 
                                                                other_reader.gt_data, 
                                                                comp_method=COMP_METHOD.LENGTH,
                                                                comp_ord=order,
                                                                trim=trim_alleles)
                        ldof_out_str += f"\t{a1_ldiff}\t{a2_ldiff}"

                        # if a1_ldiff > a1_lvdiff or a2_ldiff > a2_lvdiff:
                        #    raise Exception("\nFATAL PROGRAM ERROR\nLength difference between strings greater than Levenshtein distance.")

                        # POSDIST: calculate difference in positions between vcf files
                        vcf_start_diff = reader.pos - other_reader.pos
                        vcf_end_diff = reader.end_pos - other_reader.end_pos
                        pdof_out_str += f"\t{vcf_start_diff}\t{vcf_end_diff}"

                    else:
                        lvdof_out_str += "\tNA\tNA"
                        ldof_out_str += "\tNA\tNA"
                        pdof_out_str += "\tNA\tNA"


            # write data to output files
            bdof.write(bdof_out_str + "\n")   
            pdof.write(pdof_out_str + "\n")
            lvdof.write(lvdof_out_str + "\n")   
            ldof.write(ldof_out_str + "\n")


            # read lines for all files
            for rdr in vcf_rdrs:
                rdr.VCFParse()

                if not rdr.end_state:
                    rdr.buildGtData()   

            bed.read() 

            

    # End of Program checks
    for rdr in vcf_rdrs:
        # if any vcfs are still not at their end, then they are likely out of order
        if not rdr.end_state:
            print(f"\nWARNING: BED file finished before {rdr.path}.\nPossible chromosome ordering error.")

        # if any lines were skipped in the file, print a warning
        if rdr.skip_num > 0:
            print(f"\nWARNING: {rdr.skip_num} lines skipped in {rdr.path}")


    end_time = time.perf_counter()
    comp_time = end_time - str_time


    print("\n\n---PROGRAM COMPLETE---\n")
    print(f"Comparison time: {round(comp_time, 4)}")



mainloop()