import os
import time
from contextlib import ExitStack
from helpers.readers import BEDReader
from helpers.utils import *
from helpers.constants import *



def mainloop():
    bed_path = os.path.join(DATA_DIR, "BED_files\\benchmark-catalog-v2.vamos.bed")

    '''    # vcf_list = []
    # with os.scandir(DATA_DIR) as dir_entries:
    #     for entry in dir_entries:
    #         if entry.is_file():
    #             if entry.path.endswith(".vcf") or entry.path.endswith(".vcf.gz"):
    #                 vcf_list.append([entry.path, [0,0], False])
    '''

    '''
    # or can manually set the list using this format: 
    # the file name, the offset amount (eg. [-1, 0] for 1 based inclusive), and bool for whether the vcf is position only
    # vcf_list = [ 
    #     ["fileq.vcf", [0, 0], False], 
    #     ["file2.vcf", [0, 1], True],
    #     ]
    '''

    vcf_list = [ # the file name, and the offset amount (eg. [-1, 0] for 1 based inclusive), and whether it needs special parsing parameters
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.atarva.sorted.vcf"), [0, 0], None], 
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.strdust.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.longTR.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.straglr.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.vamos.sorted.vcf"), [0,0], SPECIAL_CASE.VAMOS],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.30x.haplotagged.strkit.sorted.vcf"), [0,0], None],
        [os.path.join(DATA_DIR, "HG007.30x\\HG007.medaka_to_ref.TR.sorted.vcf"), [0,0], None]
       # ["test_cases.vcf", [0,0], False]
        ]


    with ExitStack() as stack: 
        vcf_rdrs = []

        # create bed reader and enter the file into the stack, putting it under control of the with statement
        bed = stack.enter_context(BEDReader(bed_path))
        bed.read()
        bed.skipMetaData()

        # create list of vcf reader objects
        for i, vcf_info in enumerate(vcf_list):
            vcf_rdrs.append(setupVCFReader(vcf=vcf_info[0], 
                                           offset=vcf_info[1], 
                                           sc=vcf_info[2],
                                           stk=stack))
            
            # build first line's genotype and save it      
            vcf_rdrs[i].safeRead()
            vcf_rdrs[i].buildGt()


        # open output files and put it into the exit stack
        bdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, "bed-comp.tsv"), "w")) 
        pdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, "pos-comp.tsv"), "w"))
        lvdof = stack.enter_context(open(os.path.join(OUTPUT_DIR, "lev-comp.tsv"), "w"))
        ldof = stack.enter_context(open(os.path.join(OUTPUT_DIR, "len-comp.tsv"), "w"))
        

        bdof_meta, pdof_meta, lvdof_meta, ldof_meta = setupMetadata(vcf_rdrs)
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
                    
                    if start_diff > 500 or end_diff > 500:
                        print(f"WARNING: Large Positional difference at {bed.pos} from {reader.path}")

                    bdof_out_str += f"\t{start_diff}\t{end_diff}"


                # VCF-VCF Comparisons
                for other_reader in vcf_rdrs[i+1:]:
                    # if both readers are not paused or ended
                    if stateCheck(reader) and stateCheck(other_reader):
                        # LVDIST: calculate levenshtein distance of alleles between vcf files
                        gt_lvdiff, a1_lvdiff, a2_lvdiff, order = compareGt(reader.genotype, 
                                                                    other_reader.genotype)
                        lvdof_out_str += f"\t{a1_lvdiff}\t{a2_lvdiff}"
                    
                        # LENDIST: calculate difference in allele lengths between vcf files
                        gt_ldiff, a1_ldiff, a2_ldiff, order = compareGt(reader.genotype, 
                                                                 other_reader.genotype, 
                                                                 comp_method=COMP_METHOD.LENGTH,
                                                                 comp_ord=order)
                        ldof_out_str += f"\t{a1_ldiff}\t{a2_ldiff}"

                        if a1_ldiff > a1_lvdiff or a2_ldiff > a2_lvdiff:
                            raise Exception("\nFATAL PROGRAM ERROR\nLength difference between strings greater than Levenshtein distance.")

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


    end_time = time.perf_counter()
    comp_time = end_time - str_time


    print("\n\n---PROGRAM COMPLETE---\n")
    print(f"Comparison time: {round(comp_time, 4)}")



mainloop()