import os
import Levenshtein as lv
from contextlib import ExitStack
from readers import *

# set directory variables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH


def compareString(str1, str2):
    allowed = {'A', 'T', 'C', 'G'}
    test_string = "ATCGX"

    # if the string is not null and it doesnt contain any characters but A, T, C, or G
    if (str1 and set(str1).issubset(allowed)) and (str2 and set(str2).issubset(allowed)):
        return lv.distance(str1, str2)
    else:

    # get the positional difference and calculate the lv distance based only on positions in the actual reference
    # start_diff = str2.pos_info["start"] - str1.pos_info["start"]
    # end_diff = len(str1.ref) - str2.pos_info["end"] - str1.pos_info["end"]

    # dist2 = lv.distance(str1.ref[start_diff: end_diff], str2.ref[start_diff: end_diff])
    
    # return_dict = {
    #     "dist": dist,
    #     "dist2": dist2
    # }

        return None


def NoneToZero(num):
    if num is None:
        return 0
    else:
        return num
    

def NoneToNA(input):
    if not input:
        return 'NA'
    else:
        return input


def compareGt(gt1, gt2):
    
    # get distances for both permutations of comparing alleles
    vert_dist = [compareString(gt1[0], gt2[0]), compareString(gt1[1], gt2[1])]
    cross_dist = [compareString(gt1[1], gt2[0]), compareString(gt1[0], gt2[1])]

    v_sum = NoneToZero(vert_dist[0]) + NoneToZero(vert_dist[1])
    c_sum = NoneToZero(cross_dist[0]) + NoneToZero(cross_dist[1])

    # if the allele order is mismatched between the genotypes
    if v_sum > c_sum:

        # use the cross distance since it is less,
        # which means it is more likely to be the correct comparison ordering
        return c_sum, NoneToNA(cross_dist[0]), NoneToNA(cross_dist[1])

    else:
        return v_sum, NoneToNA(cross_dist[0]), NoneToNA(cross_dist[1])


def stateCheck(rdr):
    return (not rdr.pause) and (not rdr.end_state)


def mainloop(bed_file, vcf_list):
    vcf_rdrs = []
    comps = 0
    offsets = [[-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0], [-1, 0]]


    with ExitStack() as stack: 

        # File Prep
        bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))

        # create list of vcf reader objects
        for i, vcf in enumerate(vcf_list):
                # create VCFReader object for file. 
                vcf_rdrs.append(VCFReader(os.path.join(DATA_DIR, vcf), 
                                    start_offset=offsets[i][0], 
                                    end_offset=offsets[i][1]))
                # add vcf to file stack which puts it under the with statement
                stack.enter_context(vcf_rdrs[i])
                # move past header data
                vcf_rdrs[i].skipMetaData()

                # build first line's genotype and save it
                vcf_rdrs[i].buildGt()

        
        # print("Comparisons:\n1 (default): BED-VCF & VCF-VCF\n 2 : BED-VCF\n 3 : VCF-VCF")
        # #while comps > 3 and comps < 1:
        # comps = int(input(": "))
        comps = 1

        if comps == 1 or comps == 2:
            bof = open(os.path.join(LOCAL_DATA, "bed-comp.tsv"), "w")
            bof.write("#CHROM\tSTART\tEND")

        if comps == 1 or comps == 3:
            vof = open(os.path.join(LOCAL_DATA, "vcf-comp.tsv"), "w")
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
        # loop until the bed file has reached its end
        while not bed.end_state:
        #while not all(end_states): or loop until all files are done
            bof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}"
            vof_out_str = f"{bed.pos_info["chrom"]}\t{bed.pos_info["start"]}\t{bed.pos_info["end"]}"
            gtd_sub_str = ""
            psd_sub_str = ""


            # cycle through the vcf files and perform operations on the current line
            for i, reader in enumerate(vcf_rdrs): 

                if not reader.prev_line[0].startswith("#"):
                    prev_chrom = reader.prev_line[reader.fields["CHROM"]]

                    # ensure vcf is in order
                    if reader.pos_info["chrom"] < prev_chrom:
                        sys.exit(f"\nERROR\n{reader.path} using unknown order. Ending Program.")


                # Run Alingment Checks
                chrom_match = (reader.pos_info["chrom"] == bed.pos_info["chrom"])
                # if current vcf position is ahead of bed position range, or if the vcf chrom is ahead, then pause operations
                if ((reader.pos_info["start"] > bed.pos_info["end"]) and chrom_match) or \
                    (reader.pos_info["chrom"] > bed.pos_info["chrom"]):

                    reader.pause = True # pause current vcf from being able to move to the next line or run comparisons

                # if the current vcf is not ahead or the chromosomes dont match, continue moving forward
                else: 
                    reader.pause = False

                    # if vcf position is behind the bed, or the vcf chrom is behind, loop until the vcf catches up         
                    while (((reader.pos_info["end"] < bed.pos_info["start"]) and chrom_match) or (reader.pos_info["chrom"] < bed.pos_info["chrom"])) \
                        and not reader.end_state:

                        # move the file line forward until it is no longer behind, or the end of the file is reached
                        print(f"\nWARNING: Skipping {list(reader.pos_info.values())} from {reader.path}")
                        reader.read()


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
                            print(f"WARNING-Large Positional difference at {bed.pos_info}")

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
            bed.read()
            for rdr in vcf_rdrs:
                rdr.read()

                if not rdr.end_state:
                    rdr.buildGt()     

            
    for rdr in vcf_rdrs:
        # if any vcfs are still not at their end, then they are likely out of order
        if not rdr.end_state:
            print(f"\nWARNING: BED file finished before {reader.path}.\n Check file order.")


    print("\n\n---PROGRAM COMPLETE---\n")



mainloop("test-isolated-vc-catalog.atarva.bed",
        [
        #"HG001.PAW79146.haplotagged.URfix.atarva.vcf", 
        #"HG001.strdust.vcf.gz",
        #"HG001.PAW79146.haplotagged.URfix.strkit.vcf",
        "HG001.PAW79146.haplotagged.URfix.longTR_copy.vcf",
        #"HG001.PAW79146.haplotagged.URfix.vcf",
        #"HG001.PAW79146.haplotagged.URfix.vamos.vcf",
        #"medaka_to_ref.TR.vcf"
        ])
