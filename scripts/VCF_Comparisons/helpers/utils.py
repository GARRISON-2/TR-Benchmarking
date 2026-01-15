import Levenshtein as lv
import sys
from helpers.comp_readers import COMP_VCFReader
from helpers.readers import FileIOError
from helpers.constants import *



def compareString(str1, str2, method=COMP_METHOD.LEVENSHTEIN):
    allowed = {'A', 'T', 'C', 'G'}

    # if the string is not null and it only contains the characters A, T, C, or G
    if (str1 and set(str1).issubset(allowed)) and (str2 and set(str2).issubset(allowed)):
        if method == COMP_METHOD.LEVENSHTEIN:
            return lv.distance(str1, str2)
        
        if method == COMP_METHOD.LENGTH:
            return len(str1) - len(str2)
    
    else:
        return None


def setupMetadata(vlist):
    # output strings
    header_start = f"CHROM\tSTART\tEND"
    bdof_meta = "##<BED-VCF Position Comparison>\n"
    pdof_meta = "##<VCF Position Comparison>\n"
    lvdof_meta = "##<VCF Levenshtein Comparison>\n"
    ldof_meta = "##<VCF Length Comparison>\n"


    # Write metadata to output file
    for i, rdr in enumerate(vlist):
        file_name = getFileName(vlist[i].path)

        bdof_meta += (f"##FILE_{i}=<Name: {file_name}>\n")
        pdof_meta += (f"##FILE_{i}=<Name: {file_name}>\n")
        lvdof_meta += (f"##FILE_{i}=<Name: {file_name}>\n")
        ldof_meta += (f"##FILE_{i}=<Name: {file_name}>\n")
    
    bdof_meta += (f"##INFO=<CHROM: Chromosome of BED position>\n")
    pdof_meta += (f"##INFO=<CHROM: Chromosome of BED position>\n")
    lvdof_meta += (f"##INFO=<CHROM: Chromosome of BED position>\n")
    ldof_meta += (f"##INFO=<CHROM: Chromosome of BED position>\n")

    bdof_meta += (f"##INFO=<START: Start position from BED file>\n")
    pdof_meta += (f"##INFO=<START: Start position from BED file>\n")
    lvdof_meta += (f"##INFO=<START: Start position from BED file>\n")
    ldof_meta += (f"##INFO=<START: Start position from BED file>\n")

    bdof_meta += (f"##INFO=<END: End position from BED file>\n")
    pdof_meta += (f"##INFO=<END: End position from BED file>\n")
    lvdof_meta += (f"##INFO=<END: End position from BED file>\n")
    ldof_meta += (f"##INFO=<END: End position from BED file>\n")


    # Configure column headers for output files
    bdof_meta += (header_start)
    pdof_meta += (header_start)
    lvdof_meta += (header_start)
    ldof_meta += (header_start)
    for i in range(len(vlist)):
        for j in range(i+1, len(vlist)):
            pdof_meta += f"\tPSDIST_START_{i}-{j}\tPSDIST_END_{i}-{j}"
            lvdof_meta += f"\tLVDIST_ALL1_{i}-{j}\tLVDIST_ALL2_{i}-{j}"
            ldof_meta += f"\tLNDIST_ALL1_{i}-{j}\tLNDIST_ALL2_{i}-{j}"
        bdof_meta += f"\tBDDIST_START_{i}\tBDDIST_END_{i}"
    bdof_meta += "\n"
    pdof_meta += "\n"


    return bdof_meta, pdof_meta, lvdof_meta, ldof_meta


def cleanNum(num):
    if num is None:
        return 0
    else:
        return abs(num)
    

def NoneToNA(input):
    if input is None:
        return 'NA'
    else:
        return input


def compareGt(gt1, gt2, comp_method = COMP_METHOD.LEVENSHTEIN, comp_ord = None):

    # pad genotypes to length 2 for compatibility (eg. for handling hemizygous regions)
    gt1 += [None] * (2 - len(gt1))
    gt2 += [None] * (2 - len(gt2))


    if comp_ord == COMP_ORDER.VERTICAL:
        vert_dist = [compareString(gt1[0], gt2[0], comp_method), compareString(gt1[1], gt2[1], comp_method)]
        v_sum = sum(cleanNum(dist) for dist in vert_dist)

        return v_sum, NoneToNA(vert_dist[0]), NoneToNA(vert_dist[1]), comp_ord

    if comp_ord == COMP_ORDER.CROSS:
        cross_dist = [compareString(gt1[1], gt2[0], comp_method), compareString(gt1[0], gt2[1], comp_method)]
        c_sum = sum(cleanNum(dist) for dist in cross_dist)

        return c_sum,  NoneToNA(cross_dist[0]), NoneToNA(cross_dist[1]), comp_ord

    else: 
        # get comparisons for both permutations of comparing alleles
        vert_dist = [compareString(gt1[0], gt2[0], comp_method), compareString(gt1[1], gt2[1], comp_method)]
        cross_dist = [compareString(gt1[1], gt2[0], comp_method), compareString(gt1[0], gt2[1], comp_method)]

        v_sum = sum(cleanNum(dist) for dist in vert_dist)
        c_sum = sum(cleanNum(dist) for dist in cross_dist)

        # assume the lesser comparison value is the correct comparison order
        best_sum, best_dist, best_ord = (v_sum, vert_dist, COMP_ORDER.VERTICAL) if v_sum <= c_sum else (c_sum, cross_dist, COMP_ORDER.CROSS)

        # returns the sum of the comparisons between alleles, as well as the individual comparisons
        return best_sum, NoneToNA(best_dist[0]), NoneToNA(best_dist[1]), best_ord
    


def getFileName(path_str):
    # enumerate through the reversed file path to grab
    # the index of where the file name starts
    for i, letter in enumerate(reversed(path_str)):
        if letter == "\\":
            name_start = len(path_str) - i
            break

    # slice for the file name minus the path and the format
    return path_str[name_start:]#[name_start:-4]


def stateCheck(rdr):
    return (not rdr.pause) and (not rdr.end_state)


def setupVCFReader(vcf, stk, offset=[0,0], sc = None, skip_head = True, pos_only = False):
    try:
        # create VCFReader object for file.
        #if "vamos" in vcf.lower():
        
        
        if sc == None:
            rdr = COMP_VCFReader(file_path=vcf, 
                        start_offset=offset[0], 
                        end_offset=offset[1])
        else: 
            rdr = COMP_VCFReader(file_path=vcf, 
                        start_offset=offset[0], 
                        end_offset=offset[1],
                        sc=sc,
                        pos_only=sc.pos_only)


        # add vcf to exit stack 
        stk.enter_context(rdr)

    except FileIOError as e:
        sys.exit(f"\nERROR\nExiting program due to file error: {e}")


    if skip_head:
        # move past meta data
        rdr.skipMetaData()

    return rdr