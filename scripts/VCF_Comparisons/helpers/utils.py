import Levenshtein
from contextlib import ExitStack
import sys
from helpers.comp_readers import COMP_VCFReader, alleleData
from helpers.readers import FileIOError, FileReadError, VCFFormatError, BEDFormatError
from helpers.constants import *



def cleanNum(num: int):
    """
    Docstring for cleanNum
    
    :param num: Description
    :type num: int
    """
    if num is None:
        return 0
    else:
        return abs(num)


def compareAllele(all1: alleleData, all2: alleleData, method=COMP_METHOD.LEVENSHTEIN, trim = False):
    """
    Docstring for compareAllele
    
    :param all1: Description
    :type all1: str
    :param all2: Description
    :type all2: str
    :param method: Description
    :type method: COMP_METHOD
    """
    ALLOWED = {'A', 'T', 'C', 'G'}

    if all1 and all2:
        # if the allele string is not null and it only contains the characters A, T, C, or G
        if method == COMP_METHOD.LEVENSHTEIN and \
        (all1.allele_str and set(all1.allele_str).issubset(ALLOWED)) and \
        (all2.allele_str and set(all2.allele_str).issubset(ALLOWED)):
            
            if trim:
                return Levenshtein.distance(all1.allele_str[alleleData.start_trim:alleleData.end_trim], \
                                            all2.allele_str[alleleData.start_trim:alleleData.end_trim])
            else:
                return Levenshtein.distance(all1.allele_str, all2.allele_str)
            
        if method == COMP_METHOD.LENGTH and \
        all1.length > 0 and all2.length > 0:
            return all1.length - all2.length
    
    # default to return None if checks fail
    return None
    

def compareGt(gt1: list, gt2: list, comp_method = COMP_METHOD.LEVENSHTEIN, comp_ord = None, trim = False):
    """
    Docstring for compareGt
    
    :param gt1: Description
    :type gt1: str
    :param gt2: Description
    :type gt2: str
    :param comp_method: Description
    :param comp_ord: Description
    """
    # pad genotype list lengths to length 2 for compatibility (ie. for handling hemizygous regions)
    gt1 += [None] * (2 - len(gt1))
    gt2 += [None] * (2 - len(gt2))

    # if the comparison order is given:
    if comp_ord == COMP_ORDER.VERTICAL:
        vert_dist = [compareAllele(gt1[0], gt2[0], comp_method), compareAllele(gt1[1], gt2[1], comp_method)]
        v_sum = sum(cleanNum(dist) for dist in vert_dist)

        return v_sum, NoneToNA(vert_dist[0]), NoneToNA(vert_dist[1]), comp_ord

    if comp_ord == COMP_ORDER.CROSS:
        cross_dist = [compareAllele(gt1[1], gt2[0], comp_method), compareAllele(gt1[0], gt2[1], comp_method)]
        c_sum = sum(cleanNum(dist) for dist in cross_dist)

        return c_sum,  NoneToNA(cross_dist[0]), NoneToNA(cross_dist[1]), comp_ord

    # find the correct comparison order
    else: 
        # get comparisons for both permutations of comparing alleles
        vert_dist = [compareAllele(gt1[0], gt2[0], comp_method), compareAllele(gt1[1], gt2[1], comp_method)]
        cross_dist = [compareAllele(gt1[1], gt2[0], comp_method), compareAllele(gt1[0], gt2[1], comp_method)]

        v_sum = sum(cleanNum(dist) for dist in vert_dist)
        c_sum = sum(cleanNum(dist) for dist in cross_dist)

        # assume the lesser comparison value is the correct comparison order
        best_sum, best_dist, best_ord = (v_sum, vert_dist, COMP_ORDER.VERTICAL) if v_sum <= c_sum else (c_sum, cross_dist, COMP_ORDER.CROSS)

        # returns the sum of the comparisons between alleles, as well as the individual comparisons
        return best_sum, NoneToNA(best_dist[0]), NoneToNA(best_dist[1]), best_ord


def getFileName(path_str: str):
    """
    Docstring for getFileName
    
    :param path_str: Description
    :type path_str: str
    """
    # enumerate through the reversed file path to grab
    # the index of where the file name starts
    for i, letter in enumerate(reversed(path_str)):
        if letter == "\\":
            name_start = len(path_str) - i
            break

    # slice for the file name minus the path and the format
    return path_str[name_start:]#[name_start:-4]


def NoneToNA(input):
    """
    Docstring for NoneToNA
    
    :param input: Description
    """
    if input is None:
        return 'NA'
    else:
        return input


def setupMetadata(vlist: list[COMP_VCFReader], header_only = False):
    # output strings
    header_start = f"CHROM\tSTART\tEND"
    bdof_meta = "##<BED-VCF Position Comparison>\n"
    pdof_meta = "##<VCF Position Comparison>\n"
    lvdof_meta = "##<VCF Levenshtein Comparison>\n"
    ldof_meta = "##<VCF Length Comparison>\n"

    if not header_only:
        # Write metadata to output file
        for i, rdr in enumerate(vlist):
            file_name = getFileName(vlist[i].path)

            bdof_meta += (f"##FILE_{i}=<Name: {file_name}; Start Offset: {vlist[i].start_off}; End Offset: {vlist[i].end_off}>\n")
            pdof_meta += (f"##FILE_{i}=<Name: {file_name}; Start Offset: {vlist[i].start_off}; End Offset: {vlist[i].end_off}>\n")
            lvdof_meta += (f"##FILE_{i}=<Name: {file_name}; Start Offset: {vlist[i].start_off}; End Offset: {vlist[i].end_off}>\n")
            ldof_meta += (f"##FILE_{i}=<Name: {file_name}; Start Offset: {vlist[i].start_off}; End Offset: {vlist[i].end_off}>\n")
        
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

        bdof_meta += (f"##INFO=<BDDIST_START_i: BED start - FILE_i start>\n")  
        bdof_meta += (f"##INFO=<BDDIST_END_i: BED end - FILE_i end>\n")
        pdof_meta += (f"##INFO=<PSDIST_START_i-j: FILE_i start - FILE_j start>\n")
        pdof_meta += (f"##INFO=<PSDIST_END_i-j: FILE_i end - FILE_j end>\n")
        lvdof_meta += (f"##INFO=<LVDIST_ALL1_i-j: Levenshtein distance bewteen FILE_i allele 1 and FILE_j allele 1>\n")
        lvdof_meta += (f"##INFO=<LVDIST_ALL2_i-j: Levenshtein distance bewteen FILE_i allele 2 and FILE_j allele 2>\n")
        ldof_meta += (f"##INFO=<LNDIST_ALL1_i-j: FILE_i allele 1 length - FILE_j allele 1 length>\n")
        ldof_meta += (f"##INFO=<LNDIST_ALL2_i-j: FILE_i allele 2 length - FILE_j allele 2 length>\n")

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
    lvdof_meta += "\n"
    ldof_meta += "\n"


    return bdof_meta, pdof_meta, lvdof_meta, ldof_meta


def setupVCFReader(vcf: str, stk: ExitStack, offset=[0,0], sc = None, skip_head = True, pos_only = False):
    """
    Docstring for setupVCFReader
    
    :param vcf: Description
    :type vcf: str
    :param stk: Description
    :type stk: ExitStack
    :param offset: Description
    :param sc: Description
    :param skip_head: Description
    :param pos_only: Description
    """
    
    
    try:
        # create VCFReader object for file.  
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
        rdr.skipMetaData(end_delimiter="#CHROM")

    return rdr


def stateCheck(rdr: COMP_VCFReader):
    """
    Docstring for stateCheck
    
    :param rdr: Description
    :type rdr: COMP_VCFReader
    """
    return (not rdr.pause) and (not rdr.end_state)
