import Levenshtein as lv
from readers import *
import sys
import os



class SC_VCFReader(VCFReader):
    def __init__(self, file_path, start_offset = 0, end_offset = 0, pos_only = False, pause = False, lines_skipped = 0):
        super().__init__(file_path, start_offset, end_offset, pos_only)
        self.pause = pause
        self.skip_num = lines_skipped


    def safeRead(self):
        read_state = False

        try: 
            if not self.pause:  
                read_state = self.read()     

        except (FileReadError, VCFFormatError, BEDFormatError) as e:
            sys.exit(f"\nERROR\n{e}")   

        return read_state


    def checkOrder(self, order_method="ASCII"):
        # check VCF Chromosome ordering
        if not self.prev_line[0].startswith("#"):
            prev_chrom = self.prev_line[0]

            # ensure vcf is in order
            if self.pos_info['chrom'] < prev_chrom:
                sys.exit(f"\nERROR\n{self.path} using unknown order. Ending program.")


    def syncToBed(self, bp_info):
        # check VCF chromosome ordering
        self.checkOrder()

        # Run Alingment Checks
        chrom_match = (self.pos_info["chrom"] == bp_info["chrom"])

        # if vcf position is behind the bed, or the vcf chrom is behind, loop until the vcf catches up         
        while (((self.pos_info["end"] < bp_info["start"]) and chrom_match) or (self.pos_info["chrom"] < bp_info["chrom"])) \
            and not self.end_state:

            # move the file line forward until it is no longer behind, or the end of the file is reached
            self.safeRead()
            self.checkOrder()

            self.skip_num += 1


        # if vcf position is ahead of bed position range, or if the vcf chrom is ahead, then pause operations
        if ((self.pos_info["start"] > bp_info["end"]) and chrom_match) or \
            (self.pos_info["chrom"] > bp_info["chrom"]):

            self.pause = True # pause vcf from being able to move to the next line or run comparisons

        # if the vcf is not ahead or the chromosomes dont match, continue moving forward
        else: 
            self.pause = False





def compareString(str1, str2, method="lv"):
    allowed = {'A', 'T', 'C', 'G'}

    # if the string is not null and it only contains the characters A, T, C, or G
    if (str1 and set(str1).issubset(allowed)) and (str2 and set(str2).issubset(allowed)):
        if method == "lv":
            return lv.distance(str1, str2)
    
    else:
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

    # pad genotypes to length 2 for compatibility (eg. for handling hemizygous regions)
    gt1 += [None] * (2 - len(gt1))
    gt2 += [None] * (2 - len(gt2))

    # get distances for both permutations of comparing alleles
    vert_dist = [compareString(gt1[0], gt2[0]), compareString(gt1[1], gt2[1])]
    cross_dist = [compareString(gt1[1], gt2[0]), compareString(gt1[0], gt2[1])]

    v_sum = sum(NoneToZero(dist) for dist in vert_dist)
    c_sum = sum(NoneToZero(dist) for dist in cross_dist)

    # assume the lesser distance is the correct comparison order
    best_sum, best_dist = (v_sum, vert_dist) if v_sum <= c_sum else (c_sum, cross_dist)

	# returns the sum of the differences between alleles, as well as the individual differences
    return best_sum, NoneToNA(best_dist[0]), NoneToNA(best_dist[1])

'''
def compareGCContent(str1, str2):
    str1_gc = [count for count in str1[]]
    str2_gc = [count for count in str2[]]

    return str1_gc - str2_gc
''' 


def stateCheck(rdr):
    return (not rdr.pause) and (not rdr.end_state)
