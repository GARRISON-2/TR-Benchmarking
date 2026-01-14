import os
import sys
from contextlib import ExitStack
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(root_dir))

from helpers.utils import *


# WARNING: LINES SAVED INTO MEMORY
# USE CAUTION WHEN RUNNING ON LARGE FILES


# set directory variables for easy file i/o while testing
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(PROJ_ROOT, 'local_data')
LOCAL = os.path.join(PROJ_ROOT, 'scripts\\VCF_Comparisons')
#LOCAL = os.path.join(PROJ_ROOT, 'local_data')


def sort(vcf_rdr, order, stk):
    #return_loc = vcf_rdr.cur_loc

    sof = stack.enter_context(open(os.path.join(DATA_DIR, vcf_rdr.path + ".sorted.vcf"), "w")) 

    # copy header to new file
    while vcf_rdr.raw_line.startswith("#") and not vcf_rdr.end_state:
        sof.write(vcf_rdr.raw_line)
        vcf_rdr.read() 
    # save the file position of the end of the header/metadata
    vcf_rdr.header_end = vcf_rdr.file_obj.tell() 
    

    for chrom in order:
        chrom_data_lines = []
        while not vcf_rdr.end_state:
            if vcf_rdr.chrom == chrom:
                if chrom_data_lines:
                    i = 0
                    while vcf_rdr.pos > chrom_data_lines[i][0]:
                        i += 1
                        if i == len(chrom_data_lines):
                            break
              
                    chrom_data_lines.insert(i, [vcf_rdr.pos, vcf_rdr.raw_line])
                    
                else:
                    chrom_data_lines.append([vcf_rdr.pos, vcf_rdr.raw_line])
         
            vcf_rdr.read()
        for lines in chrom_data_lines:
            sof.write(lines[1])
        vcf_rdr._setFilePosition(vcf_rdr.header_end)
        vcf_rdr.end_state = False
        vcf_rdr.read()





def getBEDOrder(bed):
    chrom_order = []
    while not bed.end_state:
        if bed.chrom not in chrom_order:
            chrom_order.append(bed.chrom)
        bed.read()

    return chrom_order





bed_file = "test-isolated-vc-catalog.atarva.bed.gz"

vcf_list = [ # the file name, and the offset amount (eg. [-1, 0] for 1 based inclusive), and bool for whether is is position only
    #os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.atarva.vcf"),
    os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.strdust.vcf"),
    #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.strkit.vcf"), 
    os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.longTR.vcf.gz"),
    #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.straglr.vcf"),
    #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.vamos.vcf"),
    os.path.join(DATA_DIR,"medaka_to_ref.TR.vcf"),
    #"test_cases.vcf"
    ]

vcf_rdrs = []
with ExitStack() as stack: 
    # create bed reader and enter the file in stack
    bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))

    order = getBEDOrder(bed)
    print(order)

    # create list of vcf reader objects and their new sorted, output file
    for i, vcf in enumerate(vcf_list):
        vcf_rdrs.append(setupVCFReader(vcf=vcf, 
                                        skip_head=False,
                                        stk=stack))
        
        sort(vcf_rdrs[i], order, stack)