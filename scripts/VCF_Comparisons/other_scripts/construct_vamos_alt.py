import os
import sys
from contextlib import ExitStack
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(root_dir))

from helpers.constants import *
from helpers.readers import VCFReader


def constructVCF(vcf):
  
    with ExitStack() as stack: 
        # create list of vcf reader objects and their new sorted, output file
        vcf_rdr = VCFReader(vcf)


    # vcf reader will already error out if the file is not valid, so no need for extensive checks
    if vcf_rdr.path.endswith(".gz"):
        sof_name = os.path.join(DATA_DIR, vcf_rdr.path[:-7] + ".ALT_FIX.vcf")
    else:
        sof_name = os.path.join(DATA_DIR, vcf_rdr.path[:-4] + ".ALT_FIX.vcf")

    sof = stack.enter_context(open(sof_name, "w")) 

    # copy header to new file
    vcf_rdr.read() 
    while vcf_rdr.raw_line.startswith("#") and not vcf_rdr.end_state:
        sof.write(vcf_rdr.raw_line)
        vcf_rdr.read() 



    while vcf_rdr.cur_line:

        vcf_rdr.read()



