import sys
from dataclasses import dataclass
from helpers.readers import *
from helpers.constants import *


@dataclass
class alleleData:
    """
    Docstring for alleleData
    """
    allele_str: str | None = None
    is_ref: bool | None = None
    length: int = 0


class COMP_VCFReader(VCFReader):
    
    def __init__(self, file_path, start_offset = 0, end_offset = 0, pos_only = False, sc = None, pause = False):
        """
        Docstring for __init__
        
        :param self: Description
        :param file_path: Description
        :param start_offset: Description
        :param end_offset: Description
        :param pos_only: Description
        :param sc: Description
        :param pause: Description
        """
        super().__init__(file_path)
        self.start_off = start_offset
        self.end_off = end_offset
        self.pause = pause
        self.pos_only = pos_only
        self.skip_num = 0
        self.end_state = False
        self.spec_case = sc
        self.STRAGLR_svlen = 0
        self.gt_data_pack = None


    def buildGtData(self, sample_col=9, ref=None, alt = None):
        """
        Docstring for buildGtData
        
        :param self: Description
        :param sample_col: Description
        :param ref: Description
        :param alt: Description
        """
        try:
            # use class parameters unless otherwise specified
            ref = self.ref if not ref else ref
            alt = self.alt if not alt else alt
            
            # split line into list for easy data grabbing
            ls = self._raw_line.strip().split("\t")

            # grab genotype indices from current line
            sample_str = ls[sample_col]
            idx_str = sample_str.split(':')[0] # splits sample column and grabs genotype index from first index
            idx_list = list(idx_str[::2]) # grab characters every 2 indices to skip gt seperator(ie. skipping over '|' and '/')

            # convert 0 to 0|0 for straglr
            if self.spec_case == SPECIAL_CASE.STRAGLR:
                if len(idx_list) == 1:
                    idx_list.append(idx_list[0]) 
            
            choices = [self.ref, *self.alt, None] # list of possible choices to use for constructing genotype
            gt_data = []

            # loop through the index string and construct genotype
            for idx in idx_list: 
                # create new object to hold allele string and other data
                ad = alleleData()

                # check indices to make sure they are ints (eg. accounting for 1|., etc.) before adding allele to gt
                gt_idx = self._checkIdx(idx) 
                ad.allele_str = choices[gt_idx] 

                # keep track of whether the allele is the same as ref or not
                if gt_idx == 0:
                    ad.is_ref = True
                elif gt_idx > 0:
                    ad.is_ref = False

                # custom handling for straglr since it does not contain sequences
                if self.spec_case == SPECIAL_CASE.STRAGLR:
                    ad.length = self._handleStraglrLen(ad.is_ref)

                elif ad.allele_str is not None:
                    ad.length = len(ad.allele_str)

                # else leave the length as 0

                # add the allele data object to the genotype 
                gt_data.append(ad)
                
            self.gt_data_pack = gt_data

            return gt_data
        
        except (IndexError, ValueError):
            raise VCFFormatError(f"Failed to construct Genotype using '{idx_str}' from sample: {sample_str}")

        except Exception as e:
            raise VCFFormatError(f"Failed to construct Genotype due to unknown error: {e}\nFrom sample: {sample_str}")


    def checkOrder(self, order_method="ASCII"):
        """
        Docstring for checkOrder
        
        :param self: Description
        :param order_method: Description
        """
        # check VCF Chromosome ordering
        if self.prev_line and not self.prev_line[0].startswith("#"):
            prev_chrom = self.prev_line[0]

            # ensure vcf is in order
            if self.chrom < prev_chrom and order_method == "ASCII":
                raise VCFFormatError(f"\n{self.path} using unknown order.")
    

    def constructAllele(self, idxs: list, mots: list):
        """
        Docstring for constructAllele
        
        :param self: Description
        :param idxs: Description
        :type idxs: list
        :param mots: Description
        :type mots: list
        """
        allele = ""
        for idx in idxs:
            allele += mots[int(idx)]
        
        return allele


    def constructAlt(self, info: str):
        """
        Docstring for constructAlt
        
        :param self: Description
        :param info: Description
        :type info: str
        """
        alt = []

        for i, str in enumerate(info): # loop over infor and grab data for constructing allele sequence        
            if "RU" in str:
                mot_list = info[i].removeprefix("RU=").split(",")

            if "ALTANNO_H1" in str:
                mot_idxs = info[i].removeprefix("ALTANNO_H1=").split("-")
                alt.append(self.constructAllele(mot_idxs, mot_list))

            if "ALTANNO_H2" in str:
                mot_idxs = info[i].removeprefix("ALTANNO_H2=").split("-")
                alt.append(self.constructAllele(mot_idxs, mot_list))

        return alt


    def posOnlyFormat(self, ls: str):
        """
        Docstring for posOnlyFormat
        
        :param self: Description
        :param ls: Description
        :type ls: str
        """
        # split line string into list of strings
        line_list = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not line_list[0].startswith("#"):  
            try:       

                pos = int(line_list[1]) 

                info_col = line_list[7].split(';') # grab the INFO column
                end_pos = int(info_col[0].removeprefix("END="))

                # set specific data to their own parameters for better accessibility
                self.chrom = line_list[0]             
                self.pos = pos + self.start_off     
                self.end_pos = end_pos + self.end_off   
                self.ref = None
                if self.spec_case == SPECIAL_CASE.VAMOS:
                    self.alt = self.constructAlt(info_col)
                else:        
                    self.alt = [None] # functions in super class expect alt to be a list                 
                self.info = info_col
                    

            except IndexError:
                raise VCFFormatError(f"Missing parameter data from line: {line_list}")
            
            except Exception as e:
                raise VCFFormatError(f"Unexpected error setting parameters from line: {line_list}\n{e}")

        return line_list


    def formatLine(self, ls: str):
        """
        Docstring for formatLine
        
        :param self: Description
        :param ls: Description
        :type ls: str
        """
        fl = super().formatLine(ls)

        # if the vals have been set already
        if self.pos and self.end_pos:
            self.pos += self.start_off
            self.end_pos += self.end_off

        return fl


    def safeRead(self):
        """
        Docstring for safeRead
        
        :param self: Description
        """
        line = None

        try: 
            if not self.pause:  
                format_method = self.posOnlyFormat if self.pos_only else self.formatLine
                line = super().read(format_method)  

                if not line: 
                    self.end_state = True
                 
        except (FileReadError, VCFFormatError, BEDFormatError) as e:
            sys.exit(f"\nERROR\n{e}")   

        return line


    def syncToBed(self, bed: BEDReader, order_method="ASCII"):
        """
        Docstring for syncToBed
        
        :param self: Description
        :param bed: Description
        :type bed: BEDReader
        :param order_method: Description
        """
        # check VCF chromosome ordering
        try:
            self.checkOrder(order_method)
        except Exception as e:
            sys.exit(f"\nERROR\n{e}\nEnding program.")

        # Run Alingment Checks
        chrom_match = (self.chrom == bed.chrom)

        # if vcf position is behind the bed, or the vcf chrom is behind, loop until the vcf catches up         
        while (((self.end_pos < bed.pos) and chrom_match) or (self.chrom < bed.chrom)) \
            and not self.end_state:

            # move the file line forward until it is no longer behind, or the end of the file is reached
            self.safeRead()
            self.checkOrder(order_method)

            self.skip_num += 1

        # build the most current line's genotype data
        if not self.end_state:
            self.buildGtData()

        # if vcf position is ahead of bed position range, or if the vcf chrom is ahead, then pause operations
        if ((self.pos > bed.end_pos) and chrom_match) or \
            (self.chrom > bed.chrom):

            self.pause = True # pause vcf from being able to move to the next line or run comparisons

        # if the vcf is not ahead or the chromosomes dont match, continue moving forward
        else: 
            self.pause = False


    def _checkIdx(self, idx):
        """
        Docstring for _checkIdx
        
        :param self: Description
        :param idx: Description
        """
        if idx == '.':
            return -1 # will set the allele to None
        else:
            return int(idx)
        

    def _handleStraglrLen(self, is_ref):

        if is_ref:
            # for i, str in enumerate(self.info): # loop over infor and grab data for constructing allele sequence        
            #     if "SVLEN" in str:
            #         svlen = int(self.info[i].removeprefix("SVLEN=").split(","))

            return self.end_pos - self.pos 
        else:
            for i, str in enumerate(self.info): # loop over info and grab number of bases   
                if "RB" in str:
                    return int(self.info[i].removeprefix("RB=").replace(',', ""))

