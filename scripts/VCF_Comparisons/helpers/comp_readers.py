from helpers.readers import *

class SC_VCFReader(VCFReader):
    def __init__(self, file_path, start_offset = 0, end_offset = 0, pos_only = False, pause = False, lines_skipped = 0):
        super().__init__(file_path, start_offset, end_offset)
        self.pause = pause
        self.pos_only = pos_only
        self.is_data = False
        self.skip_num = lines_skipped


    def safeRead(self):
        line = None

        try: 
            if not self.pause:  
                format_method = self.formatLine if not self.pos_only else self.posOnlyFormat
                line = super().read(format_method)  
                 
        except (FileReadError, VCFFormatError, BEDFormatError) as e:
            sys.exit(f"\nERROR\n{e}")   

        return line
    

    def posOnlyFormat(self, ls):
        # split line string into list of strings
        line_list = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not line_list[0].startswith("#"):  
            try:       

                pos = int(line_list[1]) 

                info_col = line_list[7].split(';') # grab the INFO column
                end_pos = int(info_col[0].strip("END="))

                # set specific data to their own parameters for better accessibility
                self.chrom = line_list[0]             
                self.pos = pos + self.start_off     
                self.end_pos = end_pos + self.end_off   
                self.ref = None       
                self.alt = [None]
                self.info = info_col
                    
            except ValueError:
                raise VCFFormatError(f"Failed to set position '{line_list[4]}' from line: {line_list}\n")
            
            except IndexError:
                raise VCFFormatError(f"Missing parameter data from line: {line_list}")
            
            except Exception as e:
                raise VCFFormatError(f"Unexpected error setting parameters from line: {line_list}\n{e}")

        return line_list


    def checkOrder(self, order_method="ASCII"):
        # check VCF Chromosome ordering
        if not self.prev_line[0].startswith("#"):
            prev_chrom = self.prev_line[0]

            # ensure vcf is in order
            if self.chrom < prev_chrom:
                sys.exit(f"\nERROR\n{self.path} using unknown order. Ending program.")


    def syncToBed(self, bed, order_method="ASCII"):
        # check VCF chromosome ordering
        self.checkOrder(order_method)

        # Run Alingment Checks
        chrom_match = (self.chrom == bed.chrom)

        # if vcf position is behind the bed, or the vcf chrom is behind, loop until the vcf catches up         
        while (((self.end_pos < bed.pos) and chrom_match) or (self.chrom < bed.chrom)) \
            and not self.end_state:

            # move the file line forward until it is no longer behind, or the end of the file is reached
            self.safeRead()
            self.checkOrder(order_method)

            self.skip_num += 1


        # if vcf position is ahead of bed position range, or if the vcf chrom is ahead, then pause operations
        if ((self.pos > bed.end_pos) and chrom_match) or \
            (self.chrom > bed.chrom):

            self.pause = True # pause vcf from being able to move to the next line or run comparisons

        # if the vcf is not ahead or the chromosomes dont match, continue moving forward
        else: 
            self.pause = False
