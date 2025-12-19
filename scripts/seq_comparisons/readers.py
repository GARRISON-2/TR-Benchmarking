import sys
import gzip


class Reader:
    def __init__(self, file_path):
        self.file_obj = None
        self.path = file_path
        self.end_state = False # bool for whether or not the end of the file has been reached
        self.cur_line = None # formatted list version of the last line string read from the file


    def close_file(self):
        if not self.file_obj.closed:
            self.file_obj.close()


    def open_file(self):
        if self.path.endswith(".gz"):
            self.file_obj = gzip.open(self.path, "rt", encoding="utf-8")
        else:
            self.file_obj = open(self.path, "r", encoding="utf-8")

        try:
            self.read() # move to the first line in the file
            return self
        except IOError as e: 
            sys.exit(f"Failed to open file: {e}") 

 
    '''
    Returns True if the Line was read, and False otherwise
    '''
    def read(self, pause = False):
        # try to move to the next file as long is it is not already
        # at the end and it paused
        if (not self.end_state) and (not pause):
            line_string = self.file_obj.readline() # readline allows the ability to save positions in the file, as opposed to read()

            # if the line is not empty (ie. the end of the file has been reached) 
            if line_string:
                self.prev_line = self.cur_line
                # then format and set the current line
                self._formatLine(line_string)
            else: 
                self.end_state = True
                self.prev_line = self.cur_line
                self.cur_line = None
                self.close_file()

            
            return True
        
        return False


    def __enter__(self):
        return self.open_file()


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close_file()




class VCFReader(Reader):
    def __init__(self, file_path, start_offset = -1, end_offset = 0):
        super().__init__(file_path)
        self.start_off = start_offset
        self.end_off = end_offset
        self.pause = False
        self.prev_line = None
        

    def buildGt(self):
    
        # check indices to make sure they are ints (eg. accounting for 1|., etc.)
        gt_idx = (self._checkIdx(self.cur_line[-1][0]), self._checkIdx(self.cur_line[-1][2]))
        choices = [self.ref, *self.alt, None]

        # build genotype based on gt_idx indices
        gt = [choices[gt_idx[0]], choices[gt_idx[1]]]
  
        self.genotype = gt

        # add to list of genotypes for comparison between VCFs
        return gt  


    def open_file(self):
        super().open_file()

        self.fields = self._parseForFields() # set the field information to keep track of the columns in the file

        return self


    def read(self):
        return super().read(self.pause)
 

    def skipMetaData(self): 
        # the header end position is already saved, 
        # so skip to that position in the file then advance to the 1st line of data
        self._setFilePosition(self.header_end)



    def _checkIdx(self, idx):
        try:
            return int(idx)
        except ValueError:
            return -1


    def _formatLine(self, ls):
        # split line string into list of strings
        self.cur_line = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not self.cur_line[0].startswith("#"):
            # calculate the end position and add it to the end of the line list
            pos = int(self.cur_line[self.fields["POS"]])
            ref_len = len(self.cur_line[self.fields["REF"]])
            end_pos = pos + ref_len - 1

            # set specific data to their own parameters for better accessibility
            self.pos_info = {                       
                "chrom": self.cur_line[self.fields["#CHROM"]],               
                "start": pos + self.start_off,      
                "end": end_pos + self.end_off}      
            self.ref = self.cur_line[self.fields["REF"]] # the refence sequence as a string            
            self.alt = self.cur_line[self.fields["ALT"]].split(",") # returns a list of all alt alleles


    def _parseForFields(self):
        
        # loop while the current line contains meta data
        while self.cur_line[0].startswith("#"):
            if self.cur_line[0].startswith("#CHROM"):
                # save the file position of the end of the header/metadata
                self.header_end = self.file_obj.tell()
                field_info = self._setFieldInfo(self.cur_line)  
                self._setFilePosition(0) # return to the top of the file
                return field_info

            self.read()  

        # if the header information is not found, raise an error 
        raise ValueError("VCF Header Information Not Found")
        

    def _setFilePosition(self, file_pos):
        self.file_obj.seek(file_pos)
        self.read()

             
    def _setFieldInfo(self, head_list):
        # loop over the header and make dictionary of the column names and their indices
        col_dict = {}
        for i, col in enumerate(head_list):
            col_dict[f"{col}"] = i
        
        return col_dict



class BEDReader(Reader):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.prev_line = None

   
    def _formatLine(self, ls):
        # split line string into list
        self.cur_line = ls.strip().split("\t")

        # calculate the end position and add it to the end of the row list
        pos = int(self.cur_line[1])
        end_pos = int(self.cur_line[2])
        self.cur_line.append(end_pos)

        # set position information parameter for easier access
        self.pos_info = {          # eg:
                "chrom": self.cur_line[0],   # CHROM1 
                "start": pos,            # 10002
                "end": end_pos}          # 10222
        self.ref = self.cur_line[3]

        return self.cur_line

