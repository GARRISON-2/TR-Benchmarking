import sys
import gzip
import io

'''
'''
class Reader:
    def __init__(self, file_path, buffer_size = io.DEFAULT_BUFFER_SIZE):
        self.file_obj = None
        self.buffer = buffer_size
        self.path = file_path
        self.end_state = False # bool for whether or not the end of the file has been reached
        self.cur_line = None


    def close_file(self):
        if not self.file_obj.closed:
            self.file_obj.close()


    def open_file(self):
        try:
            if self.path.endswith(".gz"):
                self.file_obj = gzip.open(self.path, "rt", encoding="utf-8")
            else:
                self.file_obj = open(self.path, "r", encoding="utf-8", buffering=self.buffer)
 
            self.read() # move to the first line in the file
            return self
        except (IOError, OSError) as e: 
            raise FileIOError(f"File Opening Error: {e}")
        

    '''
    Returns True if the Line was read, and False otherwise
    '''
    def read(self, format = True):
        # try to move to the next file as long is it is not already at the end or paused 
        if not self.end_state:
            try:
                line_string = self.file_obj.readline() # readline allows the ability to save positions in the file, as opposed to read()
            except gzip.BadGzipFile:
                raise FileReadError(f"Failed to read from {self.path}\nInvalid .gz")
            except UnicodeError:
                raise FileReadError(f"Failed to read from {self.path}\nContains Invalid UTF-8 Characters")

            # if the line is not empty (ie. the end of the file has not been reached) 
            if line_string:
                self.prev_line = self.cur_line
                if format:
                    # then format and set the current line
                    self.formatLine(line_string)
                else:
                    self.cur_line = line_string
            else: 
                self.end_state = True
                self.prev_line = self.cur_line
                self.cur_line = None
                self.close_file()

            return True
        
        return False
        #raise FileReadError(f"Failed to read from {self.path}\nFile has reached end state.")


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
        self.cur_line = None
        self.prev_line = None
        self.fields = None
        


    '''
    
    '''
    def buildGt(self, sample_col=9):
        try:
            # grab genotype indices from current line
            sample_str = self.cur_line[sample_col]
            idx_str = sample_str.split(':')[0] # splits sample column by ':' and grabs genotype index from start
            choices = [self.ref, *self.alt, None] # list of possible choices to use for constructing genotype
            gt = []

            # loop through the index string and construct genotype
            for idx in idx_str[::2]: # grab characters every 2 indices to skip gt seperator
                # check indices to make sure they are ints (eg. accounting for 1|., etc.) before adding allele to gt
                gt_idx = self._checkIdx(idx) 
                gt.append(choices[gt_idx]) 

            self.genotype = gt

            # add to list of genotypes for comparison between VCFs
            return gt  
        
        except (IndexError, ValueError):
            raise VCFFormatError(f"Failed to construct Genotype using '{idx_str}' from sample: {sample_str}")

        except Exception as e:
            raise VCFFormatError(f"Failed to construct Genotype due to unknown error: {e}\nFrom sample: {sample_str}")


    '''
    
    '''
    def open_file(self):
        super().open_file()

        # set the field information to keep track of vcf columns
        self.fields = self._parseForFields() 

        return self


    '''
    
    '''
    def read(self):
        read_state = False

        if not self.pause: 
            read_state = super().read()
        
        return read_state
 

    '''
    
    '''
    def skipMetaData(self): 
        # the header end position is already saved, 
        # so skip to that position in the file then advance to the 1st line of data
        self._setFilePosition(self.header_end)


    '''
    
    '''
    def _checkIdx(self, idx):
        if idx == '.':
            return -1 # will set the allele to None
        else:
            return int(idx)


    '''
    
    '''
    def formatLine(self, ls):
        # split line string into list of strings
        self.cur_line = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not self.cur_line[0].startswith("#") and self.fields:  
            try:  
                # calculate the end position and add it to the end of the line list
                
                pos = int(self.cur_line[self.fields["POS"]])                   
                ref_len = len(self.cur_line[self.fields["REF"]])
                end_pos = pos + ref_len - 1

                # set specific data to their own parameters for better accessibility
                self.pos_info = {                       
                    "chrom": self.cur_line[self.fields["CHROM"]],               
                    "start": pos + self.start_off,      
                    "end": end_pos + self.end_off}      
                self.ref = self.cur_line[self.fields["REF"]] # the refence sequence as a string            
                self.alt = self.cur_line[self.fields["ALT"]].split(",") # returns a list of all alt alleles
            
            except ValueError:
                raise VCFFormatError(f"Failed to set position '{self.cur_line[self.fields['POS']]}' from line: {self.cur_line}\n")
            
            except IndexError:
                raise VCFFormatError(f"Missing parameter data from line: {self.cur_line}")
            
            except Exception as e:
                raise VCFFormatError(f"Unexpected error setting parameters from line: {self.cur_line}\n{e}")

  
    '''
    
    '''
    def _parseForFields(self):
        field_list = []

        # loop while the current line contains meta data
        while self.cur_line[0].startswith("#"):
            if self.cur_line[0].upper().startswith("#CHR"):
                field_list = self.cur_line
                field_list[0] = field_list[0][1:] # remove '#' from start of first field
                break
            else:
                self.read()  

        # save the file position of the end of the header/metadata
        self.header_end = self.file_obj.tell()

        # return to the top of the file
        self._setFilePosition(0)

        # if the data field names were not found
        if not field_list:           
            field_list = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
            print(f"WARNING: Header not found. Using {field_list} for data field names and indices.")


        return self._setFieldInfo(field_list) 


    '''
    
    '''
    def _setFilePosition(self, file_pos):
        self.file_obj.seek(file_pos)

             
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

   
    def formatLine(self, ls):
        # split line string into list
        self.cur_line = ls.strip().split("\t")

        try:
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

        except ValueError:
            raise BEDFormatError(f"ERROR: From file: {self.path}\nFailed to set position {self.cur_line[self.fields['POS']]} from line: {self.cur_line}\n")
        except IndexError:
            raise BEDFormatError(f"ERROR: From file: {self.path}\nMissing parameter data from line: {self.cur_line}")
        except Exception as e:
            raise BEDFormatError(f"ERROR:From file: {self.path}\nUnexpected error setting parameters from line:{self.cur_line}: {e}")


        return self.cur_line



class FileIOError(Exception):
    def __init__(self, message):
        super().__init__(message)

class FileReadError(Exception):
    def __init__(self, message):
        super().__init__(message)

class BEDFormatError(Exception):
    def __init__(self, message):
        super().__init__(message)

class VCFFormatError(Exception):
    def __init__(self, message):
        super().__init__(message)