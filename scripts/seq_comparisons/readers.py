

class VCFReader:
    def __init__(self, file, offset = 0):
        self.file_obj = file 
        self.name = self._setName()
        self.offset = offset
        self.pause = False
        self.end_state = False # bool for whether or not the end of the file has been reacher
        self.cur_line = None # formatted list version of the current line read from the file
        self.read() # move to the first line in the file
        self.fields = self._parseForFields() # set the field information to keep track of the columns in the file
       

    def _formatLine(self, ls):
        # split line string into list of strings
        line_list = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not line_list[0].startswith("#"):
            # calculate the end position and add it to the end of the line list
            pos = int(line_list[1])
            ref_len = len(line_list[3])
            end_pos = pos + ref_len - 1
            line_list.append(end_pos)

            # set specific data to their own parameters for better accessibility
            self.pos_info = {          # eg:
                "chrom": line_list[0],   # CHROM1 
                "start": pos - self.offset,            # 10002
                "end": end_pos - self.offset}          # 10222
            self.ref = line_list[self.fields["REF"]] # the refence sequence as a string            
            self.alt = line_list[self.fields["ALT"]].split(",") # returns a list of all alt alleles

        return line_list 
 

    def read(self):
        # try to move to the next file as long is it is not already
        # at the end and it paused
        if (not self.end_state) & (not self.pause):
            line_string = self.file_obj.readline()

            # if the line is not empty (ie. the end of the file has been reached) 
            if line_string:
                # then format and set the current line
                self.cur_line = self._formatLine(line_string)
            else: 
                self.end_state = True
                self.cur_line = None

 
    def skipMetaData(self): 
        # the header end position is already saved, 
        # so skip to that position in the file then advance to the 1st line of data
        self._setFilePosition(self.header_end)
        

    def _setFilePosition(self, file_pos):
        self.file_obj.seek(file_pos)
        self.read()


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
        
             
    def _setFieldInfo(self, head_list):
        # loop over the header and make dictionary of the column names and their indices
        col_dict = {}
        for i, col in enumerate(head_list):
            col_dict[f"{col}"] = i
        
        return col_dict


    def _setName(self):
        path_str = self.file_obj.name
        # enumerate through the reversed file path to grab
        # the index of where the file name starts
        for i, letter in enumerate(reversed(path_str)):
            if letter == "\\":
                name_start = len(path_str) - i
                break
                
        # slice for the file name minus the path and the format
        return path_str[name_start:-4]



class BEDReader:
    def __init__(self, file):
        self.file_obj = file
        self.name = self._setName()
        self.end_state = False
        self.line_string = None   
        self.cur_line = None
        self.read()


    def read(self):
        if (not self.end_state):
            try:
                self.line_string = next(self.file_obj)
                self.cur_line = self._formatLine(self.line_string)
            # catch exception that signals end of file
            except StopIteration:
                self.end_state = True
                self.cur_line = None


    def _formatLine(self, ls):
        # split line string into list
        line_list = ls.strip().split("\t")

        # calculate the end position and add it to the end of the row list
        pos = int(line_list[1])
        ref_len = len(line_list[3])
        end_pos = pos + ref_len - 1
        line_list.append(end_pos)

        # set position information parameter for easier access
        self.pos_info = {          # eg:
                "chrom": line_list[0],   # CHROM1 
                "start": pos,            # 10002
                "end": end_pos}          # 10222
        self.ref = line_list[3]

        return line_list


    def _setName(self):
        path_str = self.file_obj.name
        # enumerate through the reversed file path to grab
        # the index of where the file name starts
        for i, letter in enumerate(reversed(path_str)):
            if letter == "\\":
                name_start = len(path_str) - i
                break
                
        # slice for the file name minus the path and the format
        return path_str[name_start:-4]
