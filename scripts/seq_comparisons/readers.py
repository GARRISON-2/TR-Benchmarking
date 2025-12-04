

class VCFReader:
    def __init__(self, file = None):
        self.file_obj = file
        self.name = self.setName()
        self.pause = False
        self.end_state = False   
        self.parsed_header = False
        self.line_string = None
        self.cur_line = None
        self.nextLine()        


    def setHeaderInfo(self, header_ls):
        line_list = header_ls.strip().split("\t")

        col_dict = {}
        for i, col in enumerate(line_list):
            col_dict[f"{col}"] = i
        
        return col_dict


    def skipMetaData(self):
        # loop while the current line contains meta data
        while self.line_string.startswith("##"):
            self.nextLine()

        if self.line_string.startswith("#"):
            self.fields = self.setHeaderInfo(self.line_string)      
        else:
            self.fields = None

        self.parsed_header = True
        self.nextLine()


    def formatLine(self, ls, head_bool):
        # split line string into list
        line_list = ls.strip().split("\t")

        # if the file not reading the header/metadata
        if head_bool:
            # calculate the end position and add it to the end of the line list
            pos = int(line_list[1])
            ref_len = len(line_list[3])
            end_pos = pos + ref_len - 1
            line_list.append(end_pos)

            # set specific data to their own parameters for better accessibility
            self.pos_info = {          # eg:
                "chrom": line_list[0],   # CHROM1 
                "start": pos,            # 10002
                "end": end_pos}          # 10222
            self.ref = line_list[self.fields["REF"]] # the refence sequence as a string            
            self.alt = line_list[self.fields["ALT"]].split(",") # returns a list of all alt alleles

        return line_list 
 

    def nextLine(self):
        # try to moce to the next file as long is it is not already
        # at the end and it paused
        if (not self.end_state) & (not self.pause):
            try:
                self.line_string = next(self.file_obj)
                self.cur_line = self.formatLine(self.line_string, self.parsed_header)
            # catch exception that signals end of file
            except StopIteration:
                self.end_state = True
                self.cur_line = None


    def setName(self):
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
        self.name = self.setName()
        self.end_state = False
        self.line_string = None   
        self.cur_line = None
        self.nextLine()


    def nextLine(self):
        if (not self.end_state):
            try:
                self.line_string = next(self.file_obj)
                self.cur_line = self.formatLine(self.line_string)
            # catch exception that signals end of file
            except StopIteration:
                self.end_state = True
                self.cur_line = None


    def formatLine(self, ls):
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


    def setName(self):
        path_str = self.file_obj.name
        # enumerate through the reversed file path to grab
        # the index of where the file name starts
        for i, letter in enumerate(reversed(path_str)):
            if letter == "\\":
                name_start = len(path_str) - i
                break
                
        # slice for the file name minus the path and the format
        return path_str[name_start:-4]
