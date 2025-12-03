import os
import itertools 
import Levenshtein as lv

# set directoy cariables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH


class FileInfo:
    def __init__(self, file):
        self.file_obj = file
        self.name = self.setName()
        self.pause = False
        self.end_state = False   
        self.nextLine() # sets self.cur_line
        self.line_info = None


    def skipMetaData(self):
        # loop while the current line contains meta data
        while self.cur_line.startswith("##"):
            self.nextLine()
        if self.cur_line.startswith("#"):
            self.header = self.cur_line
            self.nextLine()
        else:
            self.header = None


    def nextLine(self):
        if (not self.end_state) & (not self.pause):
            try:
                self.cur_line =  next(self.file_obj)
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


class LineInfo:
    def __init__(self, ls):
        self.line_list, self.pos_info = self.formatLine(ls)

    
    def formatLine(self, ls):
        # split line string into list
        line_list = ls.strip().split("\t")

        # calculate the end position and add it to the end of the row list
        pos = int(line_list[1])
        ref_len = len(line_list[3])
        end_pos = pos + ref_len - 1

        line_list.append(end_pos)

        return line_list, (line_list[0], pos, end_pos)
    
    def setVCF(self):
        # build genotype tuple
        print("\n", self.line_list[-2][0], self.line_list[-2][2])

        # check indices to make sure they are ints
        gt_idx = (checkIdx(self.line_list[-2][0]), checkIdx(self.line_list[-2][2]))

        alleles = [self.line_list[3], self.line_list[4], self.line_list[5]]

        # build genotype based on indices that were just grabbed
        gt = gt = [alleles[gt_idx[0]], alleles[gt_idx[1]]]

        print(gt)
    

def compareString(line1, line2):
    # compare the bed info and the refs (for now) 
    #print(f"{seq1.cur_line[3]},    {seq2.cur_line[3]}")
    return_dict = {}

    dist = lv.distance(line1.line_list[3], line2.line_list[3])

    # get the positional difference and calculate the lv distance based only on positions in the actual reference
    start_diff = line2.pos_info[1] - line1.pos_info[1]
    end_diff = len(line1.line_list[3]) - line2.pos_info[2] - line1.pos_info[2]

    dist2 = lv.distance(line1.line_list[3][start_diff: end_diff], line2.line_list[3][start_diff: end_diff])
    
    return_dict = {
        "dist": dist,
        "dist2": dist2
    }

    return return_dict 


def checkIdx(idx):
    try:
        return int(idx)
    except ValueError:
        return 0


def mainloop(bed_file, file1, file2):

    heads_parsed = False

    # print("Enter catalog path:")
    # print("Enter 1st VCF path:")
    # print("Enter 2nd VCF path")

    # open all input and output files
    with open(os.path.join(DATA_DIR, bed_file) , 'r') as ref, \
    open(os.path.join(DATA_DIR, file1) , 'r') as f1, \
    open(os.path.join(DATA_DIR, file2) , 'r') as f2, \
    open(os.path.join(LOCAL_DATA, "reference-comp.tsv"), "w") as of, \
    open(os.path.join(LOCAL_DATA, "comp-metadata.tsv"), "w") as mof:
        # CURRENTLY REPRESENTS FILE STACK
        file_stack = [bed_file, f1, f2] 

        # set column headers for output tsv file
        of.write("FILE\tBEDPOS\tVCFPOS\tDIST\tDIST2\n")

        # make list of open file objects for easier management
        file_list = [FileInfo(ref)]
        for i, file in enumerate(file_stack[1:]):
            file_list.append(FileInfo(file))


        # move past header infromation
        for i, fi_info in enumerate(file_list[1:]):
            fi_info.skipMetaData()


        # loop until all files have reached their final line
        #while not any([obj.end_state for obj in file_list]):
        for j in range(10):
            out_str = ""
    
            for i, fi_info in enumerate(file_list):      
                fi_info.line_info = LineInfo(fi_info.cur_line)

                # position alignment check
                # if i == 0:
                #     out_str += (f"BED Ref ({fi_info.pos_info}): {cur_line[3]}\n")

                # if not the bed file
                if i != 0:                   
                    #print(f"File[{i}]{fi_info.pos_info[1]}:{file_list[0].pos_info[1]} - {fi_info.pos_info[1] > file_list[0].pos_info[1]}") # DEBUGGING - REMOVE FOR LAUNCH
                    
                    # run vcf specific operations

                    # if current file position is ahead of bed position range
                    if (fi_info.line_info.pos_info[1] > file_list[0].line_info.pos_info[1]):                    
                        fi_info.line_info.pause = True
                    else:
                        fi_info.pause = False

                        # run comparisons
                        diff_dict = compareString(file_list[0].line_info, fi_info.line_info)

                        # out_str += (f"File[{i}] Ref diff score: {difference} | ({fi_info.pos_info}): {cur_line[3]}\n")
                        out_str += (f"{i}\t{file_list[0].line_info.pos_info}\t{fi_info.line_info.pos_info}\t{diff_dict["dist"]}\t{diff_dict["dist2"]}\n")
                
                fi_info.nextLine()

            of.write(out_str)






mainloop("test-isolated-vc-catalog.strkit.bed", "HG001.strdust.vcf", "HG001.strkit.vcf")
