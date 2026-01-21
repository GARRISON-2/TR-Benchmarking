from enum import Enum, auto
import os

# set directory variables for easy file i/o while testing
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(PROJ_ROOT, 'local_data')
OUTPUT_DIR = os.path.join(PROJ_ROOT, 'local_data')
LOCAL = os.path.join(PROJ_ROOT, 'scripts\\VCF_Comparisons')
#LOCAL = os.path.join(PROJ_ROOT, 'local_data')


class SPECIAL_CASE(Enum):
    STRAGLR = (auto(), True)
    STRAGLR_SVLEN = (auto(), True)
    VAMOS = (auto(), True)

    def __init__(self, value ,pos_only):
        self._value_ = value
        self.pos_only = pos_only

class COMP_ORDER(Enum):
    VERTICAL = auto()
    CROSS = auto()

class COMP_METHOD(Enum):
    LEVENSHTEIN = auto()
    LENGTH = auto()
