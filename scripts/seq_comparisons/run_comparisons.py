import os

# set directoy cariables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath('')))
DATA_DIR = os.path.join(PROJ_ROOT, 'catalogs')
# TESTING PURPOSES ONLY
LOCAL_DATA = os.path.join(PROJ_ROOT, 'local_data') # REMOVE FOR FINAL LAUNCH

def main():
    pass