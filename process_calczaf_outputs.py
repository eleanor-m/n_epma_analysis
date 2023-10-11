import sys
sys.path.insert(1,'..')
from src import calczaf, helper_funs
from pathlib import Path

folderpath = Path('./calczaf_files/')

helper_funs.check_calczaf_folder_exists(folderpath)
valence_file = sorted(folderpath.glob('valence*'))[0]

calczaf.process_calczaf_outputs(folderpath, valence_file)

# For detection limits

calczaf.process_calczaf_outputs(folderpath / 'detlim/', valence_file, detlim=True)