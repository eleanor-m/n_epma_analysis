# %% Import dependencies and set plot defaults
import pandas as pd
import re
from pathlib import Path
import os
import shutil
from src import correct_quant, helper_funs

pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)

def import_file(filepath):

    """ Imports a file from a subfolder and returns a list of lines
    in the file as strings

    filepath can be a string or a Path object
    """

    if type(filepath) == str:
        filepath = Path(filepath)

    with open(filepath, 'r', errors='ignore') as fp:
        imported_lines = fp.readlines()
        fp.close()

    return imported_lines


def import_jeol_wdscans(subfolder,
                        scan_filename='data001_mm.csv',
                        cnd_filename='data001.cnd',
                        comment_line_num=78,
                        crystal_line_name="$XM_WDS_CRYSTAL_NAME%0",
                        sep=',',
                        return_metadata=False):

    """ Imports a WD scan from files produced from a JEOL instrument.
    
    The data is assumed to be in the form of two columns, the first being L (mm) 
    and the second being raw counts. The accompanying .cnd file is needed to extract
    comments, beam current, accumulation and counting times so as to convert counts
    to cps and cps_per_nA.
    
    It is assumed that current is given in amps and counting time is given in ms.

    Arguments
    ---------
    subfolder : str or pathlib.Path object
        Path to the folder containing a .csv (or .txt) and .cnd file
    scan_filename : str
        Name of the wd scan file
    cnd_filename : str or Path object
        Name of the cnd file (this contains comment and current into)
    comment_line_num : int
        Line number containing the comment in the cnd file
        (counting from 1, as in first line = line 1; this can be checked in a text editor)
    sep : str
        Separator in the wd scan file. Default is ',' for csv, but in
        txt files this may be '\t' for tab-delimited.
    return_metadata : bool
        Default is False. If True, return a pandas series containing metadata from the
        cnd file (comment, dwell time, crystal, current, etc) as a third output.
        

    Notes: Defaults are set based on what ANU's JEOL probe uses.

    Returns
    -------
    comments : str
        Comment that was in the .cnd file
    data: pd.DataFrame
        DataFrame with columns L, cps, cps_per_nA
    """

    if type(subfolder) == str:
        subfolder = Path(subfolder)

    my_data_file = subfolder / scan_filename
    my_info_file = subfolder / cnd_filename
    
    # TODO: MAKE THESE CHECK WORK
    # Check that file and folder names area all correct, provide hints if not.
    if subfolder.exists() == False:
         raise FileNotFoundError('Incorrect folder name {}'.format(subfolder))
    elif my_data_file.exists() == False:
        folder_contents = list(subfolder.glob('*'))
        raise FileNotFoundError(
               'Did not find the file {} under folder {}. \nInstead found {}. \nCheck file structure.'.format(
                     scan_filename, subfolder, folder_contents))

    data = pd.read_csv(my_data_file, names=['L', 'counts'], sep=sep)

    info = import_file(my_info_file)

    # Extract info from the cnd file:
    
    comments = info[comment_line_num - 1].strip().split()[1:]
    
    if len(comments) > 1:
        comments = '_'.join(comments)
    elif len(comments) == 1:
        comments = comments[0]
    else:
        comments = 'no_comment'      

    current_line = next(l for l in info if 'CURRENT' in l)
    current_A = float(current_line.strip().split()[1])
    current_nA = current_A * 10 ** 9  # Convert to nA
    
    dwell_time_line = next(l for l in info if 'DWELL_TIME' in l)
    dwell_time_ms = float(dwell_time_line.replace('msec', '').split()[-1])
    dwell_time_s = dwell_time_ms / 1000 # convert to seconds
    
    accumulations_line = next(l for l in info if 'ACCUM' in l)
    accumulations = int(accumulations_line.split()[-1])
    
    num_pts_line = next(l for l in info if 'POINT' in l)
    num_pts = int(num_pts_line.split()[-1])
    
    total_time_mins = round((dwell_time_s * num_pts) / 60)
    
    step_line = next(l for l in info if 'STEP' in l)
    step = float(step_line.replace('um', '').split()[-1])
    
    # Convert counts to cps and cps_per_nA
    
    data['cps'] = data['counts'] / (dwell_time_s * accumulations)
    data['cps_per_nA'] = data.cps / current_nA
    
    # Find crystal name if possible:
    crystal_line = next(l for l in info if crystal_line_name in l)
    crystal = crystal_line.split()[-1]
    
    metadata = pd.Series({'folder': subfolder,
                          'comment': comments,
                          'current_nA': current_nA,
                          'dwell_time_s': dwell_time_s,
                          'accumulations': accumulations,
                          'num_points': num_pts,
                          'total_time_mins': total_time_mins,
                          'crystal': crystal,
                          'step_size_um': step})

    if return_metadata:
    
        return comments, data, metadata
        
    else:
        
        return comments, data


def get_comment(filepath):
    if type(filepath) == str:
        filepath = Path(filepath)
    wt_lines = import_file(filepath / '1.wt')
    comment_line = next(l for l in wt_lines if 'Comment :' in l)
    comment = comment_line.split('Comment :')[-1].strip()
    return comment


class WtFile:
    """ Stores data from a 1.wt file """

    def __init__(self):
        self.accV = None
        self.nA = None
        self.spot_size_um = None
        self.data = None
        self.date = None

def read_wt_file(path, bgi=False):
    """ Reads the contents of a 1.wt file from a JEOL EPMA system.

        Inputs: path the 1.wt file (string or Path object)
        Returns: an object containing all relevent information from 1.wt file

    """

    if type(path) == str:
        path = Path(path)

    wt_lines = import_file(path / '1.wt')

    idx_Element = [i for i, line in enumerate(wt_lines) if 'Element' in line]
    current = float(
        wt_lines[idx_Element[0] - 1].rsplit('Curr.(A) : ')[1].strip())

    current_nA = current * 10 ** 9

    accV = float(wt_lines[7].split(':')[1].rsplit('(kV)')[0].strip())
    spot_size = float(wt_lines[7].split(':')[2].rsplit('Scan')[0].strip())
    date = wt_lines[8].split(' ')[3]
    date = date.replace('/', '-')

    # Find the first block of information in the file
    idx_emptyline = [i for i, line in enumerate(wt_lines) if line == '\n']

    block1 = wt_lines[idx_Element[0] + 1:idx_emptyline[4]]

    # Find all the numbers in this block (negative, integer and floats
    # all included) but no other characters
    splitlines = [re.findall(r'[-+]?\d*\.\d+|\d+', line)[1:]
                  for i, line in enumerate(block1)]
    elementList = [re.findall('[a-zA-Z]+', line)[0]
                   for i, line in enumerate(block1)]

    wt_file_headers = ['pos', 'net_cps', 'lwr_cps',
                       'upr_cps', 'stdev_net_cps',
                       'dl_ppm', 'kraw_pcnt']

    # BGI wt files do not have kraw column. Set bgi=True to handle this.
    if bgi:
        wt_file_df = pd.DataFrame.from_records(splitlines,
                                           columns=wt_file_headers[:-1])
    else:
        wt_file_df = pd.DataFrame.from_records(splitlines,
                                           columns=wt_file_headers)

    wt_file_df.insert(0, 'element', pd.Series(elementList,
                                              index=wt_file_df.index))

    # Convert strings to floats
    for i in wt_file_df:
        if i != 'element':
            wt_file_df[i] = wt_file_df[i].astype('float')

    # Store data in an object of class WtData
    wtfile = WtFile()
    wtfile.nA = current_nA
    wtfile.accV = accV
    wtfile.spot_size_um = spot_size
    wtfile.data = wt_file_df
    wtfile.date = date

    return wtfile


def read_mes_file_bgi(path):
    """ Reads the contents of a 1.mes file from BGI

        Inputs: path the 1.mes file (string or Path object)
        Returns: an object containing all relevent information from file

    """
    if type(path) == str:
        path = Path(path)

    mes_lines = import_file(path / '1.mes')

    # Find the elements and x-ray lines----------------
    idx_emptyline = [i for i, line in enumerate(mes_lines) if line == '\n']
    idx_Element = [i for i, line in enumerate(mes_lines) if 'Element' in line]
    block1 = mes_lines[idx_Element[0] + 1:idx_emptyline[3]]

    block1_strip = [re.sub(' +', ' ', line.strip())
                    for i, line in enumerate(block1)]
    block1_split = [line.split(' ') for i, line in enumerate(block1_strip)]
    [line.remove('(') for i, line in enumerate(block1_split)]
    block1_split = [line[1:] for i, line in enumerate(block1_split)]

    block1_headers = ['element', 'xray', 'crystal', 'ch', 'acc_V',
                      'pk_pos_standard', 'nm', 'lwr_pos_rel', 'upr_pos_rel']

    mesfile = pd.DataFrame.from_records(block1_split,
                                        columns=block1_headers)

    # Find the counting times ---------------
    block2 = mes_lines[idx_Element[1] + 1:idx_emptyline[4]]

    block2_strip = [re.sub(' +', ' ', line.strip())
                    for i, line in enumerate(block2)]
    block2_split = [line.split(' ') for i, line in enumerate(block2_strip)]
    pk_time = [line[2:3][0] for i, line in enumerate(block2_split)]
    bg_time = [line[3:4][0] for i, line in enumerate(block2_split)]

    mesfile['pk_time'] = pk_time
    mesfile['bg_time'] = bg_time

    return mesfile

def read_cor_file_bgi(path, mesfile):
    """ Reads the contents of a 1.cor file from BGI

        Inputs: path the 1.cor file (string or Path object)
        Returns: an object containing all relevent information from file

    """
    if type(path) == str:
        path = Path(path)

    cor_lines = import_file(path / '1.cor')

    idx_emptyline = [i for i, line in enumerate(cor_lines) if line == '\n']
    idx_Element = [i for i, line in enumerate(cor_lines) if 'Element' in line]

    # Find the names of the standards used for each element------------
    # Find index of lines containing 'Standard Data'
    block1 = cor_lines[idx_Element[0] + 1:idx_emptyline[1]]
    print(block1)

    split_block1 = [re.findall('[a-zA-Z0-9_]+', line)
                    for i, line in enumerate(block1)]
    listStandards = [None] * len(mesfile.element)

    for i, entry in enumerate(split_block1):
        entry_strip = re.sub(r'\d', '', entry[1])
        entry_strip = entry_strip.replace('O', '')
        for j, el in enumerate(mesfile.element):
            if el == entry_strip:
                listStandards[j] = entry[2]

     # Find the beam currents and cps ----
    block2 = cor_lines[idx_emptyline[1] + 3:]

    listBeamCurrents = [re.findall('\d*\.\d+E-\d+', line)[0]
                        for i, line in enumerate(block2)]
    list_cps_all = [re.findall('\d*\.\d+|\d+', line)[3:7]
                    for i, line in enumerate(block2)]
    corfile = pd.DataFrame.from_records(
        list_cps_all,
        columns=['net_cps', 'lwr_cps', 'upr_cps', 'stdev_net_cps_standard'])

    corfile['nA'] = pd.Series(listBeamCurrents, index=corfile.index)
    corfile.nA = corfile.nA.astype(float)
    corfile.loc[:, 'nA'] = corfile.loc[:, 'nA'] * 10 ** 9
    corfile['standard'] = pd.Series(listStandards, index=corfile.index)

    # combine with the mesfile
    corfile = pd.concat([corfile, mesfile], axis=1)

    # Convert strings to floats ------------------
    tofloat = ['acc_V', 'pk_pos_standard', 'nm', 'lwr_pos_rel', 'upr_pos_rel',
               'nA',
               'net_cps', 'lwr_cps', 'upr_cps', 'stdev_net_cps_standard',
               'pk_time', 'bg_time']

    for i, label in enumerate(tofloat):
        corfile[label] = corfile[label].astype('float')


    return corfile


def read_cor_file(path):
    """ Reads the contents of a 1.cor file from a JEOL EPMA system.

        Inputs: path the 1.cor file (string or Path object)
        Returns: an object containing all relevent information from 1.wt file

    """
    if type(path) == str:
        path = Path(path)

    cor_lines = import_file(path / '1.cor')

    # Find the number of elements analysed
    first_entry_in_non_blank_lines = [line.strip().split()[0]
                                      for line in cor_lines
                                      if line.strip() != '']
    num_elements = max([int(n) for n in first_entry_in_non_blank_lines
                        if n.isnumeric()])

    # Find the elements and x-ray lines----------------

    idx_block1_start = next(i for i, line in enumerate(cor_lines)
                            if line.strip() == 'Measurement Condition')
    idx_block1_start += 3 # The block starts 3 lines down from the header
    block1 = cor_lines[idx_block1_start : idx_block1_start + num_elements]
    block1_strip = [re.sub(' +', ' ', line.strip())
                    for i, line in enumerate(block1)]
    block1_split = [line.split(' ') for i, line in enumerate(block1_strip)]
    [line.remove('(') for i, line in enumerate(block1_split)]
    block1_split = [line[1:] for i, line in enumerate(block1_split)]

    block1_headers = ['element', 'xray', 'crystal', 'ch', 'acc_V',
                      'pk_pos_standard', 'nm', 'lwr_pos_rel', 'upr_pos_rel']

    corfile = pd.DataFrame.from_records(block1_split,
                                        columns=block1_headers)

    # Find the counting times ---------------
    idx_block2_start = idx_block1_start + num_elements + 2
    block2 = cor_lines[idx_block2_start : idx_block2_start + num_elements]
    block2_strip = [re.sub(' +', ' ', line.strip())
                    for i, line in enumerate(block2)]
    block2_split = [line.split(' ') for i, line in enumerate(block2_strip)]
    pk_time = [line[2:3][0] for i, line in enumerate(block2_split)]
    bg_time = [line[3:4][0] for i, line in enumerate(block2_split)]

    corfile['pk_time'] = pk_time
    corfile['bg_time'] = bg_time

    # Find the beam currents and cps ----
    idx_block5_start = next(i for i, line in enumerate(cor_lines)
                            if line.strip() == 'Standard Intensity of WDS')
    idx_block5_start += 2 # The block starts 2 lines down from the header
    block5 = cor_lines[idx_block5_start : idx_block5_start + num_elements]
    listBeamCurrents = [re.findall('\d*\.\d+E-\d+', line)[0]
                        for i, line in enumerate(block5)]
    list_cps_all = [re.findall('\d*\.\d+|\d+', line)[3:7]
                    for i, line in enumerate(block5)]
    cps_all_df = pd.DataFrame.from_records(list_cps_all,
                        columns=['net_cps', 'lwr_cps', 'upr_cps',
                                 'stdev_net_cps_standard'])
    corfile['nA'] = pd.Series(listBeamCurrents, index=corfile.index)
    corfile.nA = corfile.nA.astype(float)
    corfile.loc[:, 'nA'] = corfile.loc[:, 'nA'] * 10 ** 9
    corfile = pd.concat([corfile, cps_all_df], axis=1)

    # Find the names of the standards used for each element------------
    # Find index of lines containing 'Standard Data'
    idx_block4_start = next(i for i, line in enumerate(cor_lines)
                        if line.strip() == 'Standard Data')
    idx_block4_start += 2 # block starts 2 lines down from header
    # Note: block4 is one line shorter than the other blocks, because
    # N is present twice in other blocks but only once in this block
    # hence the -1 in the next line.
    block4 = cor_lines[idx_block4_start : idx_block4_start + (num_elements - 1)]
    split_block4 = [re.findall('[a-zA-Z0-9_]+', line)
                    for i, line in enumerate(block4)]
    listStandards = [None] * len(corfile.element)

    for i, entry in enumerate(split_block4):
        entry_strip = re.sub(r'\d', '', entry[1])
        entry_strip = entry_strip.replace('O', '')
        for j, el in enumerate(corfile.element):
            if el == entry_strip:
                listStandards[j] = entry[2]

    corfile['standard'] = pd.Series(listStandards, index=corfile.index)

    # Convert strings to floats ------------------
    tofloat = ['acc_V', 'pk_pos_standard', 'nm', 'lwr_pos_rel', 'upr_pos_rel',
               'nA',
               'net_cps', 'lwr_cps', 'upr_cps', 'stdev_net_cps_standard',
               'pk_time', 'bg_time']

    for i, label in enumerate(tofloat):
        corfile[label] = corfile[label].astype('float')

    return corfile


def find_files_and_folders(samples, sample_folders, apf_file=None, wd_scan=None):
    """ Generates a list of the data files in the data_folder, extracts the
     comments from these files and finds the appropriate N wavelength scan file

    Inputs:
        - samples: a list of samples contained in the folder, e.g. ['budd'], or
         ['Edi01', 'Edi02']
        - data_folder: path to the folder containing the data
        - apf_file: path or float. Can be path to the csv file containing apf values for each sample.
            Default (None) assumes that apf correction is not known and
            apf will be set to 1 for all samples. 
        - wd_scan: path to file containing parameters to use for background correction

    """
    if apf_file is not None:
        apf_df = pd.read_csv(apf_file, header=0, index_col=0,
                          squeeze=True)
        apf = apf_df.apf.to_dict()
        apf_sd = apf_df.apf_sd.to_dict()

    if type(apf_file) == str:
        apf_file = Path(apf_file)

    full_datalist = []

    for j, sample in enumerate(samples):

        folderlist = sorted(list(sample_folders[j].glob('[!.]*')))
        folderlist = [f for f in folderlist if f.is_dir()]

        commentlist = [None] * len(folderlist)

        for i, folder in enumerate(folderlist):
            commentlist[i] = get_comment(folder)
        print('Comments found:', commentlist)

        datalist = pd.DataFrame({'folder': folderlist, 'comment': commentlist})
        datalist["sample"] = sample
        datalist["paramfile"] = pd.Series([None] * len(datalist.index))
        datalist["apf"] = pd.Series([None] * len(datalist.index))

        for i, comment in enumerate(datalist.comment):

            # Add the wd scan data
            datalist.loc[i, "paramfile"] = wd_scan

            # Add the apf information to the table
            if apf_file is not None:
                datalist.loc[i, "apf"] = apf[sample]
                datalist.loc[i, "apf_sd"] = apf_sd[sample]
            else:
                datalist.loc[i, "apf"] = 1
                datalist.loc[i, "apf_sd"] = 0

        full_datalist.append(datalist)

    return pd.concat(full_datalist, axis=0).reset_index(drop=True)


def save_comment_list(folder, append=False):

    """Saves a list of comments in 1.wt files to a text file.
    Appends to the text file, if requested """

    if type(folder) == str:
        folder = Path(folder)

    data_folders = sorted(list(folder.glob('[!.]*')))
    print(len(data_folders))
    data_folders = [f for f in data_folders if f.is_dir()]
    print(len(data_folders))
    print(data_folders)

    comments = [get_comment(f) for f in data_folders]
    comments = [c.replace(' ','_') for c in comments]  # Replace spaces with '_'

    if append:
        mode = 'a'
    else:
        mode = 'w'
        
    with open('comment_list.txt', mode, newline='') as f:
    
        for i in range(len(data_folders)):
            f.write(str(data_folders[i].parts[1]) + ': ' + comments[i] + '\n')

def rename_folders_as_comments(folder, new_folder):
    """Rename raw data folders based on the comments stored in the 1.wt file."""

    if type(folder) == str:
        folder = Path(folder)
    if type(new_folder) == str:
        new_folder = Path(new_folder)

    data_folders = sorted(list(folder.glob('[!.]*')))
    # remove anything that isn't a directory from the list
    data_folders = [Path(f) for f in data_folders if Path(f).is_dir()]
    print(data_folders)

    comments = [get_comment(f) for f in data_folders]
    comments = [c.replace(' ','_') for c in comments]  # Replace spaces with '_'

    for i in range(len(data_folders)):

            oldname = data_folders[i].parts[1] # Remove first part of path
            newname = Path.joinpath(new_folder, oldname + '_' + comments[i])
            print('old name: {} \nnew name: {}'.format(data_folders[i], newname))

            shutil.copytree(data_folders[i], newname)
            # os.rename(data_folders[i], newname)

def read_and_organise_data(info, bgi=False, save=True):
    """ Reads data from raw files and returns data organised into
    four dataframes/series. Optionally saves these as files.
        """

    wtfile = read_wt_file(info.folder, bgi=bgi)

    if bgi:
        mesfile = read_mes_file_bgi(info.folder)
        corfile = read_cor_file_bgi(
            info.folder, mesfile)

        # BGI wt file doesn't having kraw_pcnt, but we can calculate it!
        kraw_pcnt, _ = correct_quant.calc_kraw_stdev(wtfile.data.net_cps,
                                      corfile.net_cps,
                                      wtfile.data.stdev_net_cps,
                                      corfile.stdev_net_cps_standard,
                                      wtfile.nA, corfile.nA)
        wtfile.data['kraw_pcnt'] = kraw_pcnt

    else:
        corfile = read_cor_file(info.folder)


    columns_from_wt = wtfile.data.loc[:, ['element', 'lwr_cps', 'upr_cps']]
    columns_from_cor = corfile.loc[:, ['lwr_pos_rel', 'upr_pos_rel',
                                       'bg_time']]

    columns_from_cor.rename(columns={'bg_time': 'time'}, inplace=True)

    bg = pd.concat([columns_from_wt, columns_from_cor],
                        sort=False, axis=1)

    standard = corfile.loc[:, ['element', 'xray', 'crystal',
                                    'pk_pos_standard', 'nA', 'net_cps',
                                    'lwr_cps', 'upr_cps', 'standard',
                                    'stdev_net_cps_standard']]

    standard.rename(
        columns={'pk_pos_standard': 'pos',
                 'stdev_net_cps_standard': 'stdev_net_cps'}, inplace=True)


    columns_from_wt = wtfile.data.loc[:, ['element', 'pos',
                                          'net_cps', 'stdev_net_cps',
                                          'dl_ppm', 'kraw_pcnt']]

    peak = pd.concat([columns_from_wt, corfile['pk_time']],
                          sort=False, axis=1)

    peak.rename(columns={'pk_time': 'time'},
                     inplace=True)  # Rename for simplicity

    # Handle the 'fudge' background used for Rb,
    # if it's present in the analysis
    if not bgi:
        Rb_idx = [idx for idx, el in enumerate(bg.element) if
              'Rb' in el]
        if len(Rb_idx) > 0:
            bg.loc[Rb_idx, 'lwr_cps'] = bg.loc[
                Rb_idx, 'upr_cps']


    info['nA'] = wtfile.nA
    info['date'] = wtfile.date
    info['spot_size_um'] = wtfile.spot_size_um

    if save == True:

        # Check if folder exists, if not, create it
        folder_out = Path('./extracted_data/{}_{}/'.format(
                                    info['sample'], info['comment']))
        CHECK_FOLDER = os.path.isdir(folder_out)
        if not CHECK_FOLDER:
            os.makedirs(folder_out)

        info.to_csv(folder_out / 'info.csv')
        peak.to_csv(folder_out / 'peak.csv')
        bg.to_csv(folder_out / 'bg.csv')
        standard.to_csv(folder_out / 'standard.csv')

        print('Saved extracted data as csvs for ' + info['sample'] + '_' + info['comment'])

    return peak, bg, standard, info


# Find the valence file in the calczaf folder

def read_valence_file(subfolder, pattern='valence*'):

    helper_funs.check_calczaf_folder_exists(subfolder)

    valence_file = sorted(subfolder.glob(pattern))

    if len(valence_file) == 0:
        raise Exception('No valence file in calczaf_files folder. Please add one.')
    elif len(valence_file) > 1:
        raise Exception('More than one valence files in calczaf_files folder. Please pick the correct one and delete the rest.')
    else:
        valence_dict = pd.read_csv(valence_file[0], header=0, index_col=0,
                                   squeeze=True).to_dict()

    return valence_dict