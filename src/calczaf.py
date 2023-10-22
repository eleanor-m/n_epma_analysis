import pandas as pd
import csv
import periodictable
import numpy as np
from pathlib import Path

# %% Write CalcZAF input -----------------------

def write_calczaf_input(
        mySpot, fileName, valence, standard_database, accV, 
        calcMode=2, taAngle=40, Oxide_or_Element=1, 
        correct_bg=True, correct_apf=True,
        remove_elements=None,
        elementByDifference=None, elementByStoichToStoichOxygen=None,
        stoichOxygenRatio=0, elementByStoichToOtherElement=None,
        OtherElement=None, stoichElementRatio=0,
        definedElements=None, definedElementWts=None
        ):

    """ Writes a calczaf file from an input src.correct_quant.Spot object

    Arguments:
        - mySpot = src.correct_quant.Spot object
        - fileName = string giving path of file to write results to - results
            will be appended to this file, so file should be erased prior to
            calling this function if needed.
        - valence = dictionary specifying valences of each element
        - standard_database = dictionary specifying the database numbers
            for each element's standard in the CalcZAF database
        - accV = acceleranting volatage in kV (integer)
        - calcMode = Integer flag specifying calculation mode:
            0: 'Calculate k-ratio from concentration'
            1: 'Calculate concentration from unknown and standard intensities'
            2: 'Calculate concentration from raw k-ratios' (default)
            3: 'Calculate concentration from normalised k-ratios'
        - taAngle = Take-off Angle (integer, default 40)
        - Oxide_or_Element = Integer flag specifying whether to calculate as
            oxides (1, default) or elements (2). Calculating as oxides will add
            oxygens to each cation according to the valence file.
        - correct_bg = True (default) or False - to use background-corrected data
        - correct_apf = True (default) or False - to use apf-corrected data
        - remove_elements = a list of elements present in the analysis that
            should be removed before applying matrix correction: e.g. ['Mg',
            'Ca']. Defaults to None.
        - elementByDifference = string specifying element symbol (e.g. 'h'),
            or None (default)
        - elementByStoichToStoichOxygen = string specifying element symbol,
            or None (default)
        - stoichOxygenRatio = Integer (or float?), 0 if not used
        - elementByStoichToOtherElement =  string specifying element symbol,
            or None (default)
        - OtherElement =  string specifying element symbol,
            or None (default)
        - stoichElementRatio = Integer (or float?), 0 if not used (default)
        - definedElements = list of strings e.g. ['F', 'Cl'] or None (default)
        - definedElementWts = list of floats e.g. [1.2, 17] or None (default)

    Returns:
        Saves a csv file
     """


    if remove_elements is not None:
        idx = [i for i, el in enumerate(mySpot.corrected.element)
               if el in remove_elements]
        diff_idx = mySpot.corrected.index.difference(idx)
        peak = mySpot.peak.drop(mySpot.idx_N[1]).loc[diff_idx, :]
        corrected = mySpot.corrected.loc[diff_idx, :]
        standard = mySpot.standard.loc[diff_idx, :]
        montecarlo = mySpot.montecarlo.loc[diff_idx, :]

        peak.reset_index(inplace=True, drop=True)
        corrected.reset_index(inplace=True, drop=True)
        standard.reset_index(inplace=True, drop=True)
        montecarlo.reset_index(inplace=True, drop=True)
    else:
        peak = mySpot.peak.drop(mySpot.idx_N[1])
        corrected = mySpot.corrected
        standard = mySpot.standard
        montecarlo = mySpot.montecarlo

    numElements = len(corrected)

    if elementByDifference is not None:
        numElements = numElements + 1
    if elementByStoichToStoichOxygen is not None:
        numElements = numElements + 1
    if elementByStoichToOtherElement is not None:
        numElements = numElements + 1
    if definedElements is not None:
        numElements = numElements + len(definedElements)
    if Oxide_or_Element == 1:
        numElements = numElements + 1

    comment = mySpot.info.comment

    line1 = [calcMode, numElements, accV, taAngle, comment]

    line2 = [Oxide_or_Element, elementByDifference,
             elementByStoichToStoichOxygen, stoichOxygenRatio,
             elementByStoichToOtherElement, OtherElement, stoichElementRatio]

    body_lines = [None] * len(peak.element)

    # Choose the kraw values to use. For background-corrected values,
    # we used values from the 'corrected' dataframe, except for nitrogen
    # which we take from the montecarlo simulation (expect to be slighly more
    # accurate).

    if correct_bg == True:  # Take kraw from the 'corrected' df, and N kraw
        # from the 'montecarlo' df
        if correct_apf == True:
            kraw_val = corrected.kraw_apf_pcnt / 100
            kraw_nitrogen = montecarlo.kraw_apf_pcnt[
                                montecarlo.element == 'N'].values[0] / 100
        else:
            kraw_val = corrected.kraw_pcnt / 100
            kraw_nitrogen = montecarlo.kraw_pcnt[
                                montecarlo.element == 'N'].values[0] / 100
    else:  # if correct_bg == False; Take kraw from the 'peak' df
        if correct_apf == True:
            kraw_val = peak.kraw_apf_pcnt / 100
            kraw_nitrogen = kraw_val[peak.element == 'N'].values[0]
        else:
            kraw_val = peak.kraw_pcnt / 100
            kraw_nitrogen = kraw_val[peak.element == 'N'].values[0]

    for i, el in enumerate(corrected.element):
        if el == 'N':
            kraw = kraw_nitrogen
        else:
            kraw = kraw_val[i]

        body_lines[i] = [el.lower(), standard.xray[i].lower(),
                         valence['Cations'][el], valence['Oxygens'][el],
                         standard_database[standard.standard[i]], None, kraw, 0]

    output = [line1, line2] + body_lines

    if elementByDifference is not None:
        el = elementByDifference
        output = output + [[el, None, valence['Cations'][el.capitalize()],
                            valence['Oxygens'][el.capitalize()],
                            None, None, 0, 0]]

    if elementByStoichToStoichOxygen is not None:
        el = elementByStoichToStoichOxygen
        output = output + [[el, None, valence['Cations'][el.capitalize()],
                            valence['Oxygens'][el.capitalize()],
                            None, None, 0, 0]]

    if elementByStoichToOtherElement is not None:
        el = elementByStoichToOtherElement
        output = output + [[el, None, valence['Cations'][el.capitalize()],
                            valence['Oxygens'][el.capitalize()],
                            None, None, 0, 0]]

    if definedElements is not None:
        for i, el in enumerate(definedElements):
            output = output + [[el, None, valence['Cations'][el.capitalize()],
                                valence['Oxygens'][el.capitalize()],
                                None, definedElementWts[i], 0, 0]]

    if Oxide_or_Element == 1:
        output = output + [["o", None, 1, 0, None, None, 0, 0]]

    with open(fileName, 'a', newline='') as resultFile:
        # flag 'a' is append, 'w' is (over)write
        wr = csv.writer(resultFile, delimiter=',',
                        quotechar='\"', quoting=csv.QUOTE_NONNUMERIC)
        wr.writerows(output)

    print('Wrote CALCZAF file for {} to {}'.format(mySpot.info.comment, fileName))




def get_calczaf_export_files(folderpath):
    # Define sample names and data locations:
    if type(folderpath) == str:
        folderpath = Path(folderpath)

    datafiles = sorted(folderpath.glob('*Export.dat'))
    datanames = [datafile.parts[-1].split('_Export.dat')[0] for datafile in
                 datafiles]
                 
    if len(datafiles) == 0:
        raise FileNotFoundError('No files ending in _Export.dat found in the the calczaf_files folder')

    return datafiles, datanames

def read_calczaf_outputs(datafiles):
    # Specify the order in which elements are listed in the output excel file
    element_order = 'Si Al Ca Mg Fe Mn Ti K Na P V Cr Co Ni Cu Zn ' \
                    'Ga Ge Rb Sr Mo Ba Cl N H O'.split(' ')
    element_wt_order = [el + ' WT%' for el in element_order]
    element_at_order = [el + ' AT%' for el in element_order]
    index_order = element_wt_order + ['TOTAL', 'Atomic percent:'] + element_at_order

    data_in = [None]*len(datafiles)
    wtdata = [None]*len(datafiles)
    atdata = [None]*len(datafiles)

    for i in range(len(datafiles)):
        data_in[i] = pd.read_csv(datafiles[i], delimiter='\t')
        data_in[i].columns = [colname.strip()
                           for i, colname in enumerate(data_in[i].columns)]

        # Get the column names that contain the wt% and total values
        wtcols = [colname for colname in data_in[i].columns
                  if 'WT%' in colname or 'TOTAL' in colname]

        atcols = [colname for colname in data_in[i].columns
                  if 'AT%' in colname]

        wtdata[i] = data_in[i].loc[:, wtcols].transpose()
        atdata[i] = data_in[i].loc[:, atcols].multiply(100).transpose()

        original_wtdata_index = wtdata[i].index # Save the index for later checks

        # Reorder the table so elements always come in the same order
        wtdata[i] = wtdata[i].reindex(index_order)
        atdata[i] = atdata[i].reindex(index_order)

        # Drop all lines containing non-analysed elements
        wtdata[i].dropna(inplace=True)
        atdata[i].dropna(inplace=True)

        # Check that no elements have been lost from the dataset in the
        # re-indexing process:
        old_idx = set(original_wtdata_index)
        new_idx = set(wtdata[i].index)

        if old_idx != new_idx:
            raise IndexError(
                'Warning! Elements {} dropped from original index!'.format(
                    old_idx.difference(new_idx)))

        # Clean and label indexes:
        wt_idx = wtdata[i].index.to_list()
        wtdata[i].index = [entry.split(' WT%')[0] for entry in wt_idx]
        wtdata[i].index.name = 'wt% element'
        try:
            wtdata[i].astype('float')
        except ValueError:
            raise ValueError("Check the input file, {}: CalcZAF may have "
                             "appended to existing file rather than "
                             "overwriting.".format(datafiles[i]))

        at_idx = atdata[i].index.to_list()
        atdata[i].index = [entry.split(' AT%')[0] for entry in at_idx]
        atdata[i].index.name = 'at% element'
        atdata[i].astype('float')

        # Add a total column to the at% table
        atdata[i] = atdata[i].append(pd.Series(atdata[i].sum(axis=0), name='TOTAL'))

        print(wtdata)

    return wtdata, atdata

def find_the_nitrogen_oxide(oxide_list):

    starts_with_N = [ox for ox in oxide_list if ox[0] == 'N']
    other_elements_starting_with_N = 'Na, Nb, Nd, Ne, Ni, No, Np'.split(', ')
    nitrogen = [s for s in starts_with_N 
        if s[0:2] not in other_elements_starting_with_N]
    
    return nitrogen[0]

def write_tables(folderpath, wtdata, oxidedata, atdata, datanames, detlim=False):

    if type(folderpath) == str:
        folderpath = Path(folderpath)

    outputfilename = folderpath / 'calczaf_outputs.xlsx'
    open(outputfilename, 'w').close() # Erase contents of file before starting

    with pd.ExcelWriter(outputfilename, engine='xlsxwriter') as writer:

        for i in range(len(wtdata)):

            sheet_name = datanames[i]

            if len(sheet_name) > 31:  #  sheet names must be <= 31 characters
                sheet_name = sheet_name[:31]
                print('truncated sheet name to 31 characters:', sheet_name)

           # Write the data to excel

            wtdata[i] = add_summary_cols(wtdata[i])
            oxidedata[i] = add_summary_cols(oxidedata[i])
            atdata[i] = add_summary_cols(atdata[i])

            if detlim:
                # Only retain rows for N if this is a det-lim analysis
                wt_out = wtdata[i].loc[['N'], :]
                n_oxide = find_the_nitrogen_oxide(oxidedata[i].index.to_list())
                oxide_out = oxidedata[i].loc[[n_oxide], :]
                at_out = atdata[i].loc[['N'], :]
            else:
                wt_out = wtdata[i]
                oxide_out = oxidedata[i]
                at_out = atdata[i]

            wt_out.round(3).to_excel(writer, sheet_name=sheet_name)

            oxide_out.round(3).to_excel(writer, sheet_name=sheet_name,
                                               startrow=len(wt_out.index)+2)

            at_out.round(3).to_excel(writer, sheet_name=sheet_name,
                                               startrow=(len(wt_out.index)
                                                        +len(oxide_out.index)
                                                        +4))

            # Get the xlsxwriter objects for the worksheet and format the column width
            worksheet = writer.sheets[sheet_name]
            worksheet.set_column('A:A', 12)  # Set column A width to 20

    print('saved excel file')


def write_summary_tables(folderpath, wtdata, oxidedata, atdata, datanames):

    if type(folderpath) == str:
        folderpath = Path(folderpath)

    outputfilename = folderpath / 'calczaf_oxides_summary.xlsx'
    open(outputfilename, 'w').close() # Erase contents of file before starting
    avgs = pd.DataFrame([df.mean(axis=1) for df in oxidedata], index=datanames)
    stdev = pd.DataFrame([df.std(axis=1) for df in oxidedata], index=datanames)

    counts = [df.loc['TOTAL',:].count() for df in oxidedata]

    summary_table = avgs.join(stdev, rsuffix='_sd') # join the means and stdev tables
    summary_table = summary_table[sorted(summary_table.columns)] # sort columns alphabetically
    summary_table['num_analyses'] = counts
    summary_table.index.name = 'Sample'

    with pd.ExcelWriter(outputfilename, engine='xlsxwriter') as writer:
        # Write the data to excel
         avgs.round(3).to_excel(writer, sheet_name='mean')
         stdev.round(3).to_excel(writer, sheet_name='stdev')
         summary_table.round(3).to_excel(writer, sheet_name='summary')

         # Get the xlsxwriter objects for the worksheet and format the column width
         # worksheet = writer.sheets[sheet_name]
         # worksheet.set_column('A:A', 12)  # Set column A width to 20

    print('saved excel file - oxides summary')

def convert_element_wt_to_oxide_wt(wt_pcnt, element_string, valence_file):

    # Import the valence information
    valence = pd.read_csv(
        valence_file, header=0, index_col=0, squeeze=True).to_dict()

    # Get the number of cations and oxygens in the oxide formula
    num_cations = valence['Cations'][element_string]
    num_oxygens = valence['Oxygens'][element_string]

    # Get the molar masses
    M_element = periodictable.formula(element_string).mass
    M_oxygen = periodictable.formula('O').mass
    M_oxide = (M_element * num_cations) + (M_oxygen * num_oxygens)

    # Get the formula string
    oxide_string = element_string + str(num_cations) + 'O' + str(num_oxygens)
    oxide_string = oxide_string.replace('1', '') # Remove any '1's in the formula
    oxide_string = oxide_string.replace('O0', '') # Remove any 'zero oxygens'

    if M_oxide != periodictable.formula(oxide_string).mass:
        raise ValueError('Something went wrong when calculating oxide molar mass')

    # Convert from weight% element to weight% oxide
    conversion_factor = (M_oxide / M_element) / num_cations
    oxide_pcnt = wt_pcnt * conversion_factor

    return oxide_pcnt, oxide_string


def convert_oxide_wt_to_element_wt(oxide_pcnt, oxide_string):

    # Get the molar masses
    oxide = periodictable.formula(oxide_string)
    element = oxide.structure[0][1]
    num_cations = oxide.structure[0][0]

    M_element = element.mass
    M_oxide = oxide.mass

     # Convert from weight% oxide to weight% element
    conversion_factor = (M_oxide / M_element) / num_cations
    wt_pcnt = oxide_pcnt / conversion_factor

    print('converted {:.1f}% {} to {:.1f}% {}'.format(
                oxide_pcnt, oxide_string, wt_pcnt, str(element)))

    return wt_pcnt, str(element)

def wtdata_to_oxidedata(wtdata, valence_file, datanames):

    oxidedata = [None] * len(wtdata)
    for i in range(len(wtdata)): # For each dataset
        print(datanames[i])
        ox_series = [None]*len(wtdata[i].columns)

        for j, column in enumerate(wtdata[i].columns): # for each spot analysis
            elements = wtdata[i][column].index.to_list()
            wt_values = wtdata[i][column].values

            ox_values = [None] * (len(elements)-2)
            ox_strings = [None] * (len(elements)-2)

            for k, el in enumerate(elements): # For each element, convert to oxide
                if el != 'TOTAL' and el != 'O':
                    ox_values[k], ox_strings[k] = convert_element_wt_to_oxide_wt(
                                                    wt_values[k], el, valence_file)

            ox_series[j] = pd.Series(ox_values, index=ox_strings)

        oxidedata[i] = pd.concat(ox_series, axis=1)
        oxidedata[i].index.name = 'wt% oxide'
        oxidedata[i] = oxidedata[i].append(pd.Series(oxidedata[i].sum(axis=0),
                                                name='TOTAL'))

        if np.allclose(
                oxidedata[i].loc['TOTAL', :].values,
                wtdata[i].loc['TOTAL',:].values,
                atol=0.01) == False:

            print('\n PROBLEM.... \n', wtdata[i].loc['TOTAL',:])
            print('\n', oxidedata[i].loc['TOTAL',:])
            raise ValueError(
                'Trouble converting element to oxide data - totals do not match')


    return oxidedata


def add_summary_cols(df):
    avg = df.mean(axis=1)
    stdev = df.std(axis=1)
    mini = df.min(axis=1)
    maxi = df.max(axis=1)
    summary_cols = pd.DataFrame({'average': avg,
                                 'stdev': stdev,
                                 'minimum': mini,
                                 'maximum': maxi})
    df_out = pd.concat([df, summary_cols], axis=1)
    return df_out

def process_calczaf_outputs(folderpath, valence_file, detlim = False):
    datafiles, datanames = get_calczaf_export_files(folderpath)
    wtdata, atdata = read_calczaf_outputs(datafiles)
    oxidedata = wtdata_to_oxidedata(wtdata, valence_file, datanames)
    if not detlim:
        write_summary_tables(folderpath, wtdata, oxidedata, atdata, datanames)
    write_tables(folderpath, wtdata, oxidedata, atdata, datanames, detlim=detlim)

    return dict(
        wtdata=dict(zip(datanames, wtdata)),
        atdata=dict(zip(datanames, atdata)),
        oxidedata=dict(zip(datanames, oxidedata))
    )
