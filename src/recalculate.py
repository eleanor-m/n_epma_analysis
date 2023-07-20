"""
TODO: currently input dataframe must have 'Fe2O3' and 'FeO' entries, in that order.
It would be good to make that optional.
TODO: I shouldn't really need to use pyrolite here...

TODO: make garnet endmembers - optional of which to use either Fe3, Fe2 or Fetot but
don't calculate them all.

"""
import numpy as np
import pandas as pd
import periodictable
import multiprocessing
import time

# # Define functions

def ox_strings_to_el_strings(ox_string_list):

    ''' Takes a list of oxide strings,
    Returns a list of corresponding element strings'''

    els = []

    for ox in ox_string_list:
        # Get a list of atoms in the oxide (this line returns a
        # list of periodictable.elements objects)
        atoms = list(periodictable.formula(ox).atoms.keys())
        # Turn the list of objects into a list of strings
        atom_symbols = [a.symbol for a in atoms]
        cations = [a for a in atom_symbols if a != 'O']

        if len(cations) > 1:
            raise ValueError('Multiple cations! Help!')

        els.append(cations[0])

    return els


def wt_to_mol_cations(wt, oxide_string):
    if isinstance(wt, str):
        print(
            '    Expected a wt% {} value, but found a string: "{}". Value set to 0.'.format(
                oxide_string, wt))
        wt = 0
    elif wt == None:
        print('    Expected a wt% {} value, but found: "{}". Value set to 0.'.format(
            oxide_string, wt))
        wt = 0

    cation = pyrolite.geochem.ind.get_cations(oxide_string)
    [num_cat, num_oxy] = get_num_cations_oxygens(oxide_string)
    M = cation[0].mass * num_cat + 15.999 * num_oxy

    return num_cat * wt / M


def wt_to_mol_oxygens(wt, oxide_string):
    if isinstance(wt, str):
        wt = 0
    elif wt == None:
        wt = 0

    cation = pyrolite.geochem.ind.get_cations(oxide_string)
    [num_cat, num_oxy] = get_num_cations_oxygens(oxide_string)
    M = cation[0].mass * num_cat + 15.999 * num_oxy

    return num_oxy * wt / M


def get_num_cations_oxygens(oxide_string):
    # Get the cation from the oxide string
    cation = pyrolite.geochem.ind.get_cations(oxide_string)

    # Split the oxide string at the 'O'
    split_O = oxide_string.split('O')

    # What number comes after the oxygen, if any?
    try:
        num_oxygens = int(split_O[1])
    except:
        num_oxygens = 1

    # What number comes after the cation, if any?

    if split_O[0] == cation[0].symbol:
        num_cations = 1;
    else:
        num_cations = int(split_O[0].split(cation[0].symbol)[1])

    return [num_cations, num_oxygens]


def wt_to_mol_cations_df(df):
    oxides = df.index.to_list()  # Get oxide names
    mol_cation = [None] * len(oxides)

    for i, oxide in enumerate(oxides):
        mol_cation[i] = wt_to_mol_cations(df.loc[oxide], oxide)

    df = pd.concat([df, pd.Series(mol_cation, index=oxides, name='mol_cation')], axis=1)

    return df


def wt_to_mol_oxygens_df(df):
    oxides = df.index.to_list()  # Get oxide names
    mol_oxygen = [None] * len(oxides)

    for i, oxide in enumerate(oxides):
        mol_oxygen[i] = np.array(wt_to_mol_oxygens(df.loc[oxide, 'analysis_wt'], oxide))

    df['mol_oxygen'] = mol_oxygen

    return df


def norm_cation(df, num_cat_pfu, num_oxy_pfu):
    if df.analysis_wt['Fe2O3'] > 0:
        norm_cation = num_oxy_pfu * df.mol_cation / df.mol_oxygen.sum()
        print(
            "-----> There's Fe3+ analysed; you'd better check function 'norm_cation' is behaving as expected!")
    else:
        norm_cation = num_cat_pfu * df.mol_cation / df.mol_cation.sum()

    df['norm_cation'] = norm_cation

    return df


def norm_oxygen(df, column, new_column_name):
    # Column specifies the name of the column to use in the normalisation

    norm_oxygen = [None] * len(df.index.to_list())

    for i, oxide in enumerate(df.index.to_list()):
        [num_cat, num_oxy] = get_num_cations_oxygens(oxide)

        norm_oxygen[i] = df.loc[oxide, column] * num_oxy / num_cat

    df[new_column_name] = norm_oxygen

    return df


def atom_units(df, num_oxy_pfu, mineral):
    atom_units = [None] * len(df.index.to_list())
    Fe3_atom_units = calc_Fe3_functions(mineral, df, num_oxy_pfu)

    for i, oxide in enumerate(df.index.to_list()):

        if oxide == 'Fe2O3':
            atom_units[i] = Fe3_atom_units

        elif oxide == 'FeO':
            # If Fe2O3 was analysed, keep the FeO value
            if df.analysis_wt['Fe2O3'] > 0:
                atom_units[i] = df.norm_cation['FeO']

            # If Fe2O3 was not analysed, adjust the FeO value according to the calculated Fe2O3 value
            else:
                atom_units[i] = df.norm_cation['FeO'] - Fe3_atom_units

        else:
            atom_units[i] = df.norm_cation[i]

    df['atom_units'] = atom_units

    return df


def calc_Fe3_garnet(df, num_oxy_pfu):
    calc_charge = calculate_charge(df.norm_cation)

    # If Fe2O3 was analysed, keep that value
    if df.analysis_wt['Fe2O3'] > 0:
        Fe3_atom_units = df.norm_cation['Fe2O3']
        print('>>> Fe2O3 was analysed, so is not recalculated')

    # If Fe2O3 was not analysed, calculate it
    else:
        if (2 * num_oxy_pfu - calc_charge) > 0:
            Fe3_atom_units = 2 * num_oxy_pfu - calc_charge
        else:
            Fe3_atom_units = 0

    return Fe3_atom_units


def calc_Fe3_cpx(df, num_oxy_pfu):
    # If Fe2O3 was analysed, keep that value
    if df.analysis_wt['Fe2O3'] > 0:
        Fe3_atom_units = df.norm_cation['Fe2O3']
        print('>>> Fe2O3 was analysed, so is not recalculated')

    # If Fe2O3 was not analysed, calculate it
    else:
        if (num_oxy_pfu - df.norm_oxygen.sum()) > 0:

            if df.norm_oxygen['FeO'] > (2 * (num_oxy_pfu - df.norm_oxygen.sum())):
                Fe3_atom_units = (2 * (num_oxy_pfu - df.norm_oxygen.sum()))
            else:
                Fe3_atom_units = df.norm_oxygen['FeO']

        else:
            Fe3_atom_units = 0

    return Fe3_atom_units


def calc_Fe3_functions(mineral, df, num_oxy_pfu):
    if mineral == 'garnet':
        Fe3_atom_units = calc_Fe3_garnet(df, num_oxy_pfu)
    elif mineral == 'cpx':
        Fe3_atom_units = calc_Fe3_cpx(df, num_oxy_pfu)
    else:
        print('Warning: assuming no Fe and setting a.u. to 0')
        Fe3_atom_units = 0

    return Fe3_atom_units


def calculate_charge(nc):
    # Input is a pandas series 'norm_cation' from the df

    calc_charge = 4 * (nc['SiO2'] + nc['TiO2']) + 3 * (
                nc['Al2O3'] + nc['Fe2O3'] + nc['Cr2O3']) + 2 * (
                              nc['FeO'] + nc['MnO'] + nc['MgO'] + nc['CaO'])

    return calc_charge


def calc_all_endmb(au):
    fe_type = ['Fe3', 'Fetot', 'Fe2']
    endmb_series = [None] * 3

    for i, fe in enumerate(fe_type):
        endmb_series[i] = calc_endmb(au, fe)

    endmb_table = pd.concat(endmb_series, axis=0, sort=False)

    return endmb_table


def calc_cpx_endmb(au):
    wollastonite = au['CaO'] / (au['CaO'] + au['MgO'] + au['FeO']) * 100
    enstatite = au['MgO'] / (au['CaO'] + au['MgO'] + au['FeO']) * 100
    ferrosilite = au['FeO'] / (au['CaO'] + au['MgO'] + au['FeO']) * 100

    aegerine = au['Fe2O3'] / (
                au['Fe2O3'] + (au['Na2O'] - au['Fe2O3']) + au['CaO']) * 100
    jaedite = (au['Na2O'] - au['Fe2O3']) / (
                au['Fe2O3'] + (au['Na2O'] - au['Fe2O3']) + au['CaO']) * 100
    diopside = au['CaO'] / (au['Fe2O3'] + (au['Na2O'] - au['Fe2O3']) + au['CaO']) * 100

    endmb_names = ['wollastonite', 'enstatite', 'ferrosilite', 'aegerine', 'jaedite',
                   'diopside']
    endmb_values = [wollastonite, enstatite, ferrosilite, aegerine, jaedite, diopside]

    endmb_series = pd.Series(data=endmb_values, index=endmb_names)

    return endmb_series


def calc_endmb(au, fe_type):
    # Note: au is atom_units, i.e. df.atom_units
    # Fe_type is either Fe3, Fetot or Fe2

    if fe_type == 'Fe3':

        site_sum_8fold = au['FeO'] + au['MnO'] + au['MgO'] + au['CaO']
        site_sum_6fold = au['TiO2'] + au['Al2O3'] + au['Fe2O3'] + au['Cr2O3']

        almandine = 100 * (au['FeO'] / site_sum_8fold)
        pyrope = 100 * (au['MgO'] / site_sum_8fold)
        grossular = 100 * (au['Al2O3'] / site_sum_6fold) * (au['CaO'] / site_sum_8fold)
        spessartine = 100 * (au['MnO'] / site_sum_8fold)
        uvarovite = 100 * (au['Cr2O3'] / site_sum_6fold) * (au['CaO'] / site_sum_8fold)
        andradite = 100 * (au['Fe2O3'] / site_sum_6fold) * (au['CaO'] / site_sum_8fold)
        CaTi_Gt = 100 * (au['TiO2'] / site_sum_6fold) * (au['CaO'] / site_sum_8fold)


    elif fe_type == 'Fetot':

        site_sum_8fold = au['FeO'] + au['Fe2O3'] + au['MnO'] + au['MgO'] + au['CaO']

        almandine = 100 * ((au['Fe2O3'] + au['FeO']) / site_sum_8fold)
        pyrope = 100 * (au['MgO'] / site_sum_8fold)
        grossular = 100 * (au['CaO'] / site_sum_8fold)
        spessartine = 100 * (au['MnO'] / site_sum_8fold)

        uvarovite = None
        andradite = None
        CaTi_Gt = None

    elif fe_type == 'Fe2':

        site_sum_8fold = au['FeO'] + au['MnO'] + au['MgO'] + au['CaO']

        almandine = 100 * (au['FeO'] / site_sum_8fold)
        pyrope = 100 * (au['MgO'] / site_sum_8fold)
        grossular = 100 * (au['CaO'] / site_sum_8fold)
        spessartine = 100 * (au['MnO'] / site_sum_8fold)

        uvarovite = None
        andradite = None
        CaTi_Gt = None

    endmb_names = ['almandine', 'pyrope', 'grossular', 'spessartine', 'uvarovlite',
                   'andradite', 'Ca-Ti_Gt']
    endmb_names = ['{}_{}'.format(fe_type, endmb_name) for endmb_name in endmb_names]
    endmb_values = [almandine, pyrope, grossular, spessartine, uvarovite, andradite,
                    CaTi_Gt]

    endmb_series = pd.Series(data=endmb_values, index=endmb_names)

    return endmb_series


def check_garnet_site_occupancies(au):
    cubic_site_total = au['FeO'] + au['MnO'] + au['MgO'] + au['CaO']
    octahedral_site_total = au['TiO2'] + au['Al2O3'] + au['Cr2O3'] + au['Fe2O3']
    tetrahedral_site_total = au['SiO2']

    data = [tetrahedral_site_total, octahedral_site_total, cubic_site_total]

    data_names = ['tet_total', 'oct_total', 'cubic_total']

    site_occupancies = pd.Series(data=data, index=data_names)

    return site_occupancies


def check_cpx_site_occupancies(au):
    # -------- Tetrahedral site occupancies ------------
    Si_and_Ti_tet = au['SiO2'] + au['TiO2']

    if Si_and_Ti_tet < 2:  # If there's not enough Si and Ti for the site

        if au['Al2O3'] > (
                2 - Si_and_Ti_tet):  # If there's more Al than is needed to fill the site
            tetrahedral_Al = 2 - Si_and_Ti_tet  # Assign the required Al as tetrahedral

        else:  # If there's less than what's needed,
            tetrahedral_Al = au['Al2O3']  # Assign all the Al as tetrahedral

    else:  # If there's too much Si and Ti for the site
        tetrahedral_Al = 0  # Don't add any Al into it

    tetrahedral_site_total = Si_and_Ti_tet + tetrahedral_Al

    # -------- Octahedral site occupancies ------------
    Al_Fe_Cr_Mn_oct = (au['Al2O3'] - tetrahedral_Al) + au['Cr2O3'] + au['Fe2O3'] + au[
        'FeO'] + au['MnO']

    if Al_Fe_Cr_Mn_oct < 1:  # If there are not enough cations for the site

        if (1 - Al_Fe_Cr_Mn_oct) < au[
            'MgO']:  # If there's more Mg than is needed to fill the site
            octahedral_Mg = 1 - Al_Fe_Cr_Mn_oct  # Assign the required Mg to fill the site
        else:  # If there's less than what's needed,
            octahedral_Mg = au['MgO']  # Assign all the Mg to the site

    else:  # If there are too many cations for the site already
        octahedral_Mg = 0  # Don't add any Mg into it

    octahedral_site_total = Al_Fe_Cr_Mn_oct + octahedral_Mg

    # -------- Cubic site occupancies ------------

    cubic_site_total = au['Na2O'] + au['CaO'] + (au['MgO'] - octahedral_Mg)

    # Compile data
    data = [Si_and_Ti_tet, tetrahedral_Al,
            (au['Al2O3'] - tetrahedral_Al),
            (au['Cr2O3'] + au['Fe2O3'] + au['FeO'] + au['MnO']), octahedral_Mg,
            (au['MgO'] - octahedral_Mg), au['CaO'], au['Na2O'],
            tetrahedral_site_total, octahedral_site_total, cubic_site_total]

    data_names = ['tet_SiTi', 'tet_Al', 'oct_Al', 'oct_FeCrMn', 'oct_Mg', 'cubic_Mg',
                  'cubic_Ca', 'cubic_Na',
                  'tet_total', 'oct_total', 'cubic_total']

    site_occupancies = pd.Series(data=data, index=data_names)

    return site_occupancies


def check_phlog_site_occupancies(au, excess_Al_to_octahedral_site=True):
    # -------- Tetrahedral site occupancies ------------
    T_Al = au['Al2O3']  # Starting value for tetrahedral Al is all the Al

    # But we may choose to assign some to the octahedral site:
    if excess_Al_to_octahedral_site == True:
        if au['SiO2'] < 3:
            missing_Si = 3 - au['SiO2']
            if au['Al2O3'] > 1:
                T_Al = 1 + missing_Si
            else:
                T_Al = au['Al2O3']
        else:
            T_Al = au['Al2O3']
    else:
        T_Al = au['Al2O3']

    T_Si = au['SiO2']
    T_total = T_Si + T_Al

    # -------- Octahedral site occupancies ------------
    M_Mg = au['MgO']

    if excess_Al_to_octahedral_site == True:
        M_Al = au['Al2O3'] - T_Al

    M_total = M_Mg + M_Al

    # -------- Interlayer site occupancies ------------
    A_total = au['K2O'] + au['Na2O'] + au['CaO'] + au['Rb2O'] + au['N2O']

    # Compile data
    data = [T_Si, T_Al, T_total, M_Mg, M_Al, M_total, A_total]

    data_names = ['T_Si', 'T_Al', 'T_total', 'M_Mg', 'M_Al',
                  'M_total', 'A_total']

    site_occupancies = pd.Series(data=data, index=data_names)

    return site_occupancies


# ----------- Overall functions --------
def process_garnet_data(df):
    """ Input is a pandas series that should have oxides (strings) as the index
      and wt% oxides as the values. The name of the data should be analysis_wt.

            analysis_wt
      SiO2  37.39
      TiO2  0.28
      Al2O3 10.08
      Fe2O3 0
      FeO   19.35
      Cr2O3 0
      MgO   1.55
      CaO   29.63
      """

    df = wt_to_mol_cations_df(df)
    df = wt_to_mol_oxygens_df(df)
    df = norm_cation(df, num_cat_pfu=8, num_oxy_pfu=12)
    df = norm_oxygen(df, 'norm_cation', 'norm_oxygen')
    df = atom_units(df, num_oxy_pfu=12, mineral='garnet')
    df = norm_oxygen(df, 'atom_units', 'norm_oxygen2')
    endmb_table = calc_all_endmb(df.atom_units)
    site_occupancies = check_garnet_site_occupancies(df.atom_units)

    return df, endmb_table, site_occupancies


def process_cpx_data(df):
    """ Input is a pandas series that should have oxides (strings) as the index
      and wt% oxides as the values. The name of the data should be analysis_wt.

            analysis_wt
      SiO2  54.64
      TiO2  0.08
      Al2O3 0.90
      Fe2O3 0
      FeO   7.45
      Cr2O3 0.17
      MgO   18.93
      CaO   16.76

      """

    # Add elements that aren't already in the series:
    data_index = set(df.index.to_list())
    expected_oxides = {'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO',
                       'MgO', 'CaO', 'Na2O', 'K2O'}
    missing_oxides = expected_oxides.difference(data_index)

    for ox in missing_oxides:
        df[ox] = 0.0

    df = wt_to_mol_cations_df(df)
    df = wt_to_mol_oxygens_df(df)
    df = norm_cation(df, num_cat_pfu=4, num_oxy_pfu=6)
    df = norm_oxygen(df, 'norm_cation', 'norm_oxygen')
    df = atom_units(df, num_oxy_pfu=6, mineral='cpx')
    endmb_series = calc_cpx_endmb(df.atom_units)
    site_occupancies = check_cpx_site_occupancies(df.atom_units)

    return df, endmb_series, site_occupancies


def process_phlog_data(df):
    """ Input is a pandas series that should have oxides (strings) as the index
      and wt% oxides as the values. The name of the data should be analysis_wt.

            analysis_wt
        SiO2	47.777
        Al2O3	14.136
        CaO	    0.116
        MgO	    20.692
        K2O	    9.241
        Na2O	0.894
        Rb2O	0.176
        MoO3	0.006
        N2O	    0.662
        H2O	    6.298
      """
    # Add elements that aren't already in the series:
    data_index = set(df.index.to_list())
    # expected_oxides = {'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO',
    #                      'MgO', 'CaO', 'Na2O', 'K2O'}
    expected_oxides = {'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO',
                       'MgO', 'CaO', 'Na2O', 'K2O', 'H2O', 'Rb2O'}
    missing_oxides = expected_oxides.difference(data_index)
    print('adding missing oxides {} to dataset'.format(missing_oxides))

    for ox in missing_oxides:
        df[ox] = 0.0

    df = wt_to_mol_cations_df(df)
    df = wt_to_mol_oxygens_df(df)
    df = norm_cation(df, num_cat_pfu=8, num_oxy_pfu=11)
    df = norm_oxygen(df, 'norm_cation', 'norm_oxygen')
    df = atom_units(df, num_oxy_pfu=11, mineral='phlogopite')
    endmb_series = None  # I won't bother calculating the end-members for phlogopite
    # because it should already be an endmember
    site_occupancies = check_phlog_site_occupancies(df.atom_units,
                                                    excess_Al_to_octahedral_site=True)

    return df, endmb_series, site_occupancies


def recalc_w_multiprocessing(min_db, mineral, fun_to_apply, num_processes=4):
    """ Arguments:
                    minerals_df = 'garnets' or 'cpxs' dataframe. This should be
                                    the full database but with only the required
                                     minerals selected.
                    fun_to_apply = function to use, varies depending on mineral

        Returns:
            df_out - list of pd.dataframes
            endmbs - list of pd.series?
            site_occupancies - list of pd.series?

        """

    # !!! Trial - remove all columns except for "oxides_needed" - to compare with tests
    if mineral == 'garnet':
        oxides_to_extract = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO",
                             "CaO"]
    elif mineral == 'cpx':
        oxides_to_extract = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO",
                             "CaO", "Na2O", "K2O"]
    # !!! ----

    # Construct a dataframe that contains the needed oxides
    minerals_oxides = min_db[oxides_to_extract]  # Get the columns that have names
    # in the 'oxides_to_extract' list

    # Make a new column called "Fe2O3" and fill it with zeros, place prior to FeO column
    minerals_oxides.insert(0, 'Fe2O3', 0)

    # Set up a list, each entry is a row of the oxides df.
    ox_list = [minerals_oxides.loc[row, :] for row in minerals_oxides.index]

    # Change the name of the row (pd.Series) to enable it to work in the process_data function
    for row in ox_list:
        row.name = 'analysis_wt'

    # Perform the calculations with multiprocessing, and time how long it takes
    start_time = time.time()

    with multiprocessing.Pool(num_processes) as pool:
        result = pool.map(fun_to_apply, ox_list)

    end_time = time.time()
    print('Done {} samples in {:.1f} min'.format(len(ox_list),
                                                 (end_time - start_time) / 60))

    df_out = [res[0] for res in result]
    endmbs = [res[1] for res in result]
    site_occupancies = [res[2] for res in result]

    return df_out, endmbs, site_occupancies


def calc_Mg_Cr_nums(wt):
    """ Takes a pandas series 'wt' containing wt% oxide data
    Calculates mol% oxide data
    Returns Mg num and Cr nums."""

    # Get molar masses for each index string
    M_element = [periodictable.formula(element_string).mass
                 for element_string in wt.index]

    working_df = pd.DataFrame({'wt%': wt, 'M': M_element})
    working_df['num_moles'] = working_df['wt%'] / working_df['M']
    working_df['mol%'] = 100 * (
            working_df['num_moles'] / working_df['num_moles'].sum())

    mol_pcnt = working_df['mol%']

    Mg_num = 100 * mol_pcnt['MgO'] / (mol_pcnt['MgO'] + mol_pcnt['FeO'])

    Cr_num = 100 * mol_pcnt['Cr2O3'] / (mol_pcnt['Cr2O3'] + mol_pcnt['Al2O3'])

    return Mg_num, Cr_num
