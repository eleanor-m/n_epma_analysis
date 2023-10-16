# Import dependencies and setup plot defaults
import json
from lmfit.models import ConstantModel, LorentzianModel
from lmfit import minimize, Parameters, report_fit, fit_report
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from src import readfiles, helper_funs
import os
from pathlib import Path


from scipy.signal import savgol_filter

plt.rcParams['figure.figsize'] = 6, 4
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12


def select_roi(df, roi):
    
    """For a dataframe df, select rows within L ranges defined by roi (regions of interest)"""    
    idx_list = [df[(df.L > r[0]) & (df.L < r[1])].index.to_list() for r in roi]
    
    flat_idx_list = [item for sublist in idx_list for item in sublist]
    
    df_roi = df.loc[flat_idx_list, :]
    
    return df_roi

def plot_wdscan(comments, data, save_to='default', cps_per_nA=True):

    """ Plots a WD scan.

    Arguments
    ----------
    comments : str
        name to give title and file name
    data : pd.DataFrame
        data with columns L and cps
    save_to : str or None
        If None, figure is not saved to file.
        Otherwise, specify a custom path to save to including filename
    """

    fig, ax = plt.subplots()

    if cps_per_nA:
        plt.plot(data.L, data.cps_per_nA, '.k', markersize=1)
    else:
        plt.plot(data.L, data.cps, '.k', markersize=1)

    ylims = ax.get_ylim()
    plt.ylim(-0.1, ylims[1])

    pk_pos = [145, 146]

    plt.fill_between([pk_pos[0], pk_pos[1], pk_pos[1], pk_pos[0]],
                     [-0.1, -0.1, ylims[1], ylims[1]], alpha=0.2, color='grey',
                     linewidth=0)

    plt.xlabel('L (mm)')

    if cps_per_nA:
        plt.ylabel('cps/nA') #is that right?
    else:
        plt.ylabel('cps')

    plt.title(comments)
    plt.xlim(118, 182)

    # ---- Save the figure to a file ----
    if save_to is None:
        print('Figure not saved')

    else: # if path is specified
        plt.savefig(save_to)
        print('Saved figure for ' + comments)

    return(fig, ax)

def average_spectra(subfolders):
    """Average multiple WD spectra for the same sample.
    Arguments:
        subfolders: list of paths (strings or Path objects) corresponding
                    to the folders containing spectra that need to be
                    averaged."""

    comments = [None] * len(subfolders)
    data = [None] * len(subfolders)

    # Import the data to pandas dataframes
    for i in range(len(subfolders)):
        comments[i], data[i] = readfiles.import_jeol_wdscans(subfolders[i])

    # Rehape data into np arrays
    cpsdata = np.empty(shape=(len(data[0].cps), len(subfolders)))
    cps_per_nA_data = np.empty(shape=(len(data[0].cps), len(subfolders)))
    Ldata = np.empty(shape=(len(data[0].L), len(subfolders)))

    for i in range(len(subfolders)):
        cpsdata[:, i] = data[i].cps
        cps_per_nA_data[:,i] = data[i].cps_per_nA
        Ldata[:, i] = data[i].L

    # Check L-axes are the same:
    L_axis_check = np.array([Ldata[:,n] == Ldata[:,0]
                             for n in range(len(subfolders))])


    if L_axis_check.all().all() == False:
        print(Ldata[1:10, :])
        raise Exception('L values are not the same between datasets!')

    # Average the data
    mean_df = pd.DataFrame({'L': Ldata[:,0],
                            'cps': cpsdata.mean(axis=1),
                            'cps_per_nA': cps_per_nA_data.mean(axis=1)})

    return mean_df, comments

# %% Functions to fit background ----------------------

def get_bg_idx(data, bg_pos):
    bg_idx = [None] * 4
    for n in range(4):
        # Find the indices of Laxis corresponding to bg positions
        if bg_pos[n] > data.L.max():
            raise ValueError(
                f'Background position {n} of {bg_pos[n]} mm is too high. \
                    Actual data ends at {data.L.max()} mm.')
        elif bg_pos[n] < data.L.min():
            raise ValueError(
                f'Background position {n} of {bg_pos[n]} mm is too low. \
                    Actual data begins at {data.L.min()} mm.')
        else:
            bg_idx[n] = np.nonzero(data.L.to_numpy() > bg_pos[n])
            # Strip out the extra dimension that is not needed
            bg_idx[n] = bg_idx[n][0][0]

    return bg_idx


def trim_data(data, bg_idx, use_ends=False, cps_per_nA=True):

    if cps_per_nA:
        y = 'cps_per_nA'
    else:
        y = 'cps'

    data_trim = np.array([])
    Lax_trim = np.array([])

    if use_ends:
        Lax_trim = np.append(Lax_trim, data.L[:20]) # add first 20 data points
        data_trim = np.append(data_trim, data[y][:20])

    for i in range(4):  # for each background position
        Lax_trim = np.append(Lax_trim, data.L[bg_idx[i] - 10:bg_idx[i] + 10])
        data_trim = np.append(data_trim, data[y][bg_idx[i] - 10:bg_idx[i] + 10])

    if use_ends:
        Lax_trim = np.append(Lax_trim, data.L[-20:]) # add last 20 data points
        data_trim = np.append(data_trim, data[y][-20:])

    trimmed_data = pd.concat(
        [pd.Series(Lax_trim, name='L'), pd.Series(data_trim, name=y)], axis=1)

    return trimmed_data


def trim_data_from_regions(data, fit_regions): 
    '''
    Arguments
    =========
    data : pandas df : dataframe with columns L and cps_per_nA
    fit_regions : list or array : regions of the spectrum (in L units) to 
        include in the fit. Example: the array [[L1, L2], [L3, L4]] will fit
        the regions of the spectrum from L1 to L2 and L3 to L4.
    '''
    return select_roi(data, fit_regions)


def fit_bg(data, cps_per_nA=True, max_c='default'):

    '''Fit a lorentzian plus constant to the background
    Arguments
    data = pd.DataFrame with columns L, cps, cps_per_nA
    cps_per_nA = bool, if true use cps/nA in fitting
    max_c: string or numeric : specifies whether there should be a maximum 
        imposed on the "c" constant parameter in the fit. Sometimes this 
        improves the fit, sometimes not. Default ('default') is no limit. 
        'min_y' sets the  maximum to the minimum of the cps values. Otherwise 
        you could pass a numeric value in here to directly specify a maximum.
    '''
    x = data.L
    if cps_per_nA:
        y = data.cps_per_nA
    else:
        y = data.cps

    mod = LorentzianModel() + ConstantModel()

    pars = mod.make_params(amplitude=max(y) * 2, sigma=10, center=100, c=min(y) / 2)
    pars['c'].set(min=0)
    pars['amplitude'].set(min=0)
    pars['sigma'].set(min=0)
        
    if max_c == 'min_y':
        pars['c'].set(max=min(y))
    elif max_c == 'default':
        pars['c'].set(max=None)

    pars['center'].set(max=120, min=50)
    out = mod.fit(y, pars, x=x)

    return out


def plot_bg_fit(data, trimmed_data, out, comment, pk_pos_markers=False,
                save_to='default', cps_per_nA=True, print_parameters=True,
                fitted=None):

    """ fitted = either string ('default'), or an array in cases
    where you might not want the fitted curve to be calculated directly
    from the 'out' fit object.
    """

    if fitted is None:
        fitted = out.eval(x=data.L)

    fig, ax = plt.subplots(figsize=(6, 4))

    if cps_per_nA:
        y = 'cps_per_nA'
    else:
        y = 'cps'

    plt.plot(data.L, data[y], '.', color=[0.7, 0.7, 0.7], markersize=1)

    # Make a smoothing filter over the data to make it easier to visualise this
    yhat = savgol_filter(data[y], 51, 3)  # window size 51, polynomial order 3
    plt.plot(data.L, yhat, '-', color='tab:gray', linewidth=0.5)

    plt.plot(trimmed_data.L, trimmed_data[y], '.r', markersize=1.5)
    plt.plot(data.L, fitted, '-', color='xkcd:cornflower', linewidth=1.5)

    ylims = ax.get_ylim()
    xlims = ax.get_xlim()
    yrange = ylims[1] - ylims[0]
    xrange = xlims[1] - xlims[0]

    plt.ylim(-0.1, ylims[1])

    if print_parameters:

        param_names = ['amplitude', 'sigma', 'center', 'c']
        # param_names = ['amplitude', 'sigma', 'center', 'fraction', 'c']
        param_vals = [out.params[name].value for name in param_names]

        plt.text(xlims[1] - 0.35 * xrange, ylims[1] - 0.05 * yrange, 'Fit parameters: ')

        for i in range(len(param_names)):
            plt.text(xlims[1] - 0.3 * xrange, ylims[1] - (i + 2) * 0.05 * yrange,
                     '{} = {:.1f}'.format(param_names[i], param_vals[i]), ha='left',
                     va='bottom')

    if pk_pos_markers:
        fmt = ['-r', '-k']
        for i, marker in enumerate(pk_pos_markers):
            plt.plot([marker, marker], [ylims[0], ylims[1]], fmt[i],
             linewidth=0.5)

    plt.xlabel('L (mm)')
    plt.ylabel(y)

    plt.title(comment)

    plt.tight_layout()

    # ---- Save the figure to a file ----
    if save_to is None:

        print('Figure not saved')

    elif save_to == 'default':

        helper_funs.make_folder_if_it_does_not_exist(Path('./wd_scans/fits/'))

        plt.savefig(Path('./wd_scans/fits/bg_fit_' + comment + '.pdf'))
        print(f'Saved figure showing fit to {comment} in folder "./wd_scans/fits/"')

    elif isinstance(save_to, Path): # if path is specified

        plt.savefig(save_to / f'{comment}.pdf')

        print(f'Saved figure showing fit to {comment} in folder "{save_to}"')

    else:
        raise Exception(f'The save_to argument {save_to} of type {type(save_to)} is not allowed. It must be either None, "default" or a pathlib.Path object.')

    return fig, ax


# Write out the fit parameters to a text file
def write_fit_params(out, comment, save_to='default'):

    report = out.fit_report()

    par_dict = {'sigma': out.params['sigma'].value,
                 'center': out.params['center'].value}

      # ---- Save the fit parameters to a file ----
    if save_to is None:

        print('Parameters not saved')

    elif isinstance(save_to, Path) : # if path is specified
        json.dump(par_dict,
              open(save_to / f'key_params_{comment}.txt', 'w'))

        with open(save_to / f'report_{comment}.txt', 'w') as f:
            f.writelines(report)

        print('Saved fit parameters for ' + comment)
    else:
        raise Exception(f'The save_to argument {save_to} of type {type(save_to)} is not allowed. It must be either None, or a pathlib.Path object.')

    return par_dict

def calc_curve(params, i, x):

    """calculate lorentzian from parameters for data set i"""

    # Get parameters for the ith curve
    params_i = Parameters()
    params_i.add('amplitude', value=params['amplitude_%i' % (i+1)].value)
    params_i.add('center', params['center_%i' % (i+1)].value)
    params_i.add('sigma', params['sigma_%i' % (i+1)].value)
    params_i.add('c', params['c_%i' % (i+1)].value)

    # Use these parameters to calculate a curve
    mod = LorentzianModel() + ConstantModel()
    calculated_curve = mod.eval(x=x, params=params_i)

    return calculated_curve

def objective(params, x, data):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by Lorentzian functions"""
    ndata, nx = data.shape
    resid = 0.0*data[:]
    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - calc_curve(params, i, x)
    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

def fit_scans_together(data, fit_regions, path_out,
                       max_c='default'):

    """ Fit a set of WD scans to a single background lorentzian
    with comman sigma & centre parameters, varying only amplitude
    and a constant value

    Arguments
    =========
    data : pandas df : dataframe with columns L and cps_per_nA
    fit_regions : list or array : regions of the spectrum (in L units) to 
        include in the fit. Example: the array [[L1, L2], [L3, L4]] will fit
        the regions of the spectrum from L1 to L2 and L3 to L4.
    path_out: pathlib.Path : folder in which to save fit results
    max_c: string or numeric : specifies whether there should be a maximum 
        imposed on the "c" constant parameter in the fit. Sometimes this 
        improves the fit, sometimes not. Default ('default') is no limit. 
        'min_y' sets the  maximum to the minimum of the cps values. Otherwise 
        you could pass a numeric value in here to directly specify a maximum.
    """

    # Select the data for fitting and store in a new list
    trimmed_data = []
    for i in range(len(data)):
        trimmed_data.append(select_roi(data[i], fit_regions))

    # The L axis is the same across all datasets, so:
    L = trimmed_data[0].L.to_numpy()
    trimmed_array = np.array([d.cps_per_nA.to_numpy()
                              for d in trimmed_data])

    # So the new dataset is trimmed_array, with each row is a
    # dataset and each column corresponds to an "energy" in L units

    # create sets of parameters, one per data set
    fit_params = Parameters()

    for iy, y in enumerate(trimmed_array):

        if max_c == 'min_y':
            max_c = min(trimmed_array[iy])
        elif max_c == 'default':
            max_c = None

        fit_params.add( 'amplitude_%i' % (iy+1), value=max(trimmed_array[iy])*2, min=0.0)
        fit_params.add( 'center_%i' % (iy+1), value=100, min=0,  max=120)
        fit_params.add( 'sigma_%i' % (iy+1), value=0.1, min=0, max=100)
        fit_params.add( 'c_%i' % (iy+1), value=0, min=0, max=max_c)

    # but now constrain all values of sigma & center to have the same value
    # by assigning sig_2, sig_3, .. sig_5 to be equal to sig_1
    for iy in range(2,len(data) + 1):
        fit_params['sigma_%i' % iy].expr='sigma_1'
        fit_params['center_%i' % iy].expr='center_1'

    # Perform the fit
    print("Performing the fit...")
    result = minimize(objective, fit_params, args=(L, trimmed_array))

    with open(path_out / 'fit_report.txt', 'w') as f:

        f.writelines(fit_report(result))

    # Save parameters
    param_names = ['sigma', 'center']
    # Just take the first fit for the parameters
    param_vals = [result.params[name + '_1'].value for name in param_names]

    json.dump(dict(zip(param_names, param_vals)),
              open(path_out / 'key_params.txt', 'w'))
    print('Saved fit parameters to {}'.format(path_out))

    return result, trimmed_data

def plot_fits_together(data, trimmed_data, result, comments,
                       path_out=Path('./wd_scans/fits/')):

    helper_funs.make_folder_if_it_does_not_exist(path_out)

    plt.figure()

    Lrange = np.arange(120, 180, 0.1)

    for i in range(len(data)):
        plt.plot(trimmed_data[i].L, trimmed_data[i].cps_per_nA, 'x',
                 markersize=2)

    for i in range(len(data)):
        y_fit = calc_curve(result.params, i, Lrange)
        plt.plot(Lrange, y_fit, '-', linewidth=0.5)


    plt.savefig(path_out / 'all_fits.pdf')
    plt.close()

    for i in range(len(data)):

        fitted = calc_curve(result.params, i, data[i].L)
        plot_bg_fit(data[i], trimmed_data[i], result, comments[i],
                    pk_pos_markers=False, save_to=path_out,
                    cps_per_nA=True, print_parameters=False,
                    fitted=fitted)

