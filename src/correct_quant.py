""" Functions used to process Jeol EPMA data and write CalcZAF files
"""

# %% Import dependencies and set plot defaults
import pickle
import json
from lmfit.models import PseudoVoigtModel, ConstantModel, LorentzianModel
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm  # For monte carlo number generation
from src import readfiles
from pathlib import Path
import os
import pdb
import copy

plt.rcParams['figure.figsize'] = 6, 4
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12

pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)

# %% Classes

class Spot:
    def __init__(self):

        self.info = None
        self.bg = None  # bg = background
        self.peak = None
        self.standard = None
        self.corrected = None
        self.montecarlo = None
        self.wd_scan_params = None

    def add_data(self, info, bg, peak, standard):

        """Add data to spot object from pandas series/dataframes stored
        in working memory (i.e. generated by another function).
        Alternatively, use add_data_from_files()."""

        self.info = info
        self.bg = bg  # bg = background
        self.peak = peak
        self.standard = standard
        self.idx_N = [idx for idx in self.peak.index.tolist()
                  if self.peak.element[idx] == "N"]

    def add_data_from_files(self, path):

        """Add data to spot object from files stored in the folder (path)"""

        self.info = pd.read_csv(path / 'info.csv', index_col=0, squeeze=True)

        # Convert strings to floats where necessary:
        self.info.apf = float(self.info.apf)
        self.info.apf_sd = float(self.info.apf_sd)
        self.info.nA = float(self.info.nA)
        self.info.spot_size_um = float(self.info.spot_size_um)

        self.bg = pd.read_csv(path / 'bg.csv', index_col=0)
        self.peak = pd.read_csv(path / 'peak.csv', index_col=0)
        self.standard = pd.read_csv(path / 'standard.csv', index_col=0)
        self.idx_N = [idx for idx in self.peak.index.tolist()
                  if self.peak.element[idx] == "N"]

    def add_wd_scan_params_from_file(self, path):

        """ read a text file containing wd scan parameters
        sigma and position, and add these as attribute of the spot"""
        print('\n loading path: {}'.format(path))
        self.wd_scan_params = json.load(open(path))
        print(self.wd_scan_params)

    def comprehensify_data(self):

        # --------- Compute bg positions and store them ---------
        self.bg['lwr_pos'] = self.peak['pos'] - self.bg['lwr_pos_rel']
        self.bg['upr_pos'] = self.peak['pos'] + self.bg['upr_pos_rel']

        # --------- 'Uncorrect' background, i.e. find raw cps at peak by
        # undoing the bg subtraction ---------
        bg_cps_at_pk_pos = [None] * len(self.peak.index)

        for i in self.peak.index:
            bg_cps_at_pk_pos[i] = lin_bg(
                self.peak['pos'][i],
                [self.bg['lwr_cps'][i], self.bg['upr_cps'][i]],
                [self.bg['lwr_pos'][i], self.bg['upr_pos'][i]])

        self.bg['cps_at_pk_pos'] = bg_cps_at_pk_pos
        self.peak['raw_cps'] = self.peak['net_cps'] + bg_cps_at_pk_pos

        # Calculate the relative standard deviation on kraw ----------
        self.peak['kraw_stdev_pcnt'] = np.sqrt(
            (self.peak.stdev_net_cps) ** 2 +
            (self.standard.stdev_net_cps) ** 2)

        # Calculate the relative standard deviation on the raw cps at peak ---
        self.peak['stdev_raw_cps'] = (np.sqrt(
            (self.peak['raw_cps'] * self.peak['time']))
                                      / self.peak['time'])

        self.bg['stdev_lwr'] = (np.sqrt(
            (self.bg['lwr_cps'] * self.bg['time']))
                                / self.bg['time'])

        self.bg['stdev_upr'] = (np.sqrt(
            (self.bg['upr_cps'] * self.bg['time']))
                                / self.bg['time'])

    def correct_bg(self, fit_parameter_file, path_out):

        """ Takes parameters of a lorentzian curve from fit_parameter_file and
        fits a new lorentzian curve to the four background measurements.
        The only variable parameter in the fit is the amplitude and the
        constant (c) so the 'shape' of the background curve remains the same.

        The new background curve is used to determine the 'true' background cps
        value, and hence the 'true' net_cps value. This new net_cps value is
        put into a copy of the .peak dataframe, stored now as
        .corrected and with the 'dummy' nitrogen line removed.

        fit_parameter_file is a file containing a dictionary of fit
        parameters

        """
        idx = self.idx_N[0]  # index of the element to be corrected
        idx_dupl = self.idx_N[1]  # index of the additional bg measurements

        # Get background positions for N as a numpy array
        bg_pos = compute_bg_positions(self)

        # Get background cps values as a numpy array
        lwr_bg_cps = np.array(self.bg.lwr_cps[idx:idx_dupl + 1])
        upr_bg_cps = np.array(self.bg.upr_cps[idx:idx_dupl + 1])
        bg_cps = np.concatenate((lwr_bg_cps, upr_bg_cps), axis=0)

        # Get the raw cps at peak by adding a linear background to the net_cps
        # obsolete   # pk_cps_raw = pk_cps_net + lin_bg(pk_pos, bg_cps[[0, 2]],
        # bg_pos[[0,
        # 2]])

        # Fit the background according to its shape
        modelout = fit_quant_bg(bg_pos, bg_cps, self.wd_scan_params)

        # Find the corrected background cps
        corrected_bg = modelout.eval(x=self.peak.pos[idx])
        pk_cps_corrected = self.peak.raw_cps[idx] - corrected_bg

        sample_name = self.info.date + '_' + self.info.comment

        plot_background_correction(
            modelout, self.peak.pos[idx], self.peak.raw_cps[idx], bg_pos,
            bg_cps, sample_name, save_path=Path(path_out, sample_name),
            annotate_cps=True,
            figure_axis=None)

        # Remove the 'dummy' analysis line
        self.corrected = self.peak.copy(deep=True).drop(
            idx_dupl, inplace=False)
        self.corrected.loc[idx, "net_cps"] = pk_cps_corrected

        # Update kraw and stdev for nitrogen
        kraw_pcnt, _ = calc_kraw_stdev(
            pk_cps_corrected, self.standard.net_cps[idx], 10000,
            self.standard.stdev_net_cps[idx], self.info.nA,
            self.standard.nA[idx])
        # Note stdev is not calculated here... (a dummy number of 10000 is used)
        # because there is no estimate of stdev on pk_cps_corrected, without
        # doing a monte-carlo simulation.

        self.corrected.loc[idx, 'kraw_pcnt'] = kraw_pcnt

    def resample_cps(self, num_mc_sims=200):

        """ Apply a mc simulation to all element background corrections
         to prove that it gives similar stdev results to the calculated stdev

         Returns:
            A new attribute, self.montecarlo (dataframe)

         """
        num_els = len(self.peak.index)

        net_cps_mean = [None] * num_els
        net_cps_stdev_pcnt = [None] * num_els
        bg_cps_mean = [None] * num_els
        bg_cps_stdev_pcnt = [None] * num_els
        kraw_pcnt = [None] * num_els
        kraw_stdev_pcnt = [None] * num_els

        # For each element in the analysis
        for j in range(num_els):

            # -------------Resample the cps data -------------
            # Get the analysed values for peak and bg and their stdevs
            analysed_values = [self.peak.raw_cps[j],
                               self.bg.lwr_cps[j],
                               self.bg.upr_cps[j]]

            analysed_stdevs = [self.peak.stdev_raw_cps[j],
                               self.bg.stdev_lwr[j],
                               self.bg.stdev_upr[j]]

            cps_resampled = np.empty(shape=(num_mc_sims, 3))

            # Resample the analysed values based on the stdevs
            for i in range(3):
                cps_resampled[:, i] = norm.rvs(loc=analysed_values[i],
                                               scale=analysed_stdevs[i],
                                               size=num_mc_sims,
                                               random_state=None)

            # Apply a linear background correction to each mc simulation

            bg_cps_resampled = np.empty(shape=(num_mc_sims, 1))

            for i in range(num_mc_sims):
                bg_cps_resampled[i] = lin_bg(
                    self.peak.pos[j], cps_resampled[i, [1, 2]],
                    [self.bg.lwr_pos[j], self.bg.upr_pos[j]])

            net_cps_resampled = np.reshape(
                cps_resampled[:, 0], (num_mc_sims, 1)) - bg_cps_resampled

            # Calculate means and stdevs the monte-carlo simulation
            net_cps_mean[j] = np.mean(net_cps_resampled)
            net_cps_stdev = np.std(net_cps_resampled)
            net_cps_stdev_pcnt[j] = net_cps_stdev / net_cps_mean[j] * 100

            bg_cps_mean[j] = np.mean(bg_cps_resampled)
            bg_cps_stdev = np.std(bg_cps_resampled)
            bg_cps_stdev_pcnt[j] = bg_cps_stdev / bg_cps_mean[j] * 100

            # Propagate the errors into the k_raw
            kraw_pcnt[j], kraw_stdev_pcnt[j] = calc_kraw_stdev(
                net_cps_mean[j],
                self.standard.net_cps[j],
                net_cps_stdev_pcnt[j],
                self.standard.stdev_net_cps[j],
                self.info.nA, self.standard.nA[j])

        # Store the results in a data frame
        self.montecarlo = pd.DataFrame({
            'element': self.peak.element,
            'net_cps_mc_mean': net_cps_mean,
            'net_cps_mc_stdev_pcnt': net_cps_stdev_pcnt,
            'bg_cps_mc_mean': bg_cps_mean,
            'bg_cps_mc_stdev_pcnt': bg_cps_stdev_pcnt,
            'kraw_pcnt': kraw_pcnt,
            'kraw_stdev_pcnt': kraw_stdev_pcnt})

    def correct_bg_mc(self, fit_parameter_file, path_out, sample_name,
                      num_mc_sims=200):

        """ Use a monte-carlo simulation to correct for the curved background

        Arguments:
            fit_parameter_file: file containing dictionary of fit parameters

        Returns:
            - A plot showing the monte carlo background correction results
            - A new attribute to the Spot object, corrected_mc
        """

        idx = self.idx_N[0]
        idx_dupl = self.idx_N[1]

        # Create an array with rows being a mc-simulation of the peak cps
        pk_cps_raw = norm.rvs(loc=self.peak.raw_cps[idx],
                              scale=self.peak.stdev_raw_cps[idx],
                              size=num_mc_sims, random_state=None)

        # Get background positions based on relative position from peak
        bg_pos = compute_bg_positions(self)

        # Get background cps values and stdevs
        lwr_bg_cps = np.array(self.bg.lwr_cps[idx:idx_dupl + 1])
        upr_bg_cps = np.array(self.bg.upr_cps[idx:idx_dupl + 1])
        bg_cps_means = np.concatenate((lwr_bg_cps, upr_bg_cps), axis=0)

        lwr_bg_stdev = np.array(self.bg.stdev_lwr[idx:idx_dupl + 1])
        upr_bg_stdev = np.array(self.bg.stdev_upr[idx:idx_dupl + 1])
        bg_cps_stdev = np.concatenate((lwr_bg_stdev, upr_bg_stdev), axis=0)

        # Create a matrix with rows being mc-simulations of the background
        bg_cps = np.empty(shape=(num_mc_sims, 4))

        for i in range(4):
            bg_cps[:, i] = norm.rvs(loc=bg_cps_means[i],
                                    scale=bg_cps_stdev[i],
                                    size=num_mc_sims, random_state=None)
            # Note: In the norm.rvs function, the keywords:
            # loc = midpoint of normal distribution; scale = standard deviation.


        if type(fit_parameter_file) == str: # Convert to path object
            fit_parameter_file = Path(fit_parameter_file)

        pk_cps_net_corrected = np.empty(shape=(num_mc_sims, 1))
        bg_cps_corrected = np.empty(shape=(num_mc_sims, 1))

        # Perform the background fits and plot a figure at the same time
        plt.figure()
        ax = plt.axes()

        fit_params_mc = [None] * num_mc_sims


        for i in range(num_mc_sims):
            # Fit the background according to its shape

            # Print a loop counter
            if i % 50 == 0: # If i is divisible by 50
                print('monte-carlo loop {} of {}'.format(i, num_mc_sims))

            modelout = fit_quant_bg(bg_pos, bg_cps[i, :], self.wd_scan_params)

            fit_params_mc[i] = [modelout.params['amplitude'].value,
                                modelout.params['c'].value]

            # Find the corrected background cps
            bg_cps_corrected[i] = modelout.eval(x=self.peak.pos[idx])
            pk_cps_net_corrected[i] = pk_cps_raw[i] - bg_cps_corrected[i]

            # Add this to the figure
            xlims, ylims = plot_background_correction(
                modelout, self.peak.pos[idx], pk_cps_raw[i], bg_pos,
                bg_cps[i, :], self.info.comment, save_path=None,
                annotate_cps=False,
                figure_axis=ax)

        mean_from_mc = [np.mean(pk_cps_net_corrected), np.mean(bg_cps_corrected)]
        stdev_pcnt_from_mc = [np.std(pk_cps_net_corrected) * 100 / mean_from_mc[0],
                              np.std(bg_cps_corrected) * 100 / mean_from_mc[1]]

        plt.text(120, ylims[0] + 0.01 * (ylims[1] - ylims[0]),
                 'net cps at peak from mc sims = {:.1f} ± {:.1f}%'.format(
                     mean_from_mc[0], stdev_pcnt_from_mc[0]), fontsize=8)

        plt.text(120, ylims[0] + 0.10 * (ylims[1] - ylims[0]),
                 'net cps at peak from curved bg {:.1f}'.format(
                     self.corrected.net_cps[idx]), fontsize=8)

        plt.text(120, ylims[0] + 0.05 * (ylims[1] - ylims[0]),
                 'net cps at peak from linear bg {:.1f} ± {:.1f}%'.format(
                     self.peak.net_cps[idx], self.peak.stdev_net_cps[idx]),
                 fontsize=8)

        # Annotate figure with fit parameters
        param_names = ['amplitude', 'c']
        param_vals = np.mean(fit_params_mc, axis=0)
        param_stdevs = np.std(fit_params_mc, axis=0)

        plt.text(xlims[1] - 0.45 * (xlims[1] - xlims[0]),
                 ylims[1] - 0.05 * (ylims[1] - ylims[0]),
                 'Fit parameters: ')

        for i in range(len(param_names)):
            plt.text(xlims[1] - 0.45 * (xlims[1] - xlims[0]), ylims[1] - (i + 2)
                     * 0.05 * (ylims[1] - ylims[0]),
                     '{} = {:.1f} ± {:.1f}'.format(
                         param_names[i], param_vals[i], param_stdevs[i]),
                     ha='left', va='bottom')

        filename = self.info.date + '_' + self.info.comment

        plt.savefig(Path(path_out, 'bg_mc_' + filename + '.png'), dpi=600)
        plt.close()
        print('Saved montecarlo bg correction figure for ' + filename)

        # Store the corrected value in the montecarlo dataframe -----

        # remove the 'dummy' analysis line

        self.montecarlo.drop(idx_dupl, inplace=True)

        # Update the nitrogen value and stdev
        self.montecarlo.loc[idx, 'net_cps_mc_mean'] = mean_from_mc[0]
        self.montecarlo.loc[idx, 'net_cps_mc_stdev_pcnt'] = stdev_pcnt_from_mc[0]
        self.montecarlo.loc[idx, 'bg_cps_mc_mean'] = mean_from_mc[1]
        self.montecarlo.loc[idx, 'bg_cps_mc_stdev_pcnt'] = stdev_pcnt_from_mc[1]

        # Update kraw and stdev for nitrogen
        kraw_pcnt, kraw_stdev_pcnt = calc_kraw_stdev(
            mean_from_mc[0], self.standard.net_cps[idx], stdev_pcnt_from_mc[0],
            self.standard.stdev_net_cps[idx], self.info.nA,
            self.standard.nA[idx])

        self.montecarlo.loc[idx, 'kraw_pcnt'] = kraw_pcnt
        self.montecarlo.loc[idx, 'kraw_stdev_pcnt'] = kraw_stdev_pcnt


    def save_attributes(
            self, 
            filename:str # File path for excel file
            ): 

        '''Save all relevant dataframes to sheets of excel file'''

        with pd.ExcelWriter(Path(filename),
                            engine='xlsxwriter') as writer:
            self.peak.to_excel(writer, sheet_name='peak')
            self.bg.to_excel(writer, sheet_name='background')
            self.standard.to_excel(writer, sheet_name='standard')
            self.corrected.to_excel(writer, sheet_name='corrected')
            self.montecarlo.to_excel(writer, sheet_name='montecarlo')

        print('Saved attributes to excel')

    def correct_ht_area_ratio(self):

        apf = self.info.apf
        apf_sd = self.info.apf_sd

        apf_sd_pcnt = apf_sd/apf * 100

        idx = self.idx_N[0]
        # Make a new column in the montecarlo df
        self.montecarlo["kraw_apf_pcnt"] = self.montecarlo.loc[:, "kraw_pcnt"]
        # Change N entry to corrected value
        self.montecarlo.loc[idx, "kraw_apf_pcnt"] = self.montecarlo.loc[
                                                        idx, "kraw_pcnt"] / apf
        # Add standard devation column
        self.montecarlo['kraw_stdev_apf_pcnt'] = self.montecarlo.loc[:,
                                                 "kraw_stdev_pcnt"]
        # Change N entry
        self.montecarlo.loc[idx, "kraw_stdev_apf_pcnt"] = np.sqrt(
            self.montecarlo.loc[idx, 'kraw_stdev_pcnt']**2 + apf_sd_pcnt**2)

        # Make a new column in the corrected df
        self.corrected["kraw_apf_pcnt"] = self.montecarlo.loc[:, "kraw_pcnt"]
        self.corrected.loc[idx, "kraw_apf_pcnt"] = self.montecarlo.loc[
                                                       idx, "kraw_pcnt"] / apf
        # Add standard devation column
        self.corrected['kraw_stdev_apf_pcnt'] = self.corrected.loc[:,
                                                 "kraw_stdev_pcnt"]
        # Change N entry
        self.corrected.loc[idx, "kraw_stdev_apf_pcnt"] = np.sqrt(
            self.corrected.loc[idx, 'kraw_stdev_pcnt']**2 + apf_sd_pcnt**2)

        # Also do this for the non-bg corrected data

        self.peak["kraw_apf_pcnt"] = self.peak.loc[:, "kraw_pcnt"].copy(
            deep=True)
        self.peak.loc[idx, "kraw_apf_pcnt"] = self.peak.loc[
                                                  idx, "kraw_pcnt"] / apf

        # Add standard devation column
        self.peak['kraw_stdev_apf_pcnt'] = self.peak.loc[:,
                                                 "kraw_stdev_pcnt"]
        # Change N entry
        self.peak.loc[idx, "kraw_stdev_apf_pcnt"] = np.sqrt(
            self.peak.loc[idx, 'kraw_stdev_pcnt']**2 + apf_sd_pcnt**2)

# %%  Functions to correct background ---------------------------

def compute_bg_positions(Spot):
    idx = Spot.idx_N[0]
    idx_dupl = Spot.idx_N[1]

    lwr_bg = np.array(Spot.bg.lwr_pos[idx:idx_dupl + 1])
    upr_bg = np.array(Spot.bg.upr_pos[idx:idx_dupl + 1])
    bg_pos = np.concatenate((lwr_bg, upr_bg), axis=0)

    return bg_pos

def fit_quant_bg(datax, datay, par_dict, guess_params=None):

    x = datax
    y = datay

    if guess_params is None:
        guess_params = {'amplitude': max(y) * 2,
                        'sigma': par_dict['sigma'],
                        'center': par_dict['center'],
                        'c': min(y) / 2}

    mod = LorentzianModel() + ConstantModel()

    pars = mod.make_params(amplitude=guess_params['amplitude'],
                           sigma=guess_params['sigma'],
                           center=guess_params['center'],
                           c=guess_params['c'])

    ## NOTE the last line used to be c=par_dict['c'].
    # but now I'm not saving the c parameter in the par_dict file.
    # so changed to c=0. But this has potential to cause problems
    # or differences from previous fits. >> yes, it caused problems,
    # so I changed it to min(y)/2 and added option to pass in guess
    # parameters.

    pars['center'].set(vary=False)
    pars['sigma'].set(vary=False)
    pars['c'].set(min=0, max=min(y))
    # pars['c'].set(vary=False)

    fitted_curve = mod.fit(y, pars, x=x)

    return fitted_curve


def plot_background_correction(
        out, 
        pk_pos, 
        pk_cps,
        bg_pos, 
        bg_cps, 
        sample_name,
        save_path=None,
        annotate_cps=True,
        figure_axis=None,
        ):

    newx = np.arange(120.0, 180.0, 0.2)
    corrected_bg = out.eval(x=pk_pos)
    fitted = out.eval(x=newx)

    if figure_axis == None:
        plt.figure()
        ax = plt.axes()
    else:
        ax = figure_axis

     # Add the measured data
    plt.plot(bg_pos, bg_cps, 'ok')
    plt.plot(newx, fitted, '-b', linewidth=0.5)
    plt.plot(pk_pos, pk_cps, 'xk')
    plt.plot(newx, lin_bg(newx, bg_cps[[0, 2]], bg_pos[[0, 2]]), '-r',
             linewidth=0.5)

    ylims = ax.get_ylim()
    plt.ylim(0, max(np.append(bg_cps, pk_cps)) +
             0.1 * max(np.append(bg_cps, pk_cps)))

    ylims = ax.get_ylim()
    xlims = ax.get_xlim()
    yrange = ylims[1] - ylims[0]
    xrange = xlims[1] - xlims[0]

    if annotate_cps == True:
        # Annotate figure with fit parameters
        param_names = ['amplitude', 'sigma', 'center', 'c']
        # param_names = ['amplitude', 'sigma', 'center', 'fraction', 'c']
        param_vals = [out.params[name].value for name in param_names]

        plt.text(xlims[1] - 0.35 * xrange, ylims[1] - 0.05 * yrange,
                 'Fit parameters: ')

        for i in range(len(param_names)):
            plt.text(xlims[1] - 0.35 * xrange,
                     ylims[1] - (i + 2) * 0.05 * yrange,
                     '{} = {:.1f}'.format(param_names[i], param_vals[i]),
                     ha='left', va='bottom')

        plt.text(120, ylims[0] + 0.01 * yrange,
                 'net cps at peak (curved) = {:.1f}'.format(
                     pk_cps - corrected_bg))
        plt.text(120, ylims[0] + 0.05 * yrange,
                 'net cps at peak (linear) = {:.1f}'.format(
                     pk_cps - lin_bg(pk_pos, bg_cps[[0, 2]], bg_pos[[0, 2]])))

    plt.xlabel('L (mm)')
    plt.ylabel('cps')

    plt.title(sample_name)

    if save_path is not None:
        folder_exists = os.path.isdir(Path(save_path).parent)
        if not folder_exists:
            print("Can't save because the folder doesn't exist")
        plt.savefig(Path(save_path))
        plt.close()
        print('Saved bg correction figure for ' + sample_name)

    return xlims, ylims


def lin_bg(x, yvals, xvals):
    slope = (yvals[1] - yvals[0]) / (xvals[1] - xvals[0])
    intercept = yvals[0] - slope * xvals[0]
    y = slope * x + intercept
    return y


# %%  Functions for stdev --------------------------

def check_mc_behaviour(myspot):
    # Check M.C. simulation behaves as expected
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(0.39 * 15, 0.39 * 15))
    axs = axs.ravel()

    plt.sca(axs[0])
    x = myspot.peak.net_cps[:-1]
    y = myspot.montecarlo.net_cps_mc_mean
    plt.plot(x, y, '.k', linewidth=0.5)
    plt.plot(x[myspot.idx_N[0]], y[myspot.idx_N[0]], '.r')
    for i in range(len(x)):
        plt.text(x[i], y[i] + 0.5 * y[i], myspot.peak.element[i],
                 ha='right', va='bottom', fontsize=8)

    plt.title('net cps')

    plt.sca(axs[1])
    x = myspot.peak.stdev_net_cps[:-1]
    y = myspot.montecarlo.net_cps_mc_stdev_pcnt
    plt.plot(x, y, '.k', linewidth=0.5)
    plt.plot(x[myspot.idx_N[0]], y[myspot.idx_N[0]], '.r')
    for i in range(len(x)):
        plt.text(x[i], y[i] + 0.5 * y[i], myspot.peak.element[i],
                 ha='right', va='bottom', fontsize=8)
    plt.title('stdev% net cps')

    plt.sca(axs[2])
    x = myspot.peak.kraw_pcnt[:-1]
    y = myspot.montecarlo.kraw_pcnt
    plt.plot(x, y, '.k', linewidth=0.5)
    plt.plot(x[myspot.idx_N[0]], y[myspot.idx_N[0]], '.r')
    for i in range(len(x)):
        plt.text(x[i], y[i] + 0.5 * y[i], myspot.peak.element[i],
                 ha='right', va='bottom', fontsize=8)
    plt.title('kraw')

    plt.sca(axs[3])
    x = myspot.peak.kraw_stdev_pcnt[:-1]
    y = myspot.montecarlo.kraw_stdev_pcnt
    plt.plot(x, y, '.k', linewidth=0.5)
    plt.plot(x[myspot.idx_N[0]], y[myspot.idx_N[0]], '.r')
    for i in range(len(x)):
        plt.text(x[i], y[i] + 0.5 * y[i], myspot.peak.element[i],
                 ha='right', va='bottom', fontsize=8)
    plt.title('stdev% kraw')

    for ax in axs:
        plt.sca(ax)
        plt.xlabel('original value', fontsize=10)
        plt.ylabel('monte carlo', fontsize=10)
        plt.yscale('log')
        plt.xscale('log')
        xlims = ax.get_xlim()
        plt.plot(xlims, xlims, '-b', linewidth=0.5)

    plt.savefig(Path('./mc_checks/' + myspot.info.date + '_' + myspot.info.comment + '.pdf'))
    plt.close()


def calc_kraw_stdev(unk_cps_net, std_cps_net, unk_cps_stdev_pcnt,
                    std_cps_stdev_pcnt, unk_nA, std_nA):

    kraw_pcnt = 100 * ((unk_cps_net / unk_nA) / (std_cps_net / std_nA))

    kraw_stdev_pcnt = np.sqrt(unk_cps_stdev_pcnt ** 2 + std_cps_stdev_pcnt ** 2)
    # Note: because we are using only relative standard deviations here,
    # the nA does not get propagated in the stdev calculation.

    return kraw_pcnt, kraw_stdev_pcnt



def process_datasets(myspot, datalist, num_mc_sims, path_out):

    for i in range(len(myspot)):

        print('\nProcessing dataset:', i + 1, 'of', len(datalist.folder), ':',
              myspot[i].info.comment)

        print('Correcting background') # -----------------
        myspot[i].correct_bg(datalist.paramfile[i], path_out=path_out)

        print('Resample cps to check stdev method') # -----------------
        myspot[i].resample_cps(num_mc_sims=100)

        print('Montecarlo background correction') # -----------------
        myspot[i].correct_bg_mc(
            datalist.paramfile[i],
            path_out=path_out,
            num_mc_sims=num_mc_sims,
            sample_name=myspot[i].info.comment)

        # print('check montecarlo behaviour') # -----------------
        # check_mc_behaviour(myspot[i])

        print('correct height/area ratio') # -----------------
        myspot[i].correct_ht_area_ratio()

        print("original kraw: {:.2f} ± {:.2f}%".format(
            myspot[i].peak.kraw_pcnt[myspot[i].idx_N[0]],
            myspot[i].peak.kraw_stdev_pcnt[myspot[i].idx_N[0]]))
        print("corrected kraw: {:.2f} ± {:.2f}%".format(
            myspot[i].montecarlo.kraw_apf_pcnt[myspot[i].idx_N[0]],
            myspot[i].montecarlo.kraw_stdev_apf_pcnt[myspot[i].idx_N[0]]))

    return myspot


def write_summary_excel_tables(
        spot_list:list, # List of Spot objects
        filename:str, # File name to write to (should end in .xlsx)
        ):
    '''
    Write excel file for a list of related spots (each spot will be
    one sheet in the file)
    '''
    summary = [None] * len(spot_list[0].corrected.element)
    columns = ['original.kraw_pcnt', 'original.kraw_stdev_pcnt',
               'original.kraw_apf_pcnt', 'original.kraw_stdev_apf_pcnt',
               'corrected.kraw_pcnt', 'corrected.kraw_stdev_pcnt',
               'corrected.kraw_apf_pcnt', 'corrected.kraw_stdev_apf_pcnt',
               'montecarlo.kraw_pcnt', 'montecarlo.kraw_stdev_pcnt',
               'montecarlo.kraw_apf_pcnt', 'montecarlo.kraw_stdev_apf_pcnt']

    for i, el in enumerate(spot_list[0].corrected.element):
        summaries = np.zeros(shape=(len(spot_list), len(columns)))

        for j, spot in enumerate(spot_list):

            summaries[j, 0] = spot.peak.loc[i, 'kraw_pcnt'].round(2)
            summaries[j, 1] = spot.peak.loc[i, 'kraw_stdev_pcnt'].round(2)
            summaries[j, 2] = spot.peak.loc[i, 'kraw_apf_pcnt'].round(2)
            summaries[j, 3] = spot.peak.loc[i, 'kraw_stdev_apf_pcnt'].round(2)
            summaries[j, 4] = spot.corrected.loc[i, 'kraw_pcnt'].round(2)
            summaries[j, 5] = spot.corrected.loc[i, 'kraw_stdev_pcnt'].round(2)
            summaries[j, 6] = spot.corrected.loc[i, 'kraw_apf_pcnt'].round(2)
            summaries[j, 7] = spot.corrected.loc[i, 'kraw_stdev_apf_pcnt'].round(2)
            summaries[j, 8] = spot.montecarlo.loc[i, 'kraw_pcnt'].round(2)
            summaries[j, 9] = spot.montecarlo.loc[i, 'kraw_stdev_pcnt'].round(2)
            summaries[j, 10] = spot.montecarlo.loc[i, 'kraw_apf_pcnt'].round(2)
            summaries[j, 11] = spot.montecarlo.loc[i, 'kraw_stdev_apf_pcnt'].round(2)

        summary[i] = pd.DataFrame(data=summaries, columns=columns)
        summary[i]['comment'] = [spot.info.comment for spot in spot_list]

    with pd.ExcelWriter(Path(filename),
                        engine='xlsxwriter') as writer:
        for i, el in enumerate(spot_list[0].corrected.element):
            summary[i].to_excel(writer, sheet_name=el)

    return summary

def summarise_net_cps(myspot):
    """ Produce a summary of the net_cps original and after correction with
    mc simulation. Similar to write_summary_excel_tables() but doesn't write
    out to excel and summarises net_cps instead of kraw details.

    Arguments:
        myspot: a list of spot objects after corrections have been applied
    Returns:
        summary: a pandas dataframe containing info about net_cps before and
                after corrections for each spot.
    """

    columns = ['original.net_cps', 'original.stdev_net_cps',
               'montecarlo.net_cps_mc_mean', 'montecarlo.net_cps_mc_stdev_pcnt']

    for i, el in enumerate(myspot[0].corrected.element):

        if el == 'N':
            summaries = np.zeros(shape=(len(myspot), len(columns)))

            for j, spot in enumerate(myspot):
                summaries[j, 0] = spot.peak.loc[i, 'net_cps'].round(2)
                summaries[j, 1] = spot.peak.loc[i, 'stdev_net_cps'].round(2)
                summaries[j, 2] = spot.montecarlo.loc[i, 'net_cps_mc_mean'].round(2)
                summaries[j, 3] = spot.montecarlo.loc[i, 'net_cps_mc_stdev_pcnt'].round(2)

            summary = pd.DataFrame(data=summaries, columns=columns)
            summary['comment'] = [spot.info.comment for spot in myspot]

    return summary


def create_detection_limit_spot(spot):
    """ Create a new spot object for detection limit analysis.

    This function takes data stored in a spot object, calculates three times
    the background standard deviation, and uses this to calculate a new
    kraw detection limit, which once run through Calczaf will give a detection
    limit in wt%.
    
    Arguments
    =========
    spot : corrrect_quant.Spot : Spot object that has already had data added.

    Returns
    =======
    detlim_spot : correct_quant.Spot : Copy of the original Spot object with
            the only change being that the spot.montecarlo.kraw_pcnt for N is
            now a value based on 3 times the standard deviation of the
            background.
    """

    # Calculate net cps for detection limit
    old_net_cps = spot.montecarlo.net_cps_mc_mean[spot.idx_N[0]]
    bg_stdev_pcnt = spot.montecarlo.bg_cps_mc_stdev_pcnt[spot.idx_N[0]]
    bg_stdev_absolute = (bg_stdev_pcnt/100 
                        * spot.montecarlo.bg_cps_mc_mean[spot.idx_N[0]])
    new_net_cps = 3 * bg_stdev_absolute

    print('calculating a new k-raw')
    print('instead of net_cps = {}, using 3sd on bg = {}'.format(
        old_net_cps, new_net_cps))

    # Calculate kraw
    kraw, _ = calc_kraw_stdev(
        unk_cps_net=new_net_cps,
        std_cps_net=spot.standard.loc[spot.idx_N[0], 'net_cps'],
        unk_cps_stdev_pcnt=1000,
        std_cps_stdev_pcnt=spot.standard.loc[spot.idx_N[0], 'stdev_net_cps'],
        unk_nA=spot.info.nA,
        std_nA=spot.standard.loc[spot.idx_N[0], 'nA'])

    print("old kraw:", spot.montecarlo.loc[spot.idx_N[0], 'kraw_pcnt'])

    detlim_spot = copy.copy(spot)

    # Update value
    detlim_spot.montecarlo.loc[detlim_spot.idx_N[0], 'kraw_pcnt'] = kraw  
    print("new kraw:", detlim_spot.montecarlo.loc[spot.idx_N[0], 'kraw_pcnt'])

    # Propogate through apf correction
    detlim_spot.correct_ht_area_ratio() 

    return detlim_spot

#%% Overall functions

def make_spots(datalist, bgi=False):
    """ Sets up a list of Spot objects and reads data into these.
    Also runs the 'comprehensify_data' function. """

    myspot = [None] * len(datalist.folder)

    for i in range(len(datalist.folder)):

        peak, bg, standard, info = readfiles.read_and_organise_data(datalist.loc[i,:].copy(), bgi=bgi)
        myspot[i] = Spot(info, bg, peak, standard)
        print('Read dataset:', i + 1, 'of', len(datalist), ':',
              myspot[i].info.comment)

        myspot[i].comprehensify_data()

    return myspot