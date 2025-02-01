import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib.pyplot as plt
from src.wdscan import select_roi
import sys
sys.path.insert(1, "..")
from lmfit.models import LinearModel, ConstantModel, LorentzianModel, GaussianModel

np.random.seed(105)

def get_noise_via_running_mean(data, window_size=21, plot=False):
    
    """ Fits a spline to the dataset
        Calculates difference between spline values and data"""
    
    x = data.L
    y = data.cps_per_nA
    
    weights = np.ones(window_size) / window_size
    
    y_running_mean = ndimage.convolve(y, weights, mode='nearest')
    
    residuals = y - y_running_mean
    
    if plot:

        plt.plot(x, y, '.k', label='data')
        plt.plot(x, y_running_mean, '-r', label='running mean')
        plt.plot(x, residuals, '-b', label='residuals')
        plt.legend()
        plt.show()

    return residuals, y_running_mean

def simulate_spectra_montecarlo(data, window_size=15, number_of_mc_simulations=20, plot=False):

    """ Takes residuals between data and fitted peak,
    and resamples the residuals to generate many new synthetic datasets
    with as much noise as the original dataset. """
    
    # Generate the new synthetic datasets
    
    residuals, running_mean = get_noise_via_running_mean(data, window_size)
    
    # Select zone around the peak as separate to the other zones
    pk_max_L = data.loc[data.cps_per_nA.argmax(), 'L']
    
    pk_zone_bool = (data.L > pk_max_L-3) & (data.L < pk_max_L+3)
    pk_zone_idx = data.index[pk_zone_bool]
    lwr_zone_idx = data.index[0:pk_zone_idx[0]]
    upr_zone_idx = data.index[pk_zone_idx[-1]+1:]
    zone_list = [lwr_zone_idx, pk_zone_idx, upr_zone_idx]
        
    residuals_array = np.empty(shape=(len(data), number_of_mc_simulations))
    cps_per_nA_array = np.empty(shape=(len(data), number_of_mc_simulations))

    for i in range(number_of_mc_simulations):
        for zone in zone_list: 
            residuals_array[zone, i] = residuals[zone].sample(frac=1, replace=True)
            cps_per_nA_array[zone,i] = running_mean[zone] - residuals_array[zone, i]
        
    if plot:

        plt.plot(data.L, data.cps_per_nA, '.k', label='data')
        plt.plot(data.L, running_mean, '-r', label='running mean')
        plt.plot(data.L, cps_per_nA_array, '-y', label='_nolegend', zorder=0)
        plt.title('yellow = synthetic datasets')
        plt.legend()
        plt.show()
        
    synthetic_data = pd.DataFrame(cps_per_nA_array)
    synthetic_data.index = data.L.values
    synthetic_data.index.name = 'L'

    return synthetic_data

def select_roi_index(df, roi):
    
    """For a dataframe df, select rows within index ranges defined by roi (regions of interest)"""
    
    idx_list = [df[(df.index > r[0]) & (df.index < r[1])].index.to_list() for r in roi]
    flat_idx_list = [item for sublist in idx_list for item in sublist]
    df_roi = df.loc[flat_idx_list, :]
    
    return df_roi

def fit_bg(data, bg_type='linear'):
    
    """ bg_type can be either:
            - 'linear' (default)
            - 'lorentzian_plus_c' 
    """

    x = data.L
    y = data.cps_per_nA

    if bg_type == 'linear':
        
        mod = LinearModel()
        pars = mod.make_params(slope=0, intercept=-1)
        
    elif bg_type == 'lorentzian_plus_c':
        
        mod = LorentzianModel() + ConstantModel()
        pars = mod.make_params(amplitude=max(y)*10, sigma=10, center=120, c=0)
        
        pars['amplitude'].set(min=0)
        pars['center'].set(max=130)        
        pars['c'].set(max=min(y))
    
    bg_fit_result = mod.fit(y, pars, x=x)
    
    return bg_fit_result


def plot_spectrum_and_roi(df, roi, sample=None, baseline=None):
    
    """ Plots 'regions of interest' (roi) defined for fitting background to spectrum """
    
    fig, ax = plt.subplots(figsize=(8,2))

    plt.plot(df['L'], df['cps_per_nA'], lw=1, color='k', label='data')
    
    if baseline is not None:
        plt.plot(df['L'], baseline, lw=1, color='b', label='baseline')
    
    for r in roi:
        ax.axvspan(r[0], r[1], alpha=0.1, color='red', linewidth=0)
        
    df_roi = select_roi(df, roi)
    
    ymin = df_roi['cps_per_nA'].min() - df_roi['cps_per_nA'].max()*0.05
    ymax = df_roi['cps_per_nA'].max() + df_roi['cps_per_nA'].max()*0.3
    
    plt.ylim(ymin, ymax)
    plt.title(sample)
    plt.tight_layout()


def fit_baseline_and_plot(df, roi, name=None, bg_type='linear'):
    
    """ Fit the baseline, store results as new columns in df, and plot the fit and roi."""
       
    bg_fit_result = fit_bg(select_roi(df, roi), bg_type=bg_type)
    
    baseline = bg_fit_result.eval(x=df['L'].values)
    
    corrected_data = df['cps_per_nA'].values - baseline

    plot_spectrum_and_roi(df, roi, sample=name, baseline=baseline)

    df['baseline'] = baseline
    df['cps_per_nA_corrected'] = corrected_data


def fit_mc_bg(full_synthetic_data, roi, bg_type='linear', randomise_roi=False):

    """ bg_type can be either:
            - 'linear' (default)
            - 'lorentzian_plus_c' 
    """
    
    if randomise_roi:
    
        rng = np.random.default_rng(12345)
        rfloat = rng.random(size=len(roi)*2) - 1
        roi_adjustments = rfloat.reshape(roi.shape)
        
        roi = roi + roi_adjustments
    
    synthetic_data = select_roi_index(full_synthetic_data, roi)
    x = synthetic_data.index
    
    baseline_array = np.zeros(shape=full_synthetic_data.values.shape) 
    corrected_array = np.zeros(shape=full_synthetic_data.values.shape)
    
    for i, col in enumerate(synthetic_data.columns):
    
        if i % 10 == 0:
            print(f'Fitting montecarlo simulation {i+1} of {len(synthetic_data.columns)}')
    
        if col != 'L':
        
            y = synthetic_data.loc[:, col].values

            if bg_type == 'linear':
        
                mod = LinearModel()
                pars = mod.make_params(slope=0, intercept=-1)
        
            elif bg_type == 'lorentzian_plus_c':
        
                mod = LorentzianModel() + ConstantModel()
                pars = mod.make_params(amplitude=max(y)*10, sigma=10, center=120, c=0)
        
                pars['amplitude'].set(min=0)
                pars['center'].set(max=130)        
                pars['c'].set(max=min(y))
    
            bg_fit_result = mod.fit(y, pars, x=x)
        
            baseline_array[:,i] = bg_fit_result.eval(x=full_synthetic_data.index)
            corrected_array[:,i] = full_synthetic_data.loc[:, col].values - baseline_array[:, i]
        
    baseline = pd.DataFrame(baseline_array, index=full_synthetic_data.index)
    corrected_data = pd.DataFrame(corrected_array, index=full_synthetic_data.index)
    
    return baseline, corrected_data
    

def plot_mc_bg_fits(df, synthetic_data, baseline, roi, sample=None):
    
   
    fig, ax = plt.subplots(figsize=(15,4))

    plt.plot(df['L'], df['cps_per_nA'], lw=1, color='k')
    
    baseline.reset_index().plot(x='L', lw=1, color='b', ax=ax, legend=False)
    
    for r in roi:
        ax.axvspan(r[0], r[1], alpha=0.1, color='red', linewidth=0)
        
    df_roi = select_roi(df, roi)
    
    ymin = df_roi['cps_per_nA'].min() - df_roi['cps_per_nA'].max()*0.05
    ymax = df_roi['cps_per_nA'].max() + df_roi['cps_per_nA'].max()*0.3
    
    plt.ylim(ymin, ymax)
    plt.ylabel('cps/nA', fontsize=14)
    plt.xlabel('L (mm)', fontsize=14)
    plt.title(sample)
    plt.tight_layout()
    
    
def fit_multiple_gaussians(df_full, pk_params, 
                            y_column='counts_corrected', roi=[[135,160]],
                            samplename=None, plot=True, plot_only_roi=True):
    
    """ Fit gaussians to the N peak of <samplename> 
   
    """
    
    df = select_roi(df_full, roi)

    x = df.L
    y = df[y_column]

    mod = GaussianModel(prefix='a_')
    
    pars = mod.make_params(a_amplitude=pk_params[0]['amplitude'],
                           a_center=pk_params[0]['center'],
                           a_sigma=pk_params[0]['sigma'])
    
    for i in range(len(pk_params)):
        if i > 0:
            prefix = pk_params[i]['prefix']
            extra_peak = GaussianModel(prefix=prefix)
            extra_peak_pars = extra_peak.make_params()
            extra_peak_pars.add(prefix + 'amplitude', value=pk_params[i]['amplitude'])
            extra_peak_pars.add(prefix + 'center', value=pk_params[i]['center'])
            extra_peak_pars.add(prefix + 'sigma', value=pk_params[i]['sigma'])
            
            pars.update(extra_peak_pars)
            mod = mod + extra_peak
       
    for par_name in pars.keys():
        if 'amplitude' in par_name:       
            pars[par_name].set(min=0)
            
        elif 'center' in par_name:
        
            # Check if there is a constraint specified for min or max
            prefix = par_name.replace('center', '')
            
            user_set_pk_pars = next(entry for entry in pk_params if entry['prefix'] == prefix)
            try: 
                pars[par_name].set(min=user_set_pk_pars['center_min'], max=user_set_pk_pars['center_max'])
            except KeyError:
                # If they haven't been defined, set to a sensible range
                pars[par_name].set(min=140, max=155)
                
        elif 'sigma' in par_name:
            # Check if there is a constraint specified for min or max
            prefix = par_name.replace('sigma', '')
            user_set_pk_pars = next(entry for entry in pk_params if entry['prefix'] == prefix)
            try: 
                pars[par_name].set(min=user_set_pk_pars['sigma_min'], max=user_set_pk_pars['sigma_max'])
            except KeyError:
                # If they haven't been defined, set to a sensible range
                pars[par_name].set(min=0, max=10)
            
            
    out = mod.fit(y, x=x, params=pars)
    
    if plot:
    
        fig, ax = plt.subplots(1,2,figsize=(10,5))

        out.plot_fit(show_init=True, ax=ax[0])
        ax[0].set_title('')

        fitted_components = out.eval_components(x=df_full.L)
        fitted_curve = out.eval(x=df_full.L)
        ax[1].plot(df_full.L, df_full[y_column], '-', color='lightgrey')
        ax[1].plot(df_full.L, fitted_curve, label='fit')

        for p in fitted_components.keys():
            ax[1].plot(df_full.L, fitted_components[p], label=p)
            
        for r in roi:
            ax[1].axvspan(r[0], r[1], alpha=0.05, color='red', linewidth=0)
        
        if plot_only_roi:
            ax[1].set_xlim(np.array(roi).min(), np.array(roi).max())
        
        plt.suptitle(samplename)
        plt.legend()
        plt.show()
    
    return out


def fit_multiple_gaussians_mc(synthetic_data_corrected, pk_params, 
                            roi=[[135,160]],
                            samplename=None, plot=True):
    
    synthetic_data = synthetic_data_corrected.reset_index()
    fit_result = []
    fit = np.zeros(shape=synthetic_data_corrected.shape)
    
    for i in range(len(synthetic_data_corrected.columns)):
        if i % 10 == 0:
            print(f'Fitting montecarlo simulation {i+1} of {len(synthetic_data.columns)}')

        df = synthetic_data.loc[:, ['L', i]]
        out = fit_multiple_gaussians(df, pk_params, y_column=i, roi=roi,
                                     samplename=samplename, plot=False)
        
        fit_result.append(out)
        fit[:,i] = out.eval(x=synthetic_data_corrected.index.values)
        
    fit = pd.DataFrame(fit, index=synthetic_data_corrected.index.values)
    fit.index.name = 'L'
    
    if plot:
        fig, ax = plt.subplots(1,1,figsize=(8,4))

        synthetic_data.plot(x='L', lw=0.5, color='lightgrey', ax=ax, legend=False)
        fit.reset_index().plot(x='L', lw=1, color='blue', ax=ax, legend=False)
        
        
        clrs = ['tab:green', 'tab:orange', 'tab:purple', 'tab:cyan']
        for model in fit_result:
            fitted_components = model.eval_components(x=synthetic_data.L)
            
            for i,p in enumerate(fitted_components.keys()):
                ax.plot(synthetic_data.L, fitted_components[p], color=clrs[i], lw=0.5)
        
        for r in roi:
            ax.axvspan(r[0], r[1], alpha=0.05, color='red', linewidth=0)
        
        plt.title(samplename)
        plt.show()
    
    return fit_result
        