# Imports
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import multiprocessing as mp
from scipy.stats import pearsonr
from functools import partial


def MichaelisEq(co2Conc, VmaxRatio, KmMut, KmWT=0.149):
    """
    # Michaelis-Menten equation for mutant velocity relative to WT
    # KmWT is the WT CO2 Km in mM (given above as global variable)
    # co2Conc is the CO2 concentration vector in mM
    # VmaxRatio is the ratio of the mutant to wt Vmax
    # KmMut is the mutant Km in mM
    :param co2Conc:
    :param VmaxRatio:
    :param KmMut:
    :return:
    """
    return (VmaxRatio * (KmWT + co2Conc))/(KmMut + co2Conc)


# Michaelis-Menten curve fitting function
# x is the substrate concentration vector
# row is the vector of enrichments paired with the concentrations x
# returns VmaxRatio and KmMut fit values
def fitMichaelis(x, row):
    # find indices withount nans and filter
    valid_indices = np.logical_and(~np.isnan(x), ~np.isnan(row))
    x_filtered = np.array(x)[valid_indices]
    y_filtered = np.array(row)[valid_indices]
    
    # check that at least 3 points exist, otherwise return nans for the fit
    if valid_indices.sum() <= 3:
        return np.array([np.nan, np.nan])
    
    # use SciPy curve_fit to fit MichaelisEq with bounds of [0, 20] for both parameters
    popt, pcov = curve_fit(MichaelisEq, x_filtered, y_filtered, bounds=[0, 20])
    stds = np.sqrt(np.diag(pcov))

    # check if the fits are pressed against the bounds and set fits to nan if so
    if (popt[1] < 0.001) | (popt[1] > 19.999) | (popt[0] < 0.001) | (popt[0] > 19.999):
        popt[0] = np.nan
        popt[1] = np.nan
    
    return np.array([popt[1], popt[0]])


# runs Michaelis-Menten bootstrap fits and returns various distributional information on the fits
def mut_Km_bootstrap(mutname, enrich_combine, condnames0, co2_conc):

    # grab the enrichment data for a particular mutant
    mutdf = enrich_combine[enrich_combine.index.get_level_values('mutant') == mutname]
    # initialize the empty list for the fit parameters
    mutKms = []

    # iterate over the rows of the mutdf dataframe, each row representing a difference combination
    # of pseudocount and threshold parameters.
    for row in mutdf.iterrows():
        # create a set of random indices for the bootstrap
        binds = make_bootstrap_inds(condnames0, 11)
        # convert the above indices, binds, into replicate names for calling columns
        bootname_set = bootstrap_names(binds, condnames0)

        # iterate over the bootstraps, fit the data, and append the fits
        # (using a try to catch some rare edge cases that might not actually exist anymore)
        for bnames in bootname_set:
            boot_data = row[1][bnames]
            
            try:
                fitKm = fitMichaelis(boot_data, co2_conc)
            except RuntimeError:
                fitKm = [np.nan, np.nan]
                
            mutKms.append(np.array(fitKm))
    
    # convert the fit parameters into a numpy array
    mutKms = np.array(mutKms)
    # calculate the fraction of nans, basically the fraction of fits that failed for one reason or another
    nanfrac = np.sum(np.isnan(mutKms[:,0])) / mutKms.shape[0]
    # compute the percentiles of the Vmax distribution
    percentilesVmax = np.nanpercentile(mutKms[:,1], [25, 50, 75])
    # get the median Vmax
    medianVmax = percentilesVmax[1]
    # calculate the quantile-based coefficient of variation for the Vmax
    qbcovVmax = (percentilesVmax[2]-percentilesVmax[0])/medianVmax
    
    # compute the percentiles of the Km distribution
    percentilesKm = np.nanpercentile(mutKms[:,0], [25, 50, 75])
    # get the median Km
    medianKm = percentilesKm[1]
    # calculate the quantile-based coefficient of variation for the Km
    qbcovKm = (percentilesKm[2]-percentilesKm[0])/medianKm

    # calculate the Pearson correlation between the Vmax and Km for the distribution
    good_inds = np.logical_or( ~np.isnan(mutKms[:,1]), ~np.isnan(mutKms[:,0]))
    try:
        pear_corr = pearsonr(mutKms[good_inds, 1], mutKms[good_inds, 0])[0]
    except:
        pear_corr = np.nan

    # return various statistics from the fit distributions (many of which are not used)
    return [np.nanmean(mutKms[:,0]), np.nanstd(mutKms[:,0]), medianKm, qbcovKm, np.nanmean(mutKms[:,1]), np.nanstd(mutKms[:,1]), medianVmax, qbcovVmax, nanfrac, pear_corr, percentilesVmax[0], percentilesVmax[2], percentilesKm[0], percentilesKm[2]]


# create a set of numeric indices for a set of n bootstraps
def make_bootstrap_inds(condnames0, n0):
    bootinds = np.zeros((n0, len(condnames0))).astype(int)
    for cind, colind in enumerate(condnames0):
        bootinds[:, cind] = np.random.randint(0, len(colind), n0)
    return bootinds


# convert the numeric indices into the appropriate column names
def bootstrap_names(bootinds0, condnames0):
    bootnames_set = []
    for bootind0 in bootinds0:
        bootnames = []
        for cn, bi in zip(condnames0, bootind0):
            bootnames.append(cn[bi])
        bootnames_set.append(bootnames)
    return bootnames_set


def main():
    # Snakemake I/O
    # === Inputs
    # ==
    enrichment_data = str(snakemake.input.parameter_sweep_table)
    # === Outputs
    bootstrap_data_out = str(snakemake.output.bootstrap_data)

    # enrich_combine = pd.read_csv('/groups/doudna/projects/daniel_projects/rubiscodms/py/2024_04_30_enrich_slice.csv', index_col=[0, 1, 2])

    # == Constant Variable
    # reads in large multi-index dataframe with indices of pseudocount, threshold, and mutant name.
    enrich_combine = pd.read_csv(enrichment_data, index_col=[0, 1, 2])

    # dictionary with experiments belonging to each condition
    conditionsKey = {
        'WT_only': ['NP_11_64_2'],
        't0': ['NP_11_64_9', 'NP_11_64_11', 'NP_11_64_12', 'NP_11_64_13', 'NP_11_64_15'],
        '0.04percent_CO2': ['NP_11_66_41'],
        '0.2percent_CO2': ['NP_11_66_1', 'NP_11_66_2', 'NP_11_66_3'],
        '0.3percent_CO2': ['NP_11_66_7', 'NP_11_66_8', 'NP_11_66_9'],
        '0.4percent_CO2': ['NP_11_66_4', 'NP_11_66_5', 'NP_11_66_6'],
        '1percent_CO2': ['NP_11_66_10', 'NP_11_66_11', 'NP_11_66_12'],
        '5percent_CO2_20uM_IPTG': ['NP_11_66_13', 'NP_11_66_14', 'NP_11_66_15',
                                   'NP_11_66_16', 'NP_11_66_17', 'NP_11_66_18',
                                   'NP_11_66_37', 'NP_11_66_38', 'NP_11_66_39'],
        '0uM_IPTG': ['NP_11_66_19', 'NP_11_66_20', 'NP_11_66_21'],
        '5uM_IPTG': ['NP_11_66_22', 'NP_11_66_23', 'NP_11_66_24'],
        '30uM_IPTG': ['NP_11_66_25', 'NP_11_66_26', 'NP_11_66_27'],
        '100uM_IPTG': ['NP_11_66_28', 'NP_11_66_29', 'NP_11_66_30'],
        '300uM_IPTG': ['NP_11_66_31', 'NP_11_66_32', 'NP_11_66_33'],
        '1000uM_IPTG': ['NP_11_66_34', 'NP_11_66_35', 'NP_11_66_36'],
        'Xylose': ['NP_11_66_40']

    }

    # dictionary of CO2 and IPTG titration conditions
    titrationDict = {'CO2': ['0.04percent_CO2',
                             '0.2percent_CO2',
                             '0.3percent_CO2',
                             '0.4percent_CO2',
                             '1percent_CO2',
                             '5percent_CO2_20uM_IPTG'],
                     'IPTG': ['0uM_IPTG',
                              '5uM_IPTG',
                              '5percent_CO2_20uM_IPTG',
                              '30uM_IPTG',
                              '100uM_IPTG',
                              '300uM_IPTG',
                              '1000uM_IPTG']
                     }

    # Henry's law constant at 37C
    kH = 0.0249

    # CO2 concentrations in atm
    concCO2_atm = np.array([0.04, 0.2, 0.3, 0.4, 1, 5]) / 100

    # CO2 conc (mM) and unused IPTG conc (uM)
    concDict = {'CO2': concCO2_atm * kH * 1000,
                'IPTG': [0, 5, 20, 30, 100, 300, 1000]
                }

    # get a list of all mutations
    allmuts = pd.unique(enrich_combine.index.get_level_values('mutant'))

    # set of experiment replicates grouped by condition.
    condnames = [conditionsKey[cond] for cond in titrationDict['CO2'][1:]]
    co2_conc = concDict['CO2'][1:]

    partial_function = partial(mut_Km_bootstrap,
                               enrich_combine=enrich_combine,
                               condnames0=condnames,
                               co2_conc=co2_conc)

    # create a multiprocessing pool with 36 processes to parallelize the bootstrap fitting
    pool = mp.Pool(36)
    # map Km bootstraping to all mutants
    Km_bootstrap_lists = pool.map(partial_function, allmuts)
    pool.close()
    pool.join()

    # make a dictionary of the bootstrapping summary statistics for all mutants
    Km_bootstrap = dict(zip(allmuts, Km_bootstrap_lists))

    # create a data frame
    df = pd.DataFrame(Km_bootstrap).T

    # give the columns and index meaningful names
    df.columns = ['Km_mean', 'Km_std', 'Km_median', 'Km_qbcov', 'Vmax_mean','Vmax_std', 'Vmax_median', 'Vmax_qbcov', 'NaN_fraction','PearsonCorr', 'Vmax_25p', 'Vmax_75p', 'Km_25p', 'Km_75p']
    df.index.name = 'mutant'

    # write the data to a csv
    df.to_csv(bootstrap_data_out) #('2024_04_30_testout1.csv')


if __name__ == '__main__':
    main()
