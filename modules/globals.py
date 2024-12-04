# Configuration file with Global variables shared across modules
#
# Copyright (C) 2024  Dusan Kral, Vendula Filova
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# Libraries
import pandas as pd

# definition of variables shared between files
spectra_collection = {}                     # dictionary of all spectra in pandas tables
spectra_param_collection = pd.DataFrame()   # create an empty pandas dataframe

# definition of variabels shared between files
# dictionary of variables common for whole input folder
common_params = {
    'Irradiation time': '0',
    'End irr': '01/01/2010 00:00:00',
    'Beam int': '0E+0',
    'Beam int unc': '0E+0'
    }

# Libraries ---------------------------------------------------------------------------------------------------
# collection = all available files/paths for each library/detector
# selection = currently active with one library value


# all libraries are stored here:
lib_collection = {
    'Attenuation': [r'libraries\att.lib'],
    'Efficiency': [r'libraries\eff.lib'],
    'Gamma lib. (all)': [r'libraries\gamma_all.llb'],
    'Nonlinearity': [r'libraries\nonlin.lib'],
    'Beam Fluctuation': [r'libraries\beam.lib'],
    'Activity': [r'libraries\activity.lib'],
    'Gamma lib. (eff)': [r'libraries\gamma_eff.llb']
    }

# used libraries are stored here:
lib_selection = {
    'Attenuation': r'libraries\att.lib',
    'Efficiency': r'libraries\eff.lib',
    'Gamma lib. (all)': r'libraries\gamma_all.llb',
    'Nonlinearity': r'libraries\nonlin.lib',
    'Beam Fluctuation': r'libraries\beam.lib',
    'Activity': r'libraries\activity.lib',
    'Gamma lib. (eff)': r'libraries\gamma_eff.llb'
    }

# coincidence library files
coinc_lib_collection = {}   # list of paths for each detector
# used libraries are stored here
coinc_lib_selection = {}    # single value for each detector

# all libraries are stored here:
bcg_lib_collection = {}
# used libraries are stored here
bcg_lib_selection = {}

# Correction variables ----------------------------------------------------------------------------------------

# some variable from calc_corr.py for work with coincidences
coinc_track_var = []

# background correction
bcg_track_var = []

# Background correction variables -----------------------------------------------------------------------

# dictionary of all background spectra in pandas tables
bcg_spectra_lib = {}

# create an empty pandas dataframe for background correction parameters
bcg_spectra_param = pd.DataFrame()

# Efficiency calculation window --------------------------------------------------------------------------

# data for efficiency calculation
eff_param_collection = pd.DataFrame()

# gamma library data
gamma_lib_eff = {}

# activity library data
activity_lib_df = pd.DataFrame()

# calculated efficiency results
eff_results_df = pd.DataFrame(columns=['E (keV)', 'E_ref (kev)','Nuclide', 'Int (-)', 'A0 (Bq)', 'COI_corr', 'Eff (-)', 'dEff (-)'])

# final fits parameters
eff_fit_dict = {}

# Coincidence calculation window -------------------------------------------------------------------------

# dictionary of all .tcd files and paths
tcd_files_dict = {}

# coincidences input informations dictionary
coinc_dict = {
    'tc path': '',
    'tc input': ''
    }

# Nonlin calculation window ------------------------------------------------------------------------------

# dictionary for results of nonlinearity - used in nonlinearity window (Detector MENU)
nonlin_results_df = pd.DataFrame(columns=['E_lib (keV)', 'Nuclide', 'Int_lib (%)', 'E (keV)', 'Diff (keV)'])
coinc_fit_res = None

# Calculation window -------------------------------------------------------------------------------------

# [en, en_unc, inten, int_unc, half_life, half_life_unc, atomic_number]
gamma_lib = {}  # dictionary of gamma library data TODO use for efficiency correction

# matplotlib variables
canvas_id = None
fig_id = None
ax = None

# RR instead of global variable in main file
#rr_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'RR (arb.)', 'dRR (arb.)', 'COI_corr', 'Beam_corr', 'A (s^-1)', 'dA (s^-1)', 't_delay (hr)', 'A_0 (s^-1)'])
rr_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'RR (arb.)', 'dRR (arb.)', 'C_irr,g', 'COI_corr', 'Beam_corr', 'N (-)', 'dN (-)', 't_delay (hr)', 'N_0 (-)'])

# raw data for export to (dataframe from RR calculation)
#rr_all_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'RR (arb.)', 'dRR (arb.)', 'COI_corr', 'Beam_corr', 'A (s^-1)', 'dA (s^-1)', 't_delay (hr)'])
rr_all_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'RR (arb.)', 'dRR (arb.)', 'C_irr,g', 'COI_corr', 'Beam_corr', 'N (-)', 'dN (-)', 't_delay (hr)'])

# rr results for isomer correction calculation
# rr_iso_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'N_yield,m (-)', 'dN_yield,m (-)', 'COI_corr', 'Beam_corr', 'A (s^-1)', 'dA (s^-1)', 't_delay (hr)', 'N_m0 (-)', 'dN_m0 (-)'])
rr_iso_results = pd.DataFrame(columns=['Spectrum', 'Nuclide', 'E (keV)', 'N_yield,m (-)', 'dN_yield,m (-)', 'COI_corr', 'Beam_corr', 'N (-)', 'dN (-)', 't_delay (hr)', 'N_m0 (-)', 'dN_m0 (-)'])

# weighted average of N_m0 , dN_m0, N_yield,m, dN_yield,m
iso_res_dict = {}
# iso_res_dict = {
#     'N_m0 (-)': 0,
#     'dN_m0 (-)': 0,
#     'N_yield,m (-)': 0,
#     'dN_yield,m (-)': 0,
#     'nuclide': '',
#     }

# dictionary for RR calculation to avoid global variables
widgets_dict = {}
