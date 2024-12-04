# module for storing constants and settings
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

# SEP and DEP correction parameters -----------------------------------------------------------------------------------

# n is variable, that can be chosen according to statistics (n≈1)
N_EP = 1
# m is a correction factor, which takes into account the non-Gausian nature of the fitted gamma peak (m=1)
M_EP = 1
# a is a correction factor, which takes into account the non-linearity effect (a=0.15)
A_EP = 0.15

# SEP dependece function parameters
PAR_SEP = [-2.0866027851, 13.8832738535, -8.2291578426, 1.7178618884]
PAR_SEP_ERR = [0.2860230743, 1.1099975243, 1.3992857931, 0.5675750337]
# DEP dependece function parameters
PAR_DEP = [-1.3273140355, 12.8053750577, -8.9129115117, 2.6470549618]
PAR_DEP_ERR = [0.1814312750, 0.8383854495, 1.1473203273, 0.4847393276]

# Background correction parameters ------------------------------------------------------------------------------------

# background subtraction parameters
N_BCG = 2
# correction factor for non-linearity effect
A_BCG = 0.15
# correction factor, which takes into account the non-Gausian nature of the fitted gamma peak
M_BCG = 0.5

# Calculation parameters ---------------------------------------------------------------------------------------------
N_CALC = 2
A_CALC = 0.15

# ---------------------------------------------------------------------------------------------------------------------

# Tkinter constants ---------------------------------------------------------------------------------------------------

# treeview column width
COLUMN_WIDTH = {
    'Name': 120,
    'Start of measurement': 140,
    'Live time': 70,
    'Real time': 70,
    'Detector': 70,
    'Position': 70,
    'Sample name': 90,
    'Measurement number': 81,
    'Delayed time': 85,
    'Weight': 66,
    'Weight error': 80,
    'Molar mass': 77,
    'Element': 65,
    'Fraction by weight': 110,
    'Density': 70,
    'Thickness': 70,
    'Cal en_1': 60,
    'Cal en_2': 60,
    'Non-point': 70,
    'channel': 70,
    'cherr': 55,
    'energy': 70,
    'area': 50,
    'aerr%': 55,
    'fwhm': 50,
    'chisq': 50,
    'it': 40,
    'left': 40,
    'n': 30,
    'lim': 30,
    'isotope': 60,
    'cor_en': 70,
    'en_diff': 60,
    'cor_en_err': 70,
    'bcg cor area': 80,
    'bcg cor area unc': 110,
    'sep_dep_diff': 90,
    'SEP-DEP cor area': 120,
    'SEP-DEP cor area unc': 140,
    'Spectrum': 100,
    'Nuclide': 50,
    'E (keV)': 60,
    'RR (arb.)': 80,
    'dRR (arb.)': 80,
    'Beam_corr': 65,
    'A (s^-1)': 80,
    'dA (s^-1)': 80,
    't_delay (hr)': 80,
    'A_0 (s^-1)': 80,
    'E_ref (kev)': 70,
    'Int (%)': 80,
    'A0 (Bq)': 80,
    'Eff (-)': 80,
    'dEff (-)': 80,
    'COI_corr': 65,
    'Int_lib (%)': 80,
    'E_lib (keV)': 80,
    'Diff (keV)': 80,
    'N_yield,m (-)': 80,
    'dN_yield,m (-)': 80,
    'N_m0 (-)': 70,
    'dN_m0 (-)': 70,
    'N_g0 (-)': 70,
    'dN_g0 (-)': 70,
    'C_irr,g': 65,
    'N (-)': 70,
    'dN (-)': 70,
    'N_0 (-)': 70,
}

# Calculation window constants ----------------------------------------------------------------------------------------

# figure dimensions
FIG_SIZE = (6.4, 4.8) # inches default

# Efficiency calculation constants ------------------------------------------------------------------------------------
# numbers determines the range of the polynomial fit
# add 1 to the maximum value to include the last value

# single
# min and max degree of polynomial fit
MIN_S_PARAM = 5
MAX_S_PARAM = 9

# double
# first curve parameters
MIN_D1_PARAM = 4
MAX_D1_PARAM = 8
# second curve parameters
MIN_D2_PARAM = 2
MAX_D2_PARAM = 5
