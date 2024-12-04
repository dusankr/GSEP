# nonlinearity correction calculation
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
import pathlib
import numpy as np
import tkinter
from modules import globals


# checks whether calibration energies are uploaded and values are nonzero
def check_cal_en():
    if ('Cal en_1' or 'Cal en_2') in globals.spectra_param_collection:
        if (any(globals.spectra_param_collection['Cal en_1']) or
                any(globals.spectra_param_collection['Cal en_2'])) == 0:
            return 'disable'
        else:
            return 'active'
    else:
        return 'disable'


# return nonlin functions parameters from file NonLin_py.lib
def nonlin_lib():
    with open(globals.lib_selection['Nonlinearity'], 'r', encoding='utf-8') as temporary_file:
        content = temporary_file.readlines()

    nonlin_dict = {}

    for i in range(0, len(content)):
        line = content[i].strip().split()
        if len(line) > 0 and content[i][0] != '#':
            # if first element is letter and second is 0
            if line[0].isalpha() and float(line[1]) == 0:
                det = line[0]
                en_div = float(line[1])
                line = content[i+1].strip().split()
                set1 = [float(a) for a in line]
                nonlin_dict[det] = [en_div, set1, None]
            elif line[0].isalpha() and float(line[1]) > 0:
                det = line[0]
                en_div = float(line[1])
                line = content[i+1].strip().split()
                set1 = [float(a) for a in line]   # first set of parameters
                line = content[i+2].strip().split()
                set2 = [float(a) for a in line]   # second set of parameters
                nonlin_dict[det] = [en_div, set1, set2]
            
    return nonlin_dict


# return corrected energy
def nonlin_calc(en_div, a1, a2, en, cal_en_1, cal_en_2):
    
    # return result of correction function
    def nonlin_cor_eq(energy):
        if en_div == 0:
            a = a1
        elif energy <= en_div:
            a = a1
        else:
            a = a2
        
        return a[0] + a[1] * energy + a[2] * energy**2


    p_1 = (nonlin_cor_eq(cal_en_1) - nonlin_cor_eq(cal_en_2)) / (cal_en_1 - cal_en_2)
    q_1 = nonlin_cor_eq(cal_en_1) - cal_en_1 * p_1
    
    return en - (nonlin_cor_eq(en) - (p_1 * en + q_1))


# NonLin main correction function
def nonlinearity(temp_df, temp_spec_dict, lib_suffix):

    # get nonlinearity parameters from library file
    nonlin_dict = nonlin_lib()

    # loop through all spectra
    for spectrum in temp_spec_dict:
        spec = pathlib.Path(spectrum + lib_suffix)

        # get detector, calibration energies and nonlinearity parameters for each spectrum      TODO move before loop and create dictionary instead of this
        det = temp_df.loc[temp_df['Name'] == spec, 'Detector'].item()
        cal_en_1 = temp_df.loc[temp_df['Name'] == spec, 'Cal en_1'].item()
        cal_en_2 = temp_df.loc[temp_df['Name'] == spec, 'Cal en_2'].item()

        # get parameters of nonlinearity correction function(s) and divide energy if present
        try:
            en_div, a1, a2 = nonlin_dict[det]
        except KeyError:
            tkinter.messagebox.showerror('Error', 'Nonlinearity correction parameters for detector ' + det + ' are missing.')
            continue

        # calculate corrected energy
        correct_en = []
        for en in temp_spec_dict[spectrum]['energy']:
            correct_en.append(nonlin_calc(en_div, a1, a2, en, cal_en_1, cal_en_2))

        # store corrected energy into new column 'cor_en'
        temp_spec_dict[spectrum]['cor_en'] = np.round(correct_en, decimals=4)

        # calculate energy change and store it into new column 'en_diff'
        temp_spec_dict[spectrum]['en_diff'] = np.round(temp_spec_dict[spectrum]['energy'] - temp_spec_dict[spectrum]['cor_en'], decimals=4)

        # calculate energy error and store it into new column 'cor_en_err'
        # TODO this part need redefinition, this is only temporary solution for corrected energy error
        err = []
        for en in temp_spec_dict[spectrum]['en_diff']:
            if en < 1:
                err.append(np.abs(en * 0.02))    # 2% of energy
            elif en >= 1:
                err.append(np.sqrt(np.abs(en)))    # square root of energy
        
        temp_spec_dict[spectrum]['cor_en_err'] = np.round(err, decimals=4)

        # calculate cor_fwhm_kev column from original channel, corrected energy and original fwhm in channels
        # temp_spec_dict[spectrum]['cor_fwhm_keV'] = np.round(temp_spec_dict[spectrum]['fwhm'] * temp_spec_dict[spectrum]['cor_en'] / temp_spec_dict[spectrum]['channel'], decimals=3)

    print('Nonlinearity correction done.')