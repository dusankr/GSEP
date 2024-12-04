# Background radioactivity correction calculation
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
import math
import pathlib
from tkinter import messagebox
import datetime
# Local modules
from modules import globals
from modules import settings
from modules import inputs
from modules import corr_nonlin
from modules import utilities


# read background peaks from library and send it to main bcg function
def bcg_lib():

    # check if all background libraries are selected
    if not all(globals.bcg_lib_selection.values()):
        return messagebox.showerror('User input error', 'Background libraries are not selected.')
    
    # check if all background library files exist
    for det in globals.bcg_lib_selection.values():
        if not pathlib.Path(det).exists():
            return messagebox.showerror('User input error', 'Background library file does not exist: {}'.format(det))

    # delete all values from bcg_spectra_lib and bcg_spectra_param
    globals.bcg_spectra_lib = {}
    globals.bcg_spectra_param = pd.DataFrame()

    # get all lib names from globals.bcg_lib_selection
    bcg_file_names = [pathlib.Path(pathlib.Path(p).name) for p in globals.bcg_lib_selection.values()]
    detectors = [pathlib.Path(p).stem for p in globals.bcg_lib_selection.values()]

    # create dataframe with library names
    globals.bcg_spectra_param = pd.DataFrame(bcg_file_names, columns=['Name'])
   
    # add column with detector name
    globals.bcg_spectra_param['Detector'] = detectors

    # tempoarary solution TODO: change to more general solution, each might has a different path
    folder_path = pathlib.Path(list(globals.bcg_lib_selection.values())[0]).parent

    # read all background libraries, return error if more than one spectrum or parameters section is found
    spectrum = inputs.fill_data(folder_path, globals.bcg_spectra_lib, globals.bcg_spectra_param)
    if spectrum is not None:
        return messagebox.showerror('Input format error', 'The input file contains more than one spectrum or parameters section. Please check the file: {}'.format(spectrum.name))

    # take Cal en_1 and Cal en_2 from globals.spectra_param_collection and add it to bcg_spectra_param as new columns with same names based on detector column
    if 'Cal en_1' in globals.spectra_param_collection.keys():
        det = globals.spectra_param_collection['Detector'].tolist()
        cal_en_1 = globals.spectra_param_collection['Cal en_1'].tolist()
        cal_en_2 = globals.spectra_param_collection['Cal en_2'].tolist()
    
    # create dictionary with detector name as key and calibration energies as values without duplicity
    cal_en_dict = {}
    for d, c1, c2 in zip(det, cal_en_1, cal_en_2):
        if d not in cal_en_dict.keys():
            cal_en_dict[d] = [c1, c2]
    
    # add calibration energies to bcg_spectra_param - if dict key in bcg_spectra_param column and row then add calibration energies
    for d in cal_en_dict.keys():
        if d in globals.bcg_spectra_param['Detector'].tolist():
            globals.bcg_spectra_param.loc[globals.bcg_spectra_param['Detector'] == d, 'Cal en_1'] = cal_en_dict[d][0]
            globals.bcg_spectra_param.loc[globals.bcg_spectra_param['Detector'] == d, 'Cal en_2'] = cal_en_dict[d][1]
        else:
            return messagebox.showerror('User input error', 'Detector ' + d + ' is not in background library file.')


# background main correction function
def background():

    # load background files and their parameters
    bcg_lib()

    # apply non-linearity correction if available in spectra collection - check first spectrum
    if 'cor_en' in globals.spectra_collection[list(globals.spectra_collection.keys())[0]]:
        corr_nonlin.nonlinearity(globals.bcg_spectra_param, globals.bcg_spectra_lib, '.bcg')
        print('(For background spectra.)')

    # log file for background correction - list
    bcg_log = []
    
    # go through all spectra and apply background correction from background library
    for spectrum in globals.spectra_collection:
        spec = pathlib.Path(spectrum + '.prn')
        det = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Detector'].item()
        t_live = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Live time'].item()

        if det not in globals.bcg_spectra_param['Detector'].tolist():
            return messagebox.showerror('User input error', 'Background library does not include data for detector: ' + det)

        if 'cor_en' in globals.spectra_collection[spectrum] and 'cor_en' in globals.bcg_spectra_lib[det]:
            en = globals.spectra_collection[spectrum]['cor_en'].tolist()
            en_unc = globals.spectra_collection[spectrum]['cor_en_err'].tolist()

            bcg_en = globals.bcg_spectra_lib[det]['cor_en'].tolist()
            bcg_en_unc = globals.bcg_spectra_lib[det]['cor_en_err'].tolist()
        else:
            en = globals.spectra_collection[spectrum]['energy'].tolist()
            en_unc = [e * 0.001 for e in en]     # TODO: temporary solution, change to more general solution based on channel error

            bcg_en = globals.bcg_spectra_lib[det]['energy'].tolist()
            bcg_en_unc = [e * 0.001 for e in bcg_en]     # TODO: temporary solution, change to more general solution based on channel error

        area = globals.spectra_collection[spectrum]['area'].tolist()
        # convert area error from % to number of counts
        area_unc = [ar * err / 100 for err, ar in zip(globals.spectra_collection[spectrum]['aerr%'].tolist(), area)]
        fwhm = globals.spectra_collection[spectrum]['fwhm'].tolist()
        chn = globals.spectra_collection[spectrum]['channel'].tolist()

        bcg_area = globals.bcg_spectra_lib[det]['area'].tolist()
        # convert area error from % to number of counts
        bcg_area_unc = [ar * err / 100 for err, ar in zip(globals.bcg_spectra_lib[det]['aerr%'].tolist(), bcg_area)]
        bcg_t_live = globals.bcg_spectra_param.loc[globals.bcg_spectra_param['Detector'] == det, 'Live time'].item()

        bcg_log.append('\nSpectrum: ' + spectrum + ' correction started.')

        # original areas for comparison
        area_orig = area.copy()

        # background correction
        for i in range(0, len(en)):
            for j in range(0, len(bcg_en)):
                if (abs(en[i] - bcg_en[j]) < settings.N_BCG * math.sqrt(settings.A_BCG**2 + (en_unc[i]/en[i])**2 + (bcg_en_unc[j]/bcg_en[j])**2)) or (abs(en[i] - bcg_en[j]) < settings.M_BCG * fwhm[i] * en[i] / chn[i]):
                    bcg_log.append('Peak difference: ' + str(abs(en[i] - bcg_en[j])))
                    
                    area[i] -= bcg_area[j] * t_live / bcg_t_live
                    area_unc[i] = math.sqrt((area_unc[i] / area[i])**2 + (bcg_area_unc[j] / bcg_area[j])**2) * area[i]
                    
                    # create log in file
                    bcg_log.append('Peak at {} keV corrected by background peak at {} keV with area {:.1f} counts.'.format(en[i], bcg_en[j], area_orig[i] - area[i]))
                    
                    if area[i] < 0:
                        bcg_log.append('Peak at ' + str(en[i]) + ' keV has negative area after background correction. In the next step will be removed.')
                    
        globals.spectra_collection[spectrum]['bcg_diff'] = [round(a - b, 1) for a, b in zip(area_orig, area)]
        globals.spectra_collection[spectrum]['bcg cor area'] = [round(x, 1) for x in area]
        globals.spectra_collection[spectrum]['bcg cor area unc'] = [round(x, 1) for x in area_unc]

        # remove peaks with negative or zero area
        globals.spectra_collection[spectrum] = globals.spectra_collection[spectrum][globals.spectra_collection[spectrum]['bcg cor area'] > 0]
    
    # check if the output folder exists
    utilities.check_folders()

    # save log file
    timestamp = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-'))
    try:
        bcg_log_path = globals.common_params['Output path'] / pathlib.Path(timestamp + '_bcg.log')
        with open(bcg_log_path, 'w') as f:
            for line in bcg_log:
                f.write(line + '\n')
    except Exception as e:
        return messagebox.showerror('Error', 'Error while saving log file: ' + str(e))

    print('Background correction done.')