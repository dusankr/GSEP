# single and double escape peak correction calculation
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
from tkinter import messagebox
import math
import pathlib
import datetime
# Local libraries
from modules import globals
from modules import settings
from modules import utilities

# sepdep correction function
def sep_dep_func(a, en):
    cor = 0
    for i in range(0,len(a)):
        cor += a[i] * math.log(en/1000)**i
    return math.exp(cor)/1000   # specific for SEP/DEP correction


# sepdep main correction function
def sepdep():
    
    # empty list for logging prints instead printing into terminal - save to output folder as sepdep.log
    sepdep_log = []

    for spectrum in globals.spectra_collection:
        if 'cor_en' in globals.spectra_collection[spectrum]:
            en = globals.spectra_collection[spectrum]['cor_en'].tolist()
            en_unc = globals.spectra_collection[spectrum]['cor_en_err'].tolist()
        else:
            en = globals.spectra_collection[spectrum]['energy'].tolist()
            en_unc = [e * 0.001 for e in en]     # TODO: temporary solution, change to more general solution based on channel error

            # return messagebox.showerror('User input error', 'Escape peaks correction need non-linearity correction.')

        if 'bcg cor area' in globals.spectra_collection[spectrum]:
            area = globals.spectra_collection[spectrum]['bcg cor area'].tolist()
            area_unc = globals.spectra_collection[spectrum]['bcg cor area unc'].tolist()
            # print('SEP/DEP correction works with BCG corrected area for {}.'.format(spectrum))
            sepdep_log.append('\nSEP/DEP correction works with BCG corrected area for {}.'.format(spectrum))
        else:
            area = globals.spectra_collection[spectrum]['area'].tolist()
            area_unc = [x * y /100 for x, y in zip(area, globals.spectra_collection[spectrum]['aerr%'].tolist())]
            
            # print('SEP/DEP correction works with spectrum original for {}.'.format(spectrum))
            sepdep_log.append('\nSEP/DEP correction works with original spectrum area for {}.'.format(spectrum))
        
        # independent on other corrections
        chn = globals.spectra_collection[spectrum]['channel'].tolist()
        fwhm = globals.spectra_collection[spectrum]['fwhm'].tolist()
        
        # SEP/DEP correction
        temp_corr_areas = {}    # temporary dictionary for corrections
        
        for i in range(len(en)-1, -1, -1): # loop over all peaks
            if en[i] > 1022:
                #for j in range(0, len(en)):
                for j in range(i-1, -1, -1):
                    if en[j] >= (en[i] - 500):   # if en[j] is higher than en[i] - 500 keV it does not make sense to continue (small overlap due to FWHM)
                        continue
                    elif (abs(en[j] - (en[i] - 511)) < settings.N_EP * math.sqrt(settings.A_EP**2 + (en_unc[i]/en[i])**2 + (en_unc[j]/en[j])**2)) or (abs(en[j] - (en[i] - 511)) < settings.M_EP * fwhm[i] * en[i] / chn[i]):
                        # if the condition is met, the peak j is a SEP of peak i
                        # from dependence function is calculated the correction area for peak j from peak i and added to the temporary dictionary
                        # check if the peak is already in the dictionary
                        if en[j] in temp_corr_areas and temp_corr_areas[en[j]][2] != 'DEP':
                            # if the peak is already in the dictionary, the correction area is added to the existing value
                            area_corr = area[i] * sep_dep_func(settings.PAR_SEP, en[i])
                            area_unc_corr_rel = sep_dep_func(settings.PAR_SEP_ERR, en[i])

                            temp_corr_areas[en[j]][0] = temp_corr_areas[en[j]][0] + area_corr
                            temp_corr_areas[en[j]][1] = temp_corr_areas[en[j]][0] * math.sqrt((temp_corr_areas[en[j]][1] / temp_corr_areas[en[j]][0])**2 + area_unc_corr_rel**2)
                            temp_corr_areas[en[j]][2] = 'SEP'
                        
                            # make a log
                            sepdep_log.append('Peak at {} keV has SEP at {} keV. Again.'.format(en[i], en[j]))
                        elif en[j] not in temp_corr_areas:
                            # if the peak is not in the dictionary, the correction area is added to the dictionary
                            area_corr = area[i] * sep_dep_func(settings.PAR_SEP, en[i])
                            area_unc_corr = area[i] * sep_dep_func(settings.PAR_SEP_ERR, en[i])     # relative error

                            # add the peak to the temporary dictionary
                            temp_corr_areas[en[j]] = [area_corr, area_unc_corr, 'SEP']
                            # make a log
                            sepdep_log.append('Peak at {} keV has SEP at {} keV.'.format(en[i], en[j]))
                        elif en[j] in temp_corr_areas and temp_corr_areas[en[j]][2] == 'DEP':
                            # if the peak is already in the dictionary as DEP, the SEP will not be added
                            sepdep_log.append('Peak at {} keV has SEP at {} keV. SEP will not be added. Already identified as DEP.'.format(en[i], en[j]))
                    
                    elif (abs(en[j] - (en[i] - 1022)) < settings.N_EP * math.sqrt(settings.A_EP**2 + (en_unc[i]/en[i])**2 + (en_unc[j]/en[j])**2)) or (abs(en[j] - (en[i] - 1022)) < settings.M_EP * fwhm[i] * en[i] / chn[i]):
                        # if the condition is met, the peak j is a DEP of peak i
                        # from dependence function is calculated the correction area for peak j from peak i and added to the temporary dictionary
                        # check if the peak is already in the dictionary
                        if en[j] in temp_corr_areas:
                            # if the peak is already in the dictionary, the correction area is added to the existing value
                            area_corr = area[i] * sep_dep_func(settings.PAR_DEP, en[i])
                            area_unc_corr_rel = sep_dep_func(settings.PAR_DEP_ERR, en[i])

                            temp_corr_areas[en[j]][0] = temp_corr_areas[en[j]][0] + area_corr
                            temp_corr_areas[en[j]][1] = temp_corr_areas[en[j]][0] * math.sqrt((temp_corr_areas[en[j]][1] / temp_corr_areas[en[j]][0])**2 + area_unc_corr_rel**2)
                            temp_corr_areas[en[j]][2] = 'DEP'
                            
                            # make a log
                            sepdep_log.append('Peak at {} keV has DEP at {} keV. Again.'.format(en[i], en[j]))
                        else:
                            # if the peak is not in the dictionary, the correction area is added to the dictionary
                            area_corr = area[i] * sep_dep_func(settings.PAR_DEP, en[i])
                            area_unc_corr = area[i] * sep_dep_func(settings.PAR_DEP_ERR, en[i])
                            
                            # add the peak to the temporary dictionary
                            temp_corr_areas[en[j]] = [area_corr, area_unc_corr, 'DEP']
                            # make a log
                            sepdep_log.append('Peak at {} keV has DEP at {} keV.'.format(en[i], en[j]))

        # make a copy of the original areas for correction
        area_orig = area.copy()
        
        # use the temporary dictionary to correct the areas (subtract the correction area from the peak area)
        for i in range(0, len(en)):
            if en[i] in temp_corr_areas:
                area[i] = area[i] - temp_corr_areas[en[i]][0]
                area_unc[i] = math.sqrt(area_unc[i]**2 + temp_corr_areas[en[i]][1]**2)
                
                #print('Peak at', en[i], 'keV corrected by', round(temp_corr_areas[en[i]][0], 1), '±', round(temp_corr_areas[en[i]][1], 1), 'cps.')
                # delete the peak from the temporary dictionary
                del temp_corr_areas[en[i]]

        # create new column with SEP-DEP corrections
        globals.spectra_collection[spectrum]['sep_dep_diff'] = [round(a - b, 1) for a, b in zip(area_orig, area)]

        # create new column with SEP-DEP corrected areas
        globals.spectra_collection[spectrum]['SEP-DEP cor area'] = [round(x, 1) for x in area]
        globals.spectra_collection[spectrum]['SEP-DEP cor area unc'] = [round(x, 1) for x in area_unc]

        # print if area is negative
        for i, ar in enumerate(area):
            if ar < 0:
                #print('Negative area after SEP/DEP correction in', spectrum, ' for energy', en[i], 'keV. Will be deleted.')
                sepdep_log.append('Negative area after SEP/DEP correction in {} for energy {} keV. Will be deleted.'.format(spectrum, en[i]))
        
        # delete peaks which have negative area after correction
        globals.spectra_collection[spectrum] = globals.spectra_collection[spectrum][ globals.spectra_collection[spectrum]['SEP-DEP cor area'] > 0 ]

    # check if the output folder exists
    utilities.check_folders()

    # save log to file
    timestamp = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-'))
    try:
        sepdep_log_path = globals.common_params['Output path'] / pathlib.Path(timestamp + '_sepdep.log')
        with open(sepdep_log_path, 'w') as f:
            for line in sepdep_log:
                f.write(line + '\n')
    except Exception as e:
        return messagebox.showerror('Error', 'Error while saving log file: ' + str(e))    

    print('SEP/DEP correction done.')