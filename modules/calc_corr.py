# functions defined: switch, time conversion, upload of gamma library, call and clear all corrections, weights upload
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

# TODO if compounds are in parameters and elements or fraction strings contain spaces, the uploaded data are incorrect

# Libraries
import os
import sys
import pathlib
from datetime import datetime
import tkinter as tk
# Local modules
from modules import globals
from modules import corr_nonlin
from modules import corr_eff
from modules import corr_bcg
from modules import corr_sepdep
from modules import corr_att
from modules import utilities
from modules import calc_rr

# Directory
os.chdir(os.path.dirname(sys.argv[0]))


# Functions that call corrections
def corrections(root, nonlin_cor, bcg_cor, sd_cor, att_cor, eff_cor, non_point_cor):
    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')

    # call nonlinearity correction
    if nonlin_cor:
        if 'Cal en_1' and 'Cal en_2' not in globals.spectra_param_collection:
            return tk.messagebox.showerror('User input error', 'Calibration energies were not uploaded.')
        else:
            corr_nonlin.nonlinearity(globals.spectra_param_collection, globals.spectra_collection, '.prn')

    # call background correction
    if bcg_cor:
        corr_bcg.background()

    # call SEP/DEP correction
    if sd_cor:
        corr_sepdep.sepdep()

    # call attenuation correction
    if att_cor:
        if 'Thickness' and 'Element' not in globals.spectra_param_collection:
            return tk.messagebox.showerror('User input error', 'Information about samples material and thickness were not uploaded.')
        else:
            corr_att.attenuation()

    # call efficiency correction
    if eff_cor:
        corr_eff.efficiency()

    # call non-point correction
    if non_point_cor:
        nonpoint_cor()
        
    # Ask if the user wants to open the calculation window directly
    if tk.messagebox.askyesno('Open Calculation Window', 'Corrections were successfully calculated. Do you want to open the calculation window?'):
        calc_rr.result_win(root)


def corrections_win(root):
    # check for input files
    try:
        if not len(globals.spectra_param_collection['Start of measurement']) > 0:
            raise Exception('No prn. file(s) was(were) loaded')
    except Exception as e:
        print(e)
        tk.messagebox.showerror('User input error', 'Please first choose an input file(s) directory')
        return False

    # new window
    n_win = tk.Toplevel(root)
    n_win.title('User input')
    n_win.columnconfigure(0, weight=1)

    # block a background window
    n_win.grab_set()
    n_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, n_win))

    # turn off resizing
    n_win.resizable(False, False)

    # tkinter variables and their default values
    nonlin_cor = tk.BooleanVar(value=True)
    bcg_cor = tk.BooleanVar(value=True)
    sd_cor = tk.BooleanVar(value=True)
    att_cor = tk.BooleanVar(value=True)
    eff_cor = tk.BooleanVar(value=True)
    non_point_cor = tk.BooleanVar(value=False)

    # checkbuttons frame
    check_frame = tk.ttk.LabelFrame(n_win, relief='groove', text='Select corrections')
    check_frame.grid(column=0, row=0, sticky='nswe', padx=5, pady=5)

    chk_nonlin = tk.ttk.Checkbutton(check_frame, text='Nonlinearity', var=nonlin_cor, state=corr_nonlin.check_cal_en())
    chk_nonlin.grid(column=0, row=1, sticky='nswe', padx=2, pady=2)
    chk_bcgcor = tk.ttk.Checkbutton(check_frame, text='Background', var=bcg_cor)
    chk_bcgcor.grid(column=0, row=2, sticky='nswe', padx=2, pady=2)
    chk_sepdep = tk.ttk.Checkbutton(check_frame, text='Escape peaks', var=sd_cor)
    chk_sepdep.grid(column=0, row=3, sticky='nswe', padx=2, pady=2)
    chk_attcor = tk.ttk.Checkbutton(check_frame, text='Attenuation', var=att_cor, state=corr_att.check_att())
    chk_attcor.grid(column=0, row=4, sticky='nswe', padx=2, pady=2)
    chk_effcor = tk.ttk.Checkbutton(check_frame, text='Efficiency', var=eff_cor)
    chk_effcor.grid(column=0, row=5, sticky='nswe', padx=2, pady=2)
    chk_nonpoi = tk.ttk.Checkbutton(check_frame, text='Non-point source', var=non_point_cor, state=check_nonpoint())
    chk_nonpoi.grid(column=0, row=6, sticky='nswe', padx=2, pady=2)

    # buttons frame
    button_frame = tk.ttk.Frame(n_win)
    button_frame.grid(sticky='wesn', column=0, row=1, padx=5, pady=5)
    button_frame.columnconfigure(0, weight=1)
    button_frame.columnconfigure(1, weight=1)
    
    button_corrections = tk.ttk.Button(button_frame, text='Calculate', command=lambda: [
        clear_corrections(),
        corrections(root, nonlin_cor.get(), bcg_cor.get(), sd_cor.get(), att_cor.get(), eff_cor.get(), non_point_cor.get()), utilities.close_win(root, n_win)])
    button_corrections.grid(column=0, row=0, sticky='nswe', padx=2, pady=2)

    button_quit = tk.ttk.Button(button_frame, text='Quit', command=lambda: utilities.close_win(root, n_win))
    button_quit.grid(column=1, row=0, sticky='nswe', padx=2, pady=2)


# non-point correction main function - simply take values from spectra_param_collection and put them into spectra_collection
def nonpoint_cor():
    names = list(globals.spectra_collection.keys())
    for i in range(0, len(names)):
        if 'non-point corr' not in globals.spectra_collection[names[i]].keys():
            globals.spectra_collection[names[i]]['non-point corr'] = ''
        
        for j in range(0, len(globals.spectra_param_collection['Name'])):
            if names[i] + '.prn' == str(globals.spectra_param_collection['Name'][j]):
                for k in range(0, len(globals.spectra_collection[names[i]]['energy'])):
                    globals.spectra_collection[names[i]].loc[k, 'non-point corr']  = globals.spectra_param_collection.loc[j, 'Non-point']                        


# checks whether calibration energies are uploaded and values are nonzero
def check_nonpoint():
    if 'Non-point' in globals.spectra_param_collection:
        return 'active'
    else:
        return 'disable'


# clears previously made corrections
def clear_corrections():
    for spectrum in globals.spectra_collection:
        if 'cor_en' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['cor_en']
        if 'en_diff' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['en_diff']
        if 'cor_en_err' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['cor_en_err']
        if 'eff fep' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['eff fep']
        if 'att cor' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['att cor']
        if 'bcg cor area' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['bcg cor area']
        if 'bcg cor area unc' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['bcg cor area unc']
        if 'SEP-DEP cor area' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['SEP-DEP cor area']
        if 'SEP-DEP correction' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['SEP-DEP correction']
        if 'non-point corr' in globals.spectra_collection[spectrum]:
            del globals.spectra_collection[spectrum]['non-point corr']


# read/write functions for parameters.txt file --------------------------------------------------------------------------------------------

# reads parameters.txt and adds values to spectra_param_collection and common_params (can read all parameters)
def upload_info_from_file():
    # check for input files
    try:
        if not len(globals.spectra_param_collection['Start of measurement']) > 0:
            raise Exception('No prn. file(s) was(were) loaded')
    except Exception:
        tk.messagebox.showerror('User input error', 'Please first choose an input file(s) directory')
        return False

    globals.spectra_param_collection = globals.spectra_param_collection.sort_values(
        by=['Sample name', 'Start of measurement'])

    with open((globals.common_params['Folder path'] / pathlib.Path('parameters.txt')), 'a+', encoding='utf-8') as temp_file:
        temp_file.seek(0)

        # check whether is parameters file empty or not
        if os.stat(globals.common_params['Folder path'] / pathlib.Path('parameters.txt')).st_size == 0:
            return

        content = temp_file.readlines()

        if content:
            weight_line = None
            nonpoint_line = None
            detect_line = None

            # upload MEASUREMENT INFO and get lines with WEIGHT INFO, DETECTOR INFO and NON-POINT INFO if are in the file
            # for i in range(0, len(content)):
            for i, line in enumerate(content):
                if 'MEASUREMENT INFO' in line:
                    globals.common_params['End irr'] = content[i+1].split()[3] + ' ' + content[i+1].split()[4]

                    globals.common_params['Irradiation time'] = float(content[i+2].strip().split()[3])
                    globals.common_params['Beam int'] = float(content[i+3].split()[5])
                    globals.common_params['Beam int unc'] = float(content[i+3].split()[7])

                    # add delayed time to the spectra_param_collection
                    d_time = []
                    end_irr = datetime.strptime(globals.common_params['End irr'], '%d/%m/%Y %H:%M:%S')
                    for t in globals.spectra_param_collection['Start of measurement']:
                        d_time.append((t - end_irr).total_seconds())
                    globals.spectra_param_collection['Delayed time'] = d_time

                elif 'SAMPLE INFO' in line:
                    weight_line = i
                elif 'DETECTOR INFO' in line:
                    detect_line = i
                elif 'NON-POINT INFO' in line:
                    nonpoint_line = i

            # upload WEIGHT INFO
            name = []
            weight = []
            weight_err = []
            molar = []
            material = []
            fraction = []
            density = []
            thickness = []

            if weight_line is not None:
                for i in range(weight_line + 2, len(content)):
                    if len(content[i].strip()) == 0 or content[i+1].startswith('#'):
                        break
                    else:
                        line = content[i].split()
                        name.append(pathlib.Path(line[0]))
                        weight.append(float(line[2]))
                        weight_err.append(float(line[3]))
                        molar.append(float(line[4]))
                        material.append(line[5])
                        fraction.append(line[6])
                        density.append(line[7])
                        thickness.append(float(line[8]))
                    

                globals.spectra_param_collection['Weight'] = globals.spectra_param_collection['Name'].map(dict(zip(name, weight)))
                globals.spectra_param_collection['Weight error'] = globals.spectra_param_collection['Name'].map(dict(zip(name, weight_err)))
                globals.spectra_param_collection['Molar mass'] = globals.spectra_param_collection['Name'].map(dict(zip(name, molar)))
                globals.spectra_param_collection['Element'] = globals.spectra_param_collection['Name'].map(dict(zip(name, material)))
                globals.spectra_param_collection['Fraction by weight'] = globals.spectra_param_collection['Name'].map(dict(zip(name, fraction)))
                globals.spectra_param_collection['Density'] = globals.spectra_param_collection['Name'].map(dict(zip(name, density)))
                globals.spectra_param_collection['Thickness'] = globals.spectra_param_collection['Name'].map(dict(zip(name, thickness)))

            # upload DETECTOR INFO
            detector = []
            cal_1 = []
            cal_2 = []
            if detect_line is not None:
                for i in range(detect_line + 2, len(content)):
                    line = content[i].split()
                    detector.append(line[0])
                    cal_1.append(float(line[1]))
                    cal_2.append(float(line[2]))
                    # if empty line or comment line is found, break the loop
                    if len(content[i+1].strip()) == 0 or content[i+1].startswith('#'):
                        break

                globals.spectra_param_collection['Cal en_1'] = globals.spectra_param_collection['Detector'].map(dict(zip(detector, cal_1)))
                globals.spectra_param_collection['Cal en_2'] = globals.spectra_param_collection['Detector'].map(dict(zip(detector, cal_2)))

            # upload NON-POINT INFO
            name = []
            non_point = []

            if nonpoint_line is not None:
                for i in range(nonpoint_line + 2, len(content)):
                    if len(content[i].strip()) == 0 or content[i].startswith('#'):
                        break
                    else:
                        line = content[i].split()
                        name.append(pathlib.Path(line[0]))
                        non_point.append(float(line[1]))

                globals.spectra_param_collection['Non-point'] = globals.spectra_param_collection['Name'].map(dict(zip(name, non_point)))


# a writes all information into the parameters file
def write_all_info():
    # check if paramaters file exists and if not, create it + check if is not empty in case of existing file
    if not (globals.common_params['Folder path'] / pathlib.Path('parameters.txt')).exists():
        print('parameters.txt file does not exist and will be created.')

    # new_content is a list of strings which will be written into the file
    new_content = []

    # write MEASUREMENT INFO
    new_content.append('# MEASUREMENT INFO\n')
    new_content.append('End of irradiation: {}\n'.format(globals.common_params['End irr']))
    new_content.append('Duration of irradiation: {} min\n'.format(globals.common_params['Irradiation time']))
    new_content.append('Integral beam intensity and uncertainty: {0:1.2E} +- {1:1.2E} \n'.format(
        float(globals.common_params['Beam int']),
        float(globals.common_params['Beam int unc'])
    ))
    new_content.append('\n\n')

    # temporary solution how to fix problem with itertuples column names
    original_columns = globals.spectra_param_collection.columns
    # replace columns in spectra parameter with new columns without spaces and dashes
    globals.spectra_param_collection.columns = [col.replace(' ', '_').replace('-', '_') for col in globals.spectra_param_collection.columns]

    # write SAMPLE INFO
    # check if some sample information (weight, weight err, etc.) are in spectra_param_collection
    if 'Weight' in globals.spectra_param_collection:
        new_content.append('# SAMPLE INFO\n')
        new_content.append('# Name\t\t\tDet\tWeight (g)\tWeight err (g)\tMolar mass\tMaterial\tFraction wt.\tDensity (g/cm^3)\tThickness (mm)\n')
        for row in globals.spectra_param_collection.itertuples(index=False):
            new_content.append(' {}\t{}\t{:.4f}\t{:.4f}\t\t\t{:.4f}\t{}\t\t{}\t\t\t\t{}\t\t\t\t\t{:.4f}\n'.format(row.Name, row.Detector, row.Weight, row.Weight_error, row.Molar_mass, row.Element, row.Fraction_by_weight, row.Density, row.Thickness))
        new_content.append('\n\n')
    
    # write DETECTOR INFO
    # check if some detector information (cal en_1, cal en_2) are in spectra_param_collection
    if 'Cal_en_1' in globals.spectra_param_collection:
        new_content.append('# DETECTOR INFO\n')
        new_content.append('# Detector\tCal en_1 (keV)\tCal en_2 (keV)\n')

        # new version of writing detector info
        temp_content = []
        for row in globals.spectra_param_collection.itertuples(index=False):
            temp_content.append(' {}\t\t{:.4f}\t\t\t\t{:.4f}\n'.format(row.Detector, row.Cal_en_1, row.Cal_en_2))
        # delete duplicates and sort the list
        temp_content = sorted(set(temp_content))
        
        new_content.extend(temp_content)
        new_content.append('\n\n')
   
    # write NON-POINT INFO
    # check if some non-point information are in spectra_param_collection
    if 'Non_point' in globals.spectra_param_collection:
        new_content.append('# NON-POINT INFO\n')
        new_content.append('# Name\t\t\tNon-point correction\n')
        for row in globals.spectra_param_collection.itertuples(index=False):
            new_content.append(' {}\t{:.4f}\n'.format(row.Name, row.Non_point))
        new_content.append('\n\n')

    # return original columns to the spectra_param_collection
    globals.spectra_param_collection.columns = original_columns

    # write data into the file
    with open(globals.common_params['Folder path'] / pathlib.Path('parameters.txt'), 'w', encoding='utf-8') as temp_file:
        for line in new_content:
            temp_file.writelines(line)


# begining of the library functions ------------------------------------------------------------------------------------------------

# opens the dialogue for the new library, appends a library to the list of libraries
def ask_for_file_name(label_lib, menu, sel):
    if label_lib == 'Attenuation':
        a = tk.filedialog.askopenfilename(filetypes=(("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Efficiency':
        a = tk.filedialog.askopenfilename(filetypes=(("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Nonlinearity':
        a = tk.filedialog.askopenfilename(filetypes=(("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Beam Fluctuation':
        a = tk.filedialog.askopenfilename(filetypes=(("Text file", "*.txt"), ("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Gamma lib. (all)':
        a = tk.filedialog.askopenfilename(filetypes=(("Gamma file", "*.llb"), ("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Activity':
        a = tk.filedialog.askopenfilename(filetypes=(("Gamma file", "*.llb"), ("Library file", "*.lib"), ("All Files", "*.*")))
    elif label_lib == 'Gamma lib. (eff)':
        a = tk.filedialog.askopenfilename(filetypes=(("Gamma file", "*.llb"), ("Library file", "*.lib"), ("All Files", "*.*")))
    
    globals.lib_collection[label_lib].append(a)
    # 
    original_library_selection(label_lib, menu, sel)


# selects the original library from libraries folder, adds  the option for the upload of the new file
# adds all libraries as radiobuttons of the libraries menu

def original_library_selection(label_lib, menu, sel):
    menu.delete(0, len(globals.lib_collection[label_lib]))

    # sets last inserted library as selected one
    if sel.get() == '':
        globals.lib_selection[label_lib] = globals.lib_collection[label_lib][-1]
        sel.set(globals.lib_collection[label_lib][-1])
    else:
        globals.lib_selection[label_lib] = sel.get()
        sel.set(globals.lib_collection[label_lib][-1])
    for i in range(0, len(globals.lib_collection[label_lib])):
        menu.add_radiobutton(label=globals.lib_collection[label_lib][i], value=globals.lib_collection[label_lib][i], variable=sel)

    menu.add_command(label='Select file...', command=lambda: ask_for_file_name(label_lib, menu, sel))
    sel.trace_info()
    globals.lib_selection[label_lib] = sel.get()


# creates StringVar variable for every detector in input folder
# needed for tracking the user selection
def track_variable_create():
    if 'Name' in globals.spectra_param_collection:
        detectors = tuple(set(globals.spectra_param_collection['Detector'].tolist()))
        tracking_variables = [tk.StringVar() for _ in range(0, len(detectors))]
        globals.coinc_track_var = tracking_variables


# opens window which asks for library file
def select_coinc_library(menu, detector):
    a = tk.filedialog.askopenfilename(filetypes=(("Coincidence file", "*.clib"), ("All Files", "*.*")))
    globals.coinc_lib_collection[detector].append(a)
    coinc_lib_radiobutton(menu, detector)


# creates radiobutton options with files for specified detector
def coinc_lib_radiobutton(menu, detector):
    if detector in globals.coinc_lib_collection:
        menu.delete(0, len(globals.coinc_lib_collection[detector]))
        globals.coinc_lib_selection[detector] = globals.coinc_lib_collection[detector][-1]
        for i in range(0, len(globals.coinc_lib_collection[detector])):
            globals.coinc_track_var[i].set(globals.coinc_lib_selection[detector])
            menu.add_radiobutton(label=globals.coinc_lib_collection[detector][i], value=
            globals.coinc_lib_collection[detector][i], variable=globals.coinc_track_var[i])
            globals.coinc_track_var[i].trace_info()
            globals.coinc_lib_selection[detector] = globals.coinc_track_var[i].get()
    else:
        globals.coinc_lib_collection[detector] = []
    menu.add_command(label='Select file...', command=lambda: select_coinc_library(menu, detector))


# creates the optionmenu for every detector from uploaded spectra
def coinc_library_selection(menu):
    detectors = []
    if 'Name' in globals.spectra_param_collection:
        detectors = sorted(tuple(set(globals.spectra_param_collection['Detector'].tolist())))
        for i in range(0, len(detectors)):
            sub_menu = (tk.Menu(menu, tearoff=0))
            menu.add_cascade(label='Detector ' + detectors[i], menu=sub_menu,
                             command=coinc_lib_radiobutton(sub_menu, detectors[i]))


# for background
# needed for tracking the user selection
def bcg_track_variable_create():
    if 'Name' in globals.spectra_param_collection:
        detectors = tuple(set(globals.spectra_param_collection['Detector'].tolist()))
        tracking_variables = [tk.StringVar() for _ in range(0, len(detectors))]
        globals.bcg_track_var = tracking_variables


# opens window which asks for library file
def select_bcg_library(menu, detector):
    a = tk.filedialog.askopenfilename(filetypes=(("Background file", "*.bcg"), ("All Files", "*.*")))
    globals.bcg_lib_collection[detector].append(a)
    bcg_lib_radiobutton(menu, detector)


# creates radiobutton options with files for specified detector
def bcg_lib_radiobutton(menu, detector):
    if detector in globals.bcg_lib_collection:
        menu.delete(0, len(globals.bcg_lib_collection[detector]))

        if len(globals.bcg_lib_collection[detector]) > 0:
            globals.bcg_lib_selection[detector] = globals.bcg_lib_collection[detector][-1]
        else:
            globals.bcg_lib_selection[detector] = ''

        for i in range(0, len(globals.bcg_lib_collection[detector])):
            globals.bcg_track_var[i].set(globals.bcg_lib_selection[detector])
            menu.add_radiobutton(label=globals.bcg_lib_collection[detector][i], value=globals.bcg_lib_collection[detector][i], variable=globals.bcg_track_var[i])
            globals.bcg_track_var[i].trace_info()
            globals.bcg_lib_selection[detector] = globals.bcg_track_var[i].get()
    else:
        globals.bcg_lib_collection[detector] = []
    
    menu.add_command(label='Select file...', command=lambda: select_bcg_library(menu, detector))


# creates the optionmenu for every detector from uploaded spectra
def bcg_library_selection(menu):
    # clear menu
    menu.delete(0, len(globals.bcg_lib_collection))
    
    # update menu
    detectors = []
    if 'Name' in globals.spectra_param_collection:
        detectors = sorted(tuple(set(globals.spectra_param_collection['Detector'].tolist())))
        for i in range(0, len(detectors)):
            sub_menu = (tk.Menu(menu, tearoff=0))
            menu.add_cascade(label='Detector ' + detectors[i], menu=sub_menu, command=bcg_lib_radiobutton(sub_menu, detectors[i]))
  

# end of library functions --------------------------------------------------------------------------------------------
