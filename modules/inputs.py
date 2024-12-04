# Functionality: Functions for GUI windows and data upload from input files
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
import os
import sys
import locale
import pandas as pd
import tkinter as tk
import pathlib
from datetime import datetime
import numpy as np
# Local modules
from modules import globals
from modules import settings
from modules import utilities
from modules import calc_corr

os.chdir(os.path.dirname(sys.argv[0]))


# open and read folder with input files
def open_folder(coinc_menu, bcg_menu, treeview, main_menu):
    folder_path = pathlib.Path(tk.filedialog.askdirectory(title='Choose directory with input files', initialdir=pathlib.Path.cwd()))

    if not folder_path.exists():
        return tk.messagebox.showerror('User input error', 'No directory was selected.')

    files_rr = []

    bcg_files = []

    # make a list of spectra
    for file in list(folder_path.iterdir()):
        if file.is_file():
            if file.suffix == '.prn':
                files_rr.append(pathlib.Path(file.name))
            elif file.suffix == '.bcg':
                bcg_files.append(pathlib.Path(file.stem))

    # add folder path to common parameters
    globals.common_params['Folder path'] = pathlib.Path(folder_path)

    # get relative path to working directory (relative to main python file)
    globals.common_params['Relative path'] = pathlib.Path(os.path.relpath(folder_path, os.path.dirname(sys.argv[0])))
    globals.common_params['Output path'] = globals.common_params['Relative path'] / pathlib.Path('output')
    globals.common_params['Spectra path'] = globals.common_params['Output path'] / pathlib.Path('spectra')
    globals.common_params['Results path'] = globals.common_params['Output path'] / pathlib.Path('results')
    
    print('\nRelative path to working directory: {}'.format(globals.common_params['Relative path']))
    print('Output path: {}'.format(str(globals.common_params['Output path'])))
    print('Processed spectra path: {}'.format(str(globals.common_params['Spectra path'])))
    print('Results path: {}\n'.format(str(globals.common_params['Results path'])))

    # read and load all spectra from input files folder + read live, real time and start of measurement
    if len(files_rr) > 0:
        globals.spectra_param_collection = pd.DataFrame(files_rr, columns=['Name'])
        spectrum = fill_data(folder_path, globals.spectra_collection, globals.spectra_param_collection)
        
        if spectrum is not None:
            return tk.messagebox.showerror('Input format error', 'The input file contains more than one spectrum or parameters section. Please check the file: {}'.format(spectrum.name))
        
        # extract information from file names
        data_from_name(coinc_menu, bcg_menu, treeview, bcg_files)
    else:
        return tk.messagebox.showerror('User input error', 'No spectrum files were found in the selected directory.')
       
    # change state of Menu: Open selected spectrum
    main_menu.entryconfig(2, state='normal')


# open and read all spectra from input files
def fill_data(folder_path, dict_spectra, df_parameters):      # made like this because of the need to read standard and efficiency files
    
    meas_start = []
    l_time = []
    r_time = []

    # read and load all spectra from input files folder + read live, real time and start of measurement
    for file in df_parameters['Name']:

        with open(folder_path / file, 'r', encoding='utf-8') as temp_file:
            content = temp_file.readlines()

        # check content if there are any data twice or more
        count_param = 0
        count_data = 0
        for line in content:
            if len(line.strip()) != 0:
                if 'Filename:' in line:
                    count_param += 1
                elif 'chisq' in line or 'cherr' in line:
                    count_data += 1
        if count_data > 1 or count_param > 1:
            return file

        for row in content:
            line = row.split()
            default = locale.getdefaultlocale(locale.LC_TIME)   # solve problem with local time and date format + abreveations
            if len(line) > 0 and line[0] == 'Measurement':
                # date upload (english or czech)
                try:
                    try:
                        meas_start.append(datetime.strptime(line[2] + ' ' + line[3], '%d/%b/%y %H:%M'))
                    except Exception as e:
                        print(e)
                        locale.setlocale(locale.LC_TIME, 'cs_CZ')
                        meas_start.append(datetime.strptime(line[2] + ' ' + line[3], '%d/%b/%y %H:%M'))
                        locale.setlocale(locale.LC_TIME, '')
                except Exception as e:
                    tk.messagebox.showerror('Input format error', 'Date and time format is not correct. Please check the input file.')
                    print(e)
                    locale.setlocale(locale.LC_TIME, '')
                    return
            if len(line) > 0 and line[0] == 'Live':
                l_time.append(float(line[2]))
                r_time.append(float(line[6]))
                break
        
        # find beginning of the spectrum data
        for j in range(0, len(content)):
            line = content[j].split()
            next_line = content[j+1].split()
            if len(line) > 0 and len(next_line) > 0 and line[0].startswith('===='):
                d = j + 1
                break
        
        # find end of the spectrum data
        i = 0
        for row in content[d:]:
            if len(row) > 0:
                i += 1
            else:
                break

    # read spectra to dictionary where spectrum name is key and pandas dataframe is value
        dict_spectra.update({file.stem: pd.read_csv(folder_path / file, sep=r'\s+', skiprows=d, nrows=i, header=None, engine='python',
                                                  names=['channel', 'cherr', 'energy', 'area', 'aerr%', 'fwhm', 'chisq', 'it', 'left', 'n', 'lim', 'isotope'])})
        # drop columns 'it', 'left', 'n', 'lim'
        dict_spectra[file.stem] = dict_spectra[file.stem].drop(columns=['it', 'left', 'n', 'lim'])
    
        # recalculate FWHM from channel to energy
        # dict_spectra[file.stem]['fwhm_keV'] = dict_spectra[file.stem]['fwhm'] * dict_spectra[file.stem]['energy'] / dict_spectra[file.stem]['channel']

    # fill parameters to parameters dataframe
    df_parameters['Start of measurement'] = meas_start
    df_parameters['Live time'] = l_time
    df_parameters['Real time'] = r_time


# read and fill spectra parameters from input files
def data_from_name(coinc_menu, bcg_menu, treeview, bcg_files):

    det = []
    pos = []
    name = []
    meas_num = []

    # read and load all spectra from input files folder + read live, real time and start of measurement
    for file in globals.spectra_param_collection['Name']:
        name_divided = measurement_name(file.stem)

        det.append(name_divided[0])
        pos.append(name_divided[1])
        name.append(name_divided[2])
        meas_num.append(name_divided[3])

    globals.spectra_param_collection['Detector'] = det
    globals.spectra_param_collection['Position'] = pos
    globals.spectra_param_collection['Sample name'] = name
    globals.spectra_param_collection['Measurement number'] = meas_num
    
    # sort spectra_param_collection by Sample name and Start of measurement
    globals.spectra_param_collection = globals.spectra_param_collection.sort_values(by=['Sample name', 'Start of measurement'])

    # read parameters.txt file
    calc_corr.upload_info_from_file()

    # FILL libraries options in menu

    # create track variables for coincidence and background correction
    calc_corr.track_variable_create()
    # create Menus for each library
    calc_corr.coinc_library_selection(coinc_menu)
    # track library selection from menu
    calc_corr.bcg_track_variable_create()
    # create Menus for each library
    calc_corr.bcg_library_selection(bcg_menu)
    
    # fill bcg lib collection with paths from input files
    for i in range(0, len(bcg_files)):
        bcg_files[i] = str(bcg_files[i])
    
    for file in bcg_files:
        for key in globals.bcg_lib_collection.keys():
            if file == key:
                globals.bcg_lib_collection[file] = [str(globals.common_params['Folder path'] / pathlib.Path(file + '.bcg'))]
    
    # update menu with new values
    calc_corr.bcg_library_selection(bcg_menu)
    
    # update treeview with new values
    treeview_fill(treeview)


# window with information from user input
def measur_info_win(root):

    # check for input files
    try:
        if not len(globals.spectra_param_collection['Start of measurement']) > 0:
            raise Exception('No prn. file(s) was(were) loaded')
    except Exception as e:
        print(e)
        tk.messagebox.showerror('User input error', 'Please first choose an input file(s) directory')
        return False

    # create a new window with information about input (experiment)
    exp_info_win = tk.Toplevel(root)
    exp_info_win.title('User input - irradiation info')

    # block a background window
    exp_info_win.grab_set()
    exp_info_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, exp_info_win))
    exp_info_win.resizable(False, False)

    directory_label = tk.ttk.Label(master=exp_info_win, text = 'Work directory: ' + str(globals.common_params['Folder path']))
    label_irr_time = tk.ttk.Label(master=exp_info_win, text = 'Duration of irradiation: ' + str(globals.common_params['Irradiation time']) + ' min')
    label_end_irr = tk.ttk.Label(master=exp_info_win, text = 'End of irradiation: ' + str(globals.common_params['End irr']))
    label_beam = tk.ttk.Label(master=exp_info_win, text = 'Integral beam intensity: {0:1.2E} \u00B1 {1:1.2E}'.format(float(globals.common_params['Beam int']), float(globals.common_params['Beam int unc'])))

    directory_label.grid(column=0, row=1, sticky='nswe', padx=5, pady=5)
    label_irr_time.grid(column=0, row=2, sticky='nswe', padx=5, pady=5)
    label_end_irr.grid(column=0, row=3, sticky='nswe', padx=5, pady=5)
    label_beam.grid(column=0, row=4, sticky='nswe', padx=5, pady=5)


# window for user information about measurement
def measurement_input_info_window(root, treeview):
    # check for input files
    try:
        if not len(globals.spectra_param_collection['Start of measurement']) > 0:
            raise Exception('No prn. file(s) was(were) loaded')
    except Exception as e:
        print(e)
        tk.messagebox.showerror('User input error', 'Please first choose an input file(s) directory')
        return False

    # new window
    exp_input_win = tk.Toplevel(root)
    exp_input_win.title('User input - irradiation parameters')

    # block a background window
    exp_input_win.grab_set()
    exp_input_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, exp_input_win))
    exp_input_win.resizable(False, False)
    exp_input_win.grid_propagate(True)

    # irradiation end
    end_irr_frame = tk.ttk.Frame(exp_input_win)
    end_irr_frame.grid(column=0, row=0, sticky='nswe', padx=5, pady=5)
    end_irr_frame.columnconfigure(0, weight=1)
    end_irr_label = tk.ttk.Label(end_irr_frame, text='Date and time of irradiation end:')
    end_irr_label.grid(column=0, columnspan=2 , row=0, sticky='nswe', padx=5)
    ent_end_irr = tk.ttk.Entry(end_irr_frame, width=20)
    ent_end_irr.grid(column=0, row=1, sticky='nswe', padx=5)
    end_irr_label_2 = tk.ttk.Label(end_irr_frame, text='(DD/MM/YYYY HH:MM:SS)', width=30)
    end_irr_label_2.grid(column=1, row=1, sticky='nswe', padx=5)

    ent_end_irr.insert(0, globals.common_params['End irr'])

    # irradiation time
    irr_time_frame = tk.ttk.Frame(exp_input_win)
    irr_time_frame.grid(column=0, row=1, sticky='nswe', padx=5, pady=5)
    irr_time_label = tk.ttk.Label(irr_time_frame, text='Duration of irradiation:')
    irr_time_label.grid(column=0, columnspan=2, row=0, sticky='nswe', padx=5)
    ent_time_irr = tk.ttk.Entry(irr_time_frame, width = 20)
    ent_time_irr.grid(column=0, row=1, sticky='nswe', padx=5)
    irr_time_label_2 = tk.ttk.Label(irr_time_frame, text='minutes')
    irr_time_label_2.grid(column=1, row=1, sticky='nswe', padx=5)

    ent_time_irr.insert(0, globals.common_params['Irradiation time'])

    # beam intensity
    beam_frame = tk.ttk.Frame(exp_input_win)
    beam_frame.grid(column=0, row=2, sticky='nswe', padx=5, pady=5)
    beam_label_1 = tk.ttk.Label(beam_frame, text='Integral number of particles and uncertainty:')
    beam_label_1.grid(column=0, row=0,  columnspan=3, sticky='nswe', padx=5)
    beam_int_ent = tk.ttk.Entry(beam_frame, width=15)
    beam_int_ent.grid(column=0, row=1, sticky='nswe', padx=5)
    beam_label_2 = tk.ttk.Label(beam_frame, text='\u00B1', width=1)
    beam_label_2.grid(column=1, row=1, sticky='nswe', padx=5)
    beam_int_unc_ent = tk.ttk.Entry(beam_frame, width=15)
    beam_int_unc_ent.grid(column=2, row=1, sticky='nswe', padx=5)

    beam_int_ent.insert(0, '{0:1.2E}'.format(float(globals.common_params['Beam int'])))
    beam_int_unc_ent.insert(0, '{0:1.2E}'.format(float(globals.common_params['Beam int unc'])))

    # button
    button_frame = tk.ttk.Frame(exp_input_win)
    button_frame.grid(column=0, columnspan=2, row=4, sticky='nswe', padx=5, pady=5)
    button_frame.columnconfigure(0, weight=1)
    button_frame.columnconfigure(1, weight=1)
    
    load_button = tk.ttk.Button(button_frame, text='Save', width=20, command=lambda: [
        upload_meas_info(root, exp_input_win, ent_end_irr.get(), ent_time_irr.get(), beam_int_ent.get(), beam_int_unc_ent.get()),
        treeview_fill(treeview),
        calc_corr.write_all_info()
        ])
    load_button.grid(column=0, row=0, sticky='nswe', padx=2, pady=2)

    quit_button = tk.ttk.Button(button_frame, text='Quit', width=20, command=lambda: utilities.close_win(root, exp_input_win))
    quit_button.grid(column=1, row=0, sticky='nswe', padx=2, pady=2)


# Create a new window with values of selected spectrum
def open_spectrum(sel_file, root):
    if sel_file is None:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload or no selected spectrum.')

    # new window
    spectrum_win = tk.Toplevel(root)
    spectrum_win.title('Spectrum: ' + sel_file + '.prn')
    spectrum_win.geometry('800x350')
    spectrum_win.columnconfigure(0, weight=1)
    spectrum_win.rowconfigure(0, weight=1)
       
    # block background window
    spectrum_win.grab_set()
    spectrum_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, spectrum_win))

    # frames
    tree_frame = tk.ttk.Frame(spectrum_win)
    tree_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    tree_frame.columnconfigure(0, weight=1)
    tree_frame.rowconfigure(0, weight=1)

    button_frame = tk.ttk.Frame(spectrum_win)
    button_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)

    button_close = tk.ttk.Button(button_frame, text='Close window', width=19, command=lambda: utilities.close_win(root, spectrum_win))
    button_close.grid(sticky='w', column=1, row=1)
    button_safe = tk.ttk.Button(button_frame, text='Save to XLSX', width=19, command=lambda: safe_spectrum(sel_file))
    button_safe.grid(sticky='w', column=0, row=1)

    # treeview
    spectrum_treeview = tk.ttk.Treeview(tree_frame)
    spectrum_treeview['show'] = 'headings'
    spectrum_treeview.grid(sticky='wens', column=0, columnspan=5, row=0, rowspan=5)

    # Treeview X-scrollbar
    tree_x_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal')
    tree_x_scr.grid(sticky='we', column=0, row=5, columnspan=5)
    tree_x_scr.configure(command=spectrum_treeview.xview)
    spectrum_treeview.configure(xscrollcommand=tree_x_scr.set)
    # Treeview Y-scrollbar
    tree_y_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=spectrum_treeview.yview)
    tree_y_scr.grid(sticky='ns', column=5, row=0, rowspan=5)
    spectrum_treeview.configure(yscrollcommand=tree_y_scr.set)

    spectrum_treeview['columns'] = list(globals.spectra_collection[sel_file].columns)
    
    # set columns filters
    for col in spectrum_treeview['columns']:
        spectrum_treeview.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(spectrum_treeview, _col, False))

    for col_name in list(globals.spectra_collection[sel_file].columns):
        spectrum_treeview.column(col_name, width=70, stretch=True, anchor='center')
        spectrum_treeview.heading(col_name, text=col_name)
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in spectrum_treeview['columns']:
            spectrum_treeview.column(key, width=value, stretch=True)

    for i in range(len(globals.spectra_collection[sel_file])):      # fill treeview with new values
        spectrum_treeview.insert('', index='end', values=globals.spectra_collection[sel_file].iloc[i, :].tolist())


# clear and fill treeview with new values
def treeview_fill(treeview):

    treeview['columns'] = list(globals.spectra_param_collection.columns)   # get a header from pandas dataframe

    for col_name in list(globals.spectra_param_collection.columns):
        treeview.column(col_name, width=150, stretch=False, anchor='center')
        treeview.heading(col_name, text=col_name)

    x = treeview.get_children()             # get id of all items in treeview
    for i in x:                             # delete all items
        treeview.delete(i)

    for i in range(len(globals.spectra_param_collection)):      # fill treeview with new values
        treeview.insert('', index='end', values=globals.spectra_param_collection.iloc[i, :].tolist())

    # set columns filters
    for col in treeview['columns']:
        treeview.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview, _col, False))

    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview['columns']:
            treeview.column(key, width=value, stretch=False)


# save one spectrum to .xlsx file
def safe_spectrum(spec_name):
    # check if output directory exists if not than create it
    utilities.check_folders()

    # get current date and use it in the name of the file (no time)
    now = datetime.now()
    dt_string = now.strftime('%Y-%m-%d')
    
    try:
        globals.spectra_collection[spec_name].to_excel(globals.common_params['Spectra path'] / pathlib.Path(dt_string + '_' + spec_name + '.xlsx'), index=None, header=True)
    except Exception as e:
        return tk.messagebox.showerror('User input error', 'The file was not saved. Error: ' + str(e))
    
    tk.messagebox.showinfo(title='Notification', message='Selected spectrum was saved.')


# save dataframe with all spectra to .xlsx file
def save_all():
    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload. No data is ready to save.')

    # check if output directory exists if not than create it
    utilities.check_folders()

    # get current date and use it in the name of the file (no time)
    now = datetime.now()
    dt_string = now.strftime('%Y-%m-%d')

    # save all spectra to .xlsx file
    try:
        writer = pd.ExcelWriter(globals.common_params['Spectra path'] / pathlib.Path(dt_string + '_all_spectra.xlsx'), engine='xlsxwriter')
    
        for spectrum in globals.spectra_collection:
            globals.spectra_collection[spectrum].to_excel(writer, sheet_name=spectrum, index = None, header=True)   

        writer.close()
    except Exception as e:
        return tk.messagebox.showerror('User input error', 'The file was not saved. Error: ' + str(e))

    # save parameters of all spectra to .xlsx file
    try:
        globals.spectra_param_collection.to_excel(globals.common_params['Output path'] / pathlib.Path(dt_string + '_all_spectra_parameters.xlsx'), index=None, header=True)
    except Exception as e:
        return tk.messagebox.showerror('User input error', 'The file was not saved. Error: ' + str(e))
    
    tk.messagebox.showinfo(title='Notification', message='Data were saved.')


# get neccessary values for att from user input to spectra_param_collection and upload_thickness new window
def upload_thickness(root, win, samples, dictionary):
       
    thick = []
    weight_f = []
    dens = []
    elem = []

    for s, e, w, d, t in zip(samples, dictionary['Elements'], dictionary['Weight frac'], dictionary['Density'], dictionary['Thickness']):
        # is thickness number?
        try:
            th = float(t)
        except:
            th = 0
            tk.messagebox.showerror('User input error', 'The value inserted for ' + s + ' is not a number. Zero was inserted instead')
        
        if w != '-':
            # control if compound: are fraction by weights numbers?
            try:
                control = w.split('/')
                [float(a.strip()) for a in control]
            except:
                tk.messagebox.showerror('User input error', 'The value inserted for ' + s + ' is not a number. Zero was inserted instead')
            # control if compound: is density number?
            try:
                de = float(d)
            except:
                de = 0
                tk.messagebox.showerror('User input error', 'The value inserted for ' + s + ' is not a number. Zero was inserted instead')
        else:
            de = '-'
        for i in range(0, len(globals.spectra_param_collection['Sample name'])):
            if globals.spectra_param_collection.loc[i, 'Sample name'] == s:
                elem.append(e)
                weight_f.append(w)
                dens.append(de)
                thick.append(th)

    globals.spectra_param_collection['Element'] = elem
    globals.spectra_param_collection['Fraction by weight'] = weight_f
    globals.spectra_param_collection['Density'] = dens
    globals.spectra_param_collection['Thickness'] = thick

    utilities.close_win(root, win)


# new window for information about samples calibration energies
def new_win_det(root, treeview_files):

    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')

    detectors = sorted(tuple(set(globals.spectra_param_collection['Detector'].tolist())))

    calib1 = []
    calib2 = []

    # new window
    new_win_w = tk.Toplevel(root)
    new_win_w.title('User input - en. calibration')

    # block a background window
    new_win_w.grab_set()
    new_win_w.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win_w))

    # turn off resizing
    new_win_w.resizable(False, False)

    full_frame = tk.ttk.LabelFrame(new_win_w, relief='groove', text='User input - Deimos32 calibration energy')
    full_frame.grid(sticky='wesn', column=0, row=0, padx=5, pady=5)
    #full_frame.rowconfigure(0, weight=1)
    #full_frame.columnconfigure(0, weight=1)

    cal_en_1_lab = tk.ttk.Label(full_frame, text='Cal. en. 1')
    cal_en_1_lab.grid(sticky='W', column=1, row=0)
    cal_en_2_lab = tk.ttk.Label(full_frame, text='Cal. en. 2')
    cal_en_2_lab.grid(sticky='W', column=2, row=0)

    cal_en_1__lab = tk.ttk.Label(full_frame, text='(keV)')
    cal_en_1__lab.grid(sticky='W', column=1, row=1)
    cal_en_2__lab = tk.ttk.Label(full_frame, text='(keV)')
    cal_en_2__lab.grid(sticky='W', column=2, row=1)

    # uploads user entries to the dictionary and to the file
    for i in range(0, len(detectors)):
        d = detectors[i]
        lab = tk.ttk.Label(full_frame, text=d + ' : ', width=5, anchor='e')
        lab.grid(sticky='we', column=0, row=i + 2)

        ent_cal1 = tk.ttk.Entry(full_frame, width=10)
        ent_cal1.grid(sticky='we', column=1, row=i + 2, padx=1, pady=1)

        ent_cal2 = tk.ttk.Entry(full_frame, width=10)
        ent_cal2.grid(sticky='we', column=2, row=i + 2, padx=1, pady=1)

        for j in range(0, len(globals.spectra_param_collection['Detector'])):
            if globals.spectra_param_collection.loc[j, 'Detector'] == d:
            #if globals.spectra_param_collection['Detector'][j] == d:
                if 'Cal en_1' in globals.spectra_param_collection.keys():
                    ent_cal1.insert(0, globals.spectra_param_collection.loc[j, 'Cal en_1'])
                    ent_cal2.insert(0, globals.spectra_param_collection.loc[j, 'Cal en_2'])
                else:

                    ent_cal1.insert(0, 0)
                    ent_cal2.insert(0, 0)
                break

        calib1.append(ent_cal1)
        calib2.append(ent_cal2)

    button_frame = tk.ttk.Frame(new_win_w)
    button_frame.grid(sticky='wesn', column=0, row=1, padx=5, pady=5)
    button_frame.columnconfigure(0, weight=1)
    button_frame.columnconfigure(1, weight=1)

    button_save = tk.ttk.Button(button_frame, text='Save', command=lambda: [
        upload_cal_en(root, new_win_w, detectors, calib1, calib2),
        treeview_fill(treeview_files),
        calc_corr.write_all_info()
        ])
    button_save.grid(sticky='wesn', column=0, row=0, padx=2, pady=2)

    button_save = tk.ttk.Button(button_frame, text='Quit', command=lambda: utilities.close_win(root, new_win_w))
    button_save.grid(sticky='wesn', column=1, row=0, padx=2, pady=2)


# find name of detector, name of sample, number of measurement, position of measurement
def measurement_name(file):
    file = file.replace('_', '')

    det = file[0]       # detector letter
    i = -1              # counter for position number
    for letter in reversed(file):
        if letter == 'p':
            break
        else:
            i -= 1
    
    pos = file[i:]      # position number
    name = file[1:i]    # name with number of measurement

    x = 0              # counter for measurement number
    for letter in reversed(name):
        if letter.isdigit():
            x -= 1
        else:
            break

    name_nm = int(name[x:]) # measurement number
    name = name[:x]         # clear name without measurement number

    return det, pos, name, name_nm


# takes the data from window for calibration energy from user and adds them to the dictionary
def upload_cal_en(root, win, det, cal1, cal2):

    # check if calibration energies are already in the spectra_param_collection if not add them
    if 'Cal en_1' not in globals.spectra_param_collection.keys():
        globals.spectra_param_collection['Cal en_1'] = [0.0 for x in range(0, globals.spectra_param_collection.shape[0])]

    if 'Cal en_2' not in globals.spectra_param_collection.keys():
        globals.spectra_param_collection['Cal en_2'] = [0.0 for x in range(0, globals.spectra_param_collection.shape[0])]

    for d, c1, c2 in zip(det, cal1, cal2):
        try:
            c_1 = float(c1.get())
            c_2 = float(c2.get())
        except ValueError:
            tk.messagebox.showerror('User input error',
                                 'Some of inserted values is not a number. Zero was inserted instead')
            c_1 = 0.0
            c_2 = 0.0

        for i in range(0, globals.spectra_param_collection.shape[0]):
            if globals.spectra_param_collection.loc[i, 'Detector'] == d:
                globals.spectra_param_collection.loc[i, 'Cal en_1'] = c_1
                globals.spectra_param_collection.loc[i, 'Cal en_2'] = c_2

    utilities.close_win(root, win)


# check and upload user information about measurement
def upload_meas_info(root, win, end_ir, irr_tim, beam_in, beam_int_un):
    d_time = []
    end_irr = 0

    try:
        end_irr = datetime.strptime(end_ir, '%d/%m/%Y %H:%M:%S')
    except ValueError:
        tk.messagebox.showerror('User input error', 'Wrong date of the end of measurement')
    try:
        irr_time = float(irr_tim)
    except ValueError:
        tk.messagebox.showerror('User input error', 'Wrong time of irradiation')
    try:
        beam_int = float(beam_in)
    except ValueError:
        tk.messagebox.showerror('User input error', 'Wrong beam intensity')
    try:
        beam_int_unc = float(beam_int_un)
    except ValueError:
        tk.messagebox.showerror('User input error', 'Wrong beam intensity uncertainity')

    for t in globals.spectra_param_collection['Start of measurement']:
        d_time.append((t - end_irr).total_seconds())

    globals.spectra_param_collection['Delayed time'] = d_time

    globals.common_params['Irradiation time'] = irr_time
    globals.common_params['End irr'] = end_irr.strftime('%d/%m/%Y %H:%M:%S')
    globals.common_params['Beam int'] = beam_int
    globals.common_params['Beam int unc'] = beam_int_unc

    utilities.close_win(root, win)


# takes the data from window for weights from user and adds them to the dictionary
def upload_weights(root, win, samples, dictionary):

    wght = dictionary['Weight']
    wght_err = dictionary['Weight err']
    mola = dictionary['Molar']
    wei = []
    wei_er = []
    mol = []

    for s, we, we_er, mo in zip(samples, wght, wght_err, mola):
        try:

            w = float(we)
            w_r = float(we_er)
            m = float(mo)
        except:
            w = 0
            w_r = 0
            m = 0
            tk.messagebox.showerror('User input error',
                                 'Some value inserted for ' + s + ' is not a number. Zero was inserted instead')

        for i in range(0, len(globals.spectra_param_collection['Sample name'])):
            if globals.spectra_param_collection['Sample name'][i] == s:
                wei.append(w)
                wei_er.append(w_r)
                mol.append(m)

    globals.spectra_param_collection['Weight'] = wei
    globals.spectra_param_collection['Weight error'] = wei_er
    globals.spectra_param_collection['Molar mass'] = mol

    wght.clear()
    wght_err.clear()
    mola.clear()

    utilities.close_win(root, win)


# inserts data from non-point window into the dataframe
def get_nonpoint(sam, det, pos, non_point):

    if 'Non-point' not in globals.spectra_param_collection.keys():
        globals.spectra_param_collection['Non-point'] = np.nan
    
    for s, d, p, n in zip(sam, det, pos, non_point):
        try:
            correction = float(n.get())
        except Exception as e:
            print(e)
            correction = 1

        for i in range(0, len(globals.spectra_param_collection['Sample name'])):
            # redefine the condition with .loc method
            if (globals.spectra_param_collection.loc[i, 'Sample name'] == s) and (globals.spectra_param_collection.loc[i, 'Detector'] == d) and (globals.spectra_param_collection.loc[i, 'Position'] == 'p' + str(p)):
                globals.spectra_param_collection.loc[i, 'Non-point'] = correction


# new window for information about nonpoint correction
def new_win_non_point(root, treeview):
    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')

    # Local variables
    non_point = []
    sam = []
    det = []
    pos = []

    samples = sorted(tuple(set(globals.spectra_param_collection['Sample name'].tolist())))

    # new window
    new_win = tk.Toplevel(root)
    new_win.title('User input - non-point correction')

    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))

    # turn off resizing
    new_win.resizable(False, False)

    # frame for non-point correction
    full_frame = tk.ttk.LabelFrame(new_win, relief='groove', text='Non-point correction')
    full_frame.grid(sticky='wesn', column=0, columnspan=2, row=0, padx=5, pady=5)

    for i in range(0, len(samples)):
        sample = samples[i]
        lab = tk.ttk.Label(full_frame, text=sample)
        detectors_list = []
        for j in range(0, len(globals.spectra_param_collection['Sample name'])):
            #if globals.spectra_param_collection['Sample name'][j] == sample:
            if globals.spectra_param_collection.loc[j, 'Sample name'] == sample:
                detectors_list.append(globals.spectra_param_collection['Detector'][j])
                #detectors_list.append(globals.spectra_param_collection.loc[j, 'Detector'])
        detectors = sorted(tuple(set(detectors_list)))

        for k in range(0, len(detectors)):
            detector = detectors[k]
            det_frame = tk.ttk.LabelFrame(full_frame, text='Detector ' + detector )
            det_frame.grid(sticky='W', column=k + 1, row=i + 1, padx=2, pady=2)

            positions_list = []
            for j in range(0, len(globals.spectra_param_collection['Sample name'])):
                #if (globals.spectra_param_collection['Sample name'][j] == sample) and (globals.spectra_param_collection['Detector'][j] == detector):
                if (globals.spectra_param_collection.loc[j, 'Sample name'] == sample) and (globals.spectra_param_collection.loc[j, 'Detector'] == detector):
                    positions_list.append(int(globals.spectra_param_collection['Position'][j][1:]))
                    #positions_list.append(int(globals.spectra_param_collection.loc[j, 'Position'][1:]))
            positions = sorted(tuple(set(positions_list)))
            n = 1
            for l in range(0, len(positions)):
                position = positions[l]
                n += 2
                pos_lab = tk.ttk.Label(det_frame, text='p' + str(position))
                pos_ent = tk.ttk.Entry(det_frame, width=10)
                pos_lab.grid(sticky='W', column=n, row=i + 1, padx=2, pady=2)
                pos_ent.grid(sticky='W', column=n + 1, row=i + 1, padx=2, pady=2)

                for m in range(0, len(globals.spectra_param_collection['Sample name'])):
                    #if (globals.spectra_param_collection['Sample name'][m] == sample) and (globals.spectra_param_collection['Detector'][m] == detector) and (globals.spectra_param_collection['Position'][m] == 'p' + str(position)):
                    if (globals.spectra_param_collection.loc[m, 'Sample name'] == sample) and (globals.spectra_param_collection.loc[m, 'Detector'] == detector) and (globals.spectra_param_collection.loc[m, 'Position'] == 'p' + str(position)):
                        if 'Non-point' in globals.spectra_param_collection.keys():
                            #pos_ent.insert(0, globals.spectra_param_collection['Non-point'][m])
                            pos_ent.insert(0, globals.spectra_param_collection.loc[m, 'Non-point'])
                        else:
                            pos_ent.insert(0, 1)
                        break

                sam.append(sample)
                det.append(detector)
                pos.append(position)
                non_point.append(pos_ent)
        # grid
        lab.grid(sticky='W', column=0, row=i + 1)

    button_save = tk.ttk.Button(new_win, text='Save', command=lambda: [
        get_nonpoint(sam, det, pos, non_point),
        treeview_fill(treeview),
        calc_corr.write_all_info(),
        utilities.close_win(root, new_win)
    ])

    button_save.grid(sticky='swne', column=0, row=1, padx=5, pady=5)

    button_save = tk.ttk.Button(new_win, text='Quit', command=lambda: utilities.close_win(root, new_win))
    button_save.grid(sticky='swne', column=1, row=1, padx=5, pady=5)