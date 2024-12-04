# Generates a window for calculation of coincidence correction using TrueCoinc program
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
# TODO: search for files also in home folder

# Libraries
import re
import pandas as pd
import tkinter as tk
import pathlib
import time
import numpy as np
import datetime
# Local modules
from modules import globals
from modules import utilities
from modules import settings

# FIX for pywinauto
import ctypes
ctypes.windll.shcore.SetProcessDpiAwareness(0)

# window for calculation of coincidence correction
def coinc_win(root):
    
    # variables ---------------------------------------------------------------------------------------------------
    det_list = tk.StringVar()

    # widgets -----------------------------------------------------------------------------------------------------
    new_win = tk.Toplevel(root)
    new_win.title('Coincidence correction calculation')
    # new_win.minsize(50, 50)

    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))
    new_win.resizable(False, False)
    new_win.grid_propagate(True)
    
    # Define a new style for the green button
    style = tk.ttk.Style()
    style.configure('green.TButton', background='green')

    # uploads of path to TrueCoinc and shows the adress -----------------------------------------------------------
    path_frame = tk.ttk.LabelFrame(new_win, text='Choose path')
    path_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    path_frame.grid_columnconfigure(0, weight=1)
    path_frame.grid_rowconfigure(0, weight=1)
    path_frame.grid_rowconfigure(1, weight=1)

    tc_file_button = tk.ttk.Button(path_frame, text='TrueCoinc app (.exe file)', width=50, command=lambda: [
        get_tc_exe(),
        tc_file_label_i.config(text='Path: {}'.format(globals.coinc_dict['tc path']))
        ])
    tc_file_button.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    tc_file_label_i = tk.ttk.Label(path_frame, text='Path:')
    tc_file_label_i.grid(sticky='wens', column=0, row=1, padx=5, pady=5)

    # uploads of path to TrueCoinc libraries and shows the adress
    tc_inp_button = tk.ttk.Button(path_frame, text='TC folder (libs., inputs, outputs)', width=50, command=lambda: [
        get_tc_dir(detector_combobox),
        tc_inp_label_i.config(text='Path: {}'.format(globals.coinc_dict['tc input'])),
        ])
    tc_inp_button.grid(sticky='wens', column=0, row=2, padx=5, pady=5)

    tc_inp_label_i = tk.ttk.Label(path_frame, text='Path:')
    tc_inp_label_i.grid(sticky='wens', column=0, row=3, padx=5, pady=5)
    
    # truecoinc settings ------------------------------------------------------------------------------------------------

    settings_frame = tk.ttk.LabelFrame(new_win, text='TrueCoinc automation settings')
    settings_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)

    isotopes_label = tk.ttk.Label(settings_frame, text='Fill down the isotopes for calculation (including isomers if needed).')
    isotopes_label.grid(sticky='wens', column=0, columnspan=2, row=0, padx=5, pady=5)
    
    isotopes_frame = tk.ttk.Frame(settings_frame)
    isotopes_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)
    
    isotopes_text = tk.Text(isotopes_frame, width=30, height=10, bg='white')
    isotopes_text.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    yscrollbar = tk.ttk.Scrollbar(isotopes_frame, orient='vertical', command=isotopes_text.yview)
    yscrollbar.grid(sticky='wens', column=1, row=0)
    isotopes_text.configure(yscrollcommand=yscrollbar.set)

    # buttons for isotopes selection buttons
    selection_frame = tk.ttk.Frame(settings_frame)
    selection_frame.grid(sticky='wens', column=1, row=1, padx=5, pady=5)
        
    gamma_button = tk.ttk.Button(selection_frame, text='Get all nucl.\nfrom gamma lib.', width=20, command=lambda: [
        isotopes_text.delete('1.0', tk.END),
        isotopes_text.insert(tk.END, '\n'.join(select_from_lib()))
        ])
    gamma_button.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    clear_button = tk.ttk.Button(selection_frame, text='Delete all nucl.', width=20, command=lambda: isotopes_text.delete('1.0', tk.END))
    clear_button.grid(sticky='wens', column=0, row=2, padx=5, pady=5)

    det_label = tk.ttk.Label(selection_frame, text='Choose detector for calculation:')
    det_label.grid(sticky='wens', column=0, row=3, padx=5, pady=5)

    # combo box for selection of detector
    detector_combobox = tk.ttk.Combobox(selection_frame, values=det_list, state='disabled')
    detector_combobox.grid(sticky='wens', column=0, row=4, padx=5, pady=5)

    # buttons
    bottom_frame = tk.ttk.Frame(new_win)
    bottom_frame.grid(sticky='wens', column=0, row=2, padx=5, pady=5)
    bottom_frame.grid_columnconfigure(0, weight=1)
    bottom_frame.grid_columnconfigure(1, weight=1)
    bottom_frame.grid_rowconfigure(0, weight=1)
    
    calc_button = tk.ttk.Button(bottom_frame, text='Run Truecoinc', width=24, command=lambda: [
        run_truecoinc(selected_item(isotopes_text.get("1.0", tk.END)), detector_combobox.get())
    ])
    calc_button.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    quit_button = tk.ttk.Button(bottom_frame, text='Quit', width=24, command=lambda: utilities.close_win(root, new_win))
    quit_button.grid(sticky='wens', column=1, row=0, padx=5, pady=5)


# make a list of all isotopes in gamma library
def select_from_lib():
    # read gamma library and return list of isotopes
    utilities.load_gamma_en('Gamma lib. (all)')

    # check if gama library is not empty
    if len(globals.gamma_lib) == 0:
        tk.messagebox.showerror('Library error', 'Gamma library is empty. Please check library file')
    
    return globals.gamma_lib.keys()


# asks user for path to TrueCoinc.exe
def get_tc_exe():  # TODO store path into text file and upload it on menu start
    
    tc_path = tk.filedialog.askopenfilename(title='Select path to TrueCoinc', filetypes=(("Applications", "*.exe"),), initialdir=pathlib.Path.home())
    tc_path = pathlib.Path(tc_path)

    # check if the file is existing
    if not tc_path.exists():
        tk.messagebox.showerror('File not found', 'Please upload correct path to TrueCoinc.exe')
    else:
        globals.coinc_dict['tc path'] = tc_path


# asks user for path to directory with TC folders including .tcd files and ENDSF libraries
def get_tc_dir(combobox):
    tc_inp = tk.filedialog.askdirectory(title='Select path to TrueCoinc input files', initialdir=pathlib.Path.home())
    tc_inp = pathlib.Path(tc_inp)
    
    if tc_inp.exists():
        globals.coinc_dict['tc input'] = tc_inp
    else:
        tk.messagebox.showerror('Folder not found', 'Please upload correct path to TrueCoinc input files')
        return None
    
    # check if truecoinc_output exists
    if not (tc_inp / pathlib.Path('truecoinc_output')).exists() or not (tc_inp / pathlib.Path('detector_efficiency')).is_dir() or not (tc_inp / pathlib.Path('endsf_data')).is_dir():
        tk.messagebox.showerror('Folder not found', 'Please check the folder has correct structure')
        return None

    det_folder = tc_inp / pathlib.Path('detector_efficiency')
    
    # check if subfolders with .tcd files exist
    folder_list = list(det_folder.iterdir())
    # remove all files from the list, keep only folders
    temp_list = folder_list.copy()
    for item in temp_list:
        if item.is_file():
            folder_list.remove(item)
   
    # find .tcd files in subfolders
    if folder_list != []:
        pass
    else:
        tk.messagebox.showerror('Folder not found', 'Please upload correct path to TrueCoinc input files')
        return None
    
    tcd_files = []
    for folder in folder_list:
        tcd_files.extend(list(folder.iterdir()))

    temp_list = tcd_files.copy() 
    for item in temp_list:
        # remove all non .tcd files from the list
        if item.is_dir():
            tcd_files.remove(item)
        elif not item.suffix == '.tcd':
            tcd_files.remove(item)

    globals.tcd_files_dict = {}
    for file in tcd_files:
        globals.tcd_files_dict[file.stem] = file

    combobox['values'] = list(globals.tcd_files_dict.keys())
    combobox.config(state='readonly')
    combobox.current(0)


# get selected isotopes from text widget
def selected_item(text_array):
    try:
        nuclides = text_array.strip().split('\n')
    except Exception as e:
        tk.messagebox.showerror('Error', 'Please check the filled isotopes. Error: {}'.format(e))
        return None
    
    return nuclides


# run TrueCoinc program
def run_truecoinc(nuclides, detpos_choice):
    from pywinauto.application import Application

    # check if all necessary paths are uploaded
    if globals.coinc_dict['tc path'] == '' or globals.coinc_dict['tc input'] == '':
        tk.messagebox.showerror('Paths not found', 'Please upload correct paths to TrueCoinc.exe and TC files.')
        return None
    
    det = detpos_choice.split('_')[0]
    pos = detpos_choice.split('_')[1]

    # check if output sub folder exists if not create it
    output_path = globals.coinc_dict['tc input'] / pathlib.Path('truecoinc_output')
    output_det_path = output_path / pathlib.Path(det)
    output_pos_path = output_det_path / pathlib.Path(pos)

    if not output_det_path.exists():
        output_det_path.mkdir()
        output_pos_path.mkdir()
    else:
        if not output_pos_path.exists():
            output_pos_path.mkdir()
        else:
            # check if the folder is empty
            if len(list(output_pos_path.iterdir())) > 0:
                # ask user if he wants to delete the folder
                if tk.messagebox.askyesno('Folder not empty', 'The folder is not empty. Do you want to delete all files in the folder?'):
                    for file in output_pos_path.iterdir():
                        file.unlink()

    # anounce the program will take control over the computer: start the program, load libraries, calculate and save results
    tk.messagebox.showinfo('TrueCoinc', 'The program will take control over the computer. Please do not use the computer during the calculation. Don\'t move with the mouse or press any key and wait until the program finishes.')

    # wait 1 second before the program starts
    time.sleep(1)

    try:
        app = Application(backend="win32").start(str(globals.coinc_dict['tc path']))
    except FileNotFoundError:
        tk.messagebox.showerror('Apllication not found', 'Please upload correct path to TrueCoinc.exe')

    # start the program
    main_win = app[u'TrueCoinc - [TrueCoinc1]']
    main_win.wait('ready')
    #main_win.print_control_identifiers()

    # upload tcd file
    main_win.menu_select(u"File->Open...")
    window = app.Dialog     # open file dialog
    combobox_name = window.ComboBoxEx   # combobox for file name
    combobox_name.Edit.set_edit_text(str(globals.tcd_files_dict[detpos_choice]))  # insert file name
    button = window.Button
    button.click()

    main_win = app[u'TrueCoinc - ['+ str(detpos_choice) + ']']
    main_win.wait('ready')
       
    # go through all isotopes, upload endsf files, calculate and save results
    for nuclide in nuclides:
        # split the isotope name to isomeric state and isotope name
        split_isotope = re.findall(r'[A-Z][a-z]?|[0-9]+|m[0-9]*', nuclide)
        if len(split_isotope) == 2:
            nucl_numb = int(split_isotope[0])
            isotop_name = split_isotope[1]
        elif len(split_isotope) == 3:
            nucl_numb = int(split_isotope[0])
            isomeric_state = split_isotope[1]
            isotop_name = split_isotope[2]
        
        # upload endsf files
        if (nucl_numb > 0) and (nucl_numb < 10):
            file_name = 'p00' + str(nucl_numb) + isotop_name + '.enx'
        elif (nucl_numb > 9) and (nucl_numb < 100):
            file_name = 'p0' + str(nucl_numb) + isotop_name + '.enx'
        else:
            file_name = 'p' + str(nucl_numb) + isotop_name + '.enx'
        
        file_name = pathlib.Path(file_name)

        path_lib = globals.coinc_dict['tc input'] / pathlib.Path('endsf_data/pensdf2')
        ranges = [(1, 20), (21, 40), (41, 60), (61, 80), (81, 100), (101, 120), (121, 140), (141, 160), (161, 180), (181, 200), (201, 220), (221, 240), (241, 280)]
        
        for start, end in ranges:
            if start <= nucl_numb <= end:
                path = path_lib / pathlib.Path(f'p{start:03d}_{end:03d}') / file_name
                break
        else:
            print(f'Isotope {nuclide} is out of range for current library. Program will skip this isotope.')
            continue

        # open endsf file in TrueCoinc
        main_win.menu_select('Database->ENDSF/2 File')       
        win = app.Dialog
        win.Edit.set_edit_text(str(path))
        win.Open.click()

        # calculate
        main_win.wait('ready')
        main_win.menu_select("Calculation->Calculate")
        
        # save results
        main_win.wait('ready')
        main_win.menu_select("File->Write Result")
        
        win_write = app.Dialog       
        path = output_pos_path / pathlib.Path(f'{nuclide}.txt')
        win_write.Edit.set_edit_text(str(path))
        win_write.button.click()
        
    # close the program
    try:
        main_win.wait('ready')
        main_win.menu_select("File->Exit")
    except Exception as e:
        main_win.wait('ready')
        app.kill()

    # print message to user
    tk.messagebox.showinfo('TrueCoinc', 'The calculation is finished. Results are saved in the folder: ' + str(output_pos_path))

    # create library from TrueCoinc output files
    create_libraries(det, pos, output_pos_path)


# read TrueCoinc output file to pandas dataframe
def read_tc_file(file, tc_output):
    with open(file, 'r', encoding='utf-8') as temp_file:
        content = temp_file.readlines()
    
    # if the file is empty return None
    if len(content) == 0:
        return None

    for i in range(4, len(content)):
        line = content[i]
        
        if len(line) == 0 or len(line) < 38:
            continue

        # split row by column borders
        en = float(line[0:9].strip())
        intens = float(line[9:18].strip()) / 100
        coinc = float(line[27:37].strip())

        tc_output.loc[len(tc_output)] = [file.stem, en, intens, coinc]
    
    return tc_output


def create_libraries(det, pos, output_pos_path):
    # read gamma library
    utilities.load_gamma_en('Gamma lib. (all)')
    
    tc_output = pd.DataFrame(columns=['isotope', 'tc en', 'tc intens', 'tc coinc'])
    final_results = pd.DataFrame(columns=['isotope', 'en', 'coinc'])

    # read all TC files in the output folder into pandas dataframe
    for file in output_pos_path.iterdir():
        if file.is_file():
            if file.suffix == '.txt':
                tc_output = read_tc_file(file, tc_output)

    # go through all isotopes from TC output and find them in gamma library
    # iterate through rows
    for index, row in tc_output.iterrows():
        nucl_tc = row['isotope']
        en_tc = row['tc en']
        int_tc = row['tc intens']
        coinc_tc = row['tc coinc']

        # go through gamma library
        if nucl_tc in globals.gamma_lib.keys():
            gammalb = globals.gamma_lib[nucl_tc]
            en_lib = gammalb[0]
            en_unc = gammalb[1]
            int_lib = gammalb[2]
            int_unc_lib = gammalb[3]

            for i in range(0, len(en_lib)):
                # if match in energy and intensity is found, save it to final results
                if abs(en_lib[i] - en_tc) < (settings.N_CALC * np.sqrt(settings.A_CALC**2 + (en_unc[i] / en_lib[i])**2 + (int_unc_lib[i] / int_lib[i])**2) + 1):
                    if abs(int_lib[i] - int_tc) < (settings.N_CALC * np.sqrt(settings.A_CALC**2 + (en_unc[i] / en_lib[i])**2 + (int_unc_lib[i] / int_lib[i])**2) + 2):
                        final_results.loc[len(final_results)] = [nucl_tc, en_lib[i], coinc_tc]
    
    # save final results to file - print to .clib (text file)
    # go one folder back in the path
    output_det_path = output_pos_path.parent
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    final_results.to_csv(output_det_path / pathlib.Path(f'{det}_{pos}_{timestamp}.clib'), sep='\t', index=False, header=False)
    
    # add two lines to the beginning of the file
    with open(output_det_path / pathlib.Path(f'{det}_{pos}_{timestamp}.clib'), 'r') as original: data = original.read()
    with open(output_det_path / pathlib.Path(f'{det}_{pos}_{timestamp}.clib'), 'w') as modified:
        modified.write(f'# Gamma coincidence library for {det} - detector from TrueCoinc program\n# Isotope    Energy\t\t{pos}\n' + data)
