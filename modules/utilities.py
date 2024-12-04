# module for utility functions shared among other modules
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
import tkinter as tk
import re
import pathlib
from modules import globals

# shared quit function for all windows
def close_win(root, current_window):
        # close the editor window and grab the plot window
        root.grab_set()
        current_window.destroy()


# choose the right parameters in dependency of peak energy and dividing energy
def a_switch(en, en_div, a1, a2):
    if en_div != 0:
        if en <= en_div:
            return a1
        else:
            return a2
    else:
        return a1


# convert different time units to seconds
def time_convert(i, time, unit):
    if unit == 'y':
        return time * 365 * 24 * 3600
    elif unit == 'd':
        return time * 24 * 3600
    elif unit == 'h':
        return time * 3600
    elif unit == 'm':
        return time * 60
    elif unit == 's':
        return time
    else:
        return tk.messagebox.showerror('Library error', 'In the gamma library is bad unit for half life. Please check library file line ' + str(i + 1))


# read a library with gamma energy, intensity, half-life of different radionuclides
def load_gamma_en(library):
    # clear the dictionary
    globals.gamma_lib.clear()

    with open(pathlib.Path(globals.lib_selection[library]), 'r', encoding='utf-8') as temporary_file:
        content = temporary_file.readlines()
    
    # get lines with data and create list of elements
    elements = []
    for line in content:
        line = line.strip().split()
        if len(line) > 1 and not line[0].startswith('#') and not line[0].startswith('-'):
            elements.append(line[4])

    # filter unique elements
    elements = tuple(set(elements))
    
    # create dictionary with gamma energies, intensities and half-lives
    for key in elements:
        en = []
        en_unc = []
        inten = []
        int_unc = []
        half_life = []
        half_life_unc = []

        for i, line in enumerate(content):
            line = line.strip().split()
            if len(line) < 1 or line[0].startswith('#') or line[0].startswith('-'):
                continue
            elif line[4] == key and float(line[0]) not in en:   # check if energy is not already in the list
                en.append(float(line[0]))
                en_unc.append(float(line[1]))

                # check if intensity is in percentage and convert it to fraction
                if float(line[2]) > 1:
                    inten.append(float(line[2]) / 100)
                    int_unc.append(float(line[3]) / 100)
                else:
                    inten.append(float(line[2]))
                    int_unc.append(float(line[3]))
                
                half_life.append(time_convert(i, float(line[7]), line[9].lower()))
                half_life_unc.append(time_convert(i, float(line[8]), line[9].lower()))
        
        # Sort all lists based on intensity
        sorted_lists = sorted(zip(inten, en, en_unc, int_unc, half_life, half_life_unc), key=lambda x: x[0], reverse=True)

        # Unpack the sorted tuples back into separate lists
        inten, en, en_unc, int_unc, half_life, half_life_unc = map(list, zip(*sorted_lists))

        # extract atomic number from the key
        atomic_number = int(re.match(r'\d+', key).group())

        if library == 'Gamma lib. (eff)':    # if library is efficiency library
            globals.gamma_lib_eff[key] = [en, en_unc, inten, int_unc, half_life, half_life_unc, atomic_number]
        else:
            globals.gamma_lib[key] = [en, en_unc, inten, int_unc, half_life, half_life_unc, atomic_number]

    # sort dictionary by nuclide number
    if library == 'Gamma lib. (eff)':
        globals.gamma_lib_eff = dict(sorted(globals.gamma_lib_eff.items(), key=lambda item: item[1][6]))
    else:
        globals.gamma_lib = dict(sorted(globals.gamma_lib.items(), key=lambda item: item[1][6]))

    # print read done if dictionary is not empty
    if globals.gamma_lib or globals.gamma_lib_eff:
        print('\nGamma library ' + library + ' read done.')


# sort data in columns in treeview
def treeview_sort_column(tree, col, reverse):
    data = [(tree.set(child, col), child) for child in tree.get_children('')]
    data.sort(reverse=reverse)

    # rearrange items in sorted positions
    for index, (val, child) in enumerate(data):
        tree.move(child, '', index)

    # rearrange items in sorted positions
    tree.heading(col, command=lambda: treeview_sort_column(tree, col, not reverse))


# check if output directory exists if not than create them
def check_folders():
    if not globals.common_params['Output path'].exists():
        globals.common_params['Output path'].mkdir()
        print('Output directory created inside the Work Directory.')
    
    if not globals.common_params['Spectra path'].exists():
        globals.common_params['Spectra path'].mkdir()
        print('Spectra directory created inside the Work Directory.')

    if not globals.common_params['Results path'].exists():
        globals.common_params['Results path'].mkdir()
        print('Results directory created inside the Work Directory.')


# loads results to treeview
def treeview_update(treeview, df):
    # clear treeview
    treeview.delete(*treeview.get_children())
    
    # check if df is not empty
    if df.empty:
        #tk.messagebox.showerror('Error', 'No data to display.')
        return
    
    # fill treeview with new values    
    for i in range(0, len(df)):  # fill treeview with new values
        treeview.insert('', index='end', values=df.iloc[i, :].tolist(), tags=['checked'])


# TODO add column name as parameter for ouput
# TODO move to utilities
# TODO delete other versions

# return selected tallies from treeview
def selected_tally(treeview_files):

    if len(treeview_files.get_checked()) > 0:
        # send selected tallies to plot_mod function
        selection = []
        for row in treeview_files.get_checked():
            selection.append(pathlib.Path(treeview_files.item(row)['values'][0]))
    else:
        return None
    
    return selection


# selected tally specific columns
def selected_columns(treeview_files, columns):
    if len(treeview_files.get_checked()) > 0 and len(columns) > 0:
        
        all_sel = []
        for column in columns:
            selection = []    
            for row in treeview_files.get_checked():
                selection.append(treeview_files.item(row)['values'][column])
            all_sel.append(selection)
    else:
        #tk.messagebox.showerror('Input error', 'No results were selected.')
        return None
    
    for i in range(0, len(all_sel)):
        if len(all_sel[i]) == 0:
            return None
    
    return all_sel

