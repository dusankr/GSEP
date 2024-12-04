# Main Calculation module
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
import re
import scipy
import math
import tkinter as tk
import ttkwidgets
import pandas as pd
import numpy as np
import matplotlib
import datetime
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Local libraries
from modules import globals
from modules import utilities
from modules import settings
from modules import corr_coinc
from modules import corr_beam
from modules import inputs


# Create a new window with values of selected spectrum - same window as for PR or Q calculation
def result_win(root):
    
    # warning if different values are missing
    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')
    elif 'Delayed time' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'The end of irradiation is missing.')
    
    # check if efficiency correction is done in all spectra
    for spectrum in globals.spectra_collection.keys():
        if 'eff fep' not in globals.spectra_collection[spectrum]:
            return tk.messagebox.showerror('User input error', 'Efficiency correction is necessary and Non-linearity correction is highly recommended for RR caculation.')

    # VARIABLES -------------------------------------------------------------------------------------------------------

    # tkinter variables
    opmenu_sampl_var = tk.StringVar()
    opmenu_nucl_var = tk.StringVar()
    opmenu_en_var = tk.StringVar()
    opmenu_iso_var = tk.StringVar()
    
    coinc_corr_var = tk.BooleanVar(value=True)
    beamtime_corr_var = tk.BooleanVar(value=False)
    number_particles_var = tk.BooleanVar(value=False)
    number_atoms_var = tk.BooleanVar(value=False)
    ratio_var = tk.StringVar(value='1.0')
    irrtime_var = tk.BooleanVar(value=False)
    
    axis_var = tk.StringVar(value='log')

    # other variables -------------------------------------------------------------------------------------------------

    # read gamma energy, intensity and half life from library
    utilities.load_gamma_en('Gamma lib. (all)')
    
    if globals.gamma_lib == {}:
        return tk.messagebox.showerror('Library error', 'Gamma library is empty. Please check the library file.')

    # option menu options - sample list
    sample_list = globals.spectra_param_collection['Sample name'].values.tolist()
    sample_list = sorted(list(set(sample_list)))

    # option menu options - nuclide list
    nuclide_list = list(globals.gamma_lib.keys())  

    # New window ------------------------------------------------------------------------------------------------------
    new_win = tk.Toplevel(root)
    new_win.title('Calculation - final results')

    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: [clear_data(), utilities.close_win(root, new_win)])

    new_win.resizable(False, False)
    new_win.grid_propagate(True)
    new_win.columnconfigure(0, weight=1)
    new_win.columnconfigure(1, weight=1)
    new_win.rowconfigure(0, weight=1)
    new_win.rowconfigure(1, weight=1)

    # main frames
    tree_frame = tk.ttk.Frame(new_win, width=500)
    tree_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    tree_frame.columnconfigure(0, weight=1)
    tree_frame.rowconfigure(0, weight=1)

    plot_frame = tk.ttk.Frame(new_win)
    plot_frame.grid(sticky='wnes', column=1, row=0, padx=5, pady=5)
    plot_frame.columnconfigure(0, weight=1)
    plot_frame.rowconfigure(0, weight=1)
    plot_frame.rowconfigure(1, weight=0)
    
    settings_frame = tk.LabelFrame(new_win)
    settings_frame.grid(sticky='wnes', column=0, row=1, columnspan=2, padx=5, pady=5)
    #settings_frame.columnconfigure(0, weight=0)
    #settings_frame.rowconfigure(0, weight=0)

    #plot_settings_frame = tk.LabelFrame(new_win)
    #plot_settings_frame.grid(sticky='wnes', column=1, row=1, padx=5, pady=5)
    #plot_settings_frame.columnconfigure(0, weight=0)
    #plot_settings_frame.rowconfigure(0, weight=0)

    # tree frame ------------------------------------------------------------------------------------------------------

    treeview_results = ttkwidgets.CheckboxTreeview(tree_frame)
    treeview_results.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    treeview_results['columns'] = list(globals.rr_results.columns)
    treeview_results['show'] = 'headings', 'tree'
    
    # set columns filters
    for col in treeview_results['columns']:
        treeview_results.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview_results, _col, False))
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview_results['columns']:
            treeview_results.column(key, width=value, stretch=True, anchor='center')
    
    # set width of zero column
    treeview_results.column('#0', width=45, stretch=False, anchor='w')

    # Treeview X-scrollbar
    tree_x_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal', command=treeview_results.xview)
    tree_x_scr.grid(sticky='we', column=0, row=1, )
    
    # Treeview Y-scrollbar
    tree_y_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=treeview_results.yview)
    tree_y_scr.grid(sticky='ns', column=1, row=0)

    treeview_results.configure(yscrollcommand=tree_y_scr.set)
    treeview_results.configure(xscrollcommand=tree_x_scr.set)

    # open spectrum on double click
    treeview_results.bind("<Double-1>", lambda event: inputs.open_spectrum(treeview_results.item(treeview_results.focus())['values'][0].split('.')[0], root))

    # isomer treeview -------------------------------------------------------------------------------------------------

    treeview_iso = ttkwidgets.CheckboxTreeview(tree_frame)
    treeview_iso.grid(sticky='wens', column=0, row=2, padx=5, pady=5)
    
    treeview_iso['columns'] = list(globals.rr_iso_results.columns)
    treeview_iso['show'] = 'headings', 'tree'
    
    # set columns filters
    for col in treeview_iso['columns']:
        treeview_iso.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview_iso, _col, False))
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview_iso['columns']:
            treeview_iso.column(key, width=value, stretch=True, anchor='center')
    
    # set width of zero column
    treeview_iso.column('#0', width=45, stretch=False, anchor='w')

    # Treeview X-scrollbar
    tree_x_iso_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal', command=treeview_iso.xview)
    tree_x_iso_scr.grid(sticky='we', column=0, row=3, )
    
    # Treeview Y-scrollbar
    tree_y_iso_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=treeview_iso.yview)
    tree_y_iso_scr.grid(sticky='ns', column=1, row=2)

    treeview_iso.configure(yscrollcommand=tree_y_iso_scr.set)
    treeview_iso.configure(xscrollcommand=tree_x_iso_scr.set)

    # open spectrum on double click
    treeview_iso.bind("<Double-1>", lambda event: inputs.open_spectrum(treeview_iso.item(treeview_iso.focus())['values'][0].split('.')[0], root))


    # plot frame ------------------------------------------------------------------------------------------------------

    # empty figure
    globals.fig_id = matplotlib.figure.Figure(tight_layout=True, figsize=settings.FIG_SIZE)
    globals.ax = globals.fig_id.add_subplot()
    # plot settings
    matplotlib.rcParams['axes.facecolor'] = 'white'
    #globals.ax.set(xlabel='$t_{decay}$ (h)', ylabel='$A$ (s$^{-1}$)', title='Decay curve')
    globals.ax.set(xlabel='$t_{decay}$ (h)', ylabel='$N$ (-)', title='Decay curve')
    matplotlib.rcParams['axes.labelsize'] = 12

    # set major grid
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')
    # set minor grid
    globals.ax.minorticks_on()
    globals.ax.grid(which='minor', linestyle=':', linewidth='0.5')
    
    # --------------------------------------------------------------------------------------------------------------

    globals.canvas_id = FigureCanvasTkAgg(globals.fig_id, plot_frame)
    globals.canvas_id.get_tk_widget().grid(sticky='nswe', column=0, row=0, padx=5, pady=5)

    toolbar_frame = tk.ttk.Frame(plot_frame)
    toolbar_frame.grid(sticky='we', column=0, row=1, padx=5, pady=5)
    toolbar = NavigationToolbar2Tk(globals.canvas_id, toolbar_frame)
    toolbar.update()
    
    # settings frame --------------------------------------------------------------------------------------------------

    # first row frame
    first_row_frame = tk.Frame(settings_frame)
    first_row_frame.grid(sticky='wens', column=0, row=0, padx=0, pady=0)

    optionmenu_frame = tk.LabelFrame(first_row_frame, text='Calculation settings')
    optionmenu_frame.grid(sticky='wens', column=0, row=0, padx=2, pady=2)

    label_spl = tk.ttk.Label(optionmenu_frame, text='Sample:')
    label_spl.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    opm_sample = tk.ttk.Combobox(optionmenu_frame, textvariable=opmenu_sampl_var, width=10, state="readonly")
    opm_sample.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    opm_sample['values'] = sample_list
    opm_sample.current(0)

    label_nucl = tk.ttk.Label(optionmenu_frame, text='Nuclide:')    
    label_nucl.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    opm_nuclide = tk.ttk.Combobox(optionmenu_frame, textvariable=opmenu_nucl_var, width=8, state="readonly")
    opm_nuclide.grid(sticky='w', column=3, row=0, padx=5, pady=5)
    opm_nuclide['values'] = nuclide_list
    opm_nuclide.current(0)
    opm_nuclide.bind('<<ComboboxSelected>>', lambda event: update_opmenu_en(opmenu_nucl_var.get(), opm_energy))

    label_en = tk.ttk.Label(optionmenu_frame, text='Energy - Intensity:')
    label_en.grid(sticky='w', column=4, row=0, padx=5, pady=5)

    opm_energy = tk.ttk.Combobox(optionmenu_frame, textvariable = opmenu_en_var, width=20, state="readonly")
    opm_energy.grid(sticky='w', column=5, row=0, padx=5, pady=5)
    opm_energy['values'] = []
    # run for the first time
    update_opmenu_en(opmenu_nucl_var.get(), opm_energy)

    # isomer corr -------------------------------------------------------------------------------------------------

    # isomer label frame
    iso_frame = tk.LabelFrame(first_row_frame, text='Isomer correction')
    iso_frame.grid(sticky='wens', column=1, row=0, padx=2, pady=2)
    
    label_iso = tk.ttk.Label(iso_frame, text='Isomer:')
    label_iso.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    # create isomer list from gamma library, filter only isomers

    # extract keys with isomer state from gamma library and pass them to the list for option menu
    isomer_list = [key for key in globals.gamma_lib.keys() if re.search(r'm', key)]
    isomer_list = ['None'] + isomer_list
    
    opm_isomer = tk.ttk.Combobox(iso_frame, textvariable=opmenu_iso_var, width=8, state="readonly")
    opm_isomer.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    opm_isomer['values'] = isomer_list
    opm_isomer.current(0)
    
    iso_button = tk.ttk.Button(iso_frame, text='Recalc. Avg. & Replot', width=20, state='disabled', command=lambda: [
        iso_avg_result(utilities.selected_columns(treeview_iso, [0, 3])),
        plot_decay(utilities.selected_columns(treeview_iso, [0, 3])),
        ])
    iso_button.grid(sticky='w', column=2, row=0, padx=5, pady=5)
    
    opm_isomer.bind('<<ComboboxSelected>>', lambda event: [get_variables() ,update_iso_button(iso_button, treeview_iso)])

    # additional corrections ------------------------------------------------------------------------------------------
    corr_frame = tk.LabelFrame(first_row_frame, text='Additional corrections')
    corr_frame.grid(sticky='wens', column=2, row=0, padx=2, pady=2)

    check_coinc = tk.ttk.Checkbutton(corr_frame, text='Coincidence', var=coinc_corr_var)
    check_coinc.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    check_beamtime = tk.ttk.Checkbutton(corr_frame, text='Beam fluctuation', var=beamtime_corr_var, state='disabled')
    check_beamtime.grid(sticky='w', column=1, row=0, padx=5, pady=5)

    # plot settings frame ----------------------------------------------------------------------------------------------
    plot_settings_frame = tk.LabelFrame(first_row_frame, text='Plot settings')
    plot_settings_frame.grid(sticky='wns', column=3, row=0, padx=2, pady=2)

    button_plot = tk.ttk.Button(plot_settings_frame, text='Replot', width=15, command=lambda: [
        get_variables(),
        plot_decay(utilities.selected_tally(treeview_results))
    ])
    button_plot.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    axis_label = tk.ttk.Label(plot_settings_frame, text='Y axis scale:')
    axis_label.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    
    linear_button = tk.ttk.Radiobutton(plot_settings_frame, text='Lin', value='lin', variable=axis_var)
    linear_button.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    log_button = tk.ttk.Radiobutton(plot_settings_frame, text='Log', value='log', variable=axis_var)
    log_button.grid(sticky='w', column=3, row=0, padx=5, pady=5)

    # buttons frame ---------------------------------------------------------------------------------------------------
    second_row_frame = tk.Frame(settings_frame)
    second_row_frame.grid(sticky='wens', column=0, row=1, padx=0, pady=0)

    # calculation
    calculation_frame = tk.LabelFrame(second_row_frame, text='Calculation')
    calculation_frame.grid(sticky='wens', column=0, row=0, padx=2, pady=2)
    
    # sample, nuclide, energy, coinc_corr, beam_corr, isotop_ratio, iso_nucl, iso_en, particle_norm, atom_norm
    button_calc = tk.ttk.Button(calculation_frame, text='Calculate', width=20, command=lambda: [
        get_variables(),
        calculation(),
        utilities.treeview_update(treeview_results, globals.rr_results),
        isomer_calculation(treeview_iso),
        plot_decay(utilities.selected_tally(treeview_results)),
        update_iso_button(iso_button, treeview_iso),
        ])
    button_calc.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    check_irr = tk.ttk.Checkbutton(calculation_frame, text='× t_irr', var=irrtime_var)
    check_irr.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    
    check_particle = tk.ttk.Checkbutton(calculation_frame, text='×1/particle', var=number_particles_var)
    check_particle.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    check_atoms = tk.ttk.Checkbutton(calculation_frame, text='×1/(atom·ratio)', var=number_atoms_var)
    check_atoms.grid(sticky='w', column=3, row=0, padx=5, pady=5)

    ratio_entry = tk.ttk.Entry(calculation_frame, width=10, textvariable=ratio_var)
    ratio_entry.grid(sticky='w', column=4, row=0, padx=5, pady=5)   
    
    # disable options based on available data
    if globals.common_params['Beam int'] is not None or globals.common_params['Beam int unc'] is None or globals.common_params['Irradiation time'] is None:
        check_beamtime.configure(state='disabled')

    if globals.common_params['Beam int'] is None or globals.common_params['Beam int unc'] is None or globals.common_params['Irradiation time'] is None:
        check_particle.configure(state='disabled')
    
    if globals.spectra_param_collection['Weight'] is None or globals.spectra_param_collection['Molar mass'] is None:
        check_atoms.configure(state='disabled')
   
    # export frame ----------------------------------------------------------------------------------------------------
    export_frame = tk.LabelFrame(second_row_frame, text='Export (save) results')
    export_frame.grid(sticky='wens', column=1, row=0, padx=2, pady=2)

    button_add = tk.ttk.Button(export_frame, text='Add to export', width=15, command=lambda: add_to_export(utilities.selected_tally(treeview_results)))
    button_add.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    button_clear = tk.ttk.Button(export_frame, text='Clear export', width=15, command=lambda: [
        globals.rr_all_results.drop(globals.rr_all_results.index, inplace=True),
        tk.messagebox.showinfo('Export', 'Export table was cleared.')
        ])
    button_clear.grid(sticky='w', column=1, row=0, padx=5, pady=5)

    button_export = tk.ttk.Button(export_frame, text='Export window', width=15, command=lambda: export_window(root))
    button_export.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    line2 = tk.ttk.Separator(export_frame, orient='vertical')
    line2.grid(sticky='ns', column=3, row=0, padx=5, pady=5)

    button_save = tk.ttk.Button(export_frame, text='Save2XLSX (curr. nucl.)', width=20, command=lambda: save_results(utilities.selected_tally(treeview_results)))
    button_save.grid(sticky='w', column=4, row=0, padx=5, pady=5)

    # --------------------------------------------------------------------------------------------------------------

    button_close = tk.ttk.Button(second_row_frame, text='Close window', width=15, command=lambda: [clear_data(), utilities.close_win(root, new_win)])
    button_close.grid(sticky='wesn', column=2, row=0, padx=5, pady=5)

    # function for getting variables from widgets and save them to global dictionary
    def get_variables():
        globals.widgets_dict['sample'] = opmenu_sampl_var.get()
        globals.widgets_dict['nuclide'] = opmenu_nucl_var.get()
        globals.widgets_dict['iso_nuclide'] = opmenu_iso_var.get()
        globals.widgets_dict['energy_lib'] = float(opmenu_en_var.get().strip().split(':')[0])
        globals.widgets_dict['ratio'] = float(ratio_var.get())
        
        # True/False
        globals.widgets_dict['switch_coinc'] = coinc_corr_var.get()
        globals.widgets_dict['switch_beam'] = beamtime_corr_var.get()
        globals.widgets_dict['particle_norm'] = number_particles_var.get()
        globals.widgets_dict['atom_norm'] = number_atoms_var.get()
        globals.widgets_dict['irrtime_norm'] = irrtime_var.get()
        globals.widgets_dict['y_scale'] = axis_var.get()
        

# fill the energy option menu based on selected nuclide
def update_opmenu_en(nuclide, optionmenu):
    #opm_energy['values'] = []
    
    if nuclide == 'None':
        energies = ["0 : 0"]
    else:
        en = globals.gamma_lib[nuclide][0]
        intens = globals.gamma_lib[nuclide][2]
        energies = [f"{en[i]} : {intens[i]:.4f}" for i in range(0, len(en))]
    
    optionmenu['values'] = energies
    optionmenu.current(0)


# calculate RR for sample and nuclide
def calculation():

    # clear dataframes before calculation
    globals.rr_results = globals.rr_results.iloc[0:0]
    globals.rr_iso_results = globals.rr_iso_results.iloc[0:0]

    # loop through all names and if sample name for name is the same as selected sample than calculate RR    
    for spectrum, sample_name in zip(globals.spectra_param_collection['Name'].to_list(), globals.spectra_param_collection['Sample name'].to_list()):
        if sample_name == globals.widgets_dict['sample']:       
            spectrum_calc(spectrum)


# one spectrum calculation
def spectrum_calc(spectrum):
    spectrum_stem = spectrum.stem

    particle_norm = globals.widgets_dict['particle_norm']
    atom_norm = globals.widgets_dict['atom_norm']
    irrtime_norm = globals.widgets_dict['irrtime_norm']
    beam_corr = globals.widgets_dict['switch_beam']
    ratio = globals.widgets_dict['ratio']

    # select nuclide and get list of energies
    if globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict == {}:   # calculate isomer values for correction
        nuclide = globals.widgets_dict['iso_nuclide']
        en_lib_list = globals.gamma_lib[nuclide][0]
    # isomer correction from previous calculation
    elif globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict['nuclide'] == globals.widgets_dict['iso_nuclide']:
        nuclide = globals.widgets_dict['nuclide']
        en_lib_list = [globals.widgets_dict['energy_lib']]
    else:
        nuclide = globals.widgets_dict['nuclide']   # calculate standard nuclide values
        en_lib_list = [globals.widgets_dict['energy_lib']]

    # get parameters from spectrum_param_collection for selected spectrum
    t_live = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Live time'].values[0]
    t_real = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Real time'].values[0]
    t_delay = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Delayed time'].values[0]
    mass = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Weight'].values[0]
    mass_unc = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Weight error'].values[0]
    mass_molar = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Molar mass'].values[0]

    # irradiation time in seconds
    t_irr = globals.common_params['Irradiation time'] * 60

    # correction for compounds or mixtures or natural isotopic composition
    if ratio == 0 or ratio is None:
        ratio = 1
    else:
        ratio = float(ratio)
    
    # calculate normalization per one sample atom
    n_atoms = 1
    if atom_norm is True:
        n_atoms = (scipy.constants.N_A * mass * ratio) / mass_molar
        n_atoms = 1 / n_atoms

    # calculate normalization per one projectile particle
    n_particles = 1
    if particle_norm is True:
        beam_int = globals.common_params['Beam int']
        beam_int_unc = globals.common_params['Beam int unc']
        n_particles = 1 / beam_int

    # calculate normalization
    t_irr_corr = 1
    if irrtime_norm is True:
        t_irr_corr = t_irr

    # go through all energies in the list (standard nuclide = 1 en., isomer = 1 or more en.)
    for en_lib in en_lib_list:
    
        # get values from gamma library for selected nuclide and energy
        for i in range(0, len(globals.gamma_lib[nuclide][0])):
            if en_lib == globals.gamma_lib[nuclide][0][i]:
                en_unc_lib = globals.gamma_lib[nuclide][1][i]
                int_lib = globals.gamma_lib[nuclide][2][i]
                int_unc_lib = globals.gamma_lib[nuclide][3][i]
                globals.widgets_dict['halflife_lib'] = halflife_lib = globals.gamma_lib[nuclide][4][i]
                halflife_unc_lib = globals.gamma_lib[nuclide][5][i]
                lamb = np.log(2) / halflife_lib  # decay constant
        
        # get isomer half life from gamma library if isomer is selected and dictionary is not empty
        if globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict != {}:
            halflife_iso_lib = globals.gamma_lib[globals.widgets_dict['iso_nuclide']][4][0]
            lamb_iso = np.log(2) / halflife_iso_lib

        # if coincidence correction is checked than try to find value in the library
        coi_corr = 1.0
        if globals.widgets_dict['switch_coinc'] is True:
            position = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Position'].values[0]
            detector = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spectrum, 'Detector'].values[0]

            if detector in globals.coinc_lib_selection:
                coi_corr = corr_coinc.coinc_correction(nuclide, en_lib, corr_coinc.coinc_lib(detector), position)
        
        # loop through all rows (peaks) in one spectrum and get values
        for index, row in globals.spectra_collection[spectrum_stem].iterrows():

            # chose energy from spectrum
            if 'cor_en' in globals.spectra_collection[spectrum_stem]:
                en = row['cor_en']
                en_unc = row['cor_en_err']
            else:
                en = row['energy']
                # calculate energy uncertainty from channel error, energy and channel
                en_unc = row['cherr'] * row['energy'] / row['channel']

            # calculate RR if peak energy is in the gamma library
            #if abs(en - en_lib) < (settings.N_CALC * math.sqrt(settings.A_CALC**2 + (en_unc / en)**2 + (en_unc_lib / en_lib)**2) + 1):
            if abs(en - en_lib) < (settings.N_CALC * math.sqrt(settings.A_CALC**2 + (en_unc / en)**2 + (en_unc_lib / en_lib)**2)):
                pass
            else:
                continue

            # chose area from spectrum library
            if ('SEP-DEP cor area' in globals.spectra_collection[spectrum_stem]) and ('bcg cor area' in globals.spectra_collection[spectrum_stem]):
                area = row['SEP-DEP cor area']
                area_unc = row['bcg cor area unc']
            elif ('SEP-DEP cor area' in globals.spectra_collection[spectrum_stem]) and ('bcg cor area' not in globals.spectra_collection[spectrum_stem]):
                area = row['SEP-DEP cor area']
                area_unc = row['area'] * row['aerr%'] / 100
            elif ('bcg cor area' in globals.spectra_collection[spectrum_stem]) and ('SEP-DEP cor area' not in globals.spectra_collection[spectrum_stem]):
                area = row['bcg cor area']
                area_unc = row['bcg cor area unc']
            else:
                area = row['area']
                area_unc = row['area'] * row['aerr%'] / 100

            if 'eff fep' in globals.spectra_collection[spectrum_stem]:
                eff = row['eff fep']
            else:
                eff = 1.0
                tk.messagebox.showinfo('Calculation warning: ', 'Efficiency correction is missing, results are not comparable.')

            if 'att cor' in globals.spectra_collection[spectrum_stem]:
                att = row['att cor']
            else:
                att = 1.0

            if 'non-point corr' in globals.spectra_collection[spectrum_stem]:
                nonpoint = row['non-point corr']
            else:
                nonpoint = 1.0

            if beam_corr is True:
                beam_corr = corr_beam.beam_correction(halflife_lib)
            else:
                beam_corr = 1.0

            # calculate RR for one gamma energy and save results ------------------------------------------------------------------------------------------------
            # TODO implement systematic error calculation? (beam, mcnp, neutron bcg. correction, efficiency, attenuation, non-point, coincidence)

            # standard nuclide, no isomer correction
            if globals.widgets_dict['iso_nuclide'] == 'None':
                #print('Standard nuclide, no isomer correction')
                try:
                    # calculate corrections
                    corrections = beam_corr * (1 / eff) * (1 / att) * (1 / nonpoint) * (1 / coi_corr)

                    '''
                    # decay calculation
                    A_decay = corrections * (area / (int_lib * t_live)) * (t_real * lamb / (1 - np.exp(-lamb * t_real)))    # normalized to beginning of measurement
                    dA_decay = A_decay * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)
                    A_0 = A_decay * (1 / np.exp(-lamb * t_delay))
                    '''
                    # decay calculation
                    N_decay = corrections * (area / (int_lib * t_live)) * (t_real / (1 - np.exp(-lamb * t_real)))
                    dN_decay = N_decay * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)
                    N_0 = N_decay * (1 / np.exp(-lamb * t_delay))

                    # final calculation without normalization
                    Q = corrections * (area / (int_lib * t_live)) * (1 / np.exp(-lamb * t_delay)) * (t_real * lamb / (1 - np.exp(-lamb * t_real))) * (1 / (1 - np.exp(-lamb * t_irr)))
                    # final normalization
                    Q = Q * n_particles * n_atoms * t_irr_corr

                    # calculate uncertainty what should be the uncertainty?
                    if atom_norm is True:
                        dQ = Q * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2 + (mass_unc/mass)**2)
                    else:
                        dQ = Q * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)

                except Exception as e:
                    return tk.messagebox.showerror('Calculation error - no isomer option:', 'Calculation error: ' + str(e))
            
            # calculate isomer values needed for correction, correction is not calculated in this step
            elif globals.widgets_dict['iso_nuclide'] != 'None' and not globals.iso_res_dict:
                #print('Isomer nuclide - 1st step')
                try:
                    # calculate corrections
                    corrections = beam_corr * (1 / eff) * (1 / att) * (1 / nonpoint) * (1 / coi_corr)
                    
                    # decay calculation
                    N_decay = corrections * (area / (int_lib * t_live)) * (t_real / (1 - np.exp(-lamb * t_real)))
                    dN_decay = N_decay * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)
                    # N_0
                    N_m0 = N_decay * (1 / np.exp(-lamb * t_delay))
                    dN_m0 = N_m0 * dN_decay / N_decay

                    # N_yield,m
                    N_yield_m = corrections * (area / (int_lib * t_live)) * (1 / np.exp(-lamb * t_delay)) * (t_real * lamb / (1 - np.exp(-lamb * t_real))) * (1 / (1 - np.exp(-lamb * t_irr))) * t_irr

                    dN_yield_m = N_yield_m * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)

                except Exception as e:
                    return tk.messagebox.showerror('Calculation error - isomer option - 1st step:', 'Calculation error: ' + str(e))
            
            # calculate isomer correction to get ground state values
            elif globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict['nuclide'] == globals.widgets_dict['iso_nuclide']:
                #print('Isomer nuclide - 2nd step')
                try:
                    # calculate corrections
                    corrections = beam_corr * (1 / eff) * (1 / att) * (1 / nonpoint) * (1 / coi_corr)

                    # A_decay - probably incorrect, only for testing
                    N_decay = corrections * (area / (int_lib * t_live)) * (t_real / (1 - np.exp(-lamb * t_real)))    # normalized to beginning of measurement
                    dN_decay = N_decay * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)

                    N_g0 = corrections * ((area * t_real) / (int_lib * t_live)) * (1 / np.exp(-lamb * t_delay)) * (1 / (1 - np.exp(-lamb * t_real))) + ((lamb_iso * globals.iso_res_dict['N_m0 (-)']) / (lamb - lamb_iso)) * (1 - np.exp((lamb - lamb_iso) * t_delay))

                    C_irr_g = (lamb * t_irr) / (1 - np.exp(-lamb * t_irr)) - globals.iso_res_dict['N_yield,m (-)'] / (N_g0 * (1 - np.exp(-lamb * t_irr))) * (1 - (lamb_iso * np.exp(-lamb * t_irr) - lamb * np.exp(-lamb_iso * t_irr)) / (lamb_iso - lamb))
                    
                    # N_yield,g recalculated to production rate (Bq)
                    Q = N_yield_g = N_g0 * C_irr_g / t_irr

                    # normalization based on user selection
                    Q = N_yield_g * n_particles * n_atoms * t_irr_corr

                    dQ = Q * math.sqrt(0.03 * (t_real / t_live - 1)**2 + t_delay**2 * lamb**2 * (halflife_unc_lib/halflife_lib)**2 + (area_unc/area)**2 + (int_unc_lib/int_lib)**2)
                    
                except Exception as e:
                    return tk.messagebox.showerror('Calculation error', 'Calculation error - isomer option - 2nd step: ' + str(e))

            # Adding a new row to globals.rr_results or globals.rr_iso_results
            if globals.widgets_dict['iso_nuclide'] == 'None':
                globals.rr_results.loc[len(globals.rr_results)] = [
                    spectrum,
                    nuclide,
                    en_lib,
                    round(Q, 4),
                    round(dQ, 4),
                    0.0,    # C_irr_g
                    round(float(coi_corr), 4),
                    round(float(beam_corr), 4),
                    round(N_decay, 4),
                    round(dN_decay, 4),
                    round(t_delay / 3600, 4),
                    round(N_0, 4)
                    ]
            elif globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict == {}:
                globals.rr_iso_results.loc[len(globals.rr_iso_results)] = [
                    spectrum,
                    nuclide,
                    en_lib,
                    round(N_yield_m, 4),
                    round(dN_yield_m, 4),
                    round(float(coi_corr), 4),
                    round(float(beam_corr), 4),
                    round(N_decay, 4),
                    round(dN_decay, 4),
                    round(t_delay / 3600, 4),
                    round(N_m0, 4),
                    round(dN_m0, 4)
                    ]
            elif globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict['nuclide'] == globals.widgets_dict['iso_nuclide']:
                globals.rr_results.loc[len(globals.rr_results)] = [
                    spectrum,
                    nuclide,
                    en_lib,
                    round(Q, 4),
                    round(dQ, 4),
                    round(float(C_irr_g), 4),
                    round(float(coi_corr), 4),
                    round(float(beam_corr), 4),
                    round(N_decay, 4),
                    round(dN_decay, 4),
                    round(t_delay / 3600, 4),
                    round(N_g0, 4)
                    ]


def plot_iso_decay(selection):
    
    # fit function for isomer decay - lambda_g and lambda_m have to be in hours    
    def decay_fit_isomer(t, N_g0, N_m0, lamb_g, lamb_m):
        N_g = lamb_m * N_m0 / (lamb_g - lamb_m) * (np.exp(-lamb_m * t) - np.exp(-lamb_g * t)) + N_g0 * np.exp(-lamb_g * t) 
        return N_g
    
    sample = globals.widgets_dict['sample']
    nuclide = globals.widgets_dict['nuclide']
    energy = globals.widgets_dict['energy_lib']
    lib_halflife = globals.gamma_lib[nuclide][4][0]
    lib_iso_halflife = globals.gamma_lib[globals.widgets_dict['iso_nuclide']][4][0]
    lamb_g = np.log(2) / (lib_halflife / 3600)
    lamb_m = np.log(2) / (lib_iso_halflife / 3600)

    # filter dataframe by selection
    filtered_df = globals.rr_results[globals.rr_results['Spectrum'].isin(selection)].copy()    

    # clear the plot
    globals.ax.clear()

    if filtered_df.empty or len(filtered_df) < 2:
        globals.ax.text(0.2,0.5,'Number of points is less than 2.', fontsize=12)
        return globals.canvas_id.draw()
    
    x_value = filtered_df['t_delay (hr)'].to_list()
    y_value = filtered_df['N (-)'].to_list()
    y_err = filtered_df['dN (-)'].to_list()

    # plot point with error bars
    globals.ax.errorbar(x_value, y_value, yerr=y_err, fmt='o', capsize=5, ecolor='#44b7f9', color='#44b7f9', markerfacecolor='#44b7f9', markersize=5)

    # x values range for fit
    x_fit = np.linspace(min(x_value) * 0.7, max(x_value) * 1.3, 1000)
    
    p0 = [max(y_value), max(y_value), lamb_g, lamb_m]   # initial values for fit

    # decay_fit_isomer(t, N_g0, N_m0, lamb_g, lamb_m)
    try:
        coeff, pcov = scipy.optimize.curve_fit(decay_fit_isomer, x_value, y_value, p0=p0, sigma=y_err, absolute_sigma=True)
        globals.ax.plot(x_fit, decay_fit_isomer(x_fit, *coeff), color='red', linewidth=1, label='decay fit')
    except Exception as e:
        tk.messagebox.showinfo(title='Notification', message='Fit cannot be calculated. ' + str(e))
    
    # chi2 calculation
    fitted_y = [decay_fit_isomer(t, *coeff) for t in x_value]
    chi, red_chi = fit_chi2(y_value, fitted_y, y_err)

    # text for plot
    text_halftime_g = "{:.2f}".format(lib_halflife / 3600)
    text_halftime_m = "{:.2f}".format(lib_iso_halflife / 3600)
    text_chi = "{:.2f}".format(chi)
    text_chi_red = "{:.2f}".format(red_chi)
    text_N_g0 = "{:.2e}".format(coeff[0])
    text_N_m0 = "{:.2e}".format(coeff[1])
    
    text_fig = '$T_{1/2,g}$ = ' + text_halftime_g + ' h\n$T_{1/2,m}$ = ' + text_halftime_m + ' h\n$\\chi^2 = $' + text_chi + '\n$\\chi^2_{red} = $' + text_chi_red + '\n$N_{g0} = $' + text_N_g0 + '\n$N_{m0} = $' + text_N_m0

    # add text to plot        
    globals.ax.text(0.995, 0.99, text_fig,
                 transform=globals.ax.transAxes, bbox={'facecolor': 'white', 'pad': 2, 'alpha': 0.75, 'boxstyle': 'round, pad=0.3, rounding_size=0.2'},
                 horizontalalignment='right', verticalalignment='top')
    
    # plot settings ---------------------------------------------------------------------------------------------------
    globals.ax.set(xlabel='$t_{decay}$ (h)', ylabel='$N$ (-)')
    matplotlib.rcParams['axes.labelsize'] = 12
    # set major grid
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')
    
    # set axis scale
    if globals.widgets_dict['y_scale'] == 'log':
        globals.ax.set_yscale('log')
        globals.ax.grid(which='minor', linestyle=':', linewidth='0.5')
    else:
        globals.ax.set_yscale('linear')
    
    # set title, use nuclide and energy
    globals.ax.set_title('Decay curve for ' + nuclide + ' ' + str(energy) + ' keV')
    # set legend
    legend = globals.ax.legend(loc='lower left')
    legend.get_frame().set_edgecolor('black')

    # save plot
    utilities.check_folders()
    plot_path = globals.common_params['Results path'] / pathlib.Path(f"{nuclide}-iso-corr_{str(energy)}keV_{sample}.png")
    try:
        globals.fig_id.savefig(plot_path, format='png', dpi=300)
    except Exception as e:
        return tk.messagebox.showerror('Save error', 'Error while saving plot: ' + str(e))

    return globals.canvas_id.draw()
            

# creates the plot
def plot_decay(selection):
    if selection is None:
        return
    
    if globals.widgets_dict['iso_nuclide'] != 'None' and globals.iso_res_dict != {} and globals.rr_iso_results.empty:
        return plot_iso_decay(selection)

    # fit function - standard decay
    def decay_fit(t, N_0, halflife):
        func = N_0 * np.exp(-(np.log(2) / halflife) * np.array(t))
        return func

    sample = globals.widgets_dict['sample']
    
    if globals.widgets_dict['iso_nuclide'] != 'None' and len(selection) == 2:
        nuclide = globals.widgets_dict['iso_nuclide']
        energy = 'all'
        # filter dataframe by selection        
        spectrum = [pathlib.Path(x) for x in selection[0]]
        N_yield = [float(x) for x in selection[1]]

        # filter dataframe by selection
        filtered_df = globals.rr_iso_results[globals.rr_iso_results['Spectrum'].isin(spectrum) & globals.rr_iso_results['N_yield,m (-)'].isin(N_yield)].copy()

    else:
        nuclide = globals.widgets_dict['nuclide']
        energy = globals.widgets_dict['energy_lib']
        # filter dataframe by selection
        filtered_df = globals.rr_results[globals.rr_results['Spectrum'].isin(selection)].copy()    

    # clear the plot
    globals.ax.clear()

    if filtered_df.empty or len(filtered_df) < 2:
        globals.ax.text(0.2,0.5,'Number of points is less than 2.', fontsize=12)
        return globals.canvas_id.draw()
    
    x_value = filtered_df['t_delay (hr)'].to_list()
    #y_value = filtered_df['A (s^-1)'].to_list()
    #y_err = filtered_df['dA (s^-1)'].to_list()
    y_value = filtered_df['N (-)'].to_list()
    y_err = filtered_df['dN (-)'].to_list()

    # plot point with error bars
    globals.ax.errorbar(x_value, y_value, yerr=y_err, fmt='o', capsize=5, ecolor='#44b7f9', color='#44b7f9', markerfacecolor='#44b7f9', markersize=5)

    # x values range for fit
    x_fit = np.linspace(min(x_value) * 0.8, max(x_value) * 1.2, 1000)
    # T1/2 from library
    lib_halflife = globals.gamma_lib[nuclide][4][0]
    
    # initial values for fit
    p0 = [max(y_value), lib_halflife / 3600]

    # fit
    try:
        coeff, pcov = scipy.optimize.curve_fit(decay_fit, x_value, y_value, p0=p0, sigma=y_err, absolute_sigma=True)
        globals.ax.plot(x_fit, decay_fit(x_fit, *coeff), color='red', linewidth=1, label='decay fit')

        # theoretical decay (known halflife) to show correct slope of the decay
        teor_coeff = coeff[0], lib_halflife / 3600
        globals.ax.plot(x_fit, decay_fit(x_fit, *teor_coeff), color='red', linewidth=1, label='calc. decay fit', linestyle='--')
    except Exception as e:
        tk.messagebox.showinfo(title='Notification', message='Fit cannot be calculated. ' + str(e))
        text_halftime_fit = "-"
        text_chi = "-"
        text_halftime_lib = str("{:.2f}".format(lib_halflife / 3600))

    # chi2 calculation - theoretical decay
    N_decay_teor = decay_fit(x_value, *teor_coeff)    
    chi_teor, red_chi_teor = fit_chi2(y_value, N_decay_teor, y_err)

    # chi2 calc for the 1st fit
    fitted_y = [decay_fit(t, *coeff) for t in x_value]
    chi, red_chi = fit_chi2(y_value, fitted_y, y_err)

    # halflife calc
    fitted_halflife = coeff[1]

    # text for plot
    text_halftime_fit = "{:.2f}".format(fitted_halflife)
    text_halftime_lib = "{:.2f}".format(lib_halflife / 3600)    
    text_chi = "{:.2f}".format(chi)
    text_chi_red = "{:.2f}".format(red_chi)
    text_chi_teor = "{:.2f}".format(chi_teor)
    text_chi_red_teor = "{:.2f}".format(red_chi_teor)
    text_fig = '$T_{1/2, fit}$ = ' + text_halftime_fit + ' h\n$T_{1/2, tab}$ = ' + text_halftime_lib + ' h\n$\\chi^2 = $' + text_chi + '\n$\\chi^2_{red} = $' + text_chi_red + '\n$\\chi^2_{teor} = $' + text_chi_teor + '\n$\\chi^2_{red, teor} = $' + text_chi_red_teor

    if globals.widgets_dict['iso_nuclide'] != 'None':
        text_N_m0 = "{:.2e}".format(globals.iso_res_dict['N_m0 (-)'])
        text_N_yield_m = "{:.2e}".format(globals.iso_res_dict['N_yield,m (-)'])
        text_fig = text_fig + '\n$N_{0,m}$ = ' + text_N_m0 + '\n$N_{yield,m}$ = ' + text_N_yield_m

    # add text to plot        
    globals.ax.text(0.995, 0.99, text_fig,
                 transform=globals.ax.transAxes, bbox={'facecolor': 'white', 'pad': 2, 'alpha': 0.75, 'boxstyle': 'round, pad=0.3, rounding_size=0.2'},
                 horizontalalignment='right', verticalalignment='top')
    
    # plot settings ---------------------------------------------------------------------------------------------------
    globals.ax.set(xlabel='$t_{decay}$ (h)', ylabel='$N$ (-)')
    #globals.ax.set(xlabel='$t_{decay}$ (h)', ylabel='$A$ (s$^{-1}$)')
    matplotlib.rcParams['axes.labelsize'] = 12
    # set major grid
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')
    
    # set axis scale
    if globals.widgets_dict['y_scale'] == 'log':
        globals.ax.set_yscale('log')
        globals.ax.grid(which='minor', linestyle=':', linewidth='0.5')
    else:
        globals.ax.set_yscale('linear')
    
    # set title, use nuclide and energy
    globals.ax.set_title('Decay curve for ' + nuclide + ' ' + str(energy) + ' keV')
    # set legend
    legend = globals.ax.legend(loc='lower left')
    legend.get_frame().set_edgecolor('black')

    # save plot
    utilities.check_folders()
    plot_path = globals.common_params['Results path'] / pathlib.Path(f"{nuclide}_{str(energy)}keV_{sample}.png")
    try:
        globals.fig_id.savefig(plot_path, format='png', dpi=300)
    except Exception as e:
        return tk.messagebox.showerror('Save error', 'Error while saving plot: ' + str(e))

    return globals.canvas_id.draw()
    

# save selection from treeview to XLSX and calculate weighted average, internal and external error
def save_results(selection):

    if globals.rr_results.empty or selection is None:
        return tk.messagebox.showinfo(title='Notification', message='No results to save.')
       
    # save only selected results from globals.rr_results - filter dataframe by selection, then save to xlsx
    filtered_df = globals.rr_results[globals.rr_results['Spectrum'].isin(selection)].copy()
    
    # calculate result weigths from error values for each row and add it to the dataframe (needed for calculating weighted average)
    filtered_df['w_i'] = [1 / (x ** 2) for x in filtered_df['dRR (arb.)']]
    w_sum = filtered_df['w_i'].sum()

    filtered_df['omega (-)'] = [x / w_sum for x in filtered_df['w_i']]

    RR_w_avg = (filtered_df['RR (arb.)'] * filtered_df['w_i']).sum() / w_sum
    dRR_int = 1 / math.sqrt(w_sum)
    
    # calculate external error
    numerator = (filtered_df['w_i'] * [(x - RR_w_avg)**2 for x in filtered_df['RR (arb.)'].to_list()]).to_list()
    dRR_w_ext = np.sqrt(sum(numerator) / (w_sum * (len(filtered_df) - 1)))

    # change w_i column format to scientific notation
    filtered_df['w_i'] = filtered_df['w_i'].apply(lambda x: "{:.2e}".format(x))

    # add new row with '----' in each column
    filtered_df = pd.concat([filtered_df, pd.DataFrame([['-------------'] * len(filtered_df.columns)], columns=filtered_df.columns)], ignore_index=True)

    # add three new rows for results, description will be in the first column and values in the second
    filtered_df.loc[len(filtered_df), 'Spectrum'] = 'Results:'
    filtered_df.loc[len(filtered_df), 'Spectrum'] = 'RR_avg (arb.)'
    filtered_df.loc[len(filtered_df) - 1, 'Nuclide'] = RR_w_avg
    filtered_df.loc[len(filtered_df), 'Spectrum'] = 'dRR_int (arb.)'
    filtered_df.loc[len(filtered_df) - 1, 'Nuclide'] = dRR_int
    filtered_df.loc[len(filtered_df), 'Spectrum'] = 'dRR_ext (arb.)'
    filtered_df.loc[len(filtered_df) - 1, 'Nuclide'] = dRR_w_ext

    # save to xlsx
    utilities.check_folders() # check if output folder exists
    file_path = globals.common_params['Results path'] / pathlib.Path(f"{globals.widgets_dict['sample']}_{globals.widgets_dict['nuclide']}_{int(globals.widgets_dict['energy_lib'])}keV.xlsx")
    
    try:
        filtered_df.to_excel(file_path, index=False)
    except Exception as e:
        return tk.messagebox.showerror('Save error', 'Error while saving results: ' + str(e))

    tk.messagebox.showinfo(title='Notification', message=f'Results were saved to {str(file_path)}\n')


# my chisquare (not same as the pearson's from scipy.stats.chisquare() !!!!!!
def fit_chi2(measured, fitted, error):
    
    chi2 = 0
    
    # from RubyDecan
    for y_observed, y_expected, err in zip(measured, fitted, error):
        chi2 += (y_observed - y_expected)**2 / err**2
    
    # https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    # for y_observed, y_expected, err in zip(measured, fitted, error):
    #     chi2 += (y_observed - y_expected)**2
    # chi2 = chi2 / np.var(measured)

    # weird calculation of reduced chi2, Radek uses only number of points without subtracting number of parameters
    # red_chi2 = chi2 / (len(measured) - n) where n = 2 number of parameters in the model (N0/A0, lambda)?
    red_chi2 = chi2 / len(measured)

    return chi2, red_chi2


# add selected results to export
def add_to_export(selection):
    if globals.rr_results.empty or selection is None:
        return tk.messagebox.showinfo(title='Notification', message='No results to add to export.')
    
    # filter dataframe by selection
    filtered_df = globals.rr_results[globals.rr_results['Spectrum'].isin(selection)].copy()

    # TODO during adding to export check if the same spectrum is already in the export table, if yes, then skip it    
        
    # add filtered_df to rr_all_results
    globals.rr_all_results = pd.concat([globals.rr_all_results, filtered_df], ignore_index=True)
    
    tk.messagebox.showinfo(title='Notification', message=f'Selected results were added to export. {len(globals.rr_all_results)} results in total in export table.')


# export all results - open new window with checkbox treeview and load all results from export dataframe, two buttons - save and close
def export_window(root):
    if globals.rr_all_results.empty:
        return tk.messagebox.showinfo(title='Notification', message='No results to export.')
    
    # create new window
    new_win = tk.Toplevel(root)    
    new_win.title('Export results')
    #new_win.geometry('800x600')

    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))
    new_win.columnconfigure(0, weight=1)
    new_win.rowconfigure(0, weight=1)
    new_win.rowconfigure(1, weight=0)
    
    # tree frame ------------------------------------------------------------------------------------------------------
    tree_frame = tk.ttk.Frame(new_win)
    tree_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    tree_frame.columnconfigure(0, weight=1)
    tree_frame.rowconfigure(0, weight=1)

    treeview_results = ttkwidgets.CheckboxTreeview(tree_frame)
    treeview_results.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    # set columns
    treeview_results['columns'] = list(globals.rr_all_results.columns)
    treeview_results['show'] = 'headings', 'tree'

    # fill treeview with data
    utilities.treeview_update(treeview_results, globals.rr_all_results)
    
    # set columns filters
    for col in treeview_results['columns']:
        treeview_results.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview_results, _col, False))
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview_results['columns']:
            treeview_results.column(key, width=value, stretch=True, anchor='center')
    
    # set width of zero column
    treeview_results.column('#0', width=55, stretch=False, anchor='w')

    # Treeview X-scrollbar
    tree_x_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal', command=treeview_results.xview)
    tree_x_scr.grid(sticky='we', column=0, row=4, columnspan=5)
    
    # Treeview Y-scrollbar
    tree_y_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=treeview_results.yview)
    tree_y_scr.grid(sticky='ns', column=1, row=0, rowspan=5)

    treeview_results.configure(yscrollcommand=tree_y_scr.set)
    treeview_results.configure(xscrollcommand=tree_x_scr.set)

    # open spectrum on double click
    treeview_results.bind("<Double-1>", lambda event: inputs.open_spectrum(treeview_results.item(treeview_results.focus())['values'][0].split('.')[0], root))

    # button frame -----------------------------------------------------------------------------------------------------
    button_frame = tk.LabelFrame(new_win)
    button_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)

    # save button
    button_save = tk.ttk.Button(button_frame, text='Save selected to XLSX', width=20, command=lambda: save_results_export(utilities.selected_tally(treeview_results)))
    button_save.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    # drop selected button
    button_drop = tk.ttk.Button(button_frame, text='Drop selected', width=15, command=lambda: [
        drop_selected(utilities.selected_tally(treeview_results)),
        utilities.treeview_update(treeview_results, globals.rr_all_results)
        ])
    button_drop.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    
    # close button
    button_close = tk.ttk.Button(button_frame, text='Close window', width=15, command=lambda: utilities.close_win(root, new_win))
    button_close.grid(sticky='w', column=2, row=0, padx=5, pady=5)


# save selected results from export to XLSX
def save_results_export(selection):
    # check if there are any results to save
    if globals.rr_all_results.empty or selection is None:
        return tk.messagebox.showinfo(title='Notification', message='No results to save.')

    # save only selected results from globals.rr_all_results - filter dataframe by selection, then save to xlsx
    filtered_df = globals.rr_all_results[globals.rr_all_results['Spectrum'].isin(selection)].copy()

    # sort dataframe by nuclide and energy
    filtered_df = filtered_df.sort_values(by=['Nuclide', 'E (keV)'])
        
    # calculate result weigths from error values for each row and add it to the dataframe (needed for calculating weighted average)
    filtered_df['w_i'] = [1 / (x ** 2) for x in filtered_df['dRR (arb.)']]

    # temporary dataframe for energy and nuclide results
    final_results_df = pd.DataFrame(columns=['Nuclide', 'E (keV)', 'mean RR (arb.)', 'dRR_int (arb.)', 'dRR_ext (arb.)'])

    # calculate weighted average, internal and external error for each energy group
    for nuclide in filtered_df['Nuclide'].unique():
        for energy in filtered_df.loc[filtered_df['Nuclide'] == nuclide, 'E (keV)'].unique():
            # get data for one energy group -
            # filter dataframe by nuclide and energy
            filtered_df_nucl_en = filtered_df.loc[(filtered_df['Nuclide'] == nuclide) & (filtered_df['E (keV)'] == energy)].copy()

            # check if there is only one value in the group
            if len(filtered_df_nucl_en) < 2:
                RR_w_avg = filtered_df_nucl_en['RR (arb.)'].values[0]
                dRR_int = dRR_w_ext = filtered_df_nucl_en['dRR (arb.)'].values[0]
            else:
                # calculate weighted average, internal and external error for one energy group
                w_sum = filtered_df_nucl_en['w_i'].sum()

                # weights visualization
                filtered_df_nucl_en['omega_en (-)'] = [x / w_sum for x in filtered_df_nucl_en['w_i']]

                RR_w_avg = (filtered_df_nucl_en['RR (arb.)'] * filtered_df_nucl_en['w_i']).sum() / w_sum
                dRR_int = 1 / math.sqrt(w_sum)
            
                # calculate external error
                numerator = (filtered_df_nucl_en['w_i'] * [(x - RR_w_avg)**2 for x in filtered_df_nucl_en['RR (arb.)'].to_list()]).to_list()
                dRR_w_ext = np.sqrt(sum(numerator) / (w_sum * (len(filtered_df_nucl_en) - 1)))

                # add values to the final dataframe
                final_results_df.loc[len(final_results_df)] = [nuclide, energy, RR_w_avg, dRR_int, dRR_w_ext]
        
        # calculate weighted average, internal and external error for one nuclide (use en. groups internal error as dRR)
        
        # check if more than one energy for current nuclide
        if len(filtered_df.loc[filtered_df['Nuclide'] == nuclide, 'E (keV)'].unique()) < 2:
            RR_w_avg = final_results_df.loc[final_results_df['Nuclide'] == nuclide, 'mean RR (arb.)'].values[0]
            dRR_int = dRR_w_ext = final_results_df.loc[final_results_df['Nuclide'] == nuclide, 'dRR_int (arb.)'].values[0]
            final_results_df.loc[len(final_results_df)] = [nuclide, 'all', RR_w_avg, dRR_int, dRR_w_ext]
        else:
            # filter final results dataframe by nuclide
            final_results_df_nucl = final_results_df.loc[final_results_df['Nuclide'] == nuclide].copy()
            # calculate weights for each energy
            final_results_df_nucl['w_i'] = [1 / (x ** 2) for x in final_results_df_nucl['dRR_int (arb.)']]
            w_sum = final_results_df_nucl['w_i'].sum()
            final_results_df_nucl['omega_nucl (-)'] = [x / w_sum for x in final_results_df_nucl['w_i']]
            # calculate weighted average, internal and external error for one nuclide
            RR_w_avg = (final_results_df_nucl['mean RR (arb.)'] * final_results_df_nucl['w_i']).sum() / w_sum
            dRR_int = 1 / math.sqrt(w_sum)
            # calculate external error
            numerator = (final_results_df_nucl['w_i'] * [(x - RR_w_avg)**2 for x in final_results_df_nucl['mean RR (arb.)'].to_list()]).to_list()
            dRR_w_ext = np.sqrt(sum(numerator) / (w_sum * (len(final_results_df_nucl) - 1)))
            # add values to the final dataframe
            final_results_df.loc[len(final_results_df)] = [nuclide, 'all', RR_w_avg, dRR_int, dRR_w_ext]

    # save part -----------------------------------------------------------------------------------------------------
    # TODO save everything to one sheet of xlsx file

    # go through w_i columns and change format to scientific notation
    filtered_df['w_i'] = filtered_df['w_i'].apply(lambda x: "{:.2e}".format(x))
        
    # save to xlsx
    utilities.check_folders() # check if output folder exists

    # save to xlsx with datetime in the name (no seconds) as a unique identifier with general name without nuclide and energy details
    file_path = globals.common_params['Results path'] / pathlib.Path(f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')}_results_export.xlsx")
    
    # save filtered_df and final_results_df to xlsx, each to separate sheet
    try:
        with pd.ExcelWriter(file_path) as writer:
            filtered_df.to_excel(writer, sheet_name='Results', index=False)
            final_results_df.to_excel(writer, sheet_name='Final results', index=False)
    except Exception as e:
        return tk.messagebox.showerror('Save error', 'Error while saving results: ' + str(e))

    tk.messagebox.showinfo(title='Notification', message=f'Results were saved to {str(file_path)}\n')


# drop selected values from export
def drop_selected(selection):
    if selection is None:
        return tk.messagebox.showinfo(title='Notification', message='No results to drop.')
    
    # drop rows for selected values
    globals.rr_all_results.drop(globals.rr_all_results[globals.rr_all_results['Spectrum'].isin(selection)].index, inplace=True)

    tk.messagebox.showinfo(title='Notification', message='Selected results were dropped from export.')


# isomer main function
def isomer_calculation(treeview_iso):
    if globals.widgets_dict['iso_nuclide'] == 'None':
        return
    
    if globals.rr_iso_results.empty:
        return #tk.messagebox.showinfo(title='Notification', message='No results to calculate avg. N_yield,m and N_m0.')

    # upload isomer results to treeview
    utilities.treeview_update(treeview_iso, globals.rr_iso_results),
    
    # calculate isomer results from selected values
    iso_avg_result(utilities.selected_columns(treeview_iso, [0, 3]))

    # plot isomer decay curve
    plot_decay(utilities.selected_columns(treeview_iso, [0, 3])),


# change isomer button state
def update_iso_button(iso_button, treeview):
    # check if isomer option menu is selected
    if globals.widgets_dict['iso_nuclide'] == 'None' and treeview.get_children(): # and len(globals.rr_iso_results) == 0:
        iso_button.config(state='disabled')
        treeview.delete(*treeview.get_children())   # clear treeview
        # clear isomer results
        globals.iso_res_dict = {}
    elif globals.widgets_dict['iso_nuclide'] != 'None' and treeview.get_children():
        iso_button.config(state='normal')


# calculate single isomer contribution from selected values
def iso_avg_result(sel_list):
    if globals.widgets_dict['iso_nuclide'] == 'None':
        return
    
    # check if there are any results to calculate
    if globals.widgets_dict['iso_nuclide'] != 'None' and globals.rr_iso_results.empty or sel_list is None:
        return tk.messagebox.showinfo(title='Notification', message='No results to calculate weighted average.')
    
    spectrum = [pathlib.Path(x) for x in sel_list[0]]
    N_yield = [float(x) for x in sel_list[1]]

    # filter dataframe by selection
    filtered_df = globals.rr_iso_results[globals.rr_iso_results['Spectrum'].isin(spectrum) & globals.rr_iso_results['N_yield,m (-)'].isin(N_yield)].copy()

    if filtered_df.empty:
        return tk.messagebox.showinfo(title='Notification', message='No results to calculate.')

    for val, err in zip(['N_yield,m (-)', 'N_m0 (-)'], ['dN_yield,m (-)', 'dN_m0 (-)']):
        # calculate result weigths from error values for each row and add it to the dataframe (needed for calculating weighted average)
        filtered_df['w_i'] = [1 / (x ** 2) for x in filtered_df[err]]      # 'N_yield,m (-)', 'dN_yield,m (-)'
        
        w_sum = filtered_df['w_i'].sum()

        filtered_df['omega (-)'] = [x / w_sum for x in filtered_df['w_i']]
        
        # weighted average
        val_w_avg = (filtered_df[val] * filtered_df['w_i']).sum() / w_sum
        err_int = 1 / math.sqrt(w_sum)
    
        # calculate external error
        numerator = (filtered_df['w_i'] * [(x - val_w_avg)**2 for x in filtered_df[val].to_list()]).to_list()
        err_w_ext = np.sqrt(sum(numerator) / (w_sum * (len(filtered_df) - 1)))

        # change w_i column format to scientific notation
        filtered_df['w_i'] = filtered_df['w_i'].apply(lambda x: "{:.2e}".format(x))

        if val == 'N_yield,m (-)':
            # check if the result is already in the dictionary
            globals.iso_res_dict['N_yield,m (-)'] = float(val_w_avg)
            globals.iso_res_dict['dN_yield,m (-)'] = err_int
            #globals.iso_res_dict['dN_yield,m_ext (-)'] = err_w_ext
        elif val == 'N_m0 (-)':
            globals.iso_res_dict['N_m0 (-)'] = float(val_w_avg)
            globals.iso_res_dict['dN_m0 (-)'] = err_int
            #globals.iso_res_dict['dN_m0_ext (-)'] = err_w_ext
    
    globals.iso_res_dict['nuclide'] = globals.widgets_dict['iso_nuclide']
   
    x = globals.iso_res_dict['N_yield,m (-)']
    y = globals.iso_res_dict['N_m0 (-)']
    return tk.messagebox.showinfo(title='Notification', message=f'Weigthed average for N_yield,m {x:.4e} and N_m0 {y:.4e} were calculated.')


# clear global variables for rr calculations
def clear_data():
    globals.iso_res_dict = {}
    globals.rr_results = globals.rr_results.iloc[0:0]
    globals.rr_iso_results = globals.rr_iso_results.iloc[0:0]
    globals.rr_all_results = globals.rr_all_results.iloc[0:0]
    globals.widgets_dict = {}