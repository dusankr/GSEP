# efficiency calculation 
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
# MEmatplotlib.rcHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# TODO add option for overlapping eff values for better connection between curves
# TODO add auto optimization of div energy for two curves fit
# TODO work with scipy.optimize.curve_fit outputs for better control of fit: func should not be overparametrized

# Libraries
import pandas as pd
import numpy as np
import scipy.optimize as opt
from datetime import datetime
import matplotlib
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
import ttkwidgets
import pathlib
# Local modules
from modules import globals
from modules import utilities
from modules import settings
from modules import corr_coinc
from modules import corr_eff

pd.options.mode.chained_assignment = None


# window for calculation of effectivity
def efficiency_win(root):

    # warning if different values are missing
    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')
    
    # empty globals.eff_param_collection
    globals.eff_param_collection = globals.eff_param_collection.iloc[0:0]

    # filter only spectra with measurement number >= 1000
    globals.eff_param_collection = globals.spectra_param_collection[globals.spectra_param_collection['Measurement number'] >= 1000]

    if len(globals.eff_param_collection) == 0:
        return tk.messagebox.showerror('User input error', 'No valid spectra for efficiency calculation.')

    # initialize libraries - read efficiency library, activity library and gamma library

    # read gamma energy, intensity and half life from library
    utilities.load_gamma_en('Gamma lib. (eff)')
    
    if globals.gamma_lib_eff == {}:
        return tk.messagebox.showerror('Library error', 'Gamma library is empty. Please check the library file.')
    
    # read efficiency parameters from library
    # corr_eff.eff_lib(det, pos)    # return: div energy, eff parameters, unit

    # read activity library
    activity_lib()

    # detector list
    det_list = list(set(globals.eff_param_collection['Detector']))
    det_list.sort()

    # add reference date, activity and activity uncertainty to globals.eff_param_collection
    fill_eff_parameters()

    # tkinbter VARIABLES ----------------------------------------------------------------------------------------------
    opmenu_det = tk.StringVar(value=det_list[0])
    opmenu_pos = tk.StringVar()
    opmenu_sel = tk.StringVar()
    opmenu_sel2 = tk.StringVar()
    var_curve = tk.BooleanVar(value=False)
    opmenu_plotscale = tk.StringVar()
    div_en_var = tk.StringVar(value='0.0')
    var_coi_corr = tk.BooleanVar(value=False)
    # -----------------------------------------------------------------------------------------------------------------

    new_win = tk.Toplevel(root)
    new_win.title('Detector window - FEP efficiency calculation')
    
    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))
    new_win.resizable(False, False)
    new_win.grid_propagate(True)

    # main frames
    top_frame = tk.ttk.Frame(new_win)
    top_frame.grid(sticky='wens', column=0, columnspan=2, row=0, padx=5, pady=5)

    tree_frame = tk.ttk.Frame(new_win, width=500)
    tree_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)
    tree_frame.columnconfigure(0, weight=1)
    tree_frame.rowconfigure(0, weight=1)

    plot_frame = tk.ttk.Frame(new_win)
    plot_frame.grid(sticky='wnes', column=1, row=1, padx=5, pady=5)
    plot_frame.columnconfigure(0, weight=1)
    plot_frame.rowconfigure(0, weight=1)
    plot_frame.rowconfigure(1, weight=0)
    
    bottom_frame = tk.LabelFrame(new_win)
    bottom_frame.grid(sticky='wnes', column=0, row=2, columnspan=2, padx=5, pady=5)    
    
    # top frame ------------------------------------------------------------------------------------------------------
    # input settings
    input_frame = tk.LabelFrame(top_frame, text='Input settings')
    input_frame.grid(sticky='wens', column=0, row=0, padx=2, pady=2)

    label_title = tk.ttk.Label(input_frame, text='Calculation for detector')
    label_title.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    opm_det = tk.ttk.Combobox(input_frame, textvariable=opmenu_det, width=5, state="readonly")
    opm_det.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    opm_det['values'] = det_list
    opm_det.current(0)  

    label_pos = tk.ttk.Label(input_frame, text='position')
    label_pos.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    opm_pos = tk.ttk.Combobox(input_frame, textvariable=opmenu_pos, width=5, state="readonly")
    opm_pos.grid(sticky='w', column=3, row=0, padx=5, pady=5)
    opm_pos['values'] = update_opmenu_pos(det_list[0], opm_pos)

    opm_det.bind('<<ComboboxSelected>>', lambda event: update_opmenu_pos(opmenu_det.get(), opm_pos))
    
    # add checkbox for coincidence correction
    coi_check = tk.ttk.Checkbutton(input_frame, text='COI corr.', variable=var_coi_corr)
    coi_check.grid(sticky='w', column=4, row=0, padx=5, pady=5)

    calc_button = tk.ttk.Button(input_frame, text='Calculate', width=12, command=lambda: [
        eff_calc(opmenu_det.get(), opmenu_pos.get(), var_coi_corr.get()),
        utilities.treeview_update(treeview_eff, globals.eff_results_df),
        plot_eff_fig(utilities.selected_tally(treeview_eff), div_en_entry.get(), opmenu_plotscale.get(), var_curve.get(), opmenu_det.get(), opmenu_pos.get()),
        fill_export_opmenu(opm_sel, opm_sel2, var_curve.get())
        ])
    calc_button.grid(sticky='w', column=5, row=0, padx=5, pady=5)

    # bottom frame ---------------------------------------------------------------------------------------------------

    fit_frame = tk.LabelFrame(bottom_frame, text='Fit settings')
    fit_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)

    single_curve_button = tk.ttk.Radiobutton(fit_frame, text="Single curve", variable=var_curve, value=False)
    single_curve_button.grid(sticky='w', column=0, row=0, padx=5)

    two_curves_button = tk.ttk.Radiobutton(fit_frame, text="Two curves", variable=var_curve, value=True)
    two_curves_button.grid(sticky='w', column=1, row=0, padx=5)
    
    def toggle_entry_state():
        if var_curve.get():
            div_en_entry.config(state="normal")
            opm_sel2.config(state="readonly")
        else:
            div_en_entry.config(state="disabled")
            opm_sel2.config(state="disabled")

    var_curve.trace_add('write', lambda *args: toggle_entry_state())
    
    div_en_label = tk.ttk.Label(fit_frame, text='Div. energy')
    div_en_label.grid(sticky='w', column=2, row=0, padx=5)
    
    div_en_entry = tk.ttk.Entry(fit_frame, width=6, state="disabled", textvariable=div_en_var)
    div_en_entry.grid(sticky='w', column=3, row=0, padx=5)

    unit_en_label = tk.ttk.Label(fit_frame, text='keV')
    unit_en_label.grid(sticky='w', column=4, row=0, padx=5)

    scale_label = tk.ttk.Label(fit_frame, text='Ax. scale')
    scale_label.grid(sticky='w', column=5, row=0, padx=5, pady=5)

    # add option menu for plot scale
    opm_plotscale = tk.ttk.Combobox(fit_frame, textvariable=opmenu_plotscale, width=11, state="readonly")
    opm_plotscale.grid(sticky='w', column=6, row=0, padx=5, pady=5)
    opm_plotscale['values'] = ['xy-lin', 'x-log_y-lin', 'xy-log']
    opm_plotscale.current(0)

    plot_button = tk.ttk.Button(fit_frame, text='Plot selection', width=15, command=lambda: [
                                plot_eff_fig(utilities.selected_tally(treeview_eff), div_en_entry.get(), opmenu_plotscale.get(), var_curve.get(), opmenu_det.get(), opmenu_pos.get()),
                                fill_export_opmenu(opm_sel, opm_sel2, var_curve.get())
                                ])
    plot_button.grid(sticky='w', column=7, row=0, padx=5)

    # results frame --------------------------------------------------------------------------------------------------
    
    results_frame = tk.LabelFrame(bottom_frame, text='Results')
    results_frame.grid(sticky='wens', column=1, row=0, padx=5, pady=5)
    
    choose_label = tk.ttk.Label(results_frame, text='Degree')
    choose_label.grid(sticky='w', column=0, row=0, padx=5, pady=5)
    
    opm_sel = tk.ttk.Combobox(results_frame, textvariable=opmenu_sel, width=12, state="readonly")
    opm_sel.grid(sticky='w', column=1, row=0, padx=5, pady=5)

    opm_sel2 = tk.ttk.Combobox(results_frame, textvariable=opmenu_sel2, width=12, state="disabled")
    opm_sel2.grid(sticky='w', column=2, row=0, padx=5, pady=5)
    
    save_button = tk.ttk.Button(results_frame, text='Save coefficients', width=15, command=lambda: save_eff_fit(opmenu_det.get(), opmenu_pos.get(), var_curve.get(), opmenu_sel.get(), opmenu_sel2.get()))
    save_button.grid(sticky='w', column=3, row=0, padx=5, pady=5)

    # -----------------------------------------------------------------------------------------------------------------

    button_close = tk.ttk.Button(bottom_frame, text='Close window', width=15, command=lambda: utilities.close_win(root, new_win))
    button_close.grid(sticky='wsne', column=3, row=0, padx=5, pady=5)

    # treeview frame -------------------------------------------------------------------------------------------------

    treeview_eff = ttkwidgets.CheckboxTreeview(tree_frame)
    treeview_eff.grid(sticky='wens', column=0, row=0, padx=5, pady=5)

    treeview_eff['columns'] = ['E (keV)', 'E_ref (kev)','Nuclide', 'Int (%)', 'A0 (Bq)', 'COI_corr', 'Eff (-)', 'dEff (-)']
    treeview_eff['show'] = 'headings', 'tree'

    # set columns filters
    for col in treeview_eff['columns']:
        treeview_eff.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview_eff, _col, False))
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview_eff['columns']:
            treeview_eff.column(key, width=value, stretch=True, anchor='center')
    
    # set width of zero column
    treeview_eff.column('#0', width=45, stretch=False, anchor='w')
    
    # Treeview X-scrollbar
    tree_x_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal', command=treeview_eff.xview)
    tree_x_scr.grid(sticky='we', column=0, row=1, columnspan=2)
    
    # Treeview Y-scrollbar
    tree_y_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=treeview_eff.yview)
    tree_y_scr.grid(sticky='ns', column=2, row=0, rowspan=2)
    
    treeview_eff.configure(xscrollcommand=tree_x_scr.set)
    treeview_eff.configure(yscrollcommand=tree_y_scr.set)

    # open spectrum on double click
    # treeview_eff.bind("<Double-1>", lambda event: inputs.open_spectrum(treeview_eff.item(treeview_eff.focus())['values'][0].split('.')[0], root))
    # TODO doesn't work because the treeview_eff doens't contain spectra names

    # Figure ----------------------------------------------------------------------------------------------------------

    # empty figure
    globals.fig_id = matplotlib.figure.Figure(tight_layout=True, figsize=settings.FIG_SIZE)
    globals.ax = globals.fig_id.add_subplot()
    # plot settings
    matplotlib.rcParams['axes.facecolor'] = 'white'
    globals.ax.set(xlabel='$E$ (keV)', ylabel='$eff_{FEP}$ (-)', title=f'FEP efficiency calculation - {opmenu_det.get()} - {opmenu_pos.get()}')
    matplotlib.rcParams['axes.labelsize'] = 12

    # set major grid
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')

    # figure widget
    globals.canvas_id = FigureCanvasTkAgg(globals.fig_id, plot_frame)
    globals.canvas_id.get_tk_widget().grid(sticky='nswe', column=0, row=0, padx=5, pady=5)

    toolbar_frame = tk.ttk.Frame(plot_frame)
    toolbar_frame.grid(sticky='we', column=0, row=1, padx=5, pady=5)
    toolbar = NavigationToolbar2Tk(globals.canvas_id, toolbar_frame)
    toolbar.update()


# update positions for selected detector in option menu
def update_opmenu_pos(det, optionmenu):
    # filter positions for selected detector from globals.eff_param_collection into list
    position = globals.eff_param_collection[globals.eff_param_collection['Detector'] == det]['Position'].tolist()
    # remove duplicates and sort
    position = list(set(position))
    position.sort()
    
    optionmenu['values'] = position
    optionmenu.current(0)


# read library with informations about activities of samples
def activity_lib():
    with open(globals.lib_selection['Activity'], 'r', encoding='utf-8') as temporary_file:
        content = temporary_file.readlines()

    nuclide = []
    id_num = []
    date = []
    activity = []
    activity_unc = []

    for line in content:
        line = line.strip()
        if len(line) == 0:
            continue
        elif line[0] == '#' or line[0] == '-':
            continue
        
        line = line.split()
        if len(line) == 5:
            try:
                nuclide.append(line[0])
                id_num.append(int(line[1]))
                date.append(datetime.strptime(line[2] + ' 00:00:00', '%d-%m-%Y %H:%M:%S'))
                activity.append(float(line[3]))
                activity_unc.append(float(line[4]))
            except Exception as e:
                tk.messagebox.showerror('Library error', 'Error in library file. Error: {}'.format(e))
        else:
            tk.messagebox.showerror('Library error', f'Error in nuclide {line[0]} from {line[1]}.')

    globals.activity_lib_df['nucl'] = nuclide
    globals.activity_lib_df['ID'] = id_num
    globals.activity_lib_df['ref_date'] = date
    globals.activity_lib_df['A0'] = activity
    globals.activity_lib_df['dA0_%'] = activity_unc


# extend dataframe with reference date, activity and activity uncertainty
def fill_eff_parameters():
    ref_date = []
    activity = []
    activity_unc = []

    for index, row in globals.eff_param_collection.iterrows():
        nuclide = row['Sample name']
        id = row['Measurement number']

        if nuclide in globals.activity_lib_df['nucl'].values:
            for i in range(0, len(globals.activity_lib_df['nucl'])):
                if nuclide == globals.activity_lib_df['nucl'][i] and id == globals.activity_lib_df['ID'][i]:
                    ref_date.append(globals.activity_lib_df['ref_date'][i])
                    activity.append(globals.activity_lib_df['A0'][i])
                    activity_unc.append(globals.activity_lib_df['dA0_%'][i])
                    break
        else:
            ref_date.append(None)
            activity.append(None)
            activity_unc.append(None)

    globals.eff_param_collection['ref_date'] = ref_date
    globals.eff_param_collection['A0'] = activity
    globals.eff_param_collection['dA0_%'] = activity_unc
    # calculate delayed time as difference between start of measurement and reference date in seconds
    globals.eff_param_collection['delayed_time'] = (globals.eff_param_collection['Start of measurement'] - globals.eff_param_collection['ref_date']).dt.total_seconds()


# calculation FEP efficiency from all measured data - go through all spectra
def eff_calc(det, pos, coi):
    # clear globals.eff_results_df
    globals.eff_results_df = globals.eff_results_df.iloc[0:0]

    # filter only selected detector and position
    temp_eff_df = globals.eff_param_collection[(globals.eff_param_collection['Detector'] == det) & (globals.eff_param_collection['Position'] == pos)]

    # check if there are some data
    if len(temp_eff_df) == 0:
        return tk.messagebox.showerror('User input error', f'No valid spectra for selected detector {det} and position {pos}.')
    
    # go through filtred data, compare with gamma library and calculate efficiency
    for index, row in temp_eff_df.iterrows():
        spectrum = row['Name']
        spectrum_stem = spectrum.stem
        nucl = row['Sample name']
        t_live = row['Live time']
        t_real = row['Real time']
        t_delay = row['delayed_time']
        ref_act = row['A0']
        ref_act_err = row['dA0_%']

        # get values from gamma library for current nuclide
        try:
            dict_array = globals.gamma_lib_eff[nucl]
            en_lib, en_err_lib, int_lib, int_err_lib, hl_lib, hl_err_lib, nucl_num = dict_array[0], dict_array[1], dict_array[2], dict_array[3], dict_array[4], dict_array[5], dict_array[6]
        except:
            tk.messagebox.showerror('Library error', f'No data for nuclide {nucl} in gamma library. Program will continue with next nuclide.')
            continue

        # go through spectrum and locate corresponding gamma energies
        for ind, spec_row in globals.spectra_collection[spectrum_stem].iterrows():

            # chose from spectrum based on applied corrections
            if 'cor_en' in globals.spectra_collection[spectrum_stem]:
                en = spec_row['cor_en']
                en_err = spec_row['cor_en_err']
            else:
                en = spec_row['energy']
                # calculate energy uncertainty from channel error, energy and channel
                en_unc = spec_row['cherr'] * spec_row['energy'] / spec_row['channel']
            
            for i in range(0, len(en_lib)):
                # find corresponding gamma energy in library
                if abs(en - en_lib[i]) < (settings.N_CALC * np.sqrt(settings.A_CALC**2 + (en_unc / en)**2 + (en_err_lib[i] / en_lib[i])**2) + 0.5): # TODO remove +0.5
                    lamb = np.log(2) / hl_lib[i]

                    # chose area from spectrum library
                    if ('SEP-DEP cor area' in globals.spectra_collection[spectrum_stem]) and ('bcg cor area' in globals.spectra_collection[spectrum_stem]):
                        area = spec_row['SEP-DEP cor area']
                        area_unc = spec_row['bcg cor area unc']
                    elif ('SEP-DEP cor area' in globals.spectra_collection[spectrum_stem]) and ('bcg cor area' not in globals.spectra_collection[spectrum_stem]):
                        area = spec_row['SEP-DEP cor area']
                        area_unc = spec_row['area'] * spec_row['aerr%'] / 100
                    elif ('bcg cor area' in globals.spectra_collection[spectrum_stem]) and ('SEP-DEP cor area' not in globals.spectra_collection[spectrum_stem]):
                        area = spec_row['bcg cor area']
                        area_unc = spec_row['bcg cor area unc']
                    else:
                        area = spec_row['area']
                        area_unc = spec_row['area'] * spec_row['aerr%'] / 100
                    
                    # conincidence correction
                    coi_corr = 1
                    if det in globals.coinc_lib_selection and coi:
                        coi_corr = corr_coinc.coinc_correction(nucl, en_lib[i], corr_coinc.coinc_lib(det), pos)
                    
                    # calculate efficiency
                    eff = 1 / coi_corr * (area * lamb * (t_real / t_live)) / (ref_act * int_lib[i]) * (1 / (np.exp(-lamb * t_delay) * (1 - np.exp(-lamb * t_real))))
                    eff_unc = eff * np.sqrt((ref_act_err / 100)**2 + (0.03 * ((t_real / t_live) - 1))**2 + t_delay**2 * lamb**2 * (hl_err_lib[i] / hl_lib[i])**2 + (area_unc / area)**2 + (int_err_lib[i] / int_lib[i])**2)

                    # add calculated efficiency to new row of globals.eff_results_df
                    # treeview_eff['columns'] = ['E (kev)', 'E_ref (kev)','Nuclide', 'Int (%)', 'A0 (Bq)', 'COI_corr', 'Eff (-)', 'dEff (-)']
                    globals.eff_results_df.loc[len(globals.eff_results_df)] = [
                        en,                        
                        en_lib[i],
                        nucl,
                        round(int_lib[i] * 100, 3),
                        ref_act,
                        round(coi_corr, 4),
                        round(eff, 6),
                        round(eff_unc, 6)
                        ]


def eff_func(en, *a):
    eff = 0
    for i in range(0, len(a)):
        eff += a[i] * np.log(en / 1000) ** i
    return eff


def calc_chi2(data_exp, data_fit):
    residuals = np.array(data_exp) - np.array(data_fit)
    return np.sum((residuals**2 / np.std(residuals)))


# plot efficiency figure
def plot_eff_fig(selection, div_en, scale, div_en_switch, det, pos):
    # check if results are empty
    if len(globals.eff_results_df) == 0 or len(selection) == 0:
        return tk.messagebox.showerror('User input error', 'No data for plotting.')
    
    # TODO the selection function is not general enough, most of treeview needs spectrums names in pathlib format
    selection = [float(str(i)) for i in selection]

    # clear variables
    globals.ax.clear()

    # filter data based 
    filtered_eff_df = globals.eff_results_df[globals.eff_results_df['E (keV)'].isin(selection)].copy()

    if len(filtered_eff_df) == 5:
        return tk.messagebox.showerror('User input error', 'Not enough data for the fit.')

    x_data = filtered_eff_df['E_ref (kev)'].to_list()
    y_data = filtered_eff_df['Eff (-)'].to_list()
    y_data_err = filtered_eff_df['dEff (-)'].to_list()
    y_data_err_rel = [i / j for i, j in zip(y_data_err, y_data)]

    # sort data by energy
    x_data, y_data, y_data_err, y_data_err_rel = zip(*sorted(zip(x_data, y_data, y_data_err, y_data_err_rel)))

    # plot errorbars
    globals.ax.errorbar(x_data, y_data, yerr=y_data_err, capthick='0.7', fmt='o', capsize=4, ecolor='k', color='k', markersize=3)
    
    # print chi2 values in the plot text
    chi2_list = []
    # removed fits from suspection of overparametrization
    err_fit = []

    # fit - single curve
    if not div_en_switch:
        for n in range(settings.MIN_S_PARAM, settings.MAX_S_PARAM):
            p0 = np.ones(n+1, )

            try:
                coeff, pcov = opt.curve_fit(eff_func, x_data, np.log(y_data), p0=p0, sigma=y_data_err_rel)
                coeff_err = np.sqrt(np.diag(pcov))
                y_fit = [np.exp(eff_func(i, *coeff)) for i in x_data]
                chi2 = calc_chi2(y_data, y_fit)
            except Exception as e:
                tk.messagebox.showerror('FEP efficiency curve fit error', f'Fit error: {e}')
                continue
            
            # append chi2 value to the list, including the degree and the div energy in string
            chi2_list.append(f'{n}. {chi2:.3f}\n')

            # plot fitted curve
            x_fit = np.linspace(25, 3000, 600)
            y_fit = np.exp(eff_func(x_fit, *coeff))

            # check if the fit is correct
            if max(y_fit) > 1:
                err_fit.append(f'Fit {n}. deg: incorrect fit.\n')
                chi2_list = chi2_list[:-1]
                continue

            # plot fitted curve
            globals.ax.plot(x_fit, y_fit, label=f'fit {n}. deg')
            
            # add results to globals.eff_fit_dict
            globals.eff_fit_dict[f'{det}_{pos}_n{str(n)}'] = [None, coeff, coeff_err]
    
    elif div_en_switch and div_en != '0.0':
        div_en = float(div_en)

        # Split the dataset into two subsets based on the dividing energy value
        for i in range(len(x_data)):
            if x_data[i] >= div_en:
                div_ind = i
                break

        x_data_1 = x_data[:div_ind]
        y_data_1 = y_data[:div_ind]
        y_data_err_rel_1 = y_data_err_rel[:div_ind]
        
        x_data_2 = x_data[div_ind:]
        y_data_2 = y_data[div_ind:]
        y_data_err_rel_2 = y_data_err_rel[div_ind:]

        # Fit each subset separately
        for n in range(settings.MIN_D1_PARAM, settings.MAX_D1_PARAM):
            if len(x_data_1) < n:
                print(f'Not enough data for the fit {n}. deg: (<= {div_en}).')
                continue

            p0 = np.ones(n+1, )

            try:
                # Fit the first subset
                coeff_1, pcov_1 = opt.curve_fit(eff_func, x_data_1, np.log(y_data_1), p0=p0, sigma=y_data_err_rel_1)
                coeff_err_1 = np.sqrt(np.diag(pcov_1))
                y_fit_1 = [np.exp(eff_func(i, *coeff_1)) for i in x_data_1]
                chi2 = calc_chi2(y_data_1, y_fit_1)
            except Exception as e:
                tk.messagebox.showerror('FEP efficiency curve fit error', f'Fit error: {e}')
                continue
            
            # append chi2 value to the list, including the degree and the div energy in string
            chi2_list.append(f'{n}. ($\\leq$ {div_en:.0f}): {chi2:.3f}\n')

            # Plot fitted curves for the first subset
            x_fit_1 = np.linspace(25, div_en, 300)
            y_fit_1 = np.exp(eff_func(x_fit_1, *coeff_1))
            
            # check if the fit is correct
            if max(y_fit_1) > 1:
                err_fit.append(f'Fit {n}. deg: (<= {div_en}) incorrect fit.\n')
                chi2_list = chi2_list[:-1]
                continue                

            globals.ax.plot(x_fit_1, y_fit_1, label=f'fit {n}. ($\\leq$ {div_en:.0f})')

            # Add results to globals.eff_fit_dict
            globals.eff_fit_dict[f'{det}_{pos}_n{str(n)}_1'] = [div_en, coeff_1, coeff_err_1]

        # Fit the second subset separately
        for n in range(settings.MIN_D2_PARAM, settings.MAX_D2_PARAM):
            if len(x_data_2) < n:
                print(f'Not enough data for the fit {n}. deg: (> {div_en}).')
                continue

            p0 = np.ones(n+1, )

            try:
                # Fit the second subset
                coeff_2, pcov_2 = opt.curve_fit(eff_func, x_data_2, np.log(y_data_2), p0=p0, sigma=y_data_err_rel_2)
                coeff_err_2 = np.sqrt(np.diag(pcov_2))
                y_fit_2 = [np.exp(eff_func(i, *coeff_2)) for i in x_data_2]
                chi2 = calc_chi2(y_data_2, y_fit_2)
            except Exception as e:
                tk.messagebox.showerror('FEP efficiency curve fit error', f'Fit error: {e}')
                continue
            
            # append chi2 value to the list, including the degree and the div energy in string
            chi2_list.append(f'{n}. (> {div_en:.0f}): {chi2:.3f}\n')

            # Plot fitted curves for the second subset
            x_fit_2 = np.linspace(div_en, 3000, 300)
            y_fit_2 = np.exp(eff_func(x_fit_2, *coeff_2))
            
            # check if the fit is correct
            if max(y_fit_2) > 1:
                err_fit.append(f'Fit {n}. deg: (> {div_en}) incorrect fit.\n')
                chi2_list = chi2_list[:-1]
                continue

            globals.ax.plot(x_fit_2, y_fit_2, label=f'fit {n}. (> {div_en:.0f})')

            # Add results to globals.eff_fit_dict
            globals.eff_fit_dict[f'{det}_{pos}_n{str(n)}_2'] = [div_en, coeff_2, coeff_err_2]


    # plot settings
    globals.ax.set(xlabel='$E$ (keV)', ylabel='$eff_{FEP}$ (-)', title=f'FEP efficiency calculation - {det} - {pos}')
    matplotlib.rcParams['axes.labelsize'] = 12
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')

    # set scale
    globals.ax.set_ylim([None, max(y_data)*1.2])

    if scale == 'xy-lin':
        globals.ax.set_xscale('linear')
        globals.ax.set_yscale('linear')
        globals.ax.set_xlim([0, 3000])
        text_x, text_y = 0.5, 0.95
    elif scale == 'xy-log':
        globals.ax.set_xscale('log')
        globals.ax.set_xlim([10, 4000])
        globals.ax.set_yscale('log')
        globals.ax.set_ylim([min(y_data) * 0.1, max(y_data) * 1.5])
        # add minor grid
        globals.ax.grid(which='minor', linestyle=':', linewidth='0.5')
        text_x, text_y = 0.15, 0.95
    elif scale == 'x-log_y-lin':
        globals.ax.set_xscale('log')
        globals.ax.set_xlim([10, 4000])
        globals.ax.set_yscale('linear')
        # add minor grid for log scale
        globals.ax.grid(which='minor', linestyle=':', linewidth='0.5')
        text_x, text_y = 0.15, 0.95
    
    globals.ax.legend(loc='best', fontsize='9')

    # add text with chi2 values
    chi2_list[-1] = chi2_list[-1].replace('\n', '')
    chi2_list = ['$\\chi^2$ values:\n'] + chi2_list
    globals.ax.text(text_x, text_y, ''.join(chi2_list), transform=globals.ax.transAxes, fontsize=9, verticalalignment='top',horizontalalignment='center', bbox=dict(alpha=0.75, facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))
    
    return globals.canvas_id.draw()


# fill option menu in export frame with available fits
def fill_export_opmenu(opm_sel, opm_sel2, div_en_switch):
    if not div_en_switch:
        opm_sel['values'] = list(globals.eff_fit_dict.keys())
        if len(opm_sel['values']) > 0:
            opm_sel.current(0)
    else:
        opm_sel['values'] = [key for key in globals.eff_fit_dict.keys() if '_1' in key]
        opm_sel2['values'] = [key for key in globals.eff_fit_dict.keys() if '_2' in key]
        if len(opm_sel['values']) > 0 and len(opm_sel2['values']) > 0:
            opm_sel.current(0)
            opm_sel2.current(0)


def update_eff_lib(det, pos, lines):
    # open the library file and search for the selected detector and position
    with open(globals.lib_selection['Efficiency'], 'r', encoding='utf-8') as file:
        content = file.readlines()

    for i, line in enumerate(content):
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue

        line = line.split()
        if len(line) >= 5:
            if line[0] == det and line[1] == pos:
                break

    # update the library file with the new data
    content[i] = '# ' + content[i]
    content[i+1] = '# ' + content[i+1]
    if float(line[2]) > 0:
        content[i+2] = '# ' + content[i+2]
    
    for j, line in enumerate(lines):
        content.insert(i+3+j, line)
    
    content.insert(i+3+len(lines), '\n')
    
    with open(globals.lib_selection['Efficiency'], 'w', encoding='utf-8') as file:
        for row in content:
            file.write(row)
    
    return tk.messagebox.showinfo('Library update', f'Efficiency data for {det} - {pos} was successfully updated in the library.')


# helper function for saving coefficients of the fit
def eff_line_func(coeff, coeff_err):
    eff_line = ''
    err_line = '# '
    for i in range(0,9):
        if i < len(coeff):
            eff_line += f'{coeff[i]:.6e}\t'
            err_line += f'{coeff_err[i]:.6e}\t'
        else:
            eff_line += '0\t'
            err_line += '0\t'
    
    eff_line += '\n'
    err_line += '# fit error\n'
    
    return eff_line, err_line


# save coefficients of the fit
def save_eff_fit(det, pos, div_en_switch, sel, sel2):
    # check if results dictionary is empty
    if len(globals.eff_fit_dict) == 0:
        return tk.messagebox.showerror('User input error', 'No data for saving.')
    
    # check if selected fit is value from dictionary
    if not div_en_switch:
        if sel not in globals.eff_fit_dict:
            return tk.messagebox.showerror('User input error', 'Selected fit is not in the dictionary.')
    else:
        if sel not in globals.eff_fit_dict or sel2 not in globals.eff_fit_dict:
            return tk.messagebox.showerror('User input error', 'Selected fit is not in the dictionary.')

    # prepare list of lines for saving
    lines = []
    timestamp = datetime.now().strftime('%d-%m-%Y %H:%M')
    lines.append(f'# {timestamp} - eff_FEP added from GESP\n')
    
    if not div_en_switch:
        en_div, coeff, coeff_err = globals.eff_fit_dict[sel]
        # comment line
        # example: b   p1  0 keV   MeV
        lines.append(f'{det}  {pos}  0 keV   MeV\n')
        
        eff_line, err_line = eff_line_func(coeff, coeff_err)
        
        lines.append(eff_line)
        lines.append(err_line)
    else:
        en_div, coeff_1, coeff_err_1 = globals.eff_fit_dict[sel]
        en_div, coeff_2, coeff_err_2 = globals.eff_fit_dict[sel2]
        # comment line
        lines.append(f'{det}  {pos}  {en_div} keV   MeV\n')
        
        eff_line1, err_line1 = eff_line_func(coeff_1, coeff_err_1)
        eff_line2, err_line2 = eff_line_func(coeff_2, coeff_err_2)
                
        lines.append(eff_line1)
        lines.append(eff_line2)
        lines.append(err_line1)
        lines.append(err_line2)

    # save the lines to the file    
    
    # check if eff library already exists    
    if pathlib.Path(globals.lib_selection['Efficiency']).is_file():
        try:
            en_div, a1, a2, unit = corr_eff.eff_lib(det, pos)
            tk.messagebox.showinfo('Library info', f'Efficiency data for {det} - {pos} already exists in the library.\n Div. energy: {en_div} keV\n Fit 1: {a1}\n Fit 2: {a2}\n Unit: {unit}')
            
        except ValueError:
            en_div = None
            a1 = None
            a2 = None
            unit = None

        if en_div is not None and a1 is not None and a2 is not None and unit is not None:
            # locate the line with the selected detector and position, comment it and add the new data
            try:
                update_eff_lib(det, pos, lines)
            except Exception as e:
                tk.messagebox.showerror('Library error', f'Error during updating the library file. Error: {e}')
        else:
            # open the library file and add the new data
            with open(globals.lib_selection['Efficiency'], 'a', encoding='utf-8') as file:
                file.write('\n')
                for line in lines:
                    file.write(line)
                file.write('\n')
            
            tk.messagebox.showinfo('Library info', f'Efficiency data for {det} - {pos} was successfully saved in the library.')
    else:
        # create new library file and add the new data
        globals.lib_selection['Efficiency'] = pathlib.Path('libraries/eff_test.lib')
        try:
            with open(globals.lib_selection['Efficiency'], 'w', encoding='utf-8') as file:
                file.write('# det_name\tposition\tdividing_en (in case of two curves) keV\tfit energy denominator (keV/MeV/dubna)\n')
                file.write('# a0\t\t\ta1\t\t\ta2\t\t\ta3\t\t\ta4\t\t\ta5\t\t\ta6\t\t\ta7\t\t\ta8\n\n')
                for line in lines:
                    file.write(line)
                file.write('\n')
        except Exception as e:
            tk.messagebox.showerror('Library error', f'Error during creating the library file. Error: {e}')
        
        tk.messagebox.showinfo('Library info', f'Efficiency data for {det} - {pos} was successfully saved in new library file.')
