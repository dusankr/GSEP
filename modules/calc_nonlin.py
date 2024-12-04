# This module is used for calculation of nonlinearity of the detector. It is possible to calculate nonlinearity for
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
# TODO: if not enough points for plot, give infowindow, not error

# Libraries
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from ttkwidgets import CheckboxTreeview
import datetime
import pathlib
# Local modules
from modules import globals
from modules import utilities
from modules import settings


pd.options.mode.chained_assignment = None

# window for calculation of effectivity
def nonlinearity_win(root):
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

    # detector list
    det_list = list(set(globals.eff_param_collection['Detector']))
    det_list.sort()

    # Definition of tkinter variables ---------------------------------------------------------------------------
    opmenu_det = tk.StringVar()
    opmenu_pos = tk.StringVar()
    var_curve = tk.BooleanVar(value=False)
    div_en_var = tk.StringVar()

    # new window ------------------------------------------------------------------------------------------------
    new_win = tk.Toplevel(root)
    new_win.title('Nonlinearity calculation')
    
    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))
    new_win.resizable(False, False)
    new_win.grid_propagate(True)

    # frames
    top_frame = tk.ttk.Frame(new_win)
    top_frame.grid(sticky='wens', column=0, columnspan=2, row=0, padx=5, pady=5)

    tree_frame = tk.ttk.Frame(new_win, width=400)
    tree_frame.grid(sticky='wens', column=0, row=1, padx=5, pady=5)
    tree_frame.columnconfigure(0, weight=1)
    tree_frame.rowconfigure(0, weight=1)
    
    plot_frame = tk.ttk.Frame(new_win)
    plot_frame.grid(sticky='wens', column=1, row=1, padx=5, pady=5)
    plot_frame.columnconfigure(0, weight=1)
    plot_frame.rowconfigure(0, weight=1)
   
    bottom_frame = tk.ttk.Frame(new_win)
    bottom_frame.grid(sticky='wens', column=0, columnspan=2, row=2, padx=5, pady=5)   

    # creates menu of positions for selected detector
    def update_opmenu_pos(det):
        position = []
        opm_pos['values'] = []
        for i in range(0, len(globals.eff_param_collection['Detector'])):
            if globals.eff_param_collection['Detector'][i] == str(det):
                position.append(globals.eff_param_collection['Position'][i])
        pos_list = list(set(position))
        pos_list.sort()
        opm_pos['values'] = pos_list
        opm_pos.current(0)

    # selection -------------------------------------------------------------------------------------------------
    input_frame = tk.LabelFrame(top_frame, text='Input settings')
    input_frame.grid(sticky='wens', column=0, row=0, padx=2, pady=2)

    label_title = tk.ttk.Label(input_frame, text='Calculation for detector')
    label_title.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    opm_det = tk.ttk.Combobox(input_frame, textvariable=opmenu_det, width=5, state='readonly')
    opm_det.grid(sticky='w', column=1, row=0, padx=5, pady=5)
    opm_det['values'] = det_list
    opm_det.current(0)

    label_pos = tk.ttk.Label(input_frame, text='position')
    label_pos.grid(sticky='w', column=2, row=0, padx=5, pady=5)
    
    opm_pos = tk.ttk.Combobox(input_frame, textvariable=opmenu_pos, width=5, state='readonly')
    opm_pos.grid(sticky='w', column=3, row=0, padx=5, pady=5)
    opm_pos['values'] = update_opmenu_pos(det_list[0])
    
    opm_det.bind('<<ComboboxSelected>>', lambda event: update_opmenu_pos(opmenu_det.get()))
    
    calc_button = tk.ttk.Button(input_frame, text='Calculate', width=15, command=lambda: [
        nonlin_calc(opmenu_det.get(), opmenu_pos.get()),
        utilities.treeview_update(treeview_file, globals.nonlin_results_df),
        nonlin_plot(selected_tally_test(treeview_file), div_en_var.get(), var_curve.get(), opmenu_det.get(), opmenu_pos.get())
        ])
    calc_button.grid(sticky='w', column=4, row=0, padx=5, pady=5)

    # treeview ---------------------------------------------------------------------------------------------------
    treeview_file = CheckboxTreeview(tree_frame)
    treeview_file.grid(sticky='wens', column=0, row=0, padx=5, pady=5)

    treeview_file['columns'] = ['E_lib (keV)', 'Nuclide', 'Int_lib (%)', 'E (keV)', 'Diff (keV)']
    treeview_file['show'] = 'headings', 'tree'
    
    # set columns filters
    for col in treeview_file['columns']:
        treeview_file.heading(col, text=col, anchor='center', command=lambda _col=col: utilities.treeview_sort_column(treeview_file, _col, False))
    
    # set width of columns
    for key, value in settings.COLUMN_WIDTH.items():
        if key in treeview_file['columns']:
            treeview_file.column(key, width=value, stretch=True, anchor='center')
    
    # set width of zero column
    treeview_file.column('#0', width=45, stretch=False, anchor='w')

    # Treeview Y-scrollbar
    tree_y_scr = tk.ttk.Scrollbar(tree_frame, orient='vertical', command=treeview_file.yview)
    tree_y_scr.grid(sticky='ns', column=1, row=0)

    # Treeview X-scrollbar
    tree_x_scr = tk.ttk.Scrollbar(tree_frame, orient='horizontal', command=treeview_file.xview)
    tree_x_scr.grid(sticky='we', column=0, row=1)

    treeview_file.configure(yscrollcommand=tree_y_scr.set)
    treeview_file.configure(xscrollcommand=tree_x_scr.set)

    # Figure -----------------------------------------------------------------------------------------------------

    # empty figure
    globals.fig_id = matplotlib.figure.Figure(tight_layout=True, figsize=settings.FIG_SIZE)
    globals.ax = globals.fig_id.add_subplot()
    # plot settings
    matplotlib.rcParams['axes.facecolor'] = 'white'
    globals.ax.set(xlabel='$E$ (keV)', ylabel='$\Delta E$ (keV)', title=f'Non linearity calculation - {opmenu_det.get()} - {opmenu_pos.get()}')
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

    # Button menu ------------------------------------------------------------------------------------------------

    fit_frame = tk.LabelFrame(bottom_frame, text='Fit settings')
    fit_frame.grid(sticky='wens', column=0, row=0, padx=5, pady=5)
    
    single_curve_button = tk.ttk.Radiobutton(fit_frame, text="Single curve", variable=var_curve, value=False, command=lambda: [
        div_en_entry.delete(0, 'end'),
        div_en_entry.configure(state="disabled")
        ])
    single_curve_button.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    two_curves_button = tk.ttk.Radiobutton(fit_frame, text="Two curves", variable=var_curve, value=True, command=lambda: 
                                           div_en_entry.configure(state="active"))
    two_curves_button.grid(sticky='w', column=1, row=0, padx=5, pady=5)

    div_en_label = tk.ttk.Label(fit_frame, text='Div. energy')
    div_en_label.grid(sticky='w', column=2, row=0, padx=5, pady=5)

    div_en_entry = tk.ttk.Entry(fit_frame, width=6, state="disabled", textvariable=div_en_var)
    div_en_entry.grid(sticky='w', column=3, row=0, padx=5, pady=5)

    unit_en_label = tk.ttk.Label(fit_frame, text='keV')
    unit_en_label.grid(sticky='w', column=4, row=0, padx=5, pady=5)

    plot_button = tk.ttk.Button(fit_frame, text='Plot selection', width=15, command=lambda: 
        nonlin_plot(selected_tally_test(treeview_file), div_en_var.get(), var_curve.get(), opmenu_det.get(), opmenu_pos.get())
    )
    plot_button.grid(sticky='w', column=5, row=0, padx=5, pady=5)

    # results frame
    results_frame = tk.LabelFrame(bottom_frame, text='Results')
    results_frame.grid(sticky='w', column=1, row=0, padx=5, pady=5)

    save_button = tk.ttk.Button(results_frame, text='Save coefficients', width=15, command=lambda: save_results())
    save_button.grid(sticky='w', column=0, row=0, padx=5, pady=5)

    button_close = tk.ttk.Button(bottom_frame, text='Close window', width=15, command=lambda: utilities.close_win(root, new_win))
    button_close.grid(sticky='wsne', column=2, row=0, padx=5, pady=5)


# calculate nonlinearity - go through all spectra and peaks, calculate difference between library and measured energy
def nonlin_calc(det, pos):
    # clear globals.nonlin_results_df
    globals.nonlin_results_df = globals.nonlin_results_df.iloc[0:0]

    # get list of spectra for selected detector and position
    temp_df = globals.eff_param_collection[(globals.eff_param_collection['Detector'] == det) & 
                                           (globals.eff_param_collection['Position'] == pos)]    

    # check if there are any spectra
    if len(temp_df) == 0:
        return tk.messagebox.showerror('User input error', f'No valid data for non-linearity calculation of {det} - {pos}.')
    
    spectra_list = temp_df['Name'].tolist()
    nuclide_list = temp_df['Sample name'].tolist()

    # go through all spectra and all peaks and calculate difference between library and measured energy
    for spectrum, nucl in zip(spectra_list, nuclide_list):

        # get values from gamma library for current nuclide
        try:
            dict_array = globals.gamma_lib_eff[nucl]
            en_lib = dict_array[0]
            en_err_lib = dict_array[1]
            int_lib = dict_array[2]
            #int_err_lib = dict_array[3]
        except:
            tk.messagebox.showerror('Library error', f'No data for nuclide {nucl} in gamma library. Program will continue with next nuclide.')
            continue
            
        # get measured peaks from spectrum
        en = globals.spectra_collection[spectrum.stem]['energy']
        channel = globals.spectra_collection[spectrum.stem]['channel']
        chn_err = globals.spectra_collection[spectrum.stem]['cherr']
        #fwhm = globals.spectra_collection[spectrum.stem]['fwhm']
        #area_err_rel = globals.spectra_collection[spectrum.stem]['aerr%']     

        # calculate difference between library and measured energy
        for i in range(0, len(en)):
            for j in range(0, len(en_lib)):
                if abs(en_lib[j] - en[i]) < (settings.N_CALC * np.sqrt(settings.A_CALC**2 + (en_err_lib[j] / en_lib[j])**2 + (chn_err[i] / channel[i])**2 + 1)):
                    diff = en_lib[j] - en[i]
                    globals.nonlin_results_df.loc[len(globals.nonlin_results_df)] = [
                        en_lib[j],
                        nucl,
                        round(int_lib[j] * 100, 2),
                        en[i],
                        round(diff, 4)
                    ]
                    break


# creates the plot
def nonlin_plot(selection, div_en, div_en_switch, det, pos):
    # check if results are empty
    if len(globals.nonlin_results_df) == 0 or len(selection) == 0:
        return tk.messagebox.showerror('User input error', 'No data for plotting.')

    # clear figure
    globals.ax.clear()
    globals.coinc_fit_res = None

    # convert to float
    en_sel = [float(i) for i in selection[0]]
    diff_sel = [float(i) for i in selection[1]]    

    # filter results based on selection
    temp_df = globals.nonlin_results_df[globals.nonlin_results_df['E_lib (keV)'].isin(list(en_sel)) & globals.nonlin_results_df['Diff (keV)'].isin(list(diff_sel))]

    if len(temp_df) == 0:
        return tk.messagebox.showerror('User input error', 'No data for plotting.')
    
    # plot all selected points
    x_data = temp_df['E (keV)']
    y_data = temp_df['Diff (keV)']
    globals.ax.plot(x_data, y_data, 'o', label='selected', color='blue', markersize=2)

    # fit single curve
    if not div_en_switch:
        # try polynomial fit
        try:
            popt = np.polyfit(x_data, y_data, 2, cov=False)
            x_fit = np.linspace(25, 3000, 1000)
            y_fit = np.polyval(popt, x_fit)
        except Exception as e:
            return tk.messagebox.showerror('Fitting error', f'Fit procedure failed with error: {e}')

        globals.ax.plot(x_fit, y_fit, label='fit', color='red')

        # save results
        globals.coinc_fit_res = [det, None, popt, None]

    # fit two curves
    elif div_en_switch and div_en != '':
        # divide energy
        try:
            div_en = float(div_en)
        except:
            return tk.messagebox.showerror('User input error', 'Dividing energy must be a number.')

        # divide data
        x_data_1 = x_data[x_data <= div_en]
        y_data_1 = y_data[x_data <= div_en]
        x_data_2 = x_data[x_data > div_en]
        y_data_2 = y_data[x_data > div_en]

        # check if there are enough points for fitting
        if len(x_data_1) < 3 or len(x_data_2) < 3:
            return tk.messagebox.showerror('User input error', 'Not enough points for fitting.')

        # try polynomial fit
        try:
            popt_1 = np.polyfit(x_data_1, y_data_1, 2, cov=False)
            popt_2 = np.polyfit(x_data_2, y_data_2, 2, cov=False)
            x_fit_1 = np.linspace(25, div_en, 1000)
            x_fit_2 = np.linspace(div_en, 3000, 1000)
            y_fit_1 = np.polyval(popt_1, x_fit_1)
            y_fit_2 = np.polyval(popt_2, x_fit_2)
        except Exception as e:
            return tk.messagebox.showerror('Fit error', f'Fit procedure failed with error: {e}')

        globals.ax.plot(x_fit_1, y_fit_1, label=f'fit $\\leq$ {div_en}', color='red')
        globals.ax.plot(x_fit_2, y_fit_2, label=f'fit $>$ {div_en}', color='green')

        # save results
        globals.coinc_fit_res = [det, div_en, popt_1, popt_2]

    # plot settings
    globals.ax.set(xlabel='$E$ (keV)', ylabel='$\Delta E$ (keV)', title=f'Non linearity calculation - {det} - {pos}')
    matplotlib.rcParams['axes.labelsize'] = 12
    globals.ax.grid(which='major', linestyle='-', linewidth='0.5')
    globals.ax.legend(loc='best', fontsize='9')
    # range of x-axis
    globals.ax.set_xlim(0, 3000)

    return globals.canvas_id.draw()


def selected_tally_test(treeview_files):
    if len(treeview_files.get_checked()) != 0:
        # send selected tallies to plot_mod function
        en_sel = []
        diff_sel = []
        for row in treeview_files.get_checked():
            en_sel.append(treeview_files.item(row)['values'][0])
            diff_sel.append(treeview_files.item(row)['values'][4])
    else:
        tk.messagebox.showerror('Input error', 'No results were selected.')
        return None
    
    if len(en_sel) == 0 or len(diff_sel) == 0:
        return None
    else:   
        return [en_sel, diff_sel]


# save coefficients to eff library
def save_results():
    if globals.coinc_fit_res is None:
        return tk.messagebox.showerror('User input error', 'No fit results to save.')
    
    det = globals.coinc_fit_res[0]
    energy = globals.coinc_fit_res[1]
    coeff1 = globals.coinc_fit_res[2]
    coeff2 = globals.coinc_fit_res[3]

    if energy is None:
        energy = 0

    nonlin_content = []
    nonlin_content.append('# det_name	dividing_en (in case of two curves)\n')
    nonlin_content.append('# a0\t\t\ta1\t\t\ta2\n')
    nonlin_content.append(f'{det}\t{energy} keV\n')

    if energy is None or energy == 0:
        nonlin_content.append(f'{coeff1[0]:.5f}\t{coeff1[1]:.5f}\t{coeff1[2]:.5f}\n')
    else:
        nonlin_content.append(f'{coeff1[0]:.5f}\t{coeff1[1]:.5f}\t{coeff1[2]:.5f}\n')
        nonlin_content.append(f'{coeff2[0]:.5f}\t{coeff2[1]:.5f}\t{coeff2[2]:.5f}\n')
    
    # check if output directory exists
    utilities.check_folders()

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    file_name = globals.common_params['Output path'] / pathlib.Path(f'{timestamp}_nonlinearity.lib')

    with open(file_name, 'w') as file:
        for line in nonlin_content:
            file.write(line)

    return tk.messagebox.showinfo('Save results', f'Coefficients saved to {file_name}.')
