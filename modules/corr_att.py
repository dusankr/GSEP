# self-absorption correction coefficient calculation
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
# TODO: Rewrite fit from library values to fit calculated attenuation coefficients based on thickness, density and material (lib values)

# Libraries
import pathlib
import tkinter as tk
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
# Local modules
from modules import globals
from modules import calc_corr
from modules import inputs
from modules import utilities


# new window for user input of sample parameters - weight, molar mass, thickness, material, weight fraction, density
def new_win_sample(root, treeview_files):

    if 'Name' not in globals.spectra_param_collection:
        return tk.messagebox.showerror('User input error', 'Spectra were not upload.')

    # new window
    new_win = tk.Toplevel(root)
    new_win.title('User input - sample parameters')
    #new_win.columnconfigure(0, weight=1)
    #new_win.columnconfigure(0, weight=1)

    # block a background window
    new_win.grab_set()
    new_win.protocol('WM_DELETE_WINDOW', lambda: utilities.close_win(root, new_win))

    # turn off resizing
    new_win.resizable(False, False)

    full_frame = tk.ttk.LabelFrame(new_win, relief='groove', borderwidth=2, text='Sample parameters', labelanchor='nw')
    full_frame.grid(sticky='wesn', column=0, row=0, columnspan=2, padx=5, pady=5)
    sample_lab = tk.ttk.Label(full_frame, text='Sample')
    sample_lab.grid(sticky='W', column=0, row=0)
    weight_lab = tk.ttk.Label(full_frame, text='Weight')
    weight_lab.grid(sticky='W', column=1, row=0)
    weight_err_lab = tk.ttk.Label(full_frame, text='Weight err')
    weight_err_lab.grid(sticky='W', column=2, row=0)
    molar_err_lab = tk.ttk.Label(full_frame, text='Molar mass')
    molar_err_lab.grid(sticky='W', column=3, row=0)
    mat_lab = tk.ttk.Label(full_frame, text='Material')
    mat_lab.grid(sticky='W', column=4, row=0)
    wei_lab = tk.ttk.Label(full_frame, text='Frac. by weight')
    wei_lab.grid(sticky='W', column=5, row=0)
    den_lab = tk.ttk.Label(full_frame, text='Density')
    den_lab.grid(sticky='W', column=6, row=0)
    thick_lab = tk.ttk.Label(full_frame, text='Thickness')
    thick_lab.grid(sticky='W', column=7, row=0)
    weight__lab = tk.ttk.Label(full_frame, text='(g)')
    weight__lab.grid(sticky='W', column=1, row=1)
    weight__err_lab = tk.ttk.Label(full_frame, text='(g)')
    weight__err_lab.grid(sticky='W', column=2, row=1)
    molar__lab = tk.ttk.Label(full_frame, text='(-)')
    molar__lab.grid(sticky='W', column=3, row=1)
    mat__lab = tk.ttk.Label(full_frame, text='(-)')
    mat__lab.grid(sticky='W', column=4, row=1)
    wei__lab = tk.ttk.Label(full_frame, text='(-)')
    wei__lab.grid(sticky='W', column=5, row=1)
    den__lab = tk.ttk.Label(full_frame, text='(g/cm^3)')
    den__lab.grid(sticky='W', column=6, row=1)
    thick__lab = tk.ttk.Label(full_frame, text='(mm)')
    thick__lab.grid(sticky='W', column=7, row=1)

    # local variables
    weight = []
    weight_err = []
    molar = []
    thickness = []
    elements = []
    weight_frac = []
    density = []
    samples = sorted(tuple(set(globals.spectra_param_collection['Sample name'].tolist())))

    for i in range(0, len(samples)):
        sample = samples[i]

        lab = tk.ttk.Label(full_frame, text=sample)
        lab.grid(sticky='W', column=0, row=i + 2)
        ent_weight = tk.ttk.Entry(full_frame, width=10)
        ent_weight.grid(sticky='W', column=1, row=i + 2, padx=1, pady=1)
        ent_weight_err = tk.ttk.Entry(full_frame, width=10)
        ent_weight_err.grid(sticky='W', column=2, row=i + 2, padx=1, pady=1)
        ent_molar = tk.ttk.Entry(full_frame, width=10)
        ent_molar.grid(sticky='W', column=3, row=i + 2, padx=1, pady=1)
        ent_thickness = tk.ttk.Entry(full_frame, width=10)
        ent_thickness.grid(sticky='W', column=7, row=i + 2, pady=1)
        ent_weight_f = tk.ttk.Entry(full_frame, width=10)
        ent_weight_f.grid(sticky='W', column=5, row=i + 2, padx=1, pady=1)
        ent_dens = tk.ttk.Entry(full_frame, width=10)
        ent_dens.grid(sticky='W', column=6, row=i + 2, padx=1, pady=1)
        ent_element = tk.ttk.Entry(full_frame, width=10)
        ent_element.grid(sticky='W', column=4, row=i + 2, padx=1, pady=1)

        # fill entries with values from the collection
        for j in range(0, len(globals.spectra_param_collection['Sample name'])):
            if globals.spectra_param_collection.loc[j, 'Sample name'] == sample:
                try:
                    ent_weight.insert(0, globals.spectra_param_collection.loc[j, 'Weight'])
                    ent_weight_err.insert(0, globals.spectra_param_collection.loc[j, 'Weight error'])
                    ent_molar.insert(0, globals.spectra_param_collection.loc[j, 'Molar mass'])
                    ent_element.insert(0, globals.spectra_param_collection.loc[j, 'Element'])
                    ent_weight_f.insert(0, globals.spectra_param_collection.loc[j, 'Fraction by weight'])
                    ent_dens.insert(0, globals.spectra_param_collection.loc[j, 'Density'])
                    ent_thickness.insert(0, globals.spectra_param_collection.loc[j, 'Thickness'])
                except Exception as e:
                    print(e)
                    ent_weight.insert(0, 0)
                    ent_weight_err.insert(0, 0)
                    ent_molar.insert(0, 0)
                    ent_element.insert(0, '90Th')
                    ent_weight_f.insert(0, '-')
                    ent_dens.insert(0, '-')
                    ent_thickness.insert(0, '3')
                break

        weight.append(ent_weight)
        weight_err.append(ent_weight_err)
        molar.append(ent_molar)
        elements.append(ent_element)
        weight_frac.append(ent_weight_f)
        density.append(ent_dens)
        thickness.append(ent_thickness)

    sample_parameters_dict = {'Weight': [], 'Weight err': [], 'Molar': [], 'Thickness': [], 'Weight frac': [], 'Density': [], 'Elements': []}

    # reload function
    def reload():
        for k in range(0, len(samples)):
            sample_parameters_dict['Weight'].append(weight[k].get())
            sample_parameters_dict['Weight err'].append(weight_err[k].get())
            sample_parameters_dict['Molar'].append(molar[k].get())
            sample_parameters_dict['Thickness'].append(thickness[k].get())
            sample_parameters_dict['Density'].append(density[k].get())
            # check if ' ' is in the element name or fraction by weight, if yes, just delete it
            sample_parameters_dict['Weight frac'].append((weight_frac[k].get()).replace(' ', ''))
            sample_parameters_dict['Elements'].append((elements[k].get()).replace(' ', ''))
        return sample_parameters_dict

    button_load = tk.ttk.Button(new_win, text='Save', width=25, command=lambda: [
        reload(),
        inputs.upload_weights(root, new_win, samples, sample_parameters_dict),
        inputs.upload_thickness(root, new_win, samples, sample_parameters_dict),
        inputs.treeview_fill(treeview_files),
        calc_corr.write_all_info()
    ])
    button_load.grid(sticky='wesn', column=0, row=1, padx=5, pady=5)

    button_quit = tk.ttk.Button(new_win, text='Quit', width=25, command=lambda: utilities.close_win(root, new_win))
    button_quit.grid(sticky='wesn', column=1, row=1, padx=5, pady=5)


# control of checkbox for attenuation correction - change state based on thickness column in sample parameters
def check_att():
    if 'Thickness' in globals.spectra_param_collection:
        if any(globals.spectra_param_collection['Thickness']) == 0:
            return 'disable'
        else:
            return 'active'
    else:
        return 'disable'


# read attenuation coefficients from library and send it to main att function
def att_lib():
    with open(globals.lib_selection['Attenuation'], 'r', encoding='utf-8') as temporary_file:
        content = temporary_file.readlines()

    att_lib_param = {}
    att_lib_values = {}
    
    # read elements parameters from attenuation library
    for i in range(0, len(content)):
        line = content[i].strip().split()
        if len(line) > 0 and content[i][0] != '#':
            if line[0][-1] == 'C':
                # element name, density, divide energy
                att_lib_param[line[0][:-1]] = [float(line[2]), float(line[4])]
    
    # read elements energy and attenuation coefficients from attenuation library
    for key in att_lib_param:
        en = []
        mu = []
        for i in range(0, len(content)):
            line = content[i].strip().split()
            if len(line) > 0 and content[i][0] != '#' and line[0] == key:
                en.append(float(line[1]))
                mu.append(float(line[2]))
        att_lib_values[key] = [en, mu]

    return att_lib_param, att_lib_values


# function for attenuation fit - lib en. values are already in MeV
def att_func(en, *a):
    fit = 0
    for i in range(0, len(a)):
        fit += a[i] * np.log(en)**i
    return fit


# final attenuation values calculation
def att_calc(a, en):
    mu_m = 0
    for i in range(0, len(a)):
        mu_m += a[i] * np.log(en/1000)**i
    return np.exp(mu_m)


# user check - fit library values and create a figure for each element
def att_figures(div_en, att_lib_values, coeff, elem):
    name = f'att_fit_{elem}.png'

    # check if the output folder exists
    utilities.check_folders()

    # check if the output figure already exists
    figure_path = globals.common_params['Output path'] / pathlib.Path(name)
    if figure_path.exists():
        return
    
    fig, ax = plt.subplots(figsize=(20/2.54, 15/2.54))
    ax.set_xlabel('Energy (MeV)')
    ax.set_ylabel(r'$\mu$/$\rho$ (cm$^2$/g)')
    ax.grid()
    ax.set_yscale('log')
    ax.set_xscale('log')

    # plot attenuation library values
    att_err_teor = [0.1 * i for i in att_lib_values[1]]
    ax.errorbar(att_lib_values[0], att_lib_values[1], yerr=att_err_teor, fmt='.', label='Att. data for {}'.format(elem))

    if div_en == 0:
        ax.plot(np.linspace(min(att_lib_values[0]), max(att_lib_values[0]), 500), np.exp(att_func(np.linspace(min(att_lib_values[0]), max(att_lib_values[0]), 500), *coeff)), label='fit', lw=0.7)
    else:
        for i in range(0, len(att_lib_values[0])):
            if att_lib_values[0][i] > div_en:
                break
        
        ax.plot(np.linspace(min(att_lib_values[0][:i]), div_en, 500), np.exp(att_func(np.linspace(min(att_lib_values[0][:i]), div_en, 500), *coeff[0])), label='fit: en < {:.4f} MeV'.format(div_en), lw=0.7)
        # get latest color from the plot
        color = ax.get_lines()[-1].get_color()        
        ax.plot(np.linspace(div_en, max(att_lib_values[0][i:]), 500), np.exp(att_func(np.linspace(div_en, max(att_lib_values[0][i:]), 500), *coeff[1])), label='fit: en > {:.4f} MeV'.format(div_en), lw=0.7, color=color)
    
    ax.legend(loc='best', fontsize=10)
    fig.savefig(figure_path, dpi=150)
    
    return plt.close(fig)


# function for determination of fit coefficients
def att_fit(att_lib_param, att_lib_values, elem):
    div_en = att_lib_param[1]

    if div_en == 0:
        p0 = np.ones(7, )

        # 10 % error for each point from attenuation library
        #att_err_teor = [0.1 * i for i in att_lib_values[1]]

        #coeff, pcov = opt.curve_fit(att_func, att_lib_values[0], np.log(att_lib_values[1]), p0=p0, sigma=att_err_teor, absolute_sigma=True)
        coeff, pcov = opt.curve_fit(att_func, att_lib_values[0], np.log(att_lib_values[1]), p0=p0)

        att_figures(div_en, att_lib_values, coeff, elem)      

        return coeff, 0    
    elif div_en > 0:
        for i in range(0,len(att_lib_values[0])):
            if att_lib_values[0][i] > div_en:
                break
        
        # check if first or second part of dataset has 7 or more points for fit, else return error or reduce function order
        if len(att_lib_values[0][:i]) < 7 and len(att_lib_values[0][:i]) < 3:
            return tk.messagebox.showerror('User input error', 'Attenuation library does not include enough data points for fit. Problematic dataset: {}'.format(att_lib_values[0][:i]))
        elif len(att_lib_values[0][:i]) < 7 and len(att_lib_values[0][:i]) >= 3:
            p0 = np.ones(len(att_lib_values[0][:i]), )
            coeff_1, pcov_1 = opt.curve_fit(att_func, att_lib_values[0][:i], np.log(att_lib_values[1][:i]), p0=p0)
        else:
            p0 = np.ones(7, )
            coeff_1, pcov_1 = opt.curve_fit(att_func, att_lib_values[0][:i], np.log(att_lib_values[1][:i]), p0=p0)
        
        if len(att_lib_values[0][i:]) < 7 and len(att_lib_values[0][i:]) < 3:
            return tk.messagebox.showerror('User input error', 'Attenuation library does not include enough data points for fit. Problematic dataset: {}'.format(att_lib_values[0][:i]))
        elif len(att_lib_values[0][i:]) < 7 and len(att_lib_values[0][i:]) >= 3:
            p0 = np.ones(len(att_lib_values[0][i:]), )
            coeff_2, pcov_2 = opt.curve_fit(att_func, att_lib_values[0][i:], np.log(att_lib_values[1][i:]), p0=p0)
        else:
            p0 = np.ones(7, )
            coeff_2, pcov_2 = opt.curve_fit(att_func, att_lib_values[0][i:], np.log(att_lib_values[1][i:]), p0=p0)
            
        att_figures(div_en, att_lib_values, [coeff_1, coeff_2], elem)

        return coeff_1, coeff_2


# attenuation main correction function
def attenuation():
    
    # read attenuation coefficients from library + divide energy and density
    att_lib_param, att_lib_values = att_lib()

    for spectrum in globals.spectra_collection:
        spec = pathlib.Path(spectrum + '.prn')
        element = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Element'].item()
        fraction = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Fraction by weight'].item()
        density = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Density'].item()
        thickness = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Thickness'].item()

        # switch element/compound
        if '/' in element:
            # search for element in the library
            elements = [x.strip() for x in element.split('/')]
            for x in range(0, len(elements)):
                if elements[x] not in att_lib_param.keys():
                    return tk.messagebox.showerror('User input error', 'Attenuation library does not include attenuation coefficients for' + elements[x])
            attenuation_compound(spectrum, element, thickness, fraction, density)
        else:
            if element not in att_lib_param.keys():
                return tk.messagebox.showerror('User input error', 'Attenuation library does not include attenuation coefficients for'+element)
            attenuation_element(spectrum, element, thickness, att_lib_param, att_lib_values)
    
    print('Attenuation correction done.')


# attenuation correction function for elements
def attenuation_element(spec, elem, thick, att_lib_param, att_lib_values):
    # density, thickness
    density = att_lib_param[elem][0]
    en_div_kev = att_lib_param[elem][1] * 1000

    # fit coefficients for sample element
    a1, a2 = att_fit(att_lib_param[elem], att_lib_values[elem], elem)
      
    # attenuation correction calculation
    att_cor = []

    if 'cor_en' in globals.spectra_collection[spec]:
        energy = globals.spectra_collection[spec]['cor_en']
    else:
        energy = globals.spectra_collection[spec]['energy']

    for en in energy:
        a = utilities.a_switch(en, en_div_kev, a1, a2)
        mu = att_calc(a, en) * density
        att_cor.append((1 - np.exp(-1 * mu * thick / 10)) / (mu * thick / 10))  # en. in keV and att coeff. en. in MeV => eff function for rez

    globals.spectra_collection[spec]['att cor'] = [round(x, 5) for x in att_cor]



# attenuation correction function for compounds
def attenuation_compound(spec, elem, thick, frac, dens):
    att_lib_param, att_lib_values = att_lib()
    
    att_cor = []
    elements = [x.strip() for x in elem.split('/')]
    fraction = [x.strip() for x in frac.split('/')]

    for num in fraction:
        if float(num) > 1:
            return tk.messagebox.showerror('User input error', 'Fraction by weight must be in interval <0,1>.')
    
    if sum([float(x) for x in fraction]) != 1:
        return tk.messagebox.showerror('User input error', 'Sum of fractions by weight must be equal to 1.')

    list_mu_elements = []
    
    # this part takes every element from the compound
    for x in range(0, len(elements)):
        # uploads from library for the sigle element
        a1, a2 = att_fit(att_lib_param[elements[x]], att_lib_values[elements[x]], elements[x])
        en_div_kev = att_lib_param[elements[x]][1] * 1000

        parameter = []

        # uploads parameter mu/ro for single element for given energies
        if 'cor_en' in globals.spectra_collection[spec]:
            for en in globals.spectra_collection[spec]['cor_en']:
                a = utilities.a_switch(en, en_div_kev, a1, a2)
                parameter.append(att_calc(a, en))
        else:
            for en in globals.spectra_collection[spec]['energy']:
                a = utilities.a_switch(en, en_div_kev, a1, a2)
                parameter.append(att_calc(a, en))
        
        # multiplies parameter mu/ro for every energy of single element by its fraction by weight
        mu_element = [i * float(fraction[x]) for i in parameter]
        # parameters w*mu/ro for every energy are added to the list for every element in compound
        list_mu_elements.append(mu_element)
    
    # sum of parameters of every element for one energy
    for a in range(0, len(parameter)):
        b = 0
        for i in range(0, len(list_mu_elements)):
            b += list_mu_elements[i][a]
        mu = b * dens
        att_cor.append((1 - np.exp(-1 * mu * thick / 10)) / (mu * thick / 10))  # en. in keV and att coeff. en. in MeV => eff function for rez

    globals.spectra_collection[spec]['att cor'] = [round(x, 5) for x in att_cor]

