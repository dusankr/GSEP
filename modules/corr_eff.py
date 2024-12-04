# efficiency correction calculation
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
from pathlib import Path
from tkinter import messagebox
from math import exp, log
# Local modules
from modules import globals
from modules import utilities


# calculate efficiency value on dependency of det, respectively used fitting function and its parameters
def eff_calc(a, en, unit):
    eff = 0
    if unit == 'keV':           # if fit is made with keV
        for i in range(0, len(a)):
            eff += a[i] * log(en) ** i
        return exp(eff)
    elif unit == 'MeV':           # if fit is made with MeV
        for i in range(0, len(a)):
            eff += a[i] * log(en / 1000) ** i
        return exp(eff)
    elif unit == 'dubna':         # if fit was made with dubna style
        for i in range(0, len(a)):
            eff += a[i] * log(en / 1000) ** i
        return exp(eff) / 1e6


# read efficiency parameters from library
def eff_lib(det, pos):
    with open(globals.lib_selection['Efficiency'], 'r', encoding='utf-8') as temporary_file:
        content = temporary_file.readlines()
    
    for i in range(0, len(content)):
        line = content[i].strip().split()
        if len(line) > 0 and line[0][0] != '#':
            if det == line[0] and pos == line[1] and int(line[2]) == 0:
                if line[4] == 'keV':
                    unit = 'keV'
                elif line[4] == 'MeV':
                    unit = 'MeV'
                elif (line[4]).lower() == 'dubna':
                    unit = 'dubna'
                else:
                    unit = '-'
                
                # efficeincy curve parameters
                line = content[i+1].split()
                a1 = [float(a) for a in line]
                
                return 0, a1, 0, unit
            
            elif det == line[0] and pos == line[1] and int(line[2]) > 0:
                if line[4] == 'keV':
                    unit = 'keV'
                elif line[4] == 'MeV':
                    unit = 'MeV'
                elif line[4] == 'dubna':
                    unit = 'dubna'
                else:
                    unit = '-'
                en_div = int(line[2])

                line = content[i+1].split()
                a1 = [float(a) for a in line]                     # transform a list of strings to a list of floats
                line = content[i+2].split()
                a2 = [float(a) for a in line]
                
                return en_div, a1, a2, unit


# main function for efficiency correction
def efficiency():
    for spectrum in globals.spectra_collection:
        spec = Path(spectrum + '.prn')
        det = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Detector'].item()
        pos = globals.spectra_param_collection.loc[globals.spectra_param_collection['Name'] == spec, 'Position'].item()

        # try to get efficiency curve parameters from library file      TODO read library into dictionary instead of reading it every time
        try:
            en_div, a1, a2, unit = eff_lib(det, pos)
        except Exception as e:
            return messagebox.showerror('User input error', 'Efficiency curve parameters for detector {} and position {} are not in the library file. Error: {}'.format(det, pos, e))
        
        # chose the right energy column
        if 'cor_en' in globals.spectra_collection[spectrum]:
            energy = globals.spectra_collection[spectrum]['cor_en']
        else:
            energy = globals.spectra_collection[spectrum]['energy']

        # start correction
        eff_fep = []

        for en in energy:
            a = utilities.a_switch(en, en_div, a1, a2)
            eff_fep.append(eff_calc(a, en, unit))

        globals.spectra_collection[spectrum]['eff fep'] = [round(x, 5) for x in eff_fep]

    print('Efficiency correction done.')