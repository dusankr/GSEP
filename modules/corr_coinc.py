# Coincidence correction calculation
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
from modules import globals
from modules import settings
import numpy as np


# read background peaks from library and send it to main coinc function
def coinc_lib(det):
    with open(globals.coinc_lib_selection[det], 'r', encoding='utf-8') as temporary_file:
        obsah = temporary_file.readlines()

        lib_values = {}
        elements = []
        # loads elements
        for i in range(2, len(obsah)):
            line = obsah[i].split()
            elements.append(line[0])

        elements = tuple(set(elements))

        # loads energies
        for key in elements:
            en = []
            corr = []
            # loads positions
            for i in range(2, len(obsah)):
                line = obsah[i].split()
                headers = obsah[1].split()
                if line[0] == key:
                    en.append(float(line[1]))
                    pos = []
                    coinc = []
                    # loads pairs of position and coincidence correction for given nuclide, energy and position
                    # for every energy there is list of pairs (position, correction)
                    for j in range(0, len(headers)):
                        if headers[j].startswith('P'):
                            pos.append(headers[j].lower())
                            coinc.append(float(line[j-1]))
                            pair = zip(pos, coinc)
                    corr.append((list(pair)))
            lib_values[key] = [en, corr]
    return lib_values


# main correction function
def coinc_correction(nuclide, en, coinc_lib, position):
    coinc_corr = 1

    if nuclide in coinc_lib:
        coinc_corr_nuclide = coinc_lib[nuclide]
        for i, en_lib in enumerate(coinc_corr_nuclide[0]):
            if abs(en - en_lib) < (settings.N_CALC * np.sqrt(settings.A_CALC**2 + (0.01)**2 + (0.005)**2) + 0.5):
                break
    
        # get correction for given energy and position
        for j, pos in enumerate(coinc_corr_nuclide[1][i]):
            if pos[0] == position:
                coinc_corr = pos[1]
                break
    
    return coinc_corr
