# beam fluctuation correction
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
from datetime import datetime
import math
# Local modules
from modules import globals


# Function for beam fluctuation correction

# uploads data from deck file (time and intensity)
def upload_beam_info():
    with open(globals.lib_selection['Beam info'], 'r', encoding='utf-8') as temp_file:
        obsah = temp_file.readlines()

        beam_time = []
        beam_intens = []
        for i in range(0, len(obsah)):
            line = obsah[i].split()
            if not line[0]:
                continue
            if line[0] == '#':
                continue
            beam_time.append(datetime.strptime(line[0] + ' ' + line[1], '%d.%m.%Y %H:%M:%S'))
            beam_intens.append(float(line[2]))
        
        return beam_time, beam_intens


# calculates beam fluctuation correction coefficient
def beam_correction(gamma_ht):
    lamb = math.log(2) / gamma_ht
    beam_info = upload_beam_info()
    beam_times = beam_info[0]
    beam_intensities = beam_info[1]
    sum_intens = sum(beam_intensities)
    denominator = []
    for i in range(0, len(beam_times) - 1):
        duration = (beam_times[i+1] - beam_times[i]).total_seconds()
        time_to_end = (beam_times[-1] - beam_times[i+1]).total_seconds()
        weight_intens = beam_intensities[i+1] / sum_intens
        denominator.append((weight_intens / duration) * math.exp(-lamb * time_to_end) * (1 - math.exp(-lamb * duration)))
    irr_time = (beam_times[-1] - beam_times[0]).total_seconds()
    beam_corr_coef = (1 - math.exp(-lamb * irr_time)) / (irr_time * sum(denominator))
    
    return beam_corr_coef
