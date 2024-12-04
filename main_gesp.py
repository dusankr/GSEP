# Description: Main file for gamma package. It runs tkinter mainloop
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

# TODO_list:
# TODO: save bugs into log file

# Libraries
import os
import sys
import tkinter as tk
import ttkthemes
# Local libraries
from modules import corr_att
from modules import calc_coinc
from modules import calc_corr
from modules import calc_eff
from modules import calc_nonlin
from modules import calc_rr
from modules import inputs
from modules import globals


# Functions ----------------------------------------------------------------------------------------------------
# GUI exit from program
def ask_quit():
    if tk.messagebox.askokcancel('Quit', 'Do you want to quit now?'):
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent Fatal Python Error: PyEval_RestoreThread: NULL tstate


# selected spectrum from treeview
def selected_spectrum():
    try:
        spectrum = main_treeview.item(main_treeview.focus())['values'][0].split('.')[0]
        return spectrum
    except Exception as e:
        tk.messagebox.showerror('Error', 'No spectrum selected, please select one. Error: ' + str(e))
        return None

# --------------------------------------------------------------------------------------------------------------
# Work directory
os.chdir(os.path.dirname(sys.argv[0]))

# --------------------------------------------------------------------------------------------------------------

# MAIN CODE (part with all elements of GUI)
root = ttkthemes.ThemedTk(theme='breeze')
root.title('Gamma spectrometry processing package')
root.geometry('1000x300')
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

# close main window
root.protocol('WM_DELETE_WINDOW', lambda: ask_quit())

# tkinter variables
att_var = tk.StringVar()
eff_var = tk.StringVar()
nonlin_var = tk.StringVar()
beam_var = tk.StringVar()
gamma_var = tk.StringVar()
activity_var = tk.StringVar()
gamma_eff_var = tk.StringVar()

# Frames ------------------------------------------------------------------------------------------------------
main_frame = tk.ttk.Frame(root, relief='groove', borderwidth=2)
main_frame.grid(column=0, row=0, rowspan=2, sticky='nswe', padx=5, pady=5)
main_frame.columnconfigure(0, weight=1)
main_frame.rowconfigure(0, weight=1)

# main_FRAME
# treeview
main_treeview = tk.ttk.Treeview(main_frame)
main_treeview.grid(sticky='wens', column=0, row=0, rowspan=5)
main_treeview['show'] = 'headings'
# treeview on double click open window with spectrum
main_treeview.bind("<Double-1>", lambda event: inputs.open_spectrum(main_treeview.item(main_treeview.focus())['values'][0].split('.')[0], root))

# Treeview X-scrollbar
tree_x_scroll = tk.ttk.Scrollbar(root, orient='horizontal')
tree_x_scroll.grid(sticky='wens', column=0, row=5, columnspan=5, padx=5, pady=5)
tree_x_scroll.configure(command=main_treeview.xview)
main_treeview.configure(xscrollcommand=tree_x_scroll.set)

# Treeview Y-scrollbar
tree_y_scroll = tk.ttk.Scrollbar(root, orient='vertical', command=main_treeview.yview)
tree_y_scroll.grid(sticky='wens', column=5, row=0, rowspan=5, padx=5, pady=5)
main_treeview.configure(yscrollcommand=tree_y_scroll.set)

# --------------------------------------------------------------------------------------------------------------
# MENU
menubar = tk.Menu(root)

# menu for opening dictionary and close main window
main_menu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label="Main", menu=main_menu)
main_menu.add_command(label='Work Directory...', underline=0, command=lambda: inputs.open_folder(coinc_menu, bcg_menu, main_treeview, main_menu))
main_menu.add_command(label='Uploaded Information', underline=0, command=lambda: inputs.measur_info_win(root))
main_menu.add_command(label='Open selected spectrum', underline=0, command=lambda: inputs.open_spectrum(selected_spectrum(), root), state='disabled')
main_menu.add_command(label='Save all spectra', underline=0, command=lambda: inputs.save_all())
main_menu.add_separator()
main_menu.add_command(label='Quit', underline=0, command=lambda: ask_quit())

# menu for user input - files with spectra, measurement info and weights
user_input_menu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label="Input", menu=user_input_menu)
user_input_menu.add_command(label='Irradiation', underline=0, command=lambda: inputs.measurement_input_info_window(root, main_treeview))
user_input_menu.add_command(label='Sample', underline=0, command=lambda: corr_att.new_win_sample(root, main_treeview))
user_input_menu.add_command(label='Detector', underline=0, command=lambda: inputs.new_win_det(root, main_treeview))
user_input_menu.add_command(label='Non-point', underline=0, command=lambda: inputs.new_win_non_point(root, main_treeview))

# libraries menu
libraries_menu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label="Library", menu=libraries_menu)

att_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Attenuation', menu=att_menu, command=[
    calc_corr.original_library_selection('Attenuation', att_menu, att_var),
    globals.lib_selection.update({"Attenuation": att_var.get()})
    ])

bcg_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Background', menu=bcg_menu)

eff_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Efficiency', menu=eff_menu, command=[
    calc_corr.original_library_selection('Efficiency', eff_menu, eff_var),
    globals.lib_selection.update({"Efficiency": eff_var.get()})
    ])

nonlin_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Non-linearity', menu=nonlin_menu, command=[
    calc_corr.original_library_selection('Nonlinearity', nonlin_menu, nonlin_var),
    globals.lib_selection.update({"Nonlinearity": nonlin_var.get()})    
    ])

beam_corr_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Beam Fluctuation', menu=beam_corr_menu, command=[
    calc_corr.original_library_selection('Beam Fluctuation', beam_corr_menu, beam_var),
    globals.lib_selection.update({"Beam Fluctuation": beam_var.get()})    
    ])

coinc_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Coincidences', menu=coinc_menu)

libraries_menu.add_separator()

gamma_menu_all = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Gamma lib. (all)', menu=gamma_menu_all, command=[
    calc_corr.original_library_selection('Gamma lib. (all)', gamma_menu_all, gamma_var),
    globals.lib_selection.update({"Gamma lib. (all)": gamma_var.get()})  
    ])

libraries_menu.add_separator()

act_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label='Activity', menu=act_menu, command=[
    calc_corr.original_library_selection('Activity', act_menu, activity_var),
    globals.lib_selection.update({"Activity": activity_var.get()})
    ])

gamma_ref_menu = tk.Menu(libraries_menu, tearoff=0)
libraries_menu.add_cascade(label = 'Gamma lib. (eff)', menu=gamma_ref_menu, command=[
    calc_corr.original_library_selection('Gamma lib. (eff)', gamma_ref_menu, gamma_eff_var),
    globals.lib_selection.update({"Gamma lib. (eff)": gamma_eff_var.get()})
    ])

# calculation menu
calculation_menu = tk.Menu(menubar,tearoff = 0)
menubar.add_cascade(label="Calculation", menu=calculation_menu)
calculation_menu.add_command(label='Corrections', underline=0, command=lambda: calc_corr.corrections_win(root))
#calculation_menu.add_separator()
calculation_menu.add_command(label='Calculation', underline=0, command=lambda: calc_rr.result_win(root))

# efficiency menu
efficiency_menu = tk.Menu(menubar,tearoff = 0)
menubar.add_cascade(label="Detector", menu=efficiency_menu)
efficiency_menu.add_command(label='FEP efficiency', underline=0, command=lambda: calc_eff.efficiency_win(root))
efficiency_menu.add_command(label='Non-linearity', underline=0, command=lambda: calc_nonlin.nonlinearity_win(root))
efficiency_menu.add_command(label='Coincidence (Truecoinc)', underline=0, command=lambda: calc_coinc.coinc_win(root), state='disabled')

# activate Coincidence (Truecoinc) in menu if the system is Windows
if sys.platform == 'win32':    
    efficiency_menu.entryconfig('Coincidence (Truecoinc)', state='normal')

root.config(menu=menubar)

# run UI
root.mainloop()
