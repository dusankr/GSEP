# GSEP - Gamma Spectra Evaluation Package

# Introduction

This manual provides detailed guidance on the installation, usage, and simplified theory behind the Gamma Spectra Evaluation Package. The software was developed by the Nuclear Group at the Faculty of Electrical Engineering and Communication, Brno University of Technology (BUT).

A brief list of available tools and features:

- Upload processed (fitted peaks) gamma spectra in `.prn` format from the `Deimos32` [program](https://doi.org/10.1023/a:1025448800782).
- Enter experimental, sample, and detector information for further analysis.
- Upload library values for gamma corrections.
- Calculate and apply gamma spectrometry corrections.
- Calculate reaction rates, yields, and similar values from original or corrected data.
- Statistically process results and export them.
- Additional features for library preparation:
  - Gamma peak efficiency evaluation for a selected detector and its geometry.
  - Preparation of non-linearity correction libraries for a selected detector.
  - Preparation of true coincidence summation libraries via the `Truecoinc` program automation.

---

# Installation

The code is available in the [Gamma Spectra Evaluation Package](https://github.com/dusankr/gamma_spectra_processing_package) GitHub repository.

Users can install and run the program in two basic ways:

1. **Executable files**: Download compiled files from the latest release (usually available for Linux and Windows x64) and use them with a default Python ≥ 3.10 installation. All necessary packages are included in the compiled file.
2. **Source files**: Download the source files and install Python ≥ 3.10 along with the required libraries. A Python environment manager is highly recommended. Based on experience, Conda can create overly large environments without proper settings, making simpler tools like `virtualenv`, `venv`, or other managers preferable.

## Executable files

This option offers the simplest way to use the Gamma Spectra Evaluation Package. Users need to download and install Python 3.10 or newer. The software is currently developed using Python version 3.12.

## Source files

Users can run the source code if the packages listed below and their dependencies are installed:

- `Python` ≥ 3.10
- `matplotlib`
- `numpy`
- `openpyxl`
- `pandas`
- `pathlib`
- `pywin32`
- `pywinauto`
- `scipy`
- `tkinter`
- `ttkthemes`
- `ttkwidgets`
- `XlsxWriter`

The program can be started by typing the following terminal command:

```
> python main_gamma_package.py
```

```
Program folder
├── main_gsep.py
├── README.md
├── libraries
└── modules
    ├── calc_coinc.py
    ├── calc_corr.py
    ├── calc_eff.py
    ├── calc_nonlin.py
    ├── calc_rr.py
    ├── corr_att.py
    ├── corr_bcg.py
    ├── corr_beam.py
    ├── corr_coinc.py
    ├── corr_eff.py
    ├── corr_nonlin.py
    ├── corr_sepdep.py
    ├── globals.py
    ├── inputs.py
    ├── settings.py
    └── utilities.py
```