"""
This module contains scripts for analyzing the data for the FRET fitting project
updated: Dec 2014
by Justin Chen
"""

#important to set backend for running on the server
import matplotlib
matplotlib.use("Agg")

#import plot_epsilon_native
#import FRET_experiment


import plot_depsilon_native
import pair_distance_calculator
import plot_function
import JacRunModule
import quick_jac
import free_energy_plot_1d
import free_energy_plot_2d
import merge_files
import find_frames
import recipe_log_function
import gro_reader
import find_y
import histy_dmdmd
import plot_epsilons
import find_last_crossing
import make_index
import boot_strap_method
import plot_package
import tica_weighted_sampling

##temporary likely
import convert_to_ini
import free_energy_plot_titers
