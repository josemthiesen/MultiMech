# Routine to plot the tangent operators

import os

import numpy as np

import source.tool_box.file_handling_tools as file_tools

# Defines a function to run the process

def plot_operators():

    # Sets the multiscale boundary conditions for each one of the fields

    multiscale_BCsSets = [["MinimallyConstrainedFirstOrderBC", "Minima"+
    "llyConstrainedFirstOrderBC"], ["PeriodicFirstOrderBC", "PeriodicF"+
    "irstOrderBC"]]#, ["LinearFirstOrderBC", "LinearFirstOrderBC"]]

    # Sets the basic path

    base_path = (os.getcwd()+"//tests//micropolar//tangent_operators//"+
    "results")

    # Sets a list of names for each set of parameters, which will yield
    # different simulations

    simulations_names = ["simulation_11", "simulation_12", "simulation"+
    "_13", "simulation_21", "simulation_22", "simulation_23", "simulat"+
    "ion_31", "simulation_32", "simulation_33"]

    # Iterates through the multiscale boundary conditions

    for multiscale_BCs in multiscale_BCsSets:

        # Iterates through the simulations

        for i in range(len(simulations_names)):

            # Gets the folder name

            folder_name = (base_path+"//"+simulations_names[i]+"//"+
            multiscale_BCs[0]+"_"+multiscale_BCs[1])

            # Calls the function to read the derivatives of the stress
            # tensors and plot the dispersion of the components

            plot_dispersion(folder_name)

# Defines a function to read the derivatives of the stress tensors, and 
# plot the dispersion of the components

def plot_dispersion(folder_name):

    # Reads the derivative files

    dP_dGradU = file_tools.txt_toList("dP_dGradU", parent_path=
    folder_name)

    dP_dGradPhi = file_tools.txt_toList("dP_dGradPhi", parent_path=
    folder_name)

    dPcouple_dGradU = file_tools.txt_toList("dPcouple_dGradU",
    parent_path=folder_name)

    dPcouple_dGradPhi = file_tools.txt_toList("dPcouple_dGradPhi", 
    parent_path=folder_name)

    # Gets the maximum value of magnitude

    maximum_magnitude = get_maximumMagnitude(dP_dGradU)

    maximum_magnitude = get_maximumMagnitude(dP_dGradPhi, 
    maximum_magnitude=maximum_magnitude)

    maximum_magnitude = get_maximumMagnitude(dPcouple_dGradU, 
    maximum_magnitude=maximum_magnitude)

    maximum_magnitude = get_maximumMagnitude(dPcouple_dGradPhi, 
    maximum_magnitude=maximum_magnitude)

    # Normalizes the components of the tensors

    dP_dGradU = normalize_tensors(dP_dGradU, maximum_magnitude)

    dPcouple_dGradU = normalize_tensors(dPcouple_dGradU, 
    maximum_magnitude)

    dP_dGradPhi = normalize_tensors(dP_dGradPhi, maximum_magnitude)

    dPcouple_dGradPhi = normalize_tensors(dPcouple_dGradPhi, 
    maximum_magnitude)

# Defines a function to get the maximum value of magnitude of a list of
# tensors

def get_maximumMagnitude(tensors_list, maximum_magnitude=1.0):

    # Iterates through the time values

    for time in tensors_list:

        # iterates through the components of the matrix

        for i in range(len(time[1])):

            for j in range(len(time[1][i])):

                maximum_magnitude = max(maximum_magnitude, abs(time[1][i
                ][j]))

    return maximum_magnitude

# Defines a function to normalize the values of the tensors by the maxi-
# mum magnitude found

def normalize_tensors(tensors_list, maximum_magnitude):

    # Iterates through the time values

    for time in tensors_list:

        # iterates through the components of the matrix

        for i in range(len(time[1])):

            for j in range(len(time[1][i])):

                time[1][i][j] = time[1][i][j]/maximum_magnitude

    return tensors_list

# Deines a function to plot the tensor

def plot_tensor(dP_dGradU, dP_dGradPhi, dPcouple_dGradU, 
dPcouple_dGradPhi):

    pass

plot_operators()