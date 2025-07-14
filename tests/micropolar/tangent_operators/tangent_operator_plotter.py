# Routine to plot the tangent operators

import os

import numpy as np

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.plotting_tools as plotting_tools

# Defines a function to run the process

def plot_operators():

    # Sets the basic size for the component representative circle

    basic_size = 12.0

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

            plot_dispersion(folder_name, basic_size)

# Defines a function to read the derivatives of the stress tensors, and 
# plot the dispersion of the components

def plot_dispersion(folder_name, basic_size):

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

    # Iterates through the time points

    for i in range(len(dP_dGradU)):

        # Gets the time point and the file name

        time = dP_dGradU[i][0]

        file_name = ("C_voigt_t_"+str(int(np.floor(time)))+"_"+str(time-
        np.floor(time))[2:]+".pdf")

        # Plots the tensors

        plot_tensor(dP_dGradU[i][1], dP_dGradPhi[i][1], dPcouple_dGradU[
        i][1], dPcouple_dGradPhi[i][1], folder_name, file_name, time, 
        basic_size)

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
dPcouple_dGradPhi, parent_path, file_name, time, basic_size):

    # Initializes the color list and the marker size list

    color_list = []

    marker_sizeList = []

    # Initializes the plot lists

    x_data = []

    y_data = []

    # Initializes the counters for the Voigt notation indices

    row_index = 9

    column_index = 1

    # Iterates through the derivative of the first Piola-Kirchhof couple
    # stress tensor w.r.t. the displacement gradient

    for i in range(len(dPcouple_dGradU)):

        for j in range(len(dPcouple_dGradU[i])):

            x_data.append(column_index)

            y_data.append(row_index)

            color_list.append(0.5*(dPcouple_dGradU[i][j]+1))

            #marker_sizeList.append(abs(dPcouple_dGradU[i][j])*basic_size)

            marker_sizeList.append(basic_size)

            # Updates the column counter

            column_index += 1

        # Updates the counter

        column_index = 1

        row_index -= 1

    # Iterates through the derivative of the first Piola-Kirchhof couple
    # stress tensor w.r.t. the microrotation gradient

    row_index = 9

    column_index = 10

    for i in range(len(dPcouple_dGradPhi)):

        for j in range(len(dPcouple_dGradPhi[i])):

            x_data.append(column_index)

            y_data.append(row_index)

            color_list.append(0.5*(dPcouple_dGradPhi[i][j]+1))

            #marker_sizeList.append(abs(dPcouple_dGradPhi[i][j])*
            #basic_size)

            marker_sizeList.append(basic_size)

            # Updates the column counter

            column_index += 1

        # Updates the counter

        column_index = 10

        row_index -= 1

    # Iterates through the derivative of the first Piola-Kirchhof stress 
    # tensor w.r.t. the displacement gradient

    row_index = 18

    column_index = 1

    for i in range(len(dP_dGradU)):

        for j in range(len(dP_dGradU[i])):

            x_data.append(column_index)

            y_data.append(row_index)

            color_list.append(0.5*(dP_dGradU[i][j]+1))

            #marker_sizeList.append(abs(dP_dGradU[i][j])*basic_size)

            marker_sizeList.append(basic_size)

            # Updates the column counter

            column_index += 1

        # Updates the counter

        column_index = 1

        row_index -= 1

    # Iterates through the derivative of the first Piola-Kirchhof stress 
    # tensor w.r.t. the microrotation gradient

    row_index = 18

    column_index = 10

    for i in range(len(dP_dGradPhi)):

        for j in range(len(dP_dGradPhi[i])):

            x_data.append(column_index)

            y_data.append(row_index)

            color_list.append(0.5*(dP_dGradPhi[i][j]+1))

            #marker_sizeList.append(abs(dP_dGradPhi[i][j])*basic_size)

            marker_sizeList.append(basic_size)

            # Updates the column counter

            column_index += 1

        # Updates the counter

        column_index = 10

        row_index -= 1

    # Plots and saves the figure

    plotting_tools.plane_plot(parent_path+"//"+file_name, x_data=x_data, 
    y_data=y_data, element_style="s", element_size=marker_sizeList, 
    color=color_list, color_map='coolwarm', plot_type="scatter", 
    flag_grid=False, flag_noTicks=True)

def test():

    base_path = (os.getcwd()+"//tests//micropolar//tangent_operators//"+
    "results//simulation_11//LinearFirstOrderBC_LinearFirstOrderBC")

    x_data = [1,2,3]

    y_data = [1,2,3]

    marker_sizeList = [1,2,3]

    color_list = [0.0, 0.25, 1.0]

    plotting_tools.plane_plot(base_path+"//"+"t_0", x_data=x_data, 
    y_data=y_data, element_style="x", element_size=marker_sizeList, 
    color=color_list, color_map='coolwarm', plot_type="scatter")

#test()

plot_operators()