# Routine to plot the tangent operators

import os

import numpy as np

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.plotting_tools as plotting_tools

# Defines a function to run the process

def plot_operators():

    # Creates a dictionary to convert the current indices to the origi-
    # nal indices of the tensor

    voigt_conversion = {0: 0, 1: 4, 2: 8, 3: 1, 4: 5, 5: 2, 6: 3, 7: 7,
    8: 6}

    #voigt_conversion = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7,
    #8: 8}

    # Sets the basic size for the component representative circle

    basic_size = 12.0

    # Sets the multiscale boundary conditions for each one of the fields

    multiscale_BCsSets = [["MinimallyConstrainedFirstOrderBC", "Minima"+
    "llyConstrainedFirstOrderBC"], ["PeriodicFirstOrderBC", "PeriodicF"+
    "irstOrderBC"], ["LinearFirstOrderBC", "LinearFirstOrderBC"]]

    # Sets the basic path

    base_paths = [(os.getcwd()+"//tests//micropolar//tangent_operators"+
    "//results_eps_1E_6"), (os.getcwd()+"//tests//micropolar//tangent_"+
    "operators//results_eps_1E_5")]

    # Sets a list of names for each set of parameters, which will yield
    # different simulations

    simulations_names = ["simulation_11", "simulation_12", "simulation"+
    "_13", "simulation_21", "simulation_22", "simulation_23", "simulat"+
    "ion_31", "simulation_32", "simulation_33"]

    # Iterates through the different perturbation step sizes

    for base_path in base_paths:

        # Iterates through the multiscale boundary conditions

        for multiscale_BCs in multiscale_BCsSets:

            # Iterates through the simulations

            for i in range(len(simulations_names)):

                # Gets the folder name

                folder_name = (base_path+"//"+simulations_names[i]+"//"+
                multiscale_BCs[0]+"_"+multiscale_BCs[1])

                # Calls the function to read the derivatives of the 
                # stress tensors and plot the dispersion of the compo-
                # nents

                plot_dispersion(folder_name, basic_size, 
                voigt_conversion)  

# Defines a function to read the derivatives of the stress tensors, and 
# plot the dispersion of the components

def plot_dispersion(folder_name, basic_size, voigt_conversion):

    # Reads the derivative files

    dP_dGradU = file_tools.txt_toList("dP_dGradU", parent_path=
    folder_name)

    dP_dGradPhi = file_tools.txt_toList("dP_dGradPhi", parent_path=
    folder_name)

    dPcouple_dGradU = file_tools.txt_toList("dPcouple_dGradU",
    parent_path=folder_name)

    dPcouple_dGradPhi = file_tools.txt_toList("dPcouple_dGradPhi", 
    parent_path=folder_name)

    # Gets the maximum and minimum values of components

    min_component, max_component = get_maximumMagnitude(dP_dGradU)

    min_component, max_component = get_maximumMagnitude(dP_dGradPhi, 
    max_component=max_component, min_component=min_component)

    min_component, max_component = get_maximumMagnitude(dPcouple_dGradU, 
    max_component=max_component, min_component=min_component)

    min_component, max_component = get_maximumMagnitude(
    dPcouple_dGradPhi, max_component=max_component, min_component=
    min_component)

    # Iterates through the time points

    for i in range(len(dP_dGradU)):

        # Gets the time point and the file name

        time = dP_dGradU[i][0]

        file_name = ("C_voigt_t_"+str(int(np.floor(time)))+"_"+str(time-
        np.floor(time))[2:]+".pdf")

        # Plots the tensors

        plot_tensor(dP_dGradU[i][1], dP_dGradPhi[i][1], dPcouple_dGradU[
        i][1], dPcouple_dGradPhi[i][1], folder_name, file_name, time, 
        basic_size, voigt_conversion, min_component, max_component)

# Defines a function to get the maximum value of magnitude of a list of
# tensors

def get_maximumMagnitude(tensors_list, max_component=0.0, min_component=
0.0):

    # Iterates through the time values

    for time in tensors_list:

        # iterates through the components of the matrix

        for i in range(len(time[1])):

            for j in range(len(time[1][i])):

                max_component = max(max_component, time[1][i][j])

                min_component = min(min_component, time[1][i][j])

    return min_component, max_component

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
dPcouple_dGradPhi, parent_path, file_name, time, basic_size,
voigt_conversion, min_component, max_component):

    # Initializes the color list and the marker size list

    color_list = []

    marker_sizeList = []

    # Initializes the plot lists

    x_data = []

    y_data = []

    # Defines a function to scale the components of the fourth order 
    # tensor

    def scaling(x):

        return np.sign(x)*np.log10(np.abs(x)+1.0)

    # Creates a function for the conversion of indices

    def indices_function(i,j):

        x_data.append(j+1)

        y_data.append(18-i)

        marker_sizeList.append(basic_size)

        # Verifies if the indices are in the derivative of the first Pi-
        # ola-Kirchhoff stress tensor w.r.t. the displacement gradient

        if i<9 and j<9:

            color_list.append(scaling(dP_dGradU[voigt_conversion[i]][
            voigt_conversion[j]]))

        # Verifies if the indices are in the derivative of the first Pi-
        # ola-Kirchhoff stress tensor w.r.t. the microrotation gradient

        if i<9 and j>8:

            color_list.append(scaling(dP_dGradPhi[voigt_conversion[i]][
            voigt_conversion[j-9]]))

        # Verifies if the indices are in the derivative of the first Pi-
        # ola-Kirchhoff couple stress tensor w.r.t. the displacement 
        # gradient

        if i>8 and j<9:

            color_list.append(scaling(dPcouple_dGradU[voigt_conversion[i
            -9]][voigt_conversion[j]]))

        # Verifies if the indices are in the derivative of the first Pi-
        # ola-Kirchhoff couple stress tensor w.r.t. the microrotation 
        # gradient

        if i>8 and j>8:

            color_list.append(scaling(dPcouple_dGradPhi[voigt_conversion[
            i-9]][voigt_conversion[j-9]]))

    for i in range(18):

        for j in range(18):

            indices_function(i,j)

    # Gets the extreme values for the color bar

    max_magnitude = max(abs(max_component), abs(min_component))

    # Gets the limit values for the ticks of the color bar

    ticks_min = scaling(min_component)

    ticks_max = scaling(max_component)

    # Gets the initial and final value of the integer ticks

    initial_tick = int(np.ceil(ticks_min)+1)

    final_tick = int(np.floor(ticks_max))

    color_barTicks = [ticks_min]

    for tick in range(initial_tick, final_tick):

        color_barTicks.append(tick)

    color_barTicks.append(ticks_max)

    # Plots and saves the figure

    plotting_tools.plane_plot(parent_path+"//"+file_name, x_data=x_data, 
    y_data=y_data, element_style="s", element_size=marker_sizeList, 
    color=color_list, color_map='coolwarm', plot_type="scatter", 
    flag_grid=True, flag_noTicks=True, aspect_ratio='equal', x_grid=[
    3.5, 6.5, 9.5, 12.5, 15.5], y_grid=[3.5, 6.5, 9.5, 12.5, 15.5],
    color_bar=True, color_barMaximum=scaling(max_magnitude), 
    color_barMinimum=scaling(max_magnitude*np.sign(min_component)), 
    color_barTitle="$sgn\\left(\\|\\cdot\\|\\right)log\\left(\\|\\cdot"+
    "\\|+1\\right)$", color_barTicks=color_barTicks)

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