# Routine to store methods to plot

import matplotlib.pyplot as plt

########################################################################
#                          Bidimensional plots                         #
########################################################################

# Defines a function to plot curves in the XY plane

def plane_plot(file_name, data=None, x_data=None, y_data=None, label=
None, x_label=None, y_label=None, title=None, flag_grid=True, 
highlight_points=False, color='orange'):

    print("Starts plotting")

    # Checks if all variables data, x_data, and y_data are None

    if (data is None) and (x_data is None) and (y_data is None):

        raise ValueError("The list of lists of data has not been provi"+
        "ded nor the x_data list neither the y_data. Thus, no data can"+
        " be plotted")

    # If the data was given

    elif not (data is None):

        if not isinstance(data, list):

            raise TypeError("The data list is not a list. No data can "+
            "be retrieved for plotting")

        elif len(data)==0:

            raise IndexError("The data list is empty. No data can be r"+
            "etrieved for plotting")

        # Iterates through the data point

        x_data = []

        y_data = []

        for point in data:

            # Verifies if this data point is a list

            if not isinstance(point, list):

                raise TypeError("Each data point in the data list must"+
                " be a list, for plotting them")
            
            elif len(point)<2:

                raise IndexError("The data point has length less than "+
                "2, thus, a pair of data cannot be retrieved for plott"+
                "ing")

            x_data.append(point[0])

            y_data.append(point[1])

    # Checks the x and y lists

    else:

        if not isinstance(x_data, list):

            raise TypeError("The x data is not a list, thus, cannot be"+
            " used for plotting")

        elif not isinstance(y_data, list):

            raise TypeError("The y data is not a list, thus, cannot be"+
            " used for plotting")

        elif len(x_data)!=len(y_data):

            raise IndexError("The x data and y data are lists of diffe"+
            "rent sizes. Thus, cannot be used for plotting")

    # Sets the graph to be plotted in LaTeX style

    plt.rcParams.update({"text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], "axes.labelsize": 14, "fo"+
    "nt.size": 14, "legend.fontsize": 12, "xtick.labelsize": 12, "ytic"+
    "k.labelsize": 12})

    # Plots it

    if label is None:

        plt.plot(x_data, y_data, color=color)

    else:

        plt.plot(x_data, y_data, label=label, color=color)

    plt.grid(flag_grid)

    # Plots the set of points as scattered markers

    if highlight_points:

        plt.scatter(x_data, y_data, color='black', marker='x', zorder=3)

    if not (x_label is None):

        plt.xlabel(x_label)

    if not (y_label is None):

        plt.ylabel(y_label)

    if not (title is None):

        plt.title(title)

    if not (label is None):

        plt.legend()

    # Saves the plot

    plt.savefig(file_name)

    print("Finishes plotting")