# Routine to store methods to plot

import numpy as np

import matplotlib.pyplot as plt

import matplotlib.ticker as ticker

########################################################################
#                          Bidimensional plots                         #
########################################################################

# Defines a function to plot curves in the XY plane

def plane_plot(file_name, data=None, x_data=None, y_data=None, label=
None, x_label=None, y_label=None, title=None, flag_grid=True, 
highlight_points=False, color='orange', flag_scientificNotation=False,
element_style='-', element_size=1.5,  legend=None, plot_type="line", 
color_map=False, flag_noTicks=False, aspect_ratio='auto', x_grid=None,
y_grid=None, color_bar=False, color_barMaximum=None, color_barMinimum=
None, color_barTicks=None, color_barTitle=None):

    print("Starts plotting")

    # Initializes a flag to inform if mutliple curves are supplied

    multiple_curves = False

    # Checks if all variables data, x_data, and y_data are None

    if (data is None) and (x_data is None) and (y_data is None):

        raise ValueError("The list of lists of data has not been provi"+
        "ded nor the x_data list neither the y_data. Thus, no data can"+
        " be plotted")

    # If the data was given

    elif not (data is None):

        if not isinstance(data, list):

            # Verifies if does not have the to list method

            if not hasattr(data, "tolist"):

                raise TypeError("The data list is not a list nor a num"+
                "py array neither a tensorflow tensor. No data can be "+
                "retrieved for plotting")
            
            else:

                data = data.tolist()

        elif len(data)==0:

            raise IndexError("The data list is empty. No data can be r"+
            "etrieved for plotting")

        # Iterates through the data point

        x_data = []

        y_data = []

        # Verifies if each point of data has more than two values. If so,
        # this means more curves are provided

        if len(data[0])>2:

            multiple_curves = len(data[0])-1

            for i in range(len(data[0])-1):

                y_data.append([])

        # If the plot type is catter, treats each point as a separate 
        # curve

        elif plot_type=="scatter":

            multiple_curves = len(data)

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

            if multiple_curves:

                for i in range(len(point)-1):

                    y_data[i].append(point[i+1])

            else:

                y_data.append(point[1])

    # Checks the x and y lists

    else:

        if not isinstance(x_data, list):

            # Verifies if does not have the to list method

            if not hasattr(x_data, "tolist"):

                raise TypeError("The x data is not a list nor a numpy "+
                "array neither a tensorflow tensor, thus, cannot be us"+
                "ed for plotting")
            
            else:

                x_data = x_data.tolist()

        elif not isinstance(y_data, list):

            # Verifies if does not have the to list method

            if not hasattr(y_data, "tolist"):

                raise TypeError("The y data is not a list nor a numpy "+
                "array neither a tensorflow tensor, thus, cannot be us"+
                "ed for plotting")
            
            else:

                y_data = y_data.tolist()

        elif isinstance(y_data[0], list):

            multiple_curves = len(y_data)

            if len(x_data)!=len(y_data[0]):

                raise IndexError("The x data and y data are lists of d"+
                "ifferent sizes. Thus, cannot be used for plotting")

        elif len(x_data)!=len(y_data):

            raise IndexError("The x data and y data are lists of diffe"+
            "rent sizes. Thus, cannot be used for plotting")
        
        # If the plot type is catter, treats each point as a separate 
        # curve

        elif plot_type=="scatter":

            multiple_curves = len(x_data)

    # Sets the graph to be plotted in LaTeX style

    plt.rcParams.update({"text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], "axes.labelsize": 14, "fo"+
    "nt.size": 14, "legend.fontsize": 12, "xtick.labelsize": 12, "ytic"+
    "k.labelsize": 12, "text.latex.preamble": r"\usepackage{amsmath}"})

    # Gets the color map

    if color_map:

        if isinstance(color_map, str):

            try:

                color_map = plt.get_cmap(color_map)

            except Exception as error_message:

                print("Error Message:"+str(error_message)+"\nProbably "+
                "this color map does not exist")

    # Creates the figure and the subplots

    figure, subplots_tuple = plt.subplots()

    # Verifies the nature of the line styles and of the color vector

    if multiple_curves:

        # Verifies line styles

        if isinstance(element_style, str):

            element_style = [element_style for i in range(
            multiple_curves)]

        elif not isinstance(element_style, list):

            raise TypeError("Multiple curves were given to be plotted,"+
            " but the element_style variables is neither a string nor "+
            "a list")
        
        elif len(element_style)!=multiple_curves:

            raise IndexError(str(multiple_curves)+" curves were given "+
            "to be plotted but "+str(len(element_style))+" line styles"+
            " were given")

        # Verifies line thickness

        if isinstance(element_size, float) or isinstance(element_size, 
        int):

            element_size = [element_size for i in range(
            multiple_curves)]

        elif not isinstance(element_size, list):

            raise TypeError("Multiple curves were given to be plotted,"+
            " but the element_size variables is neither a float or an "+
            "integer nor a list")
        
        elif len(element_size)!=multiple_curves:

            raise IndexError(str(multiple_curves)+" curves were given "+
            "to be plotted but "+str(len(element_size))+" line thickne"+
            "sses were given")
        
        # Verifies the color vector

        if isinstance(color, str):

            color = [color for i in range(multiple_curves)]

        elif isinstance(color, float) or isinstance(color, int):

            if color_map:

                color = [color_map(color) for i in range(multiple_curves
                )]

            else:

                color = [color for i in range(multiple_curves)]

        elif not isinstance(color, list):

            raise TypeError("Multiple curves were given to be plotted,"+
            " but the color variables is neither a string nor a list")
        
        elif len(color)!=multiple_curves:

            raise IndexError(str(multiple_curves)+" curves were given "+
            "to be plotted but "+str(len(color))+" colors were given")
        
        elif color_map:

            # Gets the extrema values of the colors

            color_min = None

            color_max = None

            try:

                color_min = min(color)

                color_max = max(color)

                # If the minimum value is not given

                if color_barMinimum is None:

                    color_barMinimum = color_min*1.0

                elif color_min<color_barMinimum:

                    color_barMinimum = color_min*1.0

                # If the maximum is not given

                if color_barMaximum is None:

                    color_barMaximum = color_max*1.0

                elif color_max>color_barMaximum:

                    color_barMaximum = color_max*1.0

            except:

                pass

            # Iterates through the color values

            for i in range(len(color)):

                if isinstance(color[i], float) or isinstance(color[i], 
                int):
                    
                    color[i] = color_map((color[i]-color_barMinimum)/(
                    color_barMaximum-color_barMinimum))

            # Updates the color map variable to account for the maximum
            # and minimum values

            if (not (color_min is None)) and (not (color_max is None)):

                color_map = [color_map, color_barMinimum, 
                color_barMaximum]
            
    else:

        # Verifies the line styles
        
        if isinstance(element_style, list):

            if len(element_style)>1:

                raise IndexError(str(len(element_style))+" line styles"+
                " were given, but there is only one curve to be plotte"+
                "d")
            
            else:

                element_style = element_style[0]

        # Verifies the line thicknesses
        
        if isinstance(element_size, list):

            if len(element_size)>1:

                raise IndexError(str(len(element_size))+" line thickne"+
                "sses were given, but there is only one curve to be pl"+
                "otted")
            
            else:

                element_size = element_size[0]

        # And verifies the line colors
        
        if isinstance(color, list):

            if len(color)>1:

                raise IndexError(str(len(color))+" colors were given, "+
                "but there is only one curve to be plotted")
            
            else:

                if color_map and (isinstance(color[0], float) or (
                isinstance(color[0], int))):

                    color = color_map(color[0])

                else:

                    color = color[0]

    # Sets the aspect ratio of the plot

    subplots_tuple.set_aspect(aspect_ratio)

    # Inititalizes the plotted entities

    plotted_entities = None

    # Plots it

    if label is None:

        if multiple_curves:

            for i in range(multiple_curves):

                if plot_type=="line":

                    plotted_entities = subplots_tuple.plot(x_data, 
                    y_data[i], linestyle=element_style[i], linewidth=
                    element_size[i], color=color[i])

                elif plot_type=="scatter":

                    # If multiple curves are plotted at once

                    if isinstance(y_data[i], list):

                        plotted_entities = subplots_tuple.scatter(x_data, 
                        y_data[i], marker=element_style[i], s=
                        element_size[i]**2, color=color[i], zorder=3)

                    # If it is just the default treatment of scatter 
                    # plots

                    else:

                        plotted_entities = subplots_tuple.scatter(x_data[
                        i], y_data[i], marker=element_style[i], s=
                        element_size[i]**2, color=color[i], zorder=3)

                else:

                    raise ValueError("There are two types of plot: 'li"+
                    "ne' for data points that are orderly joined by li"+
                    "nes, and 'scatter' for points that are plotted sc"+
                    "atter fashion")

        else:

            if plot_type=="line":

                plotted_entities = subplots_tuple.plot(x_data, y_data, 
                linestyle=element_style, linewidth=element_size, color=
                color)

            elif plot_type=="scatter":

                plotted_entities = subplots_tuple.scatter(x_data, y_data, 
                marker=element_style, s=element_size**2, color=color, 
                zorder=3)

            else:

                raise ValueError("There are two types of plot: 'line' "+
                "for data points that are orderly joined by lines, and"+
                " 'scatter' for points that are plotted scatter fashio"+
                "n")

    else:

        if len(label)!=multiple_curves:

            raise IndexError("The length of the label list, "+str(len(
            label))+", is not equal to the number of curves, "+str(
            multiple_curves))

        if multiple_curves:

            for i in range(multiple_curves):
                
                if plot_type=="line":

                    plotted_entities = subplots_tuple.plot(x_data, 
                    y_data[i], linestyle=element_style[i], linewidth=
                    element_size[i], color=color[i], label=label[i])

                elif plot_type=="scatter":

                    # If multiple curves are plotted at once

                    if isinstance(y_data[i], list):

                        plotted_entities = subplots_tuple.scatter(x_data, 
                        y_data[i], marker=element_style[i], s=
                        element_size[i]**2, color=color[i], zorder=3)

                    # If it is just the default treatment of scatter 
                    # plots

                    else:

                        plotted_entities = subplots_tuple.scatter(x_data[
                        i], y_data[i], marker=element_style[i], s=
                        element_size[i]**2, color=color[i], zorder=3)

                else:

                    raise ValueError("There are two types of plot: 'li"+
                    "ne' for data points that are orderly joined by li"+
                    "nes, and 'scatter' for points that are plotted sc"+
                    "atter fashion")

        else:

            if plot_type=="line":

                plotted_entities = subplots_tuple.plot(x_data, y_data, 
                linestyle=element_style, linewidth=element_size, color=
                color, label=label)

            elif plot_type=="scatter":

                plotted_entities = subplots_tuple.scatter(x_data, y_data, 
                marker=element_style, s=element_size**2, color=color, 
                label=label, zorder=3)

            else:

                raise ValueError("There are two types of plot: 'line' "+
                "for data points that are orderly joined by lines, and"+
                " 'scatter' for points that are plotted scatter fashio"+
                "n")

    plt.grid(flag_grid)

    # Plots the set of points as scattered markers

    if highlight_points:

        if multiple_curves:

            for i in range(multiple_curves):

                subplots_tuple.scatter(x_data, y_data[i], color='black', 
                marker='x', zorder=3)

        else:

            subplots_tuple.scatter(x_data, y_data, color='black', 
            marker='x', zorder=3)

    # Verifies if a color bar is asked for

    if color_bar:
        
        # Verifies if the color map has minimum and maximum values

        if not isinstance(color_map, list):

            raise TypeError("The color map does not have minimum and m"+
            "aximum values. Either the data has just one curve and, th"+
            "us cannot possibly have a color map, or a color map has n"+
            "ot been provided")

        # Creates two artificial points to add the color gradient

        if multiple_curves:

            # Verifies the purely scattered case

            if not isinstance(y_data[0], list):

                plotted_entities = subplots_tuple.scatter(x_data[0:2], 
                [y_data[0], y_data[1]], c=[color_map[1], color_map[2]], 
                cmap=color_map[0], vmin=color_map[1], vmax=color_map[2], 
                marker='x', zorder=3, s=0.001)

            # Otherwise, plots elements of the first curve

            else:

                plotted_entities = subplots_tuple.scatter(x_data[0:2], 
                y_data[0][0:2], c=[color_map[1], color_map[2]], cmap=
                color_map[0], vmin=color_map[1], vmax=color_map[2], 
                marker='x', zorder=3, s=0.001)

        else:

            plotted_entities = subplots_tuple.scatter(x_data[0:2], 
            y_data[0:2], c=[color_map[1], color_map[2]], cmap=color_map[
            0], vmin=color_map[1], vmax=color_map[2], marker='x', 
            zorder=3, s=0.001)

        # Creates the color bar

        color_bar = plt.colorbar(plotted_entities)

        if color_barTitle:

            color_bar.set_label(color_barTitle)

        else:

            color_bar.set_label("Magnitude")

        if color_barTicks:

            color_bar.set_ticks(color_barTicks)

        else:

            color_bar.set_ticks(np.linspace(color_map[1], color_map[2], 
            5))

        # Sets a formatter object to ensure that integer ticks are shown
        # as integer

        def formatter_function(x, _):

            if abs(x-int(x)) < 1e-6:

                return f"{int(x):>3}"  # format as integer, right-aligned to 3 spaces
            else:
                return f"{x:>5.2f}"     # Right-align float, width 5

        color_barFormatter = ticker.FuncFormatter(formatter_function)

        color_bar.ax.yaxis.set_major_formatter(color_barFormatter)

    # Applies scientific notation to the ticks

    if flag_scientificNotation:

        subplots_tuple.yaxis.set_major_formatter(ticker.ScalarFormatter(
        useMathText=True))

        subplots_tuple.xaxis.set_major_formatter(ticker.ScalarFormatter(
        useMathText=True))

        subplots_tuple.ticklabel_format(style='sci', axis='both', 
        scilimits=(0,0))

    elif flag_noTicks:

        subplots_tuple.tick_params(axis='both', which='both', length=0, 
        labelbottom=False, labelleft=False)

    # Sets the grid

    if not (x_grid is None):

        subplots_tuple.set_xticks(x_grid)

    if not (y_grid is None):

        subplots_tuple.set_yticks(y_grid)

    # Verifies and uses if necessary other optional attributes

    if not (x_label is None):

        plt.xlabel(x_label)

    if not (y_label is None):

        plt.ylabel(y_label)

    if not (title is None):

        plt.title(title)

    if not (label is None):

        if legend is None:

            plt.legend()

        else:

            plt.legend(loc=legend, bbox_to_anchor=(1.0, 1.0))

    # Adjust the size of the plot to contain ticks and labels

    plt.tight_layout()

    # Saves the plot

    plt.savefig(file_name)

    print("Finishes plotting\n")