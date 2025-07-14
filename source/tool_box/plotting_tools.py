# Routine to store methods to plot

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
color_map=False, flag_noTicks=False):

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

        try:

            color_map = plt.get_cmap(color_map)

        except Exception as error_message:

            print("Error Message:"+str(error_message)+"\nProbably this"+
            " color map does not exist")

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

            for i in range(len(color)):

                if isinstance(color[i], float) or isinstance(color[i], 
                int):
                    
                    color[i] = color_map(color[i])
            
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

                if color_map and (isinstance(color, float) or isinstance(
                color, int)):

                    color = color_map(color[0])

                else:

                    color = color[0]

    # Plots it

    if label is None:

        if multiple_curves:

            for i in range(multiple_curves):

                if plot_type=="line":

                    subplots_tuple.plot(x_data, y_data[i], linestyle=
                    element_style[i], linewidth=element_size[i], color=
                    color[i])

                elif plot_type=="scatter":

                    # If multiple curves are plotted at once

                    if isinstance(y_data[i], list):

                        subplots_tuple.scatter(x_data, y_data[i], marker=
                        element_style[i], s=element_size[i]**2, color=
                        color[i], zorder=3)

                    # If it is just the default treatment of scatter 
                    # plots

                    else:

                        subplots_tuple.scatter(x_data[i], y_data[i], 
                        marker=element_style[i], s=element_size[i]**2, 
                        color=color[i], zorder=3)

                else:

                    raise ValueError("There are two types of plot: 'li"+
                    "ne' for data points that are orderly joined by li"+
                    "nes, and 'scatter' for points that are plotted sc"+
                    "atter fashion")

        else:

            if plot_type=="line":

                subplots_tuple.plot(x_data, y_data, linestyle=
                element_style, linewidth=element_size, color=color)

            elif plot_type=="scatter":

                subplots_tuple.scatter(x_data, y_data, marker=
                element_style, s=element_size**2, color=color, zorder=3)

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

                    subplots_tuple.plot(x_data, y_data[i], linestyle=
                    element_style[i], linewidth=element_size[i], color=
                    color[i], label=label[i])

                elif plot_type=="scatter":

                    # If multiple curves are plotted at once

                    if isinstance(y_data[i], list):

                        subplots_tuple.scatter(x_data, y_data[i], marker=
                        element_style[i], s=element_size[i]**2, color=
                        color[i], zorder=3)

                    # If it is just the default treatment of scatter 
                    # plots

                    else:

                        subplots_tuple.scatter(x_data[i], y_data[i], 
                        marker=element_style[i], s=element_size[i]**2, 
                        color=color[i], zorder=3)

                else:

                    raise ValueError("There are two types of plot: 'li"+
                    "ne' for data points that are orderly joined by li"+
                    "nes, and 'scatter' for points that are plotted sc"+
                    "atter fashion")

        else:

            if plot_type=="line":

                subplots_tuple.plot(x_data, y_data, linestyle=
                element_style, linewidth=element_size, color=color,
                label=label)

            elif plot_type=="scatter":

                subplots_tuple.scatter(x_data, y_data, marker=
                element_style, s=element_size**2, color=color, label=
                label, zorder=3)

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

    # Applies scientific notation to the ticks

    if flag_scientificNotation:

        subplots_tuple.yaxis.set_major_formatter(ticker.ScalarFormatter(
        useMathText=True))

        subplots_tuple.xaxis.set_major_formatter(ticker.ScalarFormatter(
        useMathText=True))

        subplots_tuple.ticklabel_format(style='sci', axis='both', 
        scilimits=(0,0))

    elif flag_noTicks:

        subplots_tuple.set_xticks([])

        subplots_tuple.set_yticks([])

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