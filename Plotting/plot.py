import matplotlib.pyplot as plt
import pandas as pd
import os
from cycler import cycler
import matplotlib
import numpy as np

plt.rcParams.update({'font.size': 12})


def plot_time(csv_files, thresholds, threshold_labels, location, figure_name):
    custom_cycler = (cycler(marker=['o', 'x', 's', 'P', 'd']) +  cycler(ls=['-', '--', ':', '-.', (0, (5, 10))]) + cycler(color=['black', 'green', 'red', 'blue', 'brown']))
    fig, ax = plt.subplots()
    ax.set_prop_cycle(custom_cycler) # setting the cycler for the current figure 

    for csv_file in csv_files:
        # Load CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        
        # Extract x and y data
        x_values = thresholds
        y_values = [None] * len(x_values)
        for i, _ in enumerate(y_values):
            for j, k in enumerate(df["probability"]):
                if k > 10**(x_values[i]):
                    y_values[i] = df["time"][j]
                    break

        index = -1
        for i, value in enumerate(y_values):
            if value == None:
                index = i
                break
        if index != -1:
            for i, _ in enumerate(x_values):
                if i>=index:
                    # x_values[i] = np.nan
                    y_values[i] =np.nan

        # Get the file name without extension
        if "dfs" in csv_file:
            plot_name = "bDFS(\u2205)"
        elif "single_species_plus" in csv_file:
            plot_name = "bDFS(1+)"
        elif "all" in csv_file:
            plot_name = "bDFS(All)"
        elif "single_species" in csv_file:
            plot_name = "bDFS(1)"
        elif "XBF" in csv_file:
            plot_name = "XBF"

    
        # Plot the bar chart
        ax.plot(x_values, y_values, label=plot_name)

    # Add labels and title
    ax.set(xlabel = 'Probability Threshold (log scale)')
    ax.set(ylabel ='Runtime (s)')
    # plt.title('Comparison of Values from Different CSV Files')

    ax.tick_params('x')
    ax.set_xticks(x_values)  # Set the positions of ticks
    ax.set_xticklabels(threshold_labels)  # Set the labels for the ticks
    # Add legend
    ax.legend(loc=location)
    
    # ax.set(xscale ="log")
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    # ax.minorticks_off()

    # Show the plot
    # plt.savefig('./test.png')
    # plt.show()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    plt.savefig(figure_name, bbox_inches='tight')

def plot_size(csv_files, thresholds, threshold_labels, location,  figure_name):
    custom_cycler = (cycler(marker=['o', 'x', 's', 'P', 'd']) +  cycler(ls=['-', '--', ':', '-.', (0, (5, 10))]) + cycler(color=['black', 'green', 'red', 'blue', 'brown']))
    fig, ax = plt.subplots()
    ax.set_prop_cycle(custom_cycler) # setting the cycler for the current figure 

    for csv_file in csv_files:
        # Load CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        
        # Extract x and y data
        x_values = thresholds
        y_values = [None] * len(x_values)
        for i, _ in enumerate(y_values):
            for j, k in enumerate(df["probability"]):
                if k > 10**(x_values[i]):
                    y_values[i] = df["transitions"][j] / 1000000
                    break

        index = -1
        for i, value in enumerate(y_values):
            if value == None:
                index = i
                break
        if index != -1:
            for i, _ in enumerate(x_values):
                if i>=index:
                    # x_values[i] = np.nan
                    y_values[i] =np.nan

        # Get the file name without extension
        if "dfs" in csv_file:
            plot_name = "bDFS(\u2205)"
        elif "single_species_plus" in csv_file:
            plot_name = "bDFS(1+)"
        elif "all" in csv_file:
            plot_name = "bDFS(All)"
        elif "single_species" in csv_file:
            plot_name = "bDFS(1)"
        elif "XBF" in csv_file:
            plot_name = "XBF"

    
        # Plot the bar chart
        ax.plot(x_values, y_values, label=plot_name)

    # Add labels and title
    ax.set(xlabel = 'Probability Threshold (log scale)')
    ax.set(ylabel ='Transitions (x 1E6)')
    # plt.title('Comparison of Values from Different CSV Files')

    ax.tick_params('x')
    ax.set_xticks(x_values)  # Set the positions of ticks
    ax.set_xticklabels(threshold_labels)  # Set the labels for the ticks
    # Add legend
    ax.legend(loc=location)
    
    # ax.set(xscale ="log")
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    # ax.minorticks_off()

    # Show the plot
    # plt.savefig('./test.png')
    # plt.show()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    plt.savefig(figure_name, bbox_inches='tight')
# List of CSV files
#
model_name = "yeast"
thresholds = [-30, -20, -15, -6]
threshold_labels = ["1E-30", "1E-20", "1E-15", "1E-6"]
#
prefix = './csv/' + model_name + '/' + model_name
csv_files = [prefix + '_dfs.csv', prefix + '_single_species.csv', prefix + '_single_species_plus.csv', prefix + '_all.csv', prefix + '_XBF.csv']  # Add your CSV file names here

# Call the function to plot bar charts
plot_size(csv_files, thresholds, threshold_labels, 'upper left', './plots/' + model_name + '_size.pdf')
plot_time(csv_files, thresholds, threshold_labels, 'upper left', './plots/' + model_name + '_time.pdf')