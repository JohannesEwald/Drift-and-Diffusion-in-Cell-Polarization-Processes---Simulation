import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
from scipy.ndimage import gaussian_filter
import random
import pandas
import data_analysis

# refer to this https://matplotlib.org/stable/tutorials/pyplot.html
# https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame
# https://seaborn.pydata.org/tutorial/introduction.html

# https://stackoverflow.com/questions/76996936/creating-a-heatmap-from-data-that-contains-counts?noredirect=1&lq=1

plot_3d_Fanofactor = False
spin_diviation_plot_parameter = 0.01


def plot_three_plots(x1, x2, x3):
    fig, axs = plt.subplots(3)
    fig.tight_layout()
    axs[0].plot(x1)
    axs[1].plot(x2)
    axs[2].plot(x3)
    plt.show()



# currently dead function
def six_visualization(title, parameter_list, plot1, plot2, plot3, extra_plot1, extra_plot2, extra_plot3):
    y_ticks = parameter_list[0]
    x_ticks = parameter_list[1]

    for k in range(len(x_ticks)):
        if k % 3 == 1 or k % 3 == 2:
            x_ticks[k], y_ticks[k] = None, None


    seee_plot = plot1 / data_atp_average
    # seee_plot = np.log(1 - (np.abs(plot1 / data_atp_average)))


    fig, axs = plt.subplots(2, 4, figsize=(12, 5.6))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95], pad=3)
    # fig.suptitle("half chain length" + str(title) + "")


    cbar_ticks_spin_diviation = {"ticks":[1, 1+spin_diviation_plot_parameter/2, 1-spin_diviation_plot_parameter/2, 1+spin_diviation_plot_parameter, 1-spin_diviation_plot_parameter]}

    sns.heatmap(plot1, cbar=True, ax=axs[0, 0], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)  # cbar_kws={"ticks":[0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25]})
    sns.heatmap(plot2, cbar=True, ax=axs[0, 1], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)
    # sns.heatmap(data_pb_skew, cbar=True, ax=axs[0, 2], cmap='viridis')
    sns.heatmap(extra_plot1, cbar=True, ax=axs[0, 2], cmap='coolwarm', xticklabels=x_ticks, yticklabels=y_ticks, center=1, vmin=1-spin_diviation_plot_parameter, vmax=1+spin_diviation_plot_parameter, cbar_kws=cbar_ticks_spin_diviation)
    sns.heatmap(plot1 / plot2, cbar=True, ax=axs[0, 3], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)


    sns.heatmap(data_atp_average, cbar=True, ax=axs[1, 0], cmap='magma', xticklabels=x_ticks, yticklabels=y_ticks)
    # sns.heatmap(data_atp_variance, cbar=True, ax=axs[1, 1], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)
    sns.heatmap(extra_plot3, cbar=True, ax=axs[1, 1], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)
    # sns.heatmap(extra_plot2, cbar=True, ax=axs[1, 2], cmap='coolwarm', xticklabels=x_ticks, yticklabels=y_ticks, center=1, vmin=1-spin_diviation_plot_parameter, vmax=1+spin_diviation_plot_parameter)
    sns.heatmap(gaussian_filter(extra_plot2, sigma=1), cbar=True, ax=axs[1, 2], cmap='coolwarm', xticklabels=x_ticks, yticklabels=y_ticks, center=1, vmin=1-spin_diviation_plot_parameter, vmax=1+spin_diviation_plot_parameter)
    # sns.heatmap(data_atp_average / data_atp_variance, cbar=True, ax=axs[1, 3], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)
    sns.heatmap(seee_plot, cbar=True, ax=axs[1, 3], cmap='viridis', xticklabels=x_ticks, yticklabels=y_ticks)
    # axs[1,3].contour(np.arange(.5, seee_plot.shape[1]), np.arange(.5, seee_plot.shape[0]), seee_plot, colors='red', alpha=0.42)

    axs[0, 0].set_title('Interface Drift')
    axs[0, 1].set_title('Interface Diffusion')
    axs[0, 2].set_title('A-Spin-Deviation')
    axs[0, 3].set_title('Interface Fanofactor')
    axs[1, 0].set_title('ATP Consumption')
    axs[1, 1].set_title('Domain Return Time')
    axs[1, 2].set_title('B-Spin_Deviation')
    axs[1, 3].set_title('Interface Direction per ATP')

    # THIS IS FOR FmG_H
    # plt.setp(axs[-1, :], xlabel='$\ln(c^B_A)$')    # bottom label
    # plt.setp(axs[:, 0], ylabel='$\ln(a^+_-)=-\ln(b^+_-)$')     # left label

    # THIS IS FOR F_G
    plt.setp(axs[-1, :], xlabel='$\ln(b^+_-)$')    # bottom label
    plt.setp(axs[:, 0], ylabel='$\ln(a^+_-)$')     # left label

    fig.text(s=f'half chain length = {half_chain_length} | base simulation duration = {base_simulation_duration} | #chains = {number_of_chains} | Time (not steps) = {time_not_steps}', x=0.5, y=0.94, horizontalalignment='center', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))
    plt.subplots_adjust(top=0.83)

    time.strftime("%Y-%m-%d %H%M%S")
    plt.savefig("Data_Pics/Figure " + time.strftime("%Y-%m-%d %H%M%S") + "---" + str(title) + ".svg")
    plt.savefig("Data_Pics/Figure " + time.strftime("%Y-%m-%d %H%M%S") + "---" + str(title) + ".png")

    plt.close()

    if plot_3d_Fanofactor:
        z = plot1 / plot2
        x, y = np.meshgrid(range(z.shape[0]), range(z.shape[1]))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z, cmap=cm.coolwarm)
        plt.title('z as 3d height map')
        html_str = mpld3.fig_to_html(fig)
        Html_file = open("3D-Plot" + time.strftime("%Y-%m-%d %H%M%S") + ".html", "w")
        Html_file.write(html_str)
        Html_file.close()

        pickle.dump(fig, open('FigureObject' + time.strftime("%Y-%m-%d %H%M%S") + '.fig.pickle', 'wb'))
        figx = pickle.load(open('FigureObject' + time.strftime("%Y-%m-%d %H%M%S") + '.fig.pickle', 'rb'))

        figx.show()  # Show the figure, edit it, etc.!
        plt.show()


def single_visualzation(data_matrix, title, x_ticks, y_ticks, xlabel="xläbel", ylabel="yläbel", color='viridis', center=None, dontSaveOnlyPlot = False):
    dataframe = pandas.DataFrame(data_matrix, index=x_ticks, columns=y_ticks)

    fig = plt.figure()
    ax = plt.axes()

    svm = sns.heatmap(dataframe, ax=ax, cbar=True, cmap=color, xticklabels=7, yticklabels=7, center=center)

    plt.title(title, fontsize =20)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    # plt.xticks(rotation=45)
    plt.yticks(rotation=0)

    plt.savefig("testing\pictures/Figure " + time.strftime("%Y-%m-%d %H%M%S") + "---" + str(title) + ".png", bbox_inches="tight")


def plot_and_data_analysis(filename):

    para_min, para_max = -2, 2
    half_chain_length = 1
    number_of_chains = 500
    base_simulation_duration = 1000

    simulation_data = np.load(filename)
    data_shape = simulation_data.shape
    print("data shape = " + str(data_shape))
    resolution = data_shape[0]
    simulation_data = simulation_data.tolist()
    print(simulation_data)

    cumulants_0 = np.zeros(shape=(resolution, resolution))
    cumulants_1 = np.zeros(shape=(resolution, resolution))
    first_passage_time = np.zeros(shape=(resolution, resolution))

    for i in range(resolution):
        print(round(i / resolution * 100, 2))
        for j in range(resolution):
            cumulants_0[i, j], cumulants_1[i, j], first_passage_time[i, j] = data_analysis.data_analysis2(simulation_data[i][j], base_simulation_duration)

    x_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    single_visualzation(cumulants_0, "0th cumulant", x_ticks, y_ticks)
    single_visualzation(cumulants_1, "1th cumulant", x_ticks, y_ticks)
    single_visualzation(first_passage_time, "fpt", x_ticks, y_ticks)




if False:
    n = 50
    test_array = np.random.rand(n, n)
    x = [n / 100 for n in range(n)]
    y = range(n)
    print(test_array)

    single_visualzation(test_array, title = "testtitle", x_ticks=x, y_ticks=y, xlabel="$a_-^+$", ylabel="$b_-^+$")
    pass