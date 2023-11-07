import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import kstat       # make sure to import this into the interpreter!!

# adjust base_simulation_duration and actual simulation_duration for return of this function? Is
# --> simulation duraction is only adjusted regarding different chain_length, right?


linear_regression = True       # linear regression or simple linear approximation by first and last data point
sqrt_first_passage_time = False
plot_cumulants = False

Z = 1


def data_analysis2(phase_barrier_history, base_simulation_duration, t_time):   # , parameters: Parameters):      # phase_barrier_history = pb[t,n]
    # t_time = [i for i in range(base_simulation_duration)]           # actually simulation duration instead of base simulation duration?
    cummulants = [[], [], []]           # average, variance, skewness
    for time_frame in list(phase_barrier_history):
        cummulants[0].append(kstat(time_frame, n=1))
        cummulants[1].append(kstat(time_frame, n=2))
        # cummulants[2].append(kstat(time_frame, n=3))      # we dont need that currently

    if plot_cumulants:
        print(cummulants)
        plot_three_plots(cummulants[0], cummulants[1], cummulants[2])
        plt.show()

    # t_time = [i for i in range(len(cummulants[0]))]     # rescale time instead of rescaling data!!!

    # WHAT DOES THIS CODE DO? First passage time calculation??
    PH_f_swapped = np.swapaxes(phase_barrier_history, 0, 1)         # pb[n,t]
    distances = []
    for chain in PH_f_swapped:

        chain_entries = set(chain)
        for position in chain_entries:
            position_indices = np.diff(np.where(chain == position)[0])
            index = np.where(position_indices == 1)
            # print(position_indices)
            # print(index)
            position_indices = np.delete(position_indices, index)
            # print(position_indices)
            # wait()

            if len(position_indices) > 0:
                distances.append(np.average(position_indices))
            else:
                distances.append(base_simulation_duration)

    distances = np.array(distances) # * base_simulation_duration / simulation_duration


    if linear_regression:
        cumulant_vector = [np.polyfit(x=t_time, y=cummulants[0], deg=1)[0] / Z, np.polyfit(x=t_time, y=cummulants[1], deg=1)[0] / Z, 1]
    else:
        cumulant_vector = [cummulants[0][-1] / base_simulation_duration, cummulants[1][-1] / base_simulation_duration, 1]

    if sqrt_first_passage_time:
        return cumulant_vector[0], cumulant_vector[1], math.sqrt( np.average(distances) / base_simulation_duration * 100 )
    else:
        return cumulant_vector[0], cumulant_vector[1], np.average(distances) / base_simulation_duration * 100




def plot_three_plots(x1, x2, x3):
    fig, axs = plt.subplots(3)
    fig.tight_layout()
    axs[0].plot(x1)
    axs[1].plot(x2)
    axs[2].plot(x3)
    plt.show()
    



