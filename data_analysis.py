import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import kstat


linear_regression = True       # linear regression or simple linear approximation by first and last data point
sqrt_first_passage_time = False
plot_cumulants = False



def data_analysis2(phase_barrier_history, base_simulation_duration, t_time, return_SE=False):       # phase_barrier_history = pb[t,n]
    cummulants = [[], [], []]           # average, variance, skewness
    for time_frame in list(phase_barrier_history):
        cummulants[0].append(kstat(time_frame, n=1))
        cummulants[1].append(kstat(time_frame, n=2))
        # cummulants[2].append(kstat(time_frame, n=3))

    if plot_cumulants:
        print(cummulants)
        plot_three_plots(cummulants[0], cummulants[1], cummulants[2])
        plt.show()

    # calculating return times (only qualitativly due to finite simulation time)
    PH_f_swapped = np.swapaxes(phase_barrier_history, 0, 1)         # pb[n,t]
    distances = []
    for chain in PH_f_swapped:

        chain_entries = set(chain)
        for position in chain_entries:
            position_indices = np.diff(np.where(chain == position)[0])
            index = np.where(position_indices == 1)
            position_indices = np.delete(position_indices, index)

            if len(position_indices) > 0:
                distances.append(np.average(position_indices))
            else:
                distances.append(base_simulation_duration)
    distances = np.array(distances)

    # Returns

    if linear_regression:
        cumulant_vector = [np.polyfit(x=t_time, y=cummulants[0], deg=1)[0] / Z, np.polyfit(x=t_time, y=cummulants[1], deg=1)[0] / Z, 1]
    else:
        cumulant_vector = [cummulants[0][-1] / base_simulation_duration, cummulants[1][-1] / base_simulation_duration, 1]

    if sqrt_first_passage_time:
        return cumulant_vector[0], cumulant_vector[1], math.sqrt( np.average(distances) / base_simulation_duration * 100 )
    else:
        return cumulant_vector[0], cumulant_vector[1], np.average(distances) / base_simulation_duration * 100






