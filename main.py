import chains as ch
import data_IO_external as dataIO
import data_analysis
import random
import numpy as np
import matplotlib.pyplot as plt
import os, psutil
from datetime import datetime
import time
import multiprocessing
import plotting

resolution = 29     # usually 25
para_min, para_max = -5, 5
half_chain_length = 1
number_of_chains = 5 * 10 ** 2     # 5 * 10 ** 3 or 2
base_simulation_duration = 5 * 10 ** 2


def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


def simulation(parameters: ch.Parameters, base_sim_dur=base_simulation_duration):
    sim_duration = base_sim_dur * 1

    chain_collection = ch.ChainCollection(number_of_chains=number_of_chains, parameters=parameters, simulation_duration=sim_duration)
    for t in range(sim_duration):
        for n in range(number_of_chains):
            chain_collection.activity_chooser(n, t)
        chain_collection.phase_Barrier_History[t] = np.copy(chain_collection.phase_Barriers)
        chain_collection.ATP_History[t] = np.copy(chain_collection.ATPs)

    return np.array(chain_collection.phase_Barrier_History), np.array(chain_collection.ATP_History)     # multiplied by parameters.Z earlier


def array_simulation_F_G():
    F_G__not__FmG_H = False
    a, b, c, H = 1.0, 1.0, 1.0, 0
    sb_workaround1 = [None for i in range(resolution)]
    sb_workaround2 = [sb_workaround1 for i in range(resolution)]
    simulation_data_phaseBarrier = sb_workaround2                   # np.empty(shape=(2, resolution, resolution))
    simulation_data_ATP = np.copy(sb_workaround2)

    print(simulation_data_phaseBarrier)

    cumulants_0 = np.zeros(shape=(resolution, resolution))
    cumulants_1 = np.zeros(shape=(resolution, resolution))
    first_passage_time = np.zeros(shape=(resolution, resolution))
    cumulants_0_ATP = np.zeros(shape=(resolution, resolution))

    x_ticks = [para_min + (para_max - para_min) / (resolution - 1) * i for i in range(resolution)]
    y_ticks = [para_min + (para_max - para_min) / (resolution - 1) * j for j in range(resolution)]

    start_time = time.time()
    for i in range(resolution):
        F = - (para_min + (para_max - para_min) / (resolution - 1) * i)
        G = -F
        for j in range(resolution):
            H = (para_min + (para_max - para_min) / (resolution - 1) * j)

            parameters = ch.Parameters(half_chain_length=half_chain_length, a=a, b=b, c=c, F=F, G=G, H=H)
            parameters.update_probability_list()
            # parameters.print()

            simulation_data_phaseBarrier[i][j], simulation_data_ATP[i][j] = simulation(parameters)
            t_time = [i / parameters.Z for i in range(base_simulation_duration)]
            cumulants_0[i, j], cumulants_1[i, j], first_passage_time[i, j] = data_analysis.data_analysis2(simulation_data_phaseBarrier[i][j], base_simulation_duration, t_time=t_time)
            cumulants_0_ATP[i, j], _, _ = data_analysis.data_analysis2(simulation_data_ATP[i][j], base_simulation_duration, t_time=t_time)


            print_time()
            percentage = round(100 * (i * resolution + (j + 1)) / (resolution ** 2), 5)
            print("{}%".format(percentage))
            time_left = (time.time() - start_time) / (percentage / 100)
            print(time_left)

    # cumulants_1 = cumulants_1 / parameters.Z
    cumulants_0 = cumulants_0 / 2
    cumulants_0_ATP = cumulants_0_ATP / 2

    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]
    # x_ticks = [-x for x in x_ticks]

    if F_G__not__FmG_H == True:
        xlabel = "$\ln(b_+^-)$"
        ylabel = "$\ln(a_+^-)$"
    else:
        xlabel = "$\ln(c^B_A)$"
        ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"

    plotting.single_visualzation(cumulants_0, "Interface Drift", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cumulants_1, "Interface Diffusion", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel)
    plotting.single_visualzation(first_passage_time, "Domain Return Time", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='YlOrBr')

    plotting.single_visualzation(cumulants_0_ATP, "ATP Consumption", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='magma')
    plotting.single_visualzation(cumulants_0 / cumulants_0_ATP, "Drift per ATP", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cumulants_0 / cumulants_1, "Fano Factor", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel)



    # dataIO.save_to_np(simulation_data_phaseBarrier, name="sim_data_only_pb")
    dataIO.save_to_np(cumulants_0, name="sim_data_cumulants_0")
    dataIO.save_to_np(cumulants_1, name="sim_data_cumulants_1")
    dataIO.save_to_np(cumulants_0_ATP, name="sim_data_cumulants_0_ATP")
    dataIO.save_to_np(first_passage_time, name="sim_data_first_passage_time")

    #plt.imshow(simulation_data[0])
    #plt.show()
    #plt.imshow(simulation_data[1])
    #plt.show()
    #plt.imshow(simulation_data[0] / simulation_data[1])
    #plt.show()
    #plt.imshow(cumulants_0)
    #plt.show()
    #plt.imshow(cumulants_1)
    #plt.show()
    #plt.imshow(first_passage_time)
    #plt.show()


def array_simulation_FmG_H():
    pass


def array_simulation_constParameter_multipleChainLengths(length_min, length_max, duration, a, b, c, F, G, H):
    simulation_data = []
    cumulants_0 = []
    cumulants_1 = []
    first_passage_time = []

    start_time = time.time()

    for n in range(length_min, length_max + 1):
        parameters = ch.Parameters(half_chain_length=n, a=a, b=b, c=c, F=F, G=G, H=H)
        parameters.update_probability_list()
        simulation_data.append(simulation(parameters, base_sim_dur=duration)[0])
        t1, t2, t3 = data_analysis.data_analysis2(simulation_data[-1], duration)
        cumulants_0.append(t1)
        cumulants_1.append(t2 / parameters.Z)
        first_passage_time.append(t3 / (parameters.Z / 2))

        print_time()
        percentage = round(100 * n / (length_max + 1 - length_min), 5)
        print("{}%".format(percentage))
        time_total_guess = (time.time() - start_time) / (percentage / 100)
        print((100 - percentage) * time_total_guess/ 100)
        print(time_total_guess)

    x = range(length_min, length_max + 1)
    print(parameters.Z)
    print(np.polyfit(x, cumulants_0, 1))
    print(np.polyfit(x, cumulants_1, 1))
    print(np.polyfit(x, first_passage_time, 1))

    plt.plot(x, cumulants_0)
    plt.xlabel("Length of Half-Chain")
    plt.ylabel("Interface Drift")
    plt.show()
    plt.plot(x, cumulants_1)
    plt.xlabel("Length of Half-Chain")
    plt.ylabel("Interface Diffusion")
    plt.show()
    plt.plot(x, first_passage_time)
    plt.xlabel("Length of Half-Chain")
    plt.ylabel("Domain Return Time")
    plt.show()


def multiprocessing_helper_constParameter_multipleChainLengths(sim_args):           # sim_args[0] = parameters, sim_args[1] = duration
    sim_args[0].update_probability_list()
    sim_data = simulation(parameters=sim_args[0], base_sim_dur=sim_args[1])
    sim_args[0].print()
    t1, t2, t3 = data_analysis.data_analysis2(sim_data, sim_args[1])
    return t1


def array_simulation_constParameter_multipleChainLengths2(length_min, length_max, duration, a, b, c, F, G, H):
    repetitions = 30
    final_results = np.zeros(length_max + 1 - length_min)

    for _ in range(repetitions):
        print("Repetition:", _)
        if __name__ == "__main__":
            # Define the range of input values for the loop
            input_values = [[ch.Parameters(half_chain_length=n, a=a, b=b, c=c, F=F, G=G, H=H), duration] for n in range(length_min, length_max + 1)]

            # Create a Pool object with the number of CPU cores
            with multiprocessing.Pool() as pool:
                results = pool.map(multiprocessing_helper_constParameter_multipleChainLengths, input_values)

        final_results = np.add(final_results, results)

    print(final_results)
    plt.plot(final_results)
    plt.show()


# array_simulation_constParameter_multipleChainLengths2(1, 3, 5 * 10 ** 2, 1, 1, 1, 1, 1, 1)
# array_simulation_constParameter_multipleChainLengths(1, 8,  1 * 10 ** 4, 1, 1, 1, 1, -1, -1)
array_simulation_F_G()