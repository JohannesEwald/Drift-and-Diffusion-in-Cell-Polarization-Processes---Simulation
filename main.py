import chains as ch
import data_IO_external
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
import math
import MF_Caclulator

resolution = 8     # usually 7 * k + 1
para_min, para_max = -3, 3
half_chain_length = 8
number_of_chains = 1 * 10 ** 3     # 5 * 10 ** 3 or 2
base_simulation_duration = 1 * 10 ** 3


def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


def simulation(parameters: ch.Parameters, base_sim_dur=base_simulation_duration):
    sim_duration = base_sim_dur * 1

    chain_collection = ch.ChainCollection(number_of_chains=number_of_chains, parameters=parameters, simulation_duration=sim_duration)
    for t in range(sim_duration):                           # this order of iteration is very inefficient on RAM, but in current interations allows to save RAW data
        for n in range(number_of_chains):
            chain_collection.activity_chooser(n, t)
        chain_collection.phase_Barrier_History[t] = np.copy(chain_collection.phase_Barriers)
        chain_collection.ATP_History[t] = np.copy(chain_collection.ATPs)

    return np.array(chain_collection.phase_Barrier_History), np.array(chain_collection.ATP_History)


def array_simulation_F_G():
    F_G__not__FmG_H = True
    a, b, c, H = 1.0, 1.0, 1.0, 0

    cumulants_0 = np.zeros(shape=(resolution, resolution))
    cumulants_1 = np.zeros(shape=(resolution, resolution))
    first_passage_time = np.zeros(shape=(resolution, resolution))
    cumulants_0_ATP = np.zeros(shape=(resolution, resolution))

    start_time = time.time()
    for i in range(resolution):
        F = - (para_min + (para_max - para_min) / (resolution - 1) * i)
        # G = -F
        for j in range(resolution):
            G = (para_min + (para_max - para_min) / (resolution - 1) * j)

            parameters = ch.Parameters(half_chain_length=half_chain_length, a=a, b=b, c=c, F=F, G=G, H=H)
            parameters.update_probability_list()
            # parameters.print()

            sim_data_pB_current, sim_data_ATP_current = simulation(parameters)
            # simulation_data_phaseBarrier[i][j] = sim_data_pB_current
            # simulation_data_ATP[i][j] = sim_data_ATP_current

            t_time = [i / parameters.Z for i in range(base_simulation_duration)]
            cumulants_0[i, j], cumulants_1[i, j], first_passage_time[i, j] = data_analysis.data_analysis2(sim_data_pB_current, base_simulation_duration, t_time=t_time)
            cumulants_0_ATP[i, j], _, _ = data_analysis.data_analysis2(sim_data_ATP_current, base_simulation_duration, t_time=t_time)


            print_time()
            percentage = round(100 * (i * resolution + (j + 1)) / (resolution ** 2), 5)
            print("{}%".format(percentage))
            time_left = (time.time() - start_time) / (percentage / 100)
            print(time_left)


    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    if F_G__not__FmG_H == True:
        xlabel = "$\ln(b_+^-)$"
        ylabel = "$\ln(a_+^-)$"
    else:
        xlabel = "$\ln(c^B_A)$"
        ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"

    # dataIO.save_to_np(simulation_data_phaseBarrier, name="sim_data_only_pb")
    dataIO.save_to_np(cumulants_0, name="sim_data_cumulants_0")
    dataIO.save_to_np(cumulants_1, name="sim_data_cumulants_1")
    dataIO.save_to_np(cumulants_0_ATP, name="sim_data_cumulants_0_ATP")
    dataIO.save_to_np(first_passage_time, name="sim_data_first_passage_time")

    plotting.single_visualzation(cumulants_0, "Interface Drift", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cumulants_1, "Interface Diffusion", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel)
    plotting.single_visualzation(first_passage_time, "Domain Return Time", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='YlOrBr')

    plotting.single_visualzation(cumulants_0_ATP, "ATP Consumption", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='magma')
    plotting.single_visualzation(cumulants_0 / cumulants_0_ATP, "Drift per ATP", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cumulants_0 / cumulants_1, "Péclet Number", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color="PiYG")




def array_simulation_FmG_H():
    pass


def array_simulation_constParameter_multipleChainLengths(length_min, length_max, duration, a, b, c, F, G, H, plot_no_data_return = True):
    simulation_data = []
    cumulants_0 = []
    cumulants_1 = []
    first_passage_time = []

    start_time = time.time()

    for n in range(length_min, length_max + 1):
        parameters = ch.Parameters(half_chain_length=n, a=a, b=b, c=c, F=F, G=G, H=H)
        parameters.update_probability_list()
        t_time = [i / parameters.Z for i in range(duration)]

        # data = simulation(parameters, base_sim_dur=duration)[0]
        simulation_data.append(simulation(parameters, base_sim_dur=duration)[0])
        t1, t2, t3 = data_analysis.data_analysis2(simulation_data[-1], duration, t_time)
        cumulants_0.append(t1)
        cumulants_1.append(t2)
        first_passage_time.append(t3)

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
    plt.title(str(a) + str(b) + str(c) + str(F) + str(G) + str(H) + "    prediction: " + str(MF_Caclulator.calculate(a,b, c, F, G, H)[0]))
    parastring = str(str(a) + str(b) + str(c) + str(F) + str(G) + str(H))
    plt.savefig("testing\singlepara/Figure " + time.strftime("%Y-%m-%d %H%M%S") + "---" + parastring + ".png",
                bbox_inches="tight")
    plt.show()
    # plt.plot(x, cumulants_1)
    # plt.xlabel("Length of Half-Chain")
    # plt.ylabel("Interface Diffusion")
    # plt.show()
    # plt.plot(x, first_passage_time)
    # plt.xlabel("Length of Half-Chain")
    # plt.ylabel("Domain Return Time")
    # plt.show()


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


def deviation_heatmap_chainLengthEffect(length_max, duration):

    cum0_map = np.zeros(shape=(resolution, resolution))
    cum0_map_division = np.zeros(shape=(resolution, resolution))
    cum1_map = np.zeros(shape=(resolution, resolution))
    cum1_map_division = np.zeros(shape=(resolution, resolution))
    cum0_atp_map = np.zeros(shape=(resolution, resolution))
    cum0_atp_map_division = np.zeros(shape=(resolution, resolution))

    cum0_trial1 = np.zeros(shape=(resolution, resolution))
    cum0_trial4 = np.zeros(shape=(resolution, resolution))


    if __name__ == "__main__":
        multiprocessing.freeze_support()  # Add this line to support Windows

        with multiprocessing.Pool(processes=8) as pool:
            for i in range(resolution):
                F = -(para_min + (para_max - para_min) / (resolution - 1) * i)
                G = -F
                for j in range(resolution):
                    H = para_min + (para_max - para_min) / (resolution - 1) * j

                    print(F, H)


                    chain_sizes = [n for n in range(1, length_max + 1)]
                    input_parameter = [[F, G, H, n, duration] for n in chain_sizes]

                    results = pool.map(worker_function, input_parameter)

                    cum0s = [element[0] for element in results]
                    cum1s = [element[1] for element in results]
                    cum0s_atp = [element[2] for element in results]

                    cum0_map[i, j] = cum0s[0] - np.average(cum0s[1:])
                    cum1_map[i, j] = cum1s[0] - np.average(cum1s[1:])
                    cum0_atp_map[i, j] = cum0s_atp[0] - np.average(cum0s_atp[1:])

                    cum0_map_division[i, j] = cum0s[0] / np.average(cum0s[1:]) if np.average(cum0s[1:]) != 0 else None
                    cum1_map_division[i, j] = cum1s[0] / np.average(cum1s[1:]) if np.average(cum1s[1:]) != 0 else None
                    cum0_atp_map_division[i, j] = cum0s_atp[0] / np.average(cum0s_atp[1:]) if np.average(cum0s_atp[1:]) != 0 else None

                    cum0_trial1[i, j] = cum0s[0]
                    cum0_trial4[i, j] = cum0s[3]

                    print(results)
                    print_time()
                    percentage = round(100 * (i * resolution + (j + 1)) / (resolution ** 2), 5)
                    print("{}%".format(percentage))


    xlabel = "$\ln(c^B_A)$"
    ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"
    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    plotting.single_visualzation(cum0_map, "Cumulant 0", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cum1_map, "Cumulant 1", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cum0_atp_map, "Cumulant 0 ATP", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    plotting.single_visualzation(cum0_map_division, "Cumulant 0 - Division", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='magma')
    plotting.single_visualzation(cum1_map_division, "Cumulant 1 - Division", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='magma')
    plotting.single_visualzation(cum0_atp_map_division, "Cumulant 0 ATP - Division", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='magma')

    plotting.single_visualzation(cum0_trial1, "Cumulant 1 - trial1", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel)
    plotting.single_visualzation(cum0_trial4, "Cumulant 1 - trial4", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel)

    # data_IO_external.save_to_np(cum0_map, "cum0")
    # data_IO_external.save_to_np(cum1_map, "cum1")
    # data_IO_external.save_to_np(cum0_atp_map, "cum0_atp")
    # data_IO_external.save_to_np(cum0_map_division, "cum0_div")
    # data_IO_external.save_to_np(cum1_map_division, "cum1_div")
    # data_IO_external.save_to_np(cum0_atp_map_division, "cum0_atp_div")

    print("FINISHED")
    return cum0_map





def worker_function(input_paras):
    a, b, c = 1, 1, 1
    F, G, H = input_paras[0], input_paras[1], input_paras[2]
    n = input_paras[3]
    duration = input_paras[4]
    parameters = ch.Parameters(half_chain_length=n, a=a, b=b, c=c, F=F, G=G, H=H)
    parameters.update_probability_list()
    t_time = [i / parameters.Z for i in range(duration)]

    pb, atp = simulation(parameters, duration)
    cum0, cum1, _ = data_analysis.data_analysis2(pb, duration, t_time)
    atp_cum0, _, _ = data_analysis.data_analysis2(atp, duration, t_time)
    return cum0, cum1, atp_cum0


def chain1_vs_chainLarge():
    a, b, c = 1, 1, 1

    cum0_map = np.zeros(shape=(resolution, resolution))
    cum1_map = np.zeros(shape=(resolution, resolution))
    cum0_atp_map = np.zeros(shape=(resolution, resolution))

    for i in range(resolution):
        F = -(para_min + (para_max - para_min) / (resolution - 1) * i)
        G = -F
        for j in range(resolution):
            H = para_min + (para_max - para_min) / (resolution - 1) * j

            _ = MF_Caclulator.calculate(a, b, c, F, G, H)
            cum0_map[i, j], cum0_atp_map[i,j], cum1_map[i,j] = _[0], _[1], _[2]

    xlabel = "$\ln(c^B_A)$"
    ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"
    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    # plotting.single_visualzation(cum0_map, "Cumulant 0 - MF", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    # plotting.single_visualzation(cum1_map, "Cumulant 1 - MF", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')
    # plotting.single_visualzation(cum0_atp_map, "Cumulant 0 ATP - MF", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color='coolwarm')

    if True:
        cum0_map = data_IO_external.load_from_np("sim_data_cumulants_0___2023-11-11 100338.npy")
        cum0_atp_map = data_IO_external.load_from_np("sim_data_cumulants_0_ATP___2023-11-11 100338.npy")
        cum1_map = data_IO_external.load_from_np("sim_data_cumulants_1___2023-11-11 100338.npy")

    center = 0
    color = 'coolwarm'
    cum0_diff = data_IO_external.load_from_np('cum0___2023-11-10 173127.npy')
    cum0_deviation = cum0_diff / cum0_map
    plotting.single_visualzation(cum0_deviation, "Cumulant 0 - relative error", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color=color, center=center)

    cum0_atp_diff = data_IO_external.load_from_np('cum0_atp___2023-11-10 173127.npy')
    cum0_atp_deviation = cum0_atp_diff / cum0_atp_map
    plotting.single_visualzation(cum0_atp_deviation, "Cumulant 0 ATP - relative error", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color=color, center=center)
    cum1_diff = data_IO_external.load_from_np('cum1___2023-11-10 173127.npy')
    cum1_deviation = cum1_diff / cum1_map
    plotting.single_visualzation(cum1_deviation, "Cumulant 1 - relative error", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color=color, center=center)


def theory_against_simulated():
    a, b, c = 1, 1, 1
    simulated_cum0 = data_IO_external.load_from_np("sim_data_cumulants_0___2023-11-08 200413.npy")
    simulated_cum0_atp = data_IO_external.load_from_np("sim_data_cumulants_0_ATP___2023-11-08 200413.npy")
    simulated_diffusion = data_IO_external.load_from_np("sim_data_cumulants_1___2023-11-08 200413.npy")
    cum0_map = np.zeros(shape=(resolution, resolution))
    cum0_atp_map = np.zeros(shape=(resolution, resolution))
    diffusion_map = np.zeros(shape=(resolution, resolution))

    for i in range(resolution):
        F = -(para_min + (para_max - para_min) / (resolution - 1) * i)
        G = -F
        for j in range(resolution):
            H = para_min + (para_max - para_min) / (resolution - 1) * j

            _ = MF_Caclulator.calculate(a, b, c, F, G, H)
            cum0_map[i, j], cum0_atp_map[i,j], diffusion_map[i,j]  = _[0], _[1], _[2]


    xlabel = "$\ln(c^B_A)$"
    ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"
    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    center = 0
    color = "coolwarm"

    diffusion_diff = (simulated_diffusion - diffusion_map) #/ diffusion_map
    plotting.single_visualzation(diffusion_diff , "T-D_diffusion AE", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel,
                                 color=color, center=0)

    cum0_diff = (simulated_cum0 - cum0_map) #/ cum0_map
    plotting.single_visualzation(cum0_diff, "T-D_cum0 AE", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color=color, center=0)
    #
    cum0_atp_diff = (simulated_cum0_atp - cum0_atp_map) #/ cum0_atp_map
    plotting.single_visualzation(cum0_atp_diff, "T-D_cum0 ATP AE", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel, color=color, center=0)


def replot():
    xlabel = "$\ln(c^B_A)$"
    ylabel = "$\ln(a_+^-) = -\ln(b_+^-)$"
    x_ticks = [- round(para_min + (para_max - para_min) / (resolution - 1) * i, 2) for i in range(resolution)]
    y_ticks = [round(para_min + (para_max - para_min) / (resolution - 1) * j, 2) for j in range(resolution)]

    plot = data_IO_external.load_from_np("cum0_atp___2023-11-10 173127.npy")
    plotting.single_visualzation(plot, "Cumulant 0 ATP", x_ticks, y_ticks, xlabel=xlabel, ylabel=ylabel,
                                 color='coolwarm', center=0)


if __name__ == "__main__":
    # array_simulation_constParameter_multipleChainLengths2(1, 3, 5 * 10 ** 2, 1, 1, 1, 1, 1, 1)
    # array_simulation_constParameter_multipleChainLengths(1, 8,  2 * 10 ** 4, 1, 1, 1, 2, -2, 2)
    array_simulation_F_G()

    # for k in range(10):
    deviation_heatmap_chainLengthEffect(8, 10 ** 3)
    # chain1_vs_chainLarge()
    # theory_against_simulated()
    # replot()

    pass





