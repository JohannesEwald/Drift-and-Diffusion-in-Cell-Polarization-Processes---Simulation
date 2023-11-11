import math
import random
import numpy as np


class Parameters:

    def __init__(self, half_chain_length=1, a=1.0, b=1.0, c=1.0, F=1.0, G=1.0, H=1.0):   #int, a: float, b: float, c: float, F: float, G: float, H: float):
        self.half_chain_length = half_chain_length
        self.a, self.b, self.c, self.F, self.G, self.H = a, b, c, F, G, H
        self.probability_list = np.array([0.0] * 6)  # [a+, a-, b+, b-, c+, c-]
        self.activity_probability = 0
        self.flip_cases_probabilities = []
        self.conversion_probability = 0
        self.A_thermalization_probability = 0
        self.B_thermalization_probability = 0
        self.Z = 0

    def print(self):
        print("(a,b,c,F,G,H) = ({},{},{},{},{},{})".format(self.a, self.b, self.c, self.F, self.G, self.H))
        print("probability list = {}".format(self.probability_list))
        print("activity_probabilities = {}".format(self.activity_probability))
        print("flip_cases_probabilities = {}".format(self.flip_cases_probabilities))
        print("conversion_probabilities = {}".format(self.conversion_probability))
        print("A_thermalization_probability = {}".format(self.A_thermalization_probability))
        print("B_thermalization_probability = {}".format(self.B_thermalization_probability))

    # creating cumulative distribution list for all reactions / activities
    def update_probability_list(self):
        base_rate = [self.half_chain_length * self.a, self.half_chain_length * self.b, self.c]    # multiply with half chain_length here as for increasing chain_length we want more flips
        exponential_rate = [self.F, self.G, self.H]

        # creative cumulative rate list
        for j in range(3):
            self.probability_list[2 * j] = self.probability_list[2 * j - 1] + base_rate[j] * math.exp(exponential_rate[j] / 2)      # 0 2 4 in array are: a+, b+, c+
            self.probability_list[2 * j + 1] = self.probability_list[2 * j] + base_rate[j] * math.exp(-exponential_rate[j] / 2)     # 1 3 5 in array are: a-, b-, c-

        self.Z = self.probability_list[-1]
        # print(self.probability_list)

        # normalize for rate -> probability
        self.probability_list = np.array(self.probability_list) / self.probability_list[-1]     # normalization
        # --- #
        self.activity_probability = self.probability_list[3]

        self.flip_cases_probabilities = self.probability_list[:4] / self.probability_list[3] if self.probability_list[3] != 0 else 0
        self.A_thermalization_probability = (self.flip_cases_probabilities[:2] / self.flip_cases_probabilities[1])[0]
        self.B_thermalization_probability = (self.flip_cases_probabilities[2:4] - self.flip_cases_probabilities[1])
        self.B_thermalization_probability = self.B_thermalization_probability[0] / self.B_thermalization_probability[-1]

        self.conversion_probability = self.probability_list[4:6] - self.probability_list[3]
        self.conversion_probability = (self.conversion_probability / self.conversion_probability[-1])[0] if self.conversion_probability[-1] != 0 else None    # renormalization
        # here we also catch the case for c=0, as we else do 0-division. We set the value arbitrary to 0 as it will not be used


class ChainCollection:
    def __init__(self, number_of_chains: int, parameters: Parameters, simulation_duration: int):
        self.parameters = parameters
        self.chains = np.random.randint(2, size=(number_of_chains, 2 * self.parameters.half_chain_length))              # state of current chain
        self.phase_Barriers = np.zeros(shape=number_of_chains, dtype=np.intc)                                          # [t]
        self.ATPs = np.zeros(shape=number_of_chains, dtype=np.intc)
        self.ATP_History = np.zeros(shape=(simulation_duration, number_of_chains), dtype=np.intc)
        # self.chains_History = [np.copy(self.chains)] * simulation_duration                                            # [t,n,pos]               # THESE ARE VERY LARGE --> TOO MUCH RAM
        self.phase_Barrier_History = np.zeros(shape=(simulation_duration, number_of_chains), dtype=np.intc)            # [t,n]                   # THESE ARE VERY LARGE --> TOO MUCH RAM

    def activity_chooser(self, chain_position: int, time: int):
        if random.random() < self.parameters.activity_probability:
            self.charge_flip(chain_position)
            pb_change = 0
        else:
            pb_change = self.phase_barrier_shift(chain_position, time)

        self.phase_Barriers[chain_position] = self.phase_Barriers[chain_position] + pb_change
        self.ATPs[chain_position] = self.ATPs[chain_position] + abs(pb_change)

    def charge_flip(self, chain_position: int):
        flip_prob = random.random()
        flip_position = random.randint(0, self.parameters.half_chain_length - 1)

        if flip_prob < self.parameters.flip_cases_probabilities[0]:
            self.chains[chain_position, flip_position] = 1
        elif flip_prob < self.parameters.flip_cases_probabilities[1]:
            self.chains[chain_position, flip_position] = 0
        elif flip_prob < self.parameters.flip_cases_probabilities[2]:
            flip_position += self.parameters.half_chain_length
            self.chains[chain_position, flip_position] = 1
        else:
            flip_position += self.parameters.half_chain_length
            self.chains[chain_position, flip_position] = 0


    def phase_barrier_shift(self, chain_position: int, time: int):
        if self.chains[chain_position, self.parameters.half_chain_length - 1] == 1:
            if self.chains[chain_position, self.parameters.half_chain_length] == 0:
                if random.random() < self.parameters.conversion_probability:
                    self.chains[chain_position] = np.roll(self.chains[chain_position], -1)
                    self.chains[chain_position, -1] = 1 if random.random() < self.parameters.B_thermalization_probability else 0                    # thermalized B gets added to the right due to shift

                    return 1
        elif self.chains[chain_position, self.parameters.half_chain_length] == 1:
            if random.random() > self.parameters.conversion_probability:
                self.chains[chain_position] = np.roll(self.chains[chain_position], 1)
                self.chains[chain_position, 0] = 1 if random.random() < self.parameters.A_thermalization_probability else 0
                return -1
        return 0








