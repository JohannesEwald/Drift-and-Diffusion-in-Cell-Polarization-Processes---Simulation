import numpy as np
import time

# custom function for having coherent naming conventions
def save_to_np(array: np.array, name: str):
    np.save("{}___".format(name) + time.strftime("%Y-%m-%d %H%M%S"), array)


def load_from_np(name: str):
    return np.load(name)

