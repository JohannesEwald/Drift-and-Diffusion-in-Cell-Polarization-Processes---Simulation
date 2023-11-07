import numpy as np
import time


def save_to_np(array: np.array, name: str):
    np.save("{}___".format(name) + time.strftime("%Y-%m-%d %H%M%S"), array)


def load_from_np(name: str):
    pass
