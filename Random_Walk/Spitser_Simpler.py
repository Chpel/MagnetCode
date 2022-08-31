import numpy as np
from numba import njit, jit, prange

@jit(parallel=True)
def create_walk():
    Pb = np.array([
        [0.36477068, 0.05275528, 0.16623702, 0.16623702, 0.25      ],
        [0.05275528, 0.36477068, 0.16623702, 0.16623702, 0.25      ],
        [0.16623702, 0.16623702, 0.36477068, 0.05275528, 0.25      ],
        [0.16623702, 0.16623702, 0.05275528, 0.36477068, 0.25      ],
        [0      , 0        , 0      , 0      , 1.      ]], dtype=np.float_)
    Xs = np.zeros(100, dtype=np.uint8)
    Xs[0] = np.random.choice(4)
    X_end = 0
    for i in range(1,100):
        Xs[i] = np.random.choice(5, None, True, Pb[Xs[i-1]])
        if Xs[i] == 4:
            return Xs[:i]
    return Xs

@njit
def linear_unique(X):
    count = np.zeros(4, dtype=np.uint8)
    for i in range(X.shape[0]):
        count[X[i]] += 1
    return (count > 0).sum()
        

@jit
def experiment():
    Xs = create_walk()
    return linear_unique(Xs)

@jit            
def complex_experiment(Ni):
    res_log = np.zeros(5, dtype=np.int_)
    for i in range(Ni):
        res = experiment()
        res_log[res - 1] += 1
    res_log[4] = Ni
    return res_log

def write_out(res):
    res_float = np.array([0, 0, 0, 0, res[-1]], dtype=np.float64)
    res_float[:-1] = res[:-1] / res[-1]
    np.savetxt('Spitser_simpler.txt', 
               [res_float], 
               delimiter=' ', 
               fmt=['%1.6f', '%1.6f', '%1.6f', '%1.6f','%d'], 
               newline=' ', 
               header='p_N1  p_N2  p_N3  p_N4  steps\n',
               comments='')

def main_func(Ni, Nstop):
    result_log = np.zeros(5, dtype=np.ulonglong)
    k = 0
    while k < Nstop:
        result_log = result_log + complex_experiment(Ni)
        write_out(result_log)
        k += 1
        
from sys import argv

script, step, stop = argv
main_func(int(step), int(stop))