import numpy as np
from numba import njit
from time import time

@njit
def experiment(phi, R, N):
    x = np.array([round(R * np.cos(phi)), round(R * np.sin(phi))], dtype=np.int_)
    kstop = R * 1000000
    neig = np.zeros(9, dtype=np.short)
    steps = np.array([[1,0],[-1,0],[0,1],[0,-1]], dtype=np.int_)
    k = 0
    while True:
        walk_directions = np.random.randint(0,4, size=N)
        for i in range(1,N):
            x = x + steps[walk_directions[i]]
            if (np.abs(x) <= 1).all():
                if (x == 0).all():
                    return neig[1] + neig[3] + neig[5] + neig[7], k
                else:
                    neig[4 + 3 * x[1] + x[0]] = 1
        k += N
        if k > kstop:
            return 0, 0
        
@njit            
def complex_experiment(R, N, Ni):
    res_log = np.zeros(6, dtype=np.int_)
    for i in range(Ni):
        phi = 2 * np.pi * np.random.rand()
        res = experiment(phi, R, N)
        res_log[res[0]] += 1
    res_log[5] = Ni
    return res_log
    
def write_out(res, R):
    res_float = np.array([R, 0, 0, 0, 0, 0, res[-1]], dtype=np.float_)
    res_float[1:-1] = res[:-1] / res[-1]
    np.savetxt('Spitser_R_' + str(R) + '.txt', 
               [res_float], 
               delimiter=' ', 
               fmt=['%d', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f','%d'], 
               newline=' ', 
               header=' R  pError  p_N1  p_N2  p_N3  p_N4  steps\n',
               comments='')

def main_func(R, Nmem, Ntime):
    result_log = np.zeros(6, dtype=np.ulonglong)
    k = 0
    while True:
        result_log = result_log + complex_experiment(R, Nmem, Ntime)
        write_out(result_log, R)
        
        
from sys import argv

script, R, mem, step = argv
main_func(int(R), int(mem), int(step))