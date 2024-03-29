import numpy as np
from numba import jit



@jit
def experiment(N: int):
    x0 = np.array([0,0])
    steps = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    walk_directions = np.random.randint(0,4, size=N)
    walk_dots = np.zeros((N, 2), dtype=np.int16)
    walk_dots[0] = x0
    for i in range(1,N):
        walk_dots[i] = walk_dots[i-1] + steps[walk_directions[i]]
        
    walk_dots_new = np.unique(walk_dots, axis=0)
    N_new = len(walk_dots_new)

    walk_neighbors = 0
    neigh_fract_0 = np.zeros((1,5), dtype=np.float_)
    for i in range(N_new):
        walk_neighbors = 0
        for step in steps:
            potential_n = walk_dots_new[i] + step
            for j in range(N_new):
                if (walk_dots_new[j] == potential_n).all():
                    walk_neighbors += 1
        neigh_fract_0[0][walk_neighbors-1] += 1
    neigh_fract_0 /= N_new
    neigh_fract_0[0][-1] = N_new / N
    return neigh_fract_0

def write_results(N, obs_mean, obs_std, steps):
    res_array = np.array([N])
    res_array = np.append(res_array, obs_mean)
    res_array = np.append(res_array, obs_std)
    res_array = np.append(res_array, steps)
    np.savetxt('Drunk_Sailor_N_' + str(N) + '.txt', 
               [res_array], 
               delimiter=' ', 
               fmt=['%d', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%d'], 
               newline=' ', 
               header='N  n1_mean  n2_mean  n3_mean  n4_mean  uni_mean  n1_std  n2_std  n3_std  n4_std  uni_std  steps\n',
               comments='')

def save_distr(N, obs, steps):
    name = 'DS_' + str(N) + '_dists.npz'
    bs = int(max(10, 5 + np.log2(steps)))
    counts_arr = []
    bins_arr = []
    for i in range(5):
        x = obs[:, i]
        counts, bins = np.histogram(x, bs)
        counts_arr.append(counts)
        bins_arr.append(bins)
    counts_arr = np.array(counts_arr)
    bins_arr = np.array(bins_arr)
    np.savez(name, hist=counts_arr, bins=bins_arr, steps=steps)
    
    
def save_history(N, m, s, step, i):
    name = 'DS_' + str(N) + '_history.npz'
    np.savez(name, means=m, stds=s, step=step, iters=i)
        

def complex_experiment(N, step_i, stop_i):
    observables = experiment(N)
    means = np.array([])
    stds = np.array([])
    
    iters = 0
    
    while True:
        for i in range(step_i):
            observables = np.append(observables, experiment(N), axis=0)
        iters += 1
        obs_mean = observables.mean(axis=0)
        obs_std = observables.std(axis=0)
        if iters == 1:
            means = np.array([obs_mean])
            stds = np.array([obs_std])
        else:
            means = np.append(means, [obs_mean], axis=0)
            stds = np.append(stds, [obs_std], axis=0)
        
        save_distr(N, observables, iters * step_i)
        write_results(N, obs_mean, obs_std, iters * step_i)
        save_history(N, means, stds, iters, step_i)
        if iters >= stop_i:
            break
        

from sys import argv

script, N, K, K_stop = argv
complex_experiment(int(N), int(K), int(K_stop))