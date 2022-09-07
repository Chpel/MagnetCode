import numpy as np
from numba import jit, njit, prange
from time import time



@njit
def create_walk(N):
    x0 = np.array([0,0])
    steps = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    walk_directions = np.random.randint(0,4, size=N)
    walk_dots = np.zeros((N, 2), dtype=np.int16)
    walk_dots[0] = x0
    for i in range(1,N):
        walk_dots[i] = walk_dots[i-1] + steps[walk_directions[i]]
    return walk_dots
    
@njit
def calc_fractions(dots, N, last_dot):
    steps = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    N_new = dots.shape[0]
    walk_neighbors = 0
    neigh_fract_0 = np.zeros(4, dtype=np.float_)
    for i in range(N_new):
        walk_neighbors = 0
        for step in steps:
            potential_n = dots[i] + step
            for j in range(N_new):
                if (dots[j] == potential_n).all():
                    walk_neighbors += 1
        neigh_fract_0[walk_neighbors-1] += 1
        if (dots[i] == last_dot).all():
            atm = 4 - walk_neighbors
    neigh_fract_0 /= N_new
    return np.append(neigh_fract_0, N_new / N), atm

@jit
def experiment(N: int):
    walk = create_walk(N)    
    unique_dots = np.unique(walk, axis=0)
    return calc_fractions(unique_dots,N,walk[-1,:])
    
@jit(parallel=True)
def complex_experiment(N, di):
    n_batch = np.zeros((5,di), dtype=np.float_)
    atm = np.zeros(di, dtype=np.short)
    for i in prange(di):
        n_batch[:,i], atm[i] = experiment(N)
    return n_batch[0, :], n_batch[1,:], n_batch[2,:], n_batch[3,:], n_batch[4,:], atm
    
@njit
def stats(*args):
    n = len(args)
    means = np.zeros(n)
    stds = np.zeros(n)
    
    i=0
    for a in args:
        means[i], stds[i] = a.mean(), a.std()
        i+=1
    return means, stds
    
@njit
def atm_bins(arr):
    atms = np.array([0,0,0,0], np.int_)
    for i in range(4):
        atms[i] = (arr == i).sum()
    return atms
    
    
def main_func(N, di, stop_i):
    n1 = np.array([])
    n2 = np.array([])
    n3 = np.array([])
    n4 = np.array([])
    nU = np.array([])
    atm = np.array([0,0,0,0], dtype=np.int_)
    means = np.array([])
    stds = np.array([])

    iters = 0
    while iters < stop_i:
        iters += 1
        n1_n, n2_n, n3_n, n4_n, nU_n, atm_n = complex_experiment(N, di)   
        n1 = np.append(n1, n1_n)
        n2 = np.append(n2, n2_n)
        n3 = np.append(n3, n3_n)
        n4 = np.append(n4, n4_n)
        nU = np.append(nU, nU_n)
        atm = atm + atm_bins(atm_n)
        obs_mean, obs_std = stats(n1,n2,n3,n4,nU)
        p_a = atm / (iters * di)
        if iters == 1:
            means = np.array([obs_mean])
            stds = np.array([obs_std])
        else:
            means = np.append(means, [obs_mean], axis=0)
            stds = np.append(stds, [obs_std], axis=0)

        write_results(N, obs_mean, obs_std, p_a, iters * di)
        
        if iters % 10 == 0:
            save_distr(N, iters * di, n1, n2, n3, n4, nU)
            save_history(N, means, stds, iters, di)
        
        
def write_results(N, obs_mean, obs_std, p, steps):
    res_array = np.array([N])
    res_array = np.append(res_array, obs_mean)
    res_array = np.append(res_array, obs_std)
    res_array = np.append(res_array, p)
    res_array = np.append(res_array, steps)
    np.savetxt('Drunk_Sailor_N_' + str(N) + '.txt', 
               [res_array], 
               delimiter=' ', 
               fmt=['%d', 
                    '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f',
                    '%1.6f', '%1.6f', '%1.6f', '%1.6f', '%1.6f',
                    '%1.6f', '%1.6f', '%1.6f', '%1.6f',
                    '%d'], 
               newline=' ', 
               header='N ' + 
               ' n1_mean  n2_mean  n3_mean  n4_mean  uni_mean ' + 
               ' n1_std  n2_std  n3_std  n4_std  uni_std ' + 
               ' p0  p1  p2  p3 ' +
               ' steps\n',
               comments='')

def save_distr(N, steps, *arrs):
    name = 'DS_' + str(N) + '_dists.npz'
    bs = int(max(10, 5 + np.log2(steps)))
    counts_arr = []
    bins_arr = []
    for x in arrs:
        counts, bins = np.histogram(x, bs)
        counts_arr.append(counts)
        bins_arr.append(bins)
    counts_arr = np.array(counts_arr)
    bins_arr = np.array(bins_arr)
    np.savez(name, hist=counts_arr, bins=bins_arr, steps=steps)
    
    
def save_history(N, m, s, step, i):
    name = 'DS_' + str(N) + '_history.npz'
    np.savez(name, means=m, stds=s, step=step, iters=i)
        


        

from sys import argv

script, N, K, K_stop = argv
main_func(int(N), int(K), int(K_stop))