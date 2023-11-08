import time
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
from scipy.special import comb

Accuracy = 1e-12
k_max = 200


def calc_XYA(ave_k, p_k, m):
    x2 = 1.0
    y2 = 1.0
    delta = 1.0
    while delta > Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1 - math.exp(- ave_k * y1)
        dxy = 0.0
        for k in p_k:
            temp_2 = sum(comb(k,s) * (1-x1)**(k-s) * (x1**s - (x1-y1)**s) for s in range(min(m,k+1)))
            dxy += p_k[k] * temp_2
        y2 = x2 - dxy
        delta = abs(x1 - x2) + abs(y1 - y2)
    A = x2
    return A, x2, y2


def gen_possion_distribution(ave_k, k_max):
    p_k = {}
    p_k[0] = math.exp(-ave_k)
    for kk in range(1, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

def calc_PINF(y2,p_k):
    p_inf =0
    for kin in range(k_max + 1):
        for kout in range(k_max + 1):
           p_2d = p_k[kin] * p_k[kout]
           p_inf += p_2d * (1-(1 - y2) ** kin)
    return p_inf


def calc_gOUT(avg_k,m):
    ## avg_k 是平均的度
    ## m 是渗流指数

    start_time = time.time()
    comb = math.comb
    pk = gen_possion_distribution(avg_k,k_max)
    A,x,y2 = calc_XYA(avg_k,pk,m)
    P_GOUT = calc_PINF(y2,pk)
    return A, P_GOUT

def test_induced():
    for mm in range(1,7):
        for average_k in np.arange(1, 10, 0.2):
            average_k = round(average_k,4)
            index_f = str(mm)+str('_')+str(average_k)+str('.txt')
            out_file_2 = open('data/GOUT_'+index_f, 'w')
            p_k = gen_possion_distribution(average_k, k_max)
            A,x2,y2 = calc_XYA(average_k, p_k, mm)
            P_GOUT = calc_PINF(y2,p_k)
            print(x2,y2)
            if P_GOUT>0.999 and average_k>10.0:
                print(P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k, P_GOUT))
                break
            else :
                print(average_k, P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k, P_GOUT))


if __name__ == "__main__":
    test_induced()