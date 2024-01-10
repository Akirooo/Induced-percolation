import numpy as np
import math
from scipy.special import comb
from decimal import Decimal, getcontext

Accuracy = 1e-12
gamma_in = 2.5
gamma_out = 3.0
N=100000

def k_max(N,gamma):
    k_max =int(N ** (1 / (gamma - 1)))
    return k_max


def calc_XYA(ave_k, p_k, m):
    x2 = 0.5
    y2 = 0.5
    delta = 1.0
    while delta > Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1 - math.exp(-ave_k * y1)
        dxy = 0.0
        for k in p_k:
            temp_2 = sum(comb(k, s) * (1 - x1) ** (k - s) * (x1 ** s - (x1 - y1) ** s) for s in range(min(m, k + 1)))
            dxy += p_k[k] * temp_2
        y2 = x2 - dxy
        delta = abs(x1 - x2) + abs(y1 - y2)
    A = x2
    return A, x2, y2



def gen_power_law_distribution(gamma, k_min, k_max):
    p_k = {}
    for k in range(k_min, k_max + 1):
        p_k[k] = k ** (-gamma)
    # 归一化
    total = sum(p_k.values())
    for k in p_k:
        p_k[k] /= total
    return p_k

def calc_PINF(y2,pk_in,pk_out,k_max_in,k_max_out):
    p_inf =0
    for kin in range(1, k_max_in + 1):
        for kout in range(1, k_max_out + 1):
           p_2d = pk_in[kin] * pk_out[kout]
           p_inf += p_2d * (1-(1 - y2) ** kin)
    return p_inf


def calc_gOUT(avg_k,m):
    ## avg_k 是平均的度
    ## m 是渗流指数
    k_max_in = k_max(N, gamma_in)
    k_max_out = k_max(N, gamma_out)
    # 生成入度的幂律分布
    pk_in = gen_power_law_distribution(gamma_in, k_max_in)
    # 生成出度的幂律分布
    pk_out = gen_power_law_distribution(gamma_out, k_max_out)

    A,x,y = calc_XYA(avg_k,pk_in,m)
    P_GOUT = calc_PINF(y,pk_in,pk_out, k_max_in, k_max_out)
    return A, P_GOUT,x,y

def test_induced():
    for mm in range(3,4):
        for average_k in np.arange(0, 10, 0.1):
            average_k = round(average_k,4)
            index_f = str(mm)+str('_')+str(average_k)+str('.txt')
            out_file_2 = open('un_SF_data_n/GOUT_'+index_f, 'w')
            k_max_in = k_max(N, gamma_in)
            k_max_out = k_max(N, gamma_out)
            # 生成入度的幂律分布
            pk_in = gen_power_law_distribution(gamma_in, k_max_in)
            # 生成出度的幂律分布
            pk_out = gen_power_law_distribution(gamma_out, k_max_out)
            A,x,y = calc_XYA(average_k, pk_in, mm)
            P_GOUT = calc_PINF(y,pk_in,pk_out,k_max_in,k_max_out)
            if P_GOUT>0.999 and average_k>10.0:
                print(P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k, P_GOUT))
                break
            else :
                print(average_k, P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k, P_GOUT))


if __name__ == "__main__":
    test_induced()