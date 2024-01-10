import math
from scipy.special import comb
import numpy as np
import os

Accuracy = 1e-12

# 泊松分布
def gen_poisson_distribution(ave_k, k_max):
    p_k = {}
    p_k[0] = math.exp(-ave_k)
    for kk in range(1, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

# 计算公式
def tilde_v(m, p_k, avg_k):
    return sum(k * p_k.get(k, 0) / avg_k for k in range(m + 1, len(p_k)))

def tilde_v_infinity(m, p_k, tilde_t_infinity, tilde_v_infinity, avg_k):
    return sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_t_infinity - tilde_v_infinity) ** (k - 1)) for k in range(m + 1, len(p_k)))

def tilde_t_infinity(m, p_k, tilde_a_infinity, tilde_y_infinity, avg_k):
    return sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_a_infinity - tilde_y_infinity) ** (k - 1)) for k in range(1, m + 1))

def tilde_y(m, p_k, tilde_v, avg_k):
    return sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_v) ** (k - 1)) for k in range(m + 1, len(p_k)))

def tilde_a_infinity(m, p_k, tilde_y, tilde_y_infinity, tilde_a_infinity, avg_k):
    return sum(
        k * p_k.get(k, 0) / avg_k * sum(
            comb(k - 1, s) * tilde_y ** s * (1 - tilde_y) ** (k - 1 - s) *
            (1 - (1 - tilde_y_infinity / tilde_y) ** s * (1 - tilde_a_infinity / (1 - tilde_y)) ** (k - 1 - s))
            for s in range(1, k)
        ) for k in range(1, m + 1)
    )

def tilde_y_infinity(m, p_k, tilde_v, tilde_v_infinity, tilde_t_infinity, avg_k):
    return sum(
        k * p_k.get(k, 0) / avg_k * sum(
            comb(k - 1, s) * tilde_v ** s * (1 - tilde_v) ** (k - 1 - s) *
            (1 - (1 - tilde_v_infinity / tilde_v) ** s * (1 - tilde_t_infinity / (1 - tilde_v)) ** (k - 1 - s))
            for s in range(1, k)
        ) for k in range(m + 1, len(p_k))
    )

# 迭代
def P_infinity(m, p_k, avg_k):
    tilde_v_val = tilde_y_val = tilde_v_infinity_val = tilde_t_infinity_val = tilde_a_infinity_val = tilde_y_infinity_val = 0.5

    for _ in range(1000):  #迭代次数
        tilde_v_val_new = tilde_v(m, p_k, avg_k)
        tilde_y_val_new = tilde_y(m, p_k, tilde_v_val, avg_k)
        tilde_v_infinity_val_new = tilde_v_infinity(m, p_k, tilde_t_infinity_val, tilde_v_infinity_val, avg_k)
        tilde_t_infinity_val_new = tilde_t_infinity(m, p_k, tilde_a_infinity_val, tilde_y_infinity_val, avg_k)
        tilde_a_infinity_val_new = tilde_a_infinity(m, p_k, tilde_y_val, tilde_y_infinity_val, tilde_a_infinity_val, avg_k)
        tilde_y_infinity_val_new = tilde_y_infinity(m, p_k, tilde_v_val, tilde_v_infinity_val, tilde_t_infinity_val, avg_k)

        if all(abs(new_val - old_val) < Accuracy for new_val, old_val in
               [(tilde_v_val_new, tilde_v_val), (tilde_y_val_new, tilde_y_val),
                (tilde_v_infinity_val_new, tilde_v_infinity_val), (tilde_t_infinity_val_new, tilde_t_infinity_val),
                (tilde_a_infinity_val_new, tilde_a_infinity_val), (tilde_y_infinity_val_new, tilde_y_infinity_val)]):
            break

        tilde_v_val, tilde_y_val, tilde_v_infinity_val, tilde_t_infinity_val, tilde_a_infinity_val, tilde_y_infinity_val = (
            tilde_v_val_new, tilde_y_val_new, tilde_v_infinity_val_new, tilde_t_infinity_val_new,
            tilde_a_infinity_val_new, tilde_y_infinity_val_new
        )

    P_inf = sum(p_k.get(k, 0) * (1 - (1 - tilde_y_val) ** k - (1 - tilde_y_infinity_val - tilde_a_infinity_val) ** k + (
                1 - tilde_y_val - tilde_a_infinity_val) ** k) for k in range(0, m + 1))
    P_inf += sum(p_k.get(k, 0) * (
                1 - (1 - tilde_v_val) ** k - (1 - tilde_t_infinity_val - tilde_v_infinity_val) ** k + (
                    1 - tilde_v_val - tilde_t_infinity_val) ** k) for k in range(m + 1, len(p_k)))

    return P_inf

def run_calculations_and_save():
    k_max = 200  # 最大度数
    output_folder = 'test'
    os.makedirs(output_folder, exist_ok=True)  

    for average_k in np.arange(0, 10, 0.1):
        average_k = round(average_k, 4)
        for m in range(1, 7):
            file_name = f'{output_folder}/{m}_{average_k}.txt'
            with open(file_name, 'w') as out_file:
                p_k = gen_poisson_distribution(average_k, k_max)
                P_inf_value = P_infinity(m, p_k, average_k)
                out_file.write('%f\t%f\n' % (average_k, P_inf_value))
                print(f'average_k: {average_k}, m: {m}, P_inf: {P_inf_value}')

# 运行计算并保存结果
if __name__ == "__main__":
    run_calculations_and_save()

