import networkx as nx
from scipy.special import comb
import numpy as np
import os

Accuracy = 1e-12
gamma = 2.5
k_max = 1000

# 幂律分布
def gen_power_law_distribution(k_max,gamma):
    p_k = {}
    # 归一化
    norm = sum([k ** (-gamma) for k in range(1, k_max + 1)])
    for k in range(1, k_max + 1):
        p_k[k] = (k ** (-gamma)) / norm
    return p_k

# 计算公式
def calculate_tilde_v(m, p_k, avg_k):
    tilde_v = sum(k * p_k.get(k, 0) / avg_k for k in range(m + 1, len(p_k)))
    return tilde_v

def calculate_tilde_v_inf(m, p_k, tilde_t_inf, tilde_v_inf, avg_k):
    tilde_v_inf = sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_t_inf - tilde_v_inf) ** (k - 1)) for k in range(m + 1, len(p_k)))
    return tilde_v_inf

def calculate_tilde_t_inf(m, p_k, tilde_a_inf, tilde_y_inf, avg_k):
    tilde_t_inf = sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_a_inf - tilde_y_inf) ** (k - 1)) for k in range(1, m + 1))
    return tilde_t_inf

def calculate_tilde_y(m, p_k, tilde_v, avg_k):
    tilde_y = sum(k * p_k.get(k, 0) / avg_k * (1 - (1 - tilde_v) ** (k - 1)) for k in range(m + 1, len(p_k)))
    return tilde_y

def calculate_tilde_a_inf(m, p_k, tilde_y, tilde_y_inf, tilde_a_inf, avg_k):
    tilde_a_inf = sum(
        k * p_k.get(k, 0) / avg_k * sum(
            comb(k - 1, s) * tilde_y ** s * (1 - tilde_y) ** (k - 1 - s) *
            (1 - (1 - tilde_y_inf / tilde_y) ** s * (1 - tilde_a_inf / (1 - tilde_y)) ** (k - 1 - s))
            for s in range(1, k)
        ) for k in range(1, m + 1))
    return tilde_a_inf

def calculate_tilde_y_inf(m, p_k, tilde_v, tilde_v_inf, tilde_t_inf, avg_k):
    tilde_y_inf = sum(
        k * p_k.get(k, 0)/ avg_k * sum(
            comb(k - 1, s) * tilde_v ** s * (1 - tilde_v) ** (k - 1 - s) *
            (1 - (1 - tilde_v_inf / tilde_v) ** s * (1 - tilde_t_inf / (1 - tilde_v)) ** (k - 1 - s))
            for s in range(1, k)
        ) for k in range(m + 1, len(p_k)))
    return tilde_y_inf


# 迭代
def calculate_P_inf(m, p_k, avg_k):
    tilde_v_val = tilde_y_val = tilde_v_inf_val = tilde_t_inf_val = tilde_a_inf_val = tilde_y_inf_val = 0.5

    while True:
    #for _ in range(10000):
        tilde_v_val_new = calculate_tilde_v(m, p_k, avg_k)
        tilde_y_val_new = calculate_tilde_y(m, p_k, tilde_v_val, avg_k)
        tilde_v_inf_val_new = calculate_tilde_v_inf(m, p_k, tilde_t_inf_val, tilde_v_inf_val, avg_k)
        tilde_t_inf_val_new = calculate_tilde_t_inf(m, p_k, tilde_a_inf_val, tilde_y_inf_val, avg_k)
        tilde_a_inf_val_new = calculate_tilde_a_inf(m, p_k, tilde_y_val, tilde_y_inf_val, tilde_a_inf_val, avg_k)
        tilde_y_inf_val_new = calculate_tilde_y_inf(m, p_k, tilde_v_val, tilde_v_inf_val, tilde_t_inf_val, avg_k)

        if all(abs(new_val - old_val) < Accuracy for new_val, old_val in
               [(tilde_v_val_new, tilde_v_val), (tilde_y_val_new, tilde_y_val),
                (tilde_v_inf_val_new, tilde_v_inf_val), (tilde_t_inf_val_new, tilde_t_inf_val),
                (tilde_a_inf_val_new, tilde_a_inf_val), (tilde_y_inf_val_new, tilde_y_inf_val)]):
            break

        tilde_v_val, tilde_y_val, tilde_v_inf_val, tilde_t_inf_val, tilde_a_inf_val, tilde_y_inf_val = (
            tilde_v_val_new, tilde_y_val_new, tilde_v_inf_val_new, tilde_t_inf_val_new,
            tilde_a_inf_val_new, tilde_y_inf_val_new
        )

        # 计算最终 P_inf
    P_inf = sum(p_k.get(k, 0) * (1 - (1 - tilde_y_val) ** k - (1 - tilde_y_inf_val - tilde_a_inf_val) ** k + (
                1 - tilde_y_val - tilde_a_inf_val) ** k) for k in range(0, m + 1))
    P_inf += sum(p_k.get(k, 0) * (
                1 - (1 - tilde_v_val) ** k - (1 - tilde_t_inf_val - tilde_v_inf_val) ** k + (
                    1 - tilde_v_val - tilde_t_inf_val) ** k) for k in range(m + 1, len(p_k)))

    return P_inf

def run_calculations_and_save():
    output_folder = 'un_SF_data_n'
    os.makedirs(output_folder, exist_ok=True)  # 确保文件夹存在

    for average_k in np.arange(0.1, 10, 0.2):
        average_k = round(average_k, 4)
        for m in range(1, 7):
            file_name = f'{output_folder}/GOUT_{m}_{average_k}.txt'
            with open(file_name, 'w') as out_file:
                p_k = gen_power_law_distribution(k_max, gamma)
                P_inf_value = calculate_P_inf(m, p_k, average_k)
                out_file.write('%f\t%f\n' % (average_k, P_inf_value))
                print(f'average_k: {average_k}, m: {m}, P_inf: {P_inf_value}')

# 运行计算并保存结果
if __name__ == "__main__":
    run_calculations_and_save()

