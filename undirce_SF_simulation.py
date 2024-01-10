import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import random


#测试网络模型
def generate_undirected_scale_free_network(N, gamma, average_degree):

    E=int(average_degree*N/2)

    # Calculate alpha
    alpha = 1 / (gamma - 1)

    # Compute in and out probabilities for each node
    nodes = np.arange(1, N + 1)
    probs = nodes ** (-alpha)
    probs /= probs.sum()

    # Generate edges and add to the graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    edges = set()
    while len(edges) < E:
        i = np.random.choice(nodes, p=probs)
        j = np.random.choice(nodes, p=probs)
        if i != j and (i, j) not in edges and (j, i) not in edges:
            G.add_edge(i, j)
            edges.add((i, j))

    return G


def m_index(g, g_stat, i):
    m = 0
    j_list = list(g.neighbors(i))
    k_list = []
    m_list = []
    for j in j_list:
        tmp_k = 0
        if g_stat[j] == 1:
            k_list = list(g.neighbors(j))
            for k in k_list:
                if k != i and g_stat[k] == 1:
                    tmp_k += 1
        m_list.append(tmp_k)
    if len(m_list) > 0:
        m = max(m_list)
    else:
        m = 0
    return int(m)


def g_active_num(g_stat):
    num = 0
    for (v, s) in g_stat.items():
        if s == 1:
            num += 1
    return num


def simulation(m, k):
    ##变量区
    n = 100000
    gamma = 2.5 # 度指数
    g = generate_undirected_scale_free_network(n, gamma, k)

    ## 状态初始化
    g_stat = {}
    for v in g.nodes():
        g_stat[v] = 1  # 设置所有节点的初始状态为1

    # 保存原始网络的副本以进行可视化
    original_network = g.copy()
    print("Original network nodes:", original_network.number_of_nodes())

    ## 诱导
    loop = 1
    while loop:
        num1 = g_active_num(g_stat)
        for v in range(1, len(g)+1):
            if g_stat[v] == 1:
                v_m = m_index(g, g_stat, v)
                if v_m < m:
                    g_stat[v] = 0
        num2 = g_active_num(g_stat)
        loop = abs(num2 - num1)

    ## 计算Active Number的数量
    A = g_active_num(g_stat)

    for v in range(1, n+1):
        if g_stat[v] == 0:
            g.remove_node(v)

    print("Induced network nodes after deletion:", g.number_of_nodes())
    print("Active Number:", A)

    g_list = sorted(list(nx.connected_components(g)), key=len, reverse=True)
    if not g_list:  # 检查 g_list 是否非空
        return 0, A

    gCC_subgraph = g.subgraph(g_list[0]).copy()  # 创建最大连通分量的子图

    GOUT = len(gCC_subgraph) / n

    return GOUT, A

# 新增一个函数用于每个进程
def run_simulation_for_m(m):
    index_f = str(m) + '.txt'
    out_file_1 = open('un_SF_data/simu_' + index_f, 'w')
    l_GOUT = []
    for average_k in np.arange(0, 10.00001, 0.1):
        l_GOUT = []
        for i in range(1):
            print("k:" + str(average_k))
            GOUT = simulation(m, average_k)[0]
            print("GOUT:" + str(GOUT))
            l_GOUT.append(GOUT)
        sum = 0
        for ii in l_GOUT:
            sum += ii
        GOUT_AVG = sum / 1
        print("average_k:" + str(average_k))
        print("GOUT:" + str(GOUT_AVG))
        if average_k < 10:
            out_file_1.write('%f\t%f\n' % (average_k, GOUT_AVG))
        else:
            out_file_1.write('%f\t%f' % (average_k, GOUT_AVG))


def test_simulation():
    processes = []
    for m in range(1, 7):
        p = multiprocessing.Process(target=run_simulation_for_m, args=(m,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()



if __name__ == "__main__":
    test_simulation()
