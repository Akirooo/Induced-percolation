import matplotlib.pyplot as plt
import numpy as np

def plot_txt_files():
    # 设置颜色和标记
    colors = ['blueviolet', 'r', 'darkorange', 'g', 'b', 'black']
    markers = ['o', 's', '^', 'v', 'p', '*']

    plt.figure(figsize=(10, 7))

    # 处理 numerical 数据（用线表示）
    for m in range(1, 7):
        average_k_values = []
        GOUT_values = []
        for average_k in np.arange(0, 10, 0.1):
            average_k = round(average_k, 4)
            filename = 'data_n/GOUT_' + str(m) + '_' + str(average_k) + '.txt'

            with open(filename, 'r') as f:
                for line in f:
                    average_k, GOUT = map(float, line.strip().split('\t'))
                    average_k_values.append(average_k)
                    GOUT_values.append(GOUT)

        plt.plot(average_k_values, GOUT_values, label=f'm={m} numerical', color=colors[m - 1])

    # 处理 simulation 数据（用点表示）
    added_legend = set()
    for m in range(1, 7):
        filename = 'data/simu_' + str(m) + '.txt'

        with open(filename, 'r') as f:
            for line in f:
                average_k, GOUT = map(float, line.strip().split('\t'))
                label = f'm={m} simulation' if m not in added_legend else ""
                plt.scatter(average_k, GOUT, color=colors[m - 1], marker=markers[m - 1], label=label)
                added_legend.add(m)

    plt.xlabel('average_k')
    plt.ylabel('GOUT')
    plt.legend()
    plt.title('GOUT vs. average_k for different m values')
    plt.grid(True)
    plt.xlim(0, max(average_k_values) + 1)
    plt.tight_layout()
    plt.savefig('plots/GOUT_vs_average_k_combine.png')
    plt.show()

if __name__ == "__main__":
    plot_txt_files()
