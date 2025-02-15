import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

def conut_GC_perBasePair(sequence, window_size=15, start=0):
    
    gc_content_list = np.zeros(len(sequence))
    # 计算每个窗口的GC含量
    for i in range(start, len(sequence), window_size):
        # 获取当前窗口的序列
        window_seq = sequence[i:i + window_size].upper()
        # 计算GC含量
        gc_content = (window_seq.count("G")+window_seq.count("C"))/len(window_seq)
        # 添加到列表
        gc_content_list[i:i+window_size] = gc_content

    
    return gc_content_list
    

def draw_gc_plot(fasta_file, window_size=20):
    
    # For sequence range, change line 32 to your target
    # Eg. sequence = record.seq[150: 320]


    # 读取FASTA文件
    sequences = SeqIO.parse(fasta_file, "fasta")

    # 设置窗口大小
    window_size = window_size
    for record in sequences:
        
        sequence = record.seq
        
        # 计算以不同起始位点开始的窗口的GC值
        gc = list()
        for i in range(window_size):
            gc_content_list = conut_GC_perBasePair(sequence, window_size=window_size, start=i)
            gc.append(gc_content_list)
        
        # 求每个碱基处GC的平均值
        matrix = np.array(gc)
        print(matrix.shape)
        normalized_gc = np.mean(matrix, axis=0)
        
        # 显示高gc区域 TODO
        for i in range(len(normalized_gc)):
            if normalized_gc[i] > 0.58:
                print(f"High gc region: {i}-{i+window_size}")
        
        # 准备color map
        cmap = plt.get_cmap('rainbow') # 创建一个颜色映射，这里使用'viridis'作为例子
        norm = plt.Normalize(min(normalized_gc), max(normalized_gc)) # 创建一个归一化对象，将数据值映射到[0, 1]区间
        colors = cmap(norm(normalized_gc)) # 将数据值映射到颜色映射的颜色空间
        
        # 绘制直方图
        
        plt.bar(range(len(sequence)), normalized_gc*100, width=2, color=colors)

        # 添加标题和轴标签
        plt.title(f'GC Content Distribution Across {record.id}')
        plt.xlabel('Basepair Number')
        plt.ylabel('GC Content %')
        
        # 更好地展示数据
        plt.ylim(0,100)
        # 显示图例
        #plt.legend([f''])

        # 显示图形
        plt.show()

if __name__ == "__main__":
    draw_gc_plot(fasta_file = 'seeGC.fasta',window_size=20)
