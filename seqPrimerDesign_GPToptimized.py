import os
import re
from typing import List, Tuple

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

"""Basic functions"""

def calculate_GC_content(sequence: str) -> float:
    """计算DNA序列的GC含量并返回百分比。"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return round(gc_content, 1)

def calculate_tm(primer: str, mode: str = 'idt_default') -> float:
    """计算引物的Tm值，根据指定的模式返回Tm值。"""
    myseq = Seq(primer.upper())
    if mode == 'idt_default':
        Tm_value = mt.Tm_NN(myseq, dnac1=200, Na=44.2, Mg=0, dNTPs=0, saltcorr=3)
    elif mode == 'KOD544':
        Tm_value = mt.Tm_NN(myseq, dnac1=300, Na=50, Mg=2, dNTPs=0.16, saltcorr=3)
    else:
        raise ValueError('Error: wrong mode')
    return round(Tm_value, 1)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq.upper()))

"""seqPrimerDesign"""

def generate_regions(sequence: str, start: int = 31, length: int = 100, interval: int = 700) -> List[str]:
    """生成给定起始位点、长度和间隔的子序列列表。"""
    region_numbers = (len(sequence) - start) // interval + 1
    regions = [sequence[start + i * interval:start + i * interval + length] for i in range(region_numbers)]
    print(f"Region numbers: {len(regions)}")
    return regions

def find_primers(region: str, end: str = 'GCgc', length_range: Tuple[int, int] = (12, 30)) -> List[str]:
    """在给定区域中找到末尾为GC的，引物长度在指定范围内的引物。"""
    primers = []
    region_length = len(region)
    for i in range(-1, -region_length, -1):
        if region[i] not in end:
            continue
        for j in range(i - length_range[0], i - length_range[1] - 1, -1):
            if abs(j + 1) <= region_length:
                primer = region[j + 1:i] + region[i]
                primers.append(primer)
    return primers

def uni_binding_site(sequence: str, primer: str, overhang: int = 8) -> bool:
    """检查引物在目标序列中的结合位点是否唯一。"""
    pattern = re.compile(re.escape(primer[-overhang:])) 
    count_sense = len(pattern.findall(sequence))
    count_antisense = len(pattern.findall(reverse_complement(sequence)))
    count = count_sense + count_antisense

    return count == 1

def primer_sorter(filtered_primer: Tuple[str, float, float]) -> int:
    """根据引物的GC含量和Tm值排序。"""
    gc, Tm = filtered_primer[1], filtered_primer[2]
    if any(para in (False, None) for para in (gc, Tm)):
        print("Wrong primer occurs, sort_value = -1")
        return -1
    return round(1000 * Tm + gc)

def parameter_filtering(primer: str, gc_range: Tuple[int, int] = (30, 60), tm_range: Tuple[int, int] = (50, 65), gc_end: int = 2) -> Tuple[bool, Tuple[str, float, float]]:
    """根据GC含量、Tm值和末端GC数量筛选引物。"""
    if all(n.lower() in 'gc' for n in primer[-gc_end - 1:]):
        return False, (primer, None, None)
    gc = calculate_GC_content(primer)
    if not (gc_range[0] <= gc <= gc_range[1]):
        return False, (primer, gc, None)
    Tm = calculate_tm(primer, mode='idt_default') # KOD544 is not supported here.
    if not (tm_range[0] <= Tm <= tm_range[1]):
        return False, (primer, gc, Tm)
    return True, (primer, gc, Tm)

def get_seqPrimers(target_sequence: str, region_para: Tuple[int, int, int] = (31, 100, 700), overhang: int = 8, find_para: Tuple[str, Tuple[int, int]] = ('GCgc', (12, 30)), gc_range: Tuple[int, int] = (53, 60), tm_range: Tuple[int, int] = (40, 65), gc_end: int = 2) -> List[Tuple[str, List[Tuple[str, float, float]]]]:
    """设计测序引物，返回符合要求的引物候选列表。"""
    s, l, i = region_para
    regions = generate_regions(target_sequence, start=s, length=l, interval=i)
    PrimerCandidates = []
    e, l_range = find_para
    for idx, region in enumerate(regions):
        primers = find_primers(region, end=e, length_range=l_range)
        print(f"{len(primers)} primers are generated in Region_{idx}.")
        filtered_primers = []
        for primer in primers:
            if uni_binding_site(target_sequence, primer, overhang) and parameter_filtering(primer, gc_range, tm_range, gc_end)[0]:
                filtered_primers.append(parameter_filtering(primer, gc_range, tm_range, gc_end)[1])
        print(f"{len(filtered_primers)} primers match requirements in Region_{idx}.")
        sorted_primers = sorted(filtered_primers, key=primer_sorter, reverse=True)
        PrimerCandidates.append((region, sorted_primers))
    return PrimerCandidates

def seqPrimer_writer(PrimerCandidates: List[Tuple[str, List[Tuple[str, float, float]]]], output: str = 'PrimerCandidates.txt') -> None:
    """将设计的引物候选写入文件。"""
    if os.path.exists(output):
        go = input("文件已存在，是否覆盖? [y/n]]：")
        if go.lower() not in ['y', 'yes']:
            print("取消")
            return
    with open(output, 'w', encoding='utf-8') as file:
        file.write("Name\tSequence: 5'-->3'\tTm/C\tGC content/%\n")
        for i, (region, primers) in enumerate(PrimerCandidates):
            file.write(f"Region_{i}:\t{region}\tlength:{len(region)} nt\t{calculate_GC_content(region)}\n")
            for j, (primer, gc, Tm) in enumerate(primers):
                file.write(f"primer_{i}{j}:\t{primer}\t{Tm} C\t{gc}\n")
            file.write('\n')
    print("文件已写入")

def check_primer(primer: str, target_sequence: str = '', overhang: int = 8, gc_range: Tuple[int, int] = (53, 60), tm_range: Tuple[int, int] = (40, 65)) -> None:
    """检查单个引物是否符合GC含量和Tm值的范围，并打印结果。"""
    if target_sequence:
        unique = uni_binding_site(target_sequence, primer, overhang)
        print(f"Binding sites: {unique}")
    gc = calculate_GC_content(primer)
    print(f"GC/%: {gc}\nIn GC range?: {gc_range[0] <= gc <= gc_range[1]}")
    Tm = calculate_tm(primer)
    print(f"Tm/C: {Tm}\nIn Tm range?: {tm_range[0] <= Tm <= tm_range[1]}")

if __name__ == "__main__":
    """
下面的参数在0507的测序中表现得很好，very nice！
也有可能是因为测序用得质粒是大提试剂盒，质量较好！
感觉间隔可以用800 bp了
"""
    # 主程序，读取目标序列并设计引物
    with open('seq_target_sequence.txt', 'r', encoding='utf-8') as f:
        Target = f.read().strip()
        seqPrimers = get_seqPrimers(Target,
                                    region_para=(0, 80, 600),
                                    overhang=6, # 6
                                    find_para=('GCgc', (18, 24)), # 金唯智参数 ('GCgc', (18, 24))
                                    gc_range=(45, 55), # 金唯智参数 (45, 55)
                                    tm_range=(50, 60), # 金唯智参数 Warning-Tm_mode is idt_default. KOD544 is unsupported.
                                    gc_end=1)   # (50, 60)
        seqPrimer_writer(seqPrimers, output='seqPrimers.txt')
