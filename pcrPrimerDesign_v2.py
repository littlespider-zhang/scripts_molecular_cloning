import os
import re
import random
from typing import List, Tuple, Union
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# Basic functions
def generate_random_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of given length.

    Args:
        length (int): Length of the sequence.

    Returns:
        str: Random DNA sequence.
    """
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return sequence

def generate_subsequence(sequence: str, length: int = 8) -> List[str]:
    """
    Generate all subsequences of a given length from a DNA sequence.

    Args:
        sequence (str): DNA sequence.
        length (int): Length of each subsequence.

    Returns:
        List[str]: List of subsequences.
    """
    return [sequence[i:i + length] for i in range(len(sequence) - length + 1)]

def calculate_GC_content(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    Args:
        sequence (str): DNA sequence.

    Returns:
        float: GC content as a percentage.
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return round(gc_content, 1)

def calculate_tm(primer: str, mode: str = 'idt_default') -> float:
    """
    Calculate the melting temperature (Tm) of a DNA primer.

    Args:
        primer (str): DNA primer sequence.
        mode (str): Mode for Tm calculation. Default is 'idt_default'.

    Returns:
        float: Melting temperature.
    """
    myseq = Seq(primer.upper())
    if mode == 'idt_default':
        Tm_value = mt.Tm_NN(myseq, dnac1=200, Na=44.2, Mg=0, dNTPs=0, saltcorr=3)
    elif mode == 'KOD544':
        Tm_value = mt.Tm_NN(myseq, dnac1=300, Na=50, Mg=2, dNTPs=0.16, saltcorr=3)
    else:
        raise ValueError('Error: wrong mode')
    return round(Tm_value, 1)

def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq.upper()))

def self_dimers_2(primers: List[str], min_pairs: int = 4, max_pairs: int = 7) -> List[Tuple[str, str, int, int, str]]:
    """
    Check for self-dimers in a list of primers.

    Args:
        primers (List[str]): List of DNA primers.
        min_pairs (int): Minimum length of self-dimer. Default is 4.
        max_pairs (int): Maximum length of self-dimer. Default is 7.

    Returns:
        List[Tuple[str, str, int, int, str]]: List of self-dimers.
    """
    dimers = [
        (primer, subseq, i, i + length, rev_subseq)
        for primer in primers
        for length in range(min_pairs, min(max_pairs + 1, len(primer)))
        for i in range(len(primer) - length + 1)
        for subseq, rev_subseq in [(primer[i:i + length], reverse_complement(primer[i:i + length]))]
        if rev_subseq in primer
    ]
    return dimers

def generate_3PrimeExtension(sequence: str, length_range: Tuple[int, int] = (3, 8)) -> List[str]:
    """
    Generate 3' end extension sequences.

    Args:
        sequence (str): DNA sequence.
        length_range (Tuple[int, int]): Range of lengths for 3' end extension.

    Returns:
        List[str]: List of 3' end extension sequences.
    """
    PrimeExtension = []
    s, e = length_range
    e = min(e, len(sequence))
    for i in range(-s, -e - 1, -1):
        PrimeExtension.append(sequence[i:])
    return PrimeExtension


# Functions in seqPrimerDesign
def uni_binding_site(template: str, primer: str, overhang: int = 8) -> Union[bool, None]:
    """
    Check if a primer binds uniquely to a template sequence.

    Args:
        template (str): Template DNA sequence.
        primer (str): Primer DNA sequence.
        overhang (int): Length of the overhang region. Default is 8.

    Returns:
        Union[bool, None]: True if unique binding, False if not, None if no binding site.
    """
    pattern = re.compile(re.escape(primer[-overhang:])) 
    count_sense = len(pattern.findall(template))
    count_antisense = len(pattern.findall(reverse_complement(template)))

    count = count_sense + count_antisense
    if count == 1:
        return True
    elif count == 0:
        print("Notice! No binding site for the primer.")
        return None
    else:
        print(f"# of Binding sites: {count}")
        return False

def parameter_filtering(primer: str, GC_range: Tuple[int, int] = (30, 60), Tm_range: Tuple[int, int] = (50, 65), Tm_mode: str = 'idt_default', GC_end: int = 2) -> Tuple[bool, Tuple[str, Union[float, None], Union[float, None]]]:
    """
    Filter a primer based on GC content and Tm range.

    Args:
        primer (str): Primer DNA sequence.
        GC_range (Tuple[int, int]): Acceptable GC content range.
        Tm_range (Tuple[int, int]): Acceptable Tm range.
        Tm_mode (str): Mode for Tm calculation. Default is 'idt_default'.
        GC_end (int): Number of bases at the end of the primer to check for GC content.

    Returns:
        Tuple[bool, Tuple[str, Union[float, None], Union[float, None]]]: Whether primer passes filtering, and primer information.
    """
    all_gc = all(n.lower() in 'gc' for n in primer[-GC_end:])
    if all_gc:
        return False, (primer, None, None)
    
    gc = calculate_GC_content(primer)
    if not (GC_range[0] <= gc <= GC_range[1]):
        return False, (primer, gc, None)
    
    Tm = calculate_tm(primer, mode=Tm_mode)
    if not (Tm_range[0] <= Tm <= Tm_range[1]):
        return False, (primer, gc, Tm)
    
    return True, (primer, gc, Tm)

def dimer_forms_PTJ(primer1: str, primer2: str, max_PTJ: int = 8, min_PTJ: int = 4) -> bool:
    """
    Check if two primers form a primer-template junction (PTJ).

    Args:
        primer1 (str): First primer sequence.
        primer2 (str): Second primer sequence.
        max_PTJ (int): Maximum length of PTJ. Default is 8.
        min_PTJ (int): Minimum length of PTJ. Default is 4.

    Returns:
        bool: True if PTJ forms, False otherwise.
    """
    if min_PTJ == max_PTJ:
        return False
    
    primer2_rc = reverse_complement(primer2)
    ends = generate_3PrimeExtension(primer1, length_range=(min_PTJ, max_PTJ))
    for e in sorted(ends, key=len, reverse=True):
        if e in primer2_rc:
            print(f"3'end anneals with {len(e)} base pairs.")
            print(f"Primer1:{primer1}-->{e}\nPrimer2:{primer2}-->{reverse_complement(e)}\n{e}")
            return True
    
    primer1_rc = reverse_complement(primer1)
    ends = generate_3PrimeExtension(primer2, length_range=(min_PTJ, max_PTJ))
    for e in sorted(ends, key=len, reverse=True):
        if e in primer1_rc:
            print(f"3'end anneals with {len(e)} base pairs.")
            print(f"Primer1:{primer1}-->{reverse_complement(e)}\nPrimer2:{primer2}-->{e}\n")
            return True
    
    return False

# pcrPrimerDesign functions
def generate_pcrPrimers(template: str, length_range: Tuple[int, int] = (17, 30), one_GorC=False) -> Tuple[List[str], List[str]]:
    """
    Generate all possible primers of given length range from a template sequence.

    Args:
        template (str): Template DNA sequence.
        length_range (Tuple[int, int]): Length range for primers.

    Returns:
        Tuple[List[str], List[str]]: Forward and reverse primers.
    """
    sense_strand = template
    antisense_strand = reverse_complement(template)
    primers_F = [sense_strand[:i] for i in range(length_range[0], length_range[1] + 1)]
    primers_R = [antisense_strand[:i] for i in range(length_range[0], length_range[1] + 1)]
    
    if one_GorC:
        primers_F_GorC = []
        for pf in primers_F:
            if pf[-1] in 'GCgc':
                primers_F_GorC.append(pf)
        
        primers_R_GorC = []
        for pr in primers_R:
            if pr[-1] in 'GCgc':
                primers_R_GorC.append(pr)
        
        return primers_F_GorC, primers_R_GorC

    return primers_F, primers_R

def primer_filter(template: str, primers: List[str], overhang: int = 6, GC_range: Tuple[int, int] = (30, 60), Tm_range: Tuple[int, int] = (50, 65), Tm_mode: str = 'KOD544', GC_end: int = 2) -> List[Tuple[str, float, float]]:
    """
    Filter primers based on binding site uniqueness, GC content, and Tm range.

    Args:
        template (str): Template DNA sequence.
        primers (List[str]): List of primer sequences.
        overhang (int): Length of the overhang region. Default is 6.
        GC_range (Tuple[int, int]): Acceptable GC content range. Default is (30, 60).
        Tm_range (Tuple[int, int]): Acceptable Tm range. Default is (50, 65).
        Tm_mode (str): Mode for Tm calculation. Default is 'KOD544'.
        GC_end (int): Number of bases at the end of the primer to check for GC content.

    Returns:
        List[Tuple[str, float, float]]: List of filtered primers with their GC content and Tm values.
    """
    filtered_primers = []
    filtered_uni = 0
    filtered_fullfil = 0
    for primer in primers:
        unique = uni_binding_site(template, primer, overhang=overhang)
        fullfil, primer_info = parameter_filtering(primer, GC_range=GC_range, Tm_range=Tm_range, Tm_mode=Tm_mode, GC_end=GC_end)
        if unique and fullfil:
            filtered_primers.append(primer_info)
        if not unique:
            filtered_uni += 1
        if not fullfil:
            filtered_fullfil += 1
    print(f"{filtered_uni} primers have multiple binding sites.\n{filtered_fullfil} primers do not fullfil requirements.")
    return filtered_primers

def primerpairs_sorter(primerpairs: Tuple[Tuple[str, float, float], Tuple[str, float, float]]) -> float:
    """
    Sort primer pairs based on the difference in their Tm values.

    Args:
        primerpairs (Tuple[Tuple[str, float, float], Tuple[str, float, float]]): Primer pair information.

    Returns:
        float: Difference in Tm values between the primer pair.
    """
    pf, pf_gc, pf_Tm = primerpairs[0]
    pr, pr_gc, pr_Tm = primerpairs[1]
    return abs(pf_Tm - pr_Tm)*1000 + min(pf_gc, pr_gc)

def get_prcPrimerPairs(sense_strand: str, backbone, length_range: Tuple[int, int] = (17, 30), mode: str = 'KOD544', overhang: int = 6, GC_range: Tuple[int, int] = (30, 60), Tm_range: Tuple[int, int] = (50, 65), GC_end: int = 2, PTJ_range: Tuple[int, int] = (4, 8), one_GorC=False) -> List[Tuple[Tuple[str, float, float], Tuple[str, float, float]]]:
    """
    Generate and filter PCR primer pairs.

    Args:
        sense_strand (str): Template DNA sequence.
        length_range (Tuple[int, int]): Length range for primers. Default is (17, 30).
        mode (str): Mode for Tm calculation. Default is 'KOD544'.
        overhang (int): Length of the overhang region. Default is 6.
        GC_range (Tuple[int, int]): Acceptable GC content range. Default is (30, 60).
        Tm_range (Tuple[int, int]): Acceptable Tm range. Default is (50, 65).
        GC_end (int): Number of bases at the end of the primer to check for GC content. Default is 2.
        PTJ_range (Tuple[int, int]): Range for primer-template junction. Default is (4, 8).

    Returns:
        List[Tuple[Tuple[str, float, float], Tuple[str, float, float]]]: List of filtered primer pairs.
    """
    forward, reverse = generate_pcrPrimers(sense_strand, length_range=length_range, one_GorC=one_GorC)
    print(f"{len(forward)} forward primer pairs are generated.")
    print(f"{len(reverse)} reverse primer pairs are generated.\n")
    antisense_strand = reverse_complement(sense_strand)
    
    print("Forward primers:")
    filtered_f = primer_filter(backbone, forward, overhang=overhang, GC_range=GC_range, Tm_range=Tm_range, Tm_mode=mode, GC_end=GC_end)
    
    print("Reverse primers:")
    filtered_r = primer_filter(backbone, reverse, overhang=overhang, GC_range=GC_range, Tm_range=Tm_range, Tm_mode=mode, GC_end=GC_end)
    
    primer_pairs = []
    min_PTJ, max_PTJ = PTJ_range
    for f_primer in filtered_f:
        for r_primer in filtered_r:
            if not dimer_forms_PTJ(f_primer[0], r_primer[0], max_PTJ=max_PTJ, min_PTJ=min_PTJ):
                primer_pairs.append((f_primer, r_primer))
    print(f"\n{len(primer_pairs)} primer pairs got.")
    
    sorted_pairs = sorted(primer_pairs, key=primerpairs_sorter)
    return sorted_pairs

def pcrPrimer_writer(PrimerPairs: List[Tuple[Tuple[str, float, float], Tuple[str, float, float]]], output: str = 'pcrPrimers.txt', mode='w') -> Union[bool, None]:
    """
    Write primer pairs to a file.

    Args:
        PrimerPairs (List[Tuple[Tuple[str, float, float], Tuple[str, float, float]]]): List of primer pairs.
        output (str): Output file name. Default is 'pcrPrimers.txt'.

    Returns:
        Union[bool, None]: True if file is written, None if canceled.
    """
    if output in os.listdir('.'):
        go = input("文件已存在，是否覆盖? [y/n]: ")
        if go not in ['y', 'yes', 'Y', 'YES']:
            print("取消")
            return None
    
    with open(output, mode, encoding='utf-8') as f:
        f.write("Pair#\tPrimer\tLength\tGC content/%\tTm/°C\tDelta-Tm\n")
        for i, pair in enumerate(PrimerPairs):
            pf, pf_gc, pf_Tm = pair[0]
            pr, pr_gc, pr_Tm = pair[1]
            f.write(f"{i}-F\t{pf}\t{len(pf)} nt\t{pf_gc*100} %\t{pf_Tm} °C\t\n")
            f.write(f"{i}-R\t{pr}\t{len(pr)} nt\t{pr_gc*100} %\t{pr_Tm} °C\t{abs(pf_Tm-pr_Tm)}\n")
            f.write('\n')
    
    print("文件已写入")
    return True


# Auto infusion primer design
def primerpairs_sorter_infusion(primerpairs: Tuple[Tuple[str, float, float], Tuple[str, float, float]]) -> float:
    """
    Sort primer pairs based on requirments
    top: Tm differences
    middle: total length
    last: gc
    """
    pf, pf_gc, pf_Tm = primerpairs[0]
    pr, pr_gc, pr_Tm = primerpairs[1]
    return abs(pf_Tm - pr_Tm)*10**5 + len(pf+pr)*10**3 + min(pf_gc, pr_gc)

def get_prcPrimerPairs_amplify(sense_strand,backbone, mode='KOD544',overhang=6):
    """默认PCR扩增序列参数"""
    primer_pairs = get_prcPrimerPairs(sense_strand, # top strand in snapgene
                                      backbone,
                                      overhang=overhang,   # check binding site
                                      length_range=(12,48),
                                      GC_range=(20,70),
                                      mode=mode,    # for Tm calculation
                                      Tm_range=(50,75),
                                      GC_end=3, # 3'end gc
                                      PTJ_range=(4,8),
                                      one_GorC=True) # check hetero dimer
    
    #pcrPrimer_writer(primer_pairs, output=output_file)
    return primer_pairs

def get_prcPrimerPairs_HRarm(sense_strand,backbone,mode='idt_default'):
    homologous_arms = get_prcPrimerPairs(sense_strand, # top strand in snapgene
                                      backbone,
                                      overhang=10,   # no need here
                                      length_range=(15,20), # C112 Kit
                                      GC_range=(20,70),
                                      mode=mode,    # for Tm calculation
                                      Tm_range=(37,75), # C112 Kit
                                      GC_end=3, # no need here
                                      PTJ_range=(4,4), # no need here
                                      one_GorC=False) # no need here
                                      
    # pcrPrimer_writer(homologous_arms, output=output_file)
    return homologous_arms

def auto_infusion(insertion='', line_vector='', mode='KOD544', output='infusion.txt',insertion_overhang=7,vector_overhang=6):
    """automatically generate primers for 1 fragment infusion"""
    start_overhang = 10
    # primers_for_fragment: specific p + HR-arm
    frag_p = ''
    start_overhang = 10
    while start_overhang:
        frag_p = get_prcPrimerPairs_amplify(insertion,mode='KOD544',overhang=start_overhang)
        if len(frag_p) <= 0:
            break
        start_overhang -= 1
    sorted_frag_pairs = sorted(frag_p, key=primerpairs_sorter_infusion)
    
    HR_arms = get_prcPrimerPairs_HRarm(sense_strand,mode='idt_default')
    HR_arm_F_R_exchange = [(j,i) for i,j in HR_arms]
    sorted_arm_pairs = sorted(HR_arm_F_R_exchange, key=primerpairs_sorter_infusion)

    # primers_for_vector: specific p
    vector_p = ''
    start_overhang = 10
    while start_overhang:
        vector_p = get_prcPrimerPairs_amplify(line_vector,mode='KOD544',overhang=start_overhang)
        if len(vector_p) <= 0:
            break
        start_overhang -= 1

    sorted_vector_pairs = sorted(primer_pairs, key=primerpairs_sorter_infusion)
    
    pcrPrimer_writer(sorted_frag_pairs,output=output, mode='w')
    pcrPrimer_writer(sorted_arm_pairs, output=output, mode='w+')
    pcrPrimer_writer(sorted_vector_pairs, output=output, mode='w+')



if __name__ == "__main__":
    target_fragment = 'insertion.txt'
    template = 'backbone.txt'
    output_file = 'pcrPrimer.txt'
    backbone_sequence = ''
    with open(template, 'r', encoding='utf-8') as f:
        backbone_sequence = f.read().strip()

    with open(target_fragment, 'r', encoding='utf-8') as f:
        sense_strand = f.read().strip()
        # print(sense_strand)
        # Main
        primer_pairs = get_prcPrimerPairs(sense_strand,# top strand in snapgene
                                      backbone_sequence,
                                      overhang=8,   # 6 check binding site
                                      length_range=(12,48), # (12,48)
                                      GC_range=(20,65), # (20,70)
                                      mode='KOD544',    # for Tm calculation 'idt_default' 'KOD544'
                                      Tm_range=(50,70), # (50, 75)
                                      GC_end=3, # 3
                                      PTJ_range=(3,8), # (4,8)
                                      one_GorC=False) # True
        print(len(primer_pairs))
        # primer_pairs = get_prcPrimerPairs_HRarm(sense_strand, backbone_sequence, mode='idt_default')
        pcrPrimer_writer(primer_pairs, output=output_file)
        
        
        
        
        
        
        
