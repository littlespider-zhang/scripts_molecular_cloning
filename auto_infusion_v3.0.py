# Modified primer_sorter logic 24/12/31
# auto_infusion() is added
import re
import os
import random
from decimal import *
from typing import List, Tuple
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq




###########################################################################
# Note: 'ATGC' means 5'-'ATGC'-3' unless otherwise noted [direction is important]
# Basic functions
def read_fasta(filename):
    """Return a list of dicts containing information about each seq"""
    with open(filename, encoding='utf-8') as f:
        reader = f.read().split('>')

    pre = [item.split('\n') for item in reader[1:]]

    # Return a list with elements representing each string.
    sequences = [{'id': seq[0], 'seq': ''.join(seq[1:])}
                 for index, seq in enumerate(pre)]

    return sequences

def generate_random_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of given length.

    Args:
        length (int): Length of the sequence.

    Returns:
        str: Random DNA sequence.
    """
    bases = ('A', 'C', 'G', 'T')
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
    gc_counts = sequence.count('G') + sequence.count('C')
    gc_content = (gc_counts / len(sequence)) * 100
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


def highlight(sequence: str, sub: str, maker='*'):
    """
    Input: "ATATGCGTAGCC",'GCG'
    Output: '5'-ATAT *GCG* CGTAGCC-3''

    Note: only first sub is highlighted.
    """
    pattern = re.compile(re.escape(sub))
    hl = ''
    start = 0  # start from position 0
    for i in pattern.finditer(sequence):
        a, b = i.span()  # splicing position
        hl += f"{sequence[start:a]} {maker}{sub}{maker} "
        start = a

    hl = "5'-" + hl + "-3'"
    return hl


def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        seq (str): DNA sequence. 5' -> 3'

    Returns:
        str: Reverse complement sequence. 5' -> 3'
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq.upper()))


def former_RC_in_latter(former: str, latter: str):
    """
    RC = reverse complement

    Return:
        # of RC site
    """
    rc_former = reverse_complement(former)
    pattern = re.compile(re.escape(rc_former))

    count = len(pattern.findall(latter))
    if count != 0:
        # print(f"Matched:\t", end='')
        # print(highlight(latter, rc_former, maker='*'))
        pass
    # print(f"RC_in_latter: {count}")
    return count

## Need change: return
def PTJ_possible(primer: str, template: str, PTJ_bp: int = 4, message=False) -> (int, int):
    """

    """
    # Format sequences
    # print(f"In PTJ_possible():")
    if PTJ_bp < 3:
        # print(f"\twarning: random binding occurs frequently at this length: {PTJ_bp}")
        # print(f"\treturn None.")
        return None

    template, primer = template.upper(), primer.upper()

    # print(f"\tsearching:\t5'-{primer[-PTJ_bp:]}-3'")

    count = former_RC_in_latter(primer[-PTJ_bp:], template)
    if count > 0:
        if message:
            print(f"In PTJ_possible():")
            print(f"\tsearching:\t5'-{primer[-PTJ_bp:]}-3'")
            print(f"\tmatched at:",end='')
            print(f"\t{highlight(template, reverse_complement(primer[-PTJ_bp:]))}")
        pass

    return count



# Dimer check and uni_PTJ_check has the same standard
def no_self_dimer(primer: str, min_PTJ_length: int = 4):
    """
    If true, no self dimer forms
    If false, self dimer forms
    """
    # print(f"Checking self dimer for: 5'-{primer}-3';\tlength: {len(primer)} nt")
    for i in range(min_PTJ_length, len(primer)+1, 1):
        check = PTJ_possible(primer, primer, PTJ_bp=i)
        if check == None:
            # print(f"\twarning: PTJ length too short for self dimer check")
            pass
        elif check:
            return False, i


    return True, min_PTJ_length


def uni_PTJ_site(primer: str, backbone: str, min_PTJ_length: int = 4) -> tuple:
    """
    If True, x; means at PTJ length x, only 1 binding site; if x may > min_PTJ_length
    If False, x; means entire primer has multiple binding sites
    """
    # print(f"Checking uni_PTJ for: 5'-{primer}-3'")
    for i in range(min_PTJ_length, len(primer) + 1):
        count = PTJ_possible(primer, backbone, PTJ_bp=i)
        if count == 0:
            print(f"In uni_PTJ_site()")
            print(f"\twarning: no binding site found for: 5'-{primer}-3'")
            print(f"\treturn None")
        elif count == 1:
            # print(f"In uni_PTJ_site()")
            # print(f"\tuni-binding-site at length {i} for: 5'-{primer}-3'")
            return True, i

    return False, len(primer)


def hetero_dimer_check(primer1, primer2, min_PTJ_hetero=5, message=False):
    """
    Check if 2 primers can form PTJ
    """
    # print(f"Checking hetero dimer between:\t5'-{primer1}-3' and 5'-{primer2}-3'")
    # Check if 3' of primer1 anneals to primer2
    check_1 = True, len(primer1)
    for i in range(min_PTJ_hetero, len(primer1) + 1):
        count = PTJ_possible(primer1, primer2, PTJ_bp=i, message=message)
        if count == None:
            # print(f"\twarning: PTJ length too short for self dimer check")
            pass
        elif count:
            check_1 = False, i


    # Check if 3' of primer2 anneals to primer1
    check_2 = True, len(primer2)
    for i in range(min_PTJ_hetero, len(primer2) + 1):
        count = PTJ_possible(primer2, primer1, PTJ_bp=i, message=message)
        if count == None:
            # print(f"\twarning: PTJ length too short for self dimer check")
            pass
        elif count:
            check_2 = False, i


    # Explaination:
    check_3 = (check_1[0] and check_2[0]), max(check_1[1],check_2[1])

    return check_3


### pcrPrimerDesign functions
def generate_pcrPrimers(template: str, length_range: Tuple[int, int] = (17, 30)) -> Tuple[
    List[str], List[str]]:
    """
    Generate all primers with the given length range from 2 ends of the template 

    Args:
        template (str): Template DNA sequence. (linear DNA)
        length_range (Tuple[int, int]): Length range for primers.

    Returns:
        Tuple[List[str], List[str]]: Forward and reverse primers.
    """
    print("\nGeneration")
    print('-' * 40)
    F_strand = template
    primers_F = [F_strand[:i] for i in range(length_range[0], length_range[1] + 1)]
    print(f"F primers: {len(primers_F)} are generated.")

    R_strand = reverse_complement(template)
    primers_R = [R_strand[:i] for i in range(length_range[0], length_range[1] + 1)]
    print(f"R primers: {len(primers_R)} are generated.\n")

    return primers_F, primers_R


def GC_Tm_checker(primer: str, evaluation = {},
                  GC_range: Tuple[int, int] = (30, 60), oneGC=0,
                  Tm_mode: str = 'idt_default', Tm_range: Tuple[int, int] = (50, 65)):
    """
    Args:
        primer (str): Primer DNA sequence.
        GC_range (Tuple[int, int]): Acceptable GC content range.
        oneGC: if primer has and only has 1 G/C at 3' end; oneGC is the value of checked length at 3' end
        Tm_mode (str): Mode for Tm calculation. Default is 'idt_default'.
        Tm_range (Tuple[int, int]): Acceptable Tm range.


    Returns:
        a dictionary with keys:
            'seq'
            'result'
            'gc_content', 'oneGC'
            'Tm_mode', 'Tm'
    """
    evaluation.update({'seq': primer, 'gc_content': 0, 'oneGC': False, 'Tm': 0, 'Tm_mode': Tm_mode, 'GC_Tm_result':False})

    check_points = 0  # monitor if primers fulfill the requirements

    # gc_content
    gc_content = calculate_GC_content(primer)
    if GC_range[0] <= gc_content <= GC_range[1]:
        evaluation['gc_content'] = gc_content
        check_points += 1

    # check if only one G or C at the end with given end length; oneGC == 0 means no restriction
    count_GC = primer.upper()[-oneGC:].count('G') + primer.upper()[-oneGC:].count('C')
    check_GC = (primer.upper()[-1] in 'GC' and (count_GC == 1))
    if  check_GC or (oneGC == 0):
        # print(f"Checking: oneGC = {oneGC} is {check_GC or (oneGC == 0)} for 5'-{primer}-3'")
        evaluation['oneGC'] = True
        check_points += 1

    # Tm
    Tm = calculate_tm(primer, mode=Tm_mode)
    if Tm_range[0] <= Tm <= Tm_range[1]:
        evaluation['Tm'] = Tm
        check_points += 1

    if check_points == 3:

        evaluation['GC_Tm_result'] = True

    return evaluation


def dimer_backbone_PTJ_checker(primer: str, backbone: str, evaluation,
                      PTJ_dimer: int = 5, PTJ_backbone: int = 6):
    """
    Check self-dimer & uni_PTJ_site

    Return:
        'self_dimer'
        'self_dimer_length'
        'uni_PTJ'
        'uni_PTJ_length'
    """
    evaluation.update({'dimer_backbone_PTJ_result': False})

    check_points = 0
    # Check self-dimer
    check_self_dimer, dimer_backbone_PTJ = no_self_dimer(primer, min_PTJ_length=PTJ_dimer)
    evaluation.update({'no self_dimer': (check_self_dimer, dimer_backbone_PTJ)})
    if check_self_dimer:
        check_points += 1

    # Check PTJ between primer and backbone
    check_PTJ, length = uni_PTJ_site(primer, backbone, min_PTJ_length=PTJ_backbone)
    evaluation.update({'uni_PTJ': (check_PTJ, length)})
    if check_PTJ:
        check_points += 1

    if check_points == 2:
        evaluation['dimer_backbone_PTJ_result'] = True

    return evaluation


def primer_checker(primer: str, backbone: str,
                   GC_range: Tuple[int, int] = (30, 60), oneGC=0,
                   Tm_mode: str = 'KOD544', Tm_range: Tuple[int, int] = (50, 65),
                   PTJ_dimer: int = 5, PTJ_backbone: int = 8):
    """
    Filter primers based on GC content, oneGC, and Tm range, Dimer, PTJ ...

    Just a combination of checkers

    Returns:
        List[Tuple[str, float, float]]: List of filtered primers with their GC content and Tm values.
    """
    evaluation_1 = GC_Tm_checker(primer,
                                 GC_range=GC_range, oneGC=oneGC,
                                 Tm_mode=Tm_mode, Tm_range=Tm_range)
    evaluation_2 = dimer_backbone_PTJ_checker(primer, backbone, evaluation=evaluation_1,
                                     PTJ_dimer=PTJ_dimer, PTJ_backbone=PTJ_backbone)


    return evaluation_2


def primerpairs_sorter(primerpairs) -> float:
    """
    Criteria: based on experiences of zyy
        higher Tm
        less PTJ formation
        only one GC at the end
    """
    pf_length, pf_gc, pf_Tm = len(primerpairs[0]["seq"]), primerpairs[0]["gc_content"], primerpairs[0]['Tm']
    pr_length, pr_gc, pr_Tm = len(primerpairs[1]["seq"]), primerpairs[1]['gc_content'], primerpairs[1]['Tm']
    pf_GCend, pr_GCend = int(primerpairs[0]['seq'][-1] in 'GCgc'), int(primerpairs[1]['seq'][-1] in 'GCgc')

    Tm_difference = 0
    if abs(pf_Tm - pr_Tm) >= 0.5:
        Tm_difference = abs(pf_Tm - pr_Tm)*1.2

    GCend = pf_GCend + pr_GCend

    pf_self, pf_uniPTJ = primerpairs[0]['no self_dimer'][1], primerpairs[0]['uni_PTJ'][1]
    pr_self, pr_uniPTJ = primerpairs[1]['no self_dimer'][1], primerpairs[1]['uni_PTJ'][1]

    GC_Tm = Tm_difference*1000 - GCend * 100 - max(pf_gc, pr_gc) + min(pf_length, pr_length)
    PTJ = max(pf_uniPTJ, pr_uniPTJ) * 10 + max(pf_self,pr_self)
    return  GC_Tm*100 - PTJ


def get_pcrPrimerPairs(template: str, backbone: str, length_range=(17, 30),
                       GC_range: Tuple[int, int] = (30, 60), oneGC=0,
                       Tm_mode: str = 'KOD544', Tm_range: Tuple[int, int] = (50, 65),
                       PTJ_dimer: int = 5, PTJ_backbone: int = 8,
                       PTJ_hetero: int = 5):
    """
    Generate and filter PCR primer pairs.
    backbone 5' -> 3'
    f aneals to reverse_complement of backbone
    """
    forward, reverse = generate_pcrPrimers(template, length_range=length_range)
    # Be really careful here ****************
    rc_backbone = reverse_complement(backbone)

    print("Filtration.")
    print('-' * 40)
    # F
    filtered_f = []
    f_check_gc = 0
    f_check_ptj = 0
    for f in forward:
        f_result = primer_checker(f, rc_backbone,
                                  GC_range=GC_range, oneGC=oneGC,
                                  Tm_mode=Tm_mode, Tm_range=Tm_range,
                                  PTJ_dimer=PTJ_dimer, PTJ_backbone=PTJ_backbone)
        # print(f_result)
        if f_result['GC_Tm_result']:
            f_check_gc += 1
        if f_result['dimer_backbone_PTJ_result']:
            f_check_ptj += 1

        if f_result['GC_Tm_result'] and f_result['dimer_backbone_PTJ_result']:
            filtered_f.append(f_result.copy())  ### really important ### interesting

    print(f"F primers: {len(filtered_f)} passed checking..."
          f"\t{f_check_gc} passed GC checking"
          f"\t{f_check_ptj} passed PTJ checking")
    # R
    filtered_r = []
    r_check_gc = 0
    r_check_ptj = 0
    for r in reverse:
        r_result = primer_checker(r, backbone,
                                  GC_range=GC_range, oneGC=oneGC,
                                  Tm_mode=Tm_mode, Tm_range=Tm_range,
                                  PTJ_dimer=PTJ_dimer, PTJ_backbone=PTJ_backbone)
        # print(r_result)
        if r_result['GC_Tm_result']:
            r_check_gc += 1
        if r_result['dimer_backbone_PTJ_result']:
            r_check_ptj += 1

        if r_result['GC_Tm_result'] and r_result['dimer_backbone_PTJ_result']:
            filtered_r.append(r_result.copy())



    print(f"R primers: {len(filtered_r)} passed checking..."
          f"\t{r_check_gc} passed GC checking"
          f"\t{r_check_ptj} passed PTJ checking")

    # Make good pairs on filtered primers
    primer_pairs = []
    for f_primer in filtered_f:
        for r_primer in filtered_r:
            # print(f"{f_primer['seq']}\t{r_primer['seq']}")
            no_hetero, PTJ_hetero = hetero_dimer_check(f_primer['seq'], r_primer['seq'], min_PTJ_hetero=PTJ_hetero, message=False)
            if no_hetero:
                f_primer['hetero_dimer_check_length'] = PTJ_hetero
                r_primer['hetero_dimer_check_length'] = PTJ_hetero
                primer_pairs.append((f_primer.copy(), r_primer.copy()))

    print(f"\nPair generation.\n"
          f"{'-' * 40}\n"
          f"F/R pairs: {len(primer_pairs)} good pairs got !\t\tTotal  {len(filtered_f) * len(filtered_r)} paris")

    sorted_pairs = sorted(primer_pairs, key=primerpairs_sorter)
    return sorted_pairs

def primer_writer(primers, output='auto_infusion.txt', mode='w', name=''):
    """
    write primer pairs information to txt
    """
    if output in os.listdir('.'):
        go = input(f"\n{output} already exits，overwrite ? [y/n]: ")
        if go not in ['y', 'yes', 'Y', 'YES']:
            print("Cancelled")
            return None

    with open(output, mode, encoding='utf-8') as f:
        f.write("Pair#\tPrimer\tLength\tGC content/%\tTm/°C\tScore & Delta-Tm/°C\n")
        for i, pair in enumerate(primers):
            pf, pr= pair
            if name:
                i = name

            f.write(f"{i}-F\t{pf['seq']}\t{len(pf['seq'])} nt\t{round(pf['gc_content'])} %\t{round(pf['Tm'])} °C\t{round(primerpairs_sorter(pair))}\n")
            f.write(f"{i}-R\t{pr['seq']}\t{len(pr['seq'])} nt\t{round(pr['gc_content'])} %\t{round(pr['Tm'])} °C\t{abs(round(pf['Tm']-pr['Tm'],1))} °C\n")
            f.write('\n')   # round(d) * 100 VS round(d * 100)


def auto_infusion(insertion, insertion_template, backbone, backbone_template, Tm_mode='KOD544'):
    """
    select the best primer pair
    """
    print(f"Insertion ...")
    insertion_pair = get_pcrPrimerPairs(insertion, insertion_template,
                                 length_range=(12, 48),  # (12, 48)
                                 GC_range=(20, 60),  # (30 ,60) **********
                                 Tm_range=(50, 65),  # (50, 65) ***********
                                 oneGC=0,  # 0 means free *****************
                                 PTJ_backbone=6,  # 6 *********************  check primer unspecific binding to template
                                 Tm_mode=Tm_mode,
                                 PTJ_dimer=4,  # 4  check in-primer unspecific binding
                                 PTJ_hetero=4)[0]  # 4 check cross-primer unspecific binding

    print(f"\nBackbone ...")
    backbone_pair = get_pcrPrimerPairs(backbone, backbone_template,
                                        length_range=(12, 48),  # (12, 48)
                                        GC_range=(20, 60),  # (30 ,60) **********
                                        Tm_range=(55, 68),  # (50, 65) ***********
                                        oneGC=0,  # 0 means free *****************
                                        PTJ_backbone=6, # 6 *********************  check primer unspecific binding to template
                                        Tm_mode=Tm_mode,
                                        PTJ_dimer=4,  # 4  check in-primer unspecific binding
                                        PTJ_hetero=4)[0]  # 4 check cross-primer unspecific binding
    primer_writer([insertion_pair], output='auto_infusion.txt', mode='w',name='insert')
    primer_writer([backbone_pair], output='auto_infusion.txt', mode='a', name='Line')



if __name__ == "__main__":
    # Function test: auto_infusion()
    """"""
    insertion = read_fasta('infusion_insertion.txt')
    backbone = read_fasta('infusion_backbone.txt')

    auto_infusion(insertion[0]['seq'], insertion[1]['seq'], backbone[0]['seq'], backbone[1]['seq'], Tm_mode='KOD544')


