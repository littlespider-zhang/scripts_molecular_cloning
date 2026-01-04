# ================================
# Primer Design Utility (Refactored, FIXED)
# ================================
# IMPORTANT:
# - Output and behavior are IDENTICAL to the original script
# - Only code organization, readability, and comments were added

import os
import re
import random
from typing import List, Tuple, Union
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


# ------------------------------------------------
# Section 1. Basic sequence utilities
# ------------------------------------------------

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    Example: "ATGC" -> "GCAT"
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[b] for b in reversed(seq.upper()))


def calculate_GC_content(sequence: str) -> float:
    """
    Calculate GC content (%) of a DNA sequence.
    Returns rounded value to 1 decimal place.
    """
    sequence = sequence.upper()
    gc = sequence.count('G') + sequence.count('C')
    return round(gc / len(sequence) * 100, 1)


def calculate_tm(primer: str, mode: str = 'idt_default') -> float:
    """
    Calculate melting temperature (Tm) of a primer using nearest-neighbor method.

    Parameters:
        primer: DNA sequence of the primer
        mode: calculation parameters ('idt_default' or 'KOD544')

    Returns:
        Tm in Celsius (rounded to 1 decimal)
    """
    seq = Seq(primer.upper())
    if mode == 'idt_default':
        tm = mt.Tm_NN(seq, dnac1=200, Na=44.2, Mg=0, dNTPs=0, saltcorr=3)
    elif mode == 'KOD544':
        tm = mt.Tm_NN(seq, dnac1=300, Na=50, Mg=2, dNTPs=0.16, saltcorr=3)
    else:
        raise ValueError('Error: wrong mode')
    return round(tm, 1)


# ------------------------------------------------
# Section 2. Random primer generation (for testing)
# ------------------------------------------------

def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def write_random_sequence(length: int = 25, numbers: int = 10) -> None:
    """
    Generate multiple random sequences and write them to 'random_primers.txt'.

    Parameters:
        length: length of each primer
        numbers: number of primers to generate
    """
    primer_list = [generate_random_sequence(length) for _ in range(numbers)]
    with open('random_primers.txt', 'w', encoding='utf-8') as f:
        for i, primer in enumerate(primer_list):
            f.write(f'{i}\t{primer}\n')


def Batch_Tm() -> None:
    """
    Read primers from 'random_primers.txt' and print their Tm and GC content.
    """
    with open('random_primers.txt', 'r', encoding='utf-8') as f:
        primers = [line.split('\t') for line in f.readlines()]
        for primer_info in primers:
            seq = primer_info[1].strip()
            print(f"{seq}\t{calculate_tm(seq)}\t{calculate_GC_content(seq)}")


# ------------------------------------------------
# Section 3. Subsequence / extension generators
# ------------------------------------------------

def generate_subsequence(sequence: str, length: int = 8) -> List[str]:
    """
    Generate all subsequences of given length from a sequence.

    Returns:
        List of substrings
    """
    return [sequence[i:i + length] for i in range(len(sequence) - length + 1)]


def generate_3PrimeExtension(sequence: str,
                             length_range: Tuple[int, int] = (3, 8)) -> List[str]:
    """
    Generate possible 3' extensions of a primer within specified length range.

    Returns:
        List of 3' end subsequences (short -> long)
    """
    PrimeExtension = []
    s, e = length_range
    e = min(e, len(sequence))
    for i in range(-s, -e - 1, -1):
        PrimeExtension.append(sequence[i:])
    return PrimeExtension


# ------------------------------------------------
# Section 4. Binding & dimer checks
# ------------------------------------------------

def uni_binding_site(template: str,
                     primer: str,
                     overhang: int = 8) -> Union[bool, None]:
    """
    Check if a primer binds uniquely to a template.

    Parameters:
        template: DNA template
        primer: primer sequence
        overhang: number of bases at 3' end used for checking

    Returns:
        True if unique binding, False if multiple binding, None if no binding
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


def dimer_forms_PTJ(primer1: str,
                    primer2: str,
                    max_PTJ: int = 8,
                    min_PTJ: int = 4) -> bool:
    """
    Check if two primers can form a primer dimer at the 3' ends (partial terminal junction).

    Returns:
        True if dimer forms, False otherwise
    """
    if min_PTJ == max_PTJ:
        return False

    # Check primer1's 3' ends against primer2's reverse complement
    primer2_rc = reverse_complement(primer2)
    ends = generate_3PrimeExtension(primer1, (min_PTJ, max_PTJ))
    for e in sorted(ends, key=len, reverse=True):
        if e in primer2_rc:
            print(f"3'end anneals with {len(e)} base pairs.")
            print(f"Primer1:{primer1}-->{e}\nPrimer2:{primer2}-->{reverse_complement(e)}\n{e}")
            return True

    # Check primer2's 3' ends against primer1's reverse complement
    primer1_rc = reverse_complement(primer1)
    ends = generate_3PrimeExtension(primer2, (min_PTJ, max_PTJ))
    for e in sorted(ends, key=len, reverse=True):
        if e in primer1_rc:
            print(f"3'end anneals with {len(e)} base pairs.")
            print(f"Primer1:{primer1}-->{reverse_complement(e)}\nPrimer2:{primer2}-->{e}\n")
            return True

    return False


# ------------------------------------------------
# Section 5. Primer parameter filtering
# ------------------------------------------------

def parameter_filtering(primer: str,
                        GC_range: Tuple[int, int] = (30, 60),
                        Tm_range: Tuple[int, int] = (50, 65),
                        Tm_mode: str = 'idt_default',
                        GC_end: int = 2) -> Tuple[bool, Tuple[str, Union[float, None], Union[float, None]]]:
    """
    Filter primers based on GC content, Tm, and 3' end GC composition.

    Returns:
        Tuple:
            - True/False if primer passes filters
            - (primer sequence, GC %, Tm)
    """
    # Check for consecutive G/C at 3' end
    all_gc = all(n.lower() in 'gc' for n in primer[-GC_end:])
    if all_gc:
        return False, (primer, None, None)

    gc = calculate_GC_content(primer)
    if not (GC_range[0] <= gc <= GC_range[1]):
        return False, (primer, gc, None)

    tm = calculate_tm(primer, mode=Tm_mode)
    if not (Tm_range[0] <= tm <= Tm_range[1]):
        return False, (primer, gc, tm)

    return True, (primer, gc, tm)


# ------------------------------------------------
# Section 6. PCR primer generation
# ------------------------------------------------

def generate_pcrPrimers(template: str,
                        length_range: Tuple[int, int] = (17, 30),
                        one_GorC: bool = False) -> Tuple[List[str], List[str]]:
    """
    Generate candidate forward and reverse primers from a template sequence.

    Returns:
        Tuple of forward primers list, reverse primers list
    """
    sense_strand = template
    antisense_strand = reverse_complement(template)

    primers_F = [sense_strand[:i] for i in range(length_range[0], length_range[1] + 1)]
    primers_R = [antisense_strand[:i] for i in range(length_range[0], length_range[1] + 1)]

    if one_GorC:
        primers_F = [p for p in primers_F if p[-1] in 'GCgc']
        primers_R = [p for p in primers_R if p[-1] in 'GCgc']

    return primers_F, primers_R


def primer_filter(template: str,
                  primers: List[str],
                  overhang: int = 6,
                  GC_range: Tuple[int, int] = (30, 60),
                  Tm_range: Tuple[int, int] = (50, 65),
                  Tm_mode: str = 'KOD544',
                  GC_end: int = 2) -> List[Tuple[str, float, float]]:
    """
    Filter a list of primers based on uniqueness, GC, and Tm.

    Returns:
        List of primers passing all criteria: [(sequence, GC%, Tm), ...]
    """
    filtered_primers = []
    filtered_uni = 0
    filtered_fullfil = 0

    for primer in primers:
        unique = uni_binding_site(template, primer, overhang=overhang)
        fullfil, primer_info = parameter_filtering(
            primer,
            GC_range=GC_range,
            Tm_range=Tm_range,
            Tm_mode=Tm_mode,
            GC_end=GC_end
        )
        if unique and fullfil:
            filtered_primers.append(primer_info)
        if not unique:
            filtered_uni += 1
        if not fullfil:
            filtered_fullfil += 1

    print(
        f"{filtered_uni} primers have multiple binding sites.\n{filtered_fullfil} primers do not fullfil requirements.")
    return filtered_primers


# ------------------------------------------------
# Section 7. Primer pair construction
# ------------------------------------------------

def primerpairs_sorter(primerpairs):
    """
    Sorting key function for primer pairs.
    Prioritizes small Tm difference, then lowest GC content.
    """
    pf, pf_gc, pf_Tm = primerpairs[0]
    pr, pr_gc, pr_Tm = primerpairs[1]
    return abs(pf_Tm - pr_Tm) * 1000 + min(pf_gc, pr_gc)


def get_prcPrimerPairs(sense_strand: str,
                       backbone: str,
                       length_range: Tuple[int, int] = (17, 30),
                       mode: str = 'KOD544',
                       overhang: int = 6,
                       GC_range: Tuple[int, int] = (30, 60),
                       Tm_range: Tuple[int, int] = (50, 65),
                       GC_end: int = 2,
                       PTJ_range: Tuple[int, int] = (4, 8),
                       one_GorC: bool = False):
    """
    Generate and filter PCR primer pairs for a template and backbone.

    Returns:
        Sorted list of primer pairs: [(forward primer info, reverse primer info), ...]
    """
    forward, reverse = generate_pcrPrimers(sense_strand, length_range, one_GorC)
    print(f"{len(forward)} forward primer pairs are generated.")
    print(f"{len(reverse)} reverse primer pairs are generated.\n")

    print("Forward primers:")
    filtered_f = primer_filter(backbone, forward, overhang, GC_range, Tm_range, mode, GC_end)

    print("Reverse primers:")
    filtered_r = primer_filter(backbone, reverse, overhang, GC_range, Tm_range, mode, GC_end)

    primer_pairs = []
    min_PTJ, max_PTJ = PTJ_range
    for f_primer in filtered_f:
        for r_primer in filtered_r:
            # Exclude primers forming dimers
            if not dimer_forms_PTJ(
                    f_primer[0],
                    r_primer[0],
                    max_PTJ=max_PTJ,
                    min_PTJ=min_PTJ
            ):
                primer_pairs.append((f_primer, r_primer))

    print(f"\n{len(primer_pairs)} primer pairs got.")
    return sorted(primer_pairs, key=primerpairs_sorter)


# ------------------------------------------------
# Section 8. Output
# ------------------------------------------------

def pcrPrimer_writer(PrimerPairs,
                     output: str = 'pcrPrimers.txt',
                     mode: str = 'w') -> Union[bool, None]:
    """
    Write primer pairs to a text file with detailed info.

    Parameters:
        PrimerPairs: list of tuples from get_prcPrimerPairs
        output: filename to write
        mode: file write mode ('w' or 'a')

    Returns:
        True if file written, None if cancelled
    """
    # Check if file exists and ask for overwrite
    if output in os.listdir('.'):
        go = input("文件已存在，是否覆盖? [y/n]: ")
        if go not in ['y', 'yes', 'Y', 'YES']:
            print("取消")
            return None

    with open(output, mode, encoding='utf-8') as f:
        f.write("Pair#\tPrimer\tLength\tGC content/%\tTm/C\tDelta-Tm\n")
        for i, pair in enumerate(PrimerPairs):
            pf, pf_gc, pf_Tm = pair[0]
            pr, pr_gc, pr_Tm = pair[1]
            f.write(f"{i}-F\t{pf}\t{len(pf)} nt\t{pf_gc} %\t{pf_Tm} C\t\n")
            f.write(f"{i}-R\t{pr}\t{len(pr)} nt\t{pr_gc} %\t{pr_Tm} C\t{round(abs(pf_Tm - pr_Tm), 1)}\n")
            f.write('\n')

    print("文件已写入")
    return True


# ------------------------------------------------
# Main
# ------------------------------------------------

if __name__ == '__main__':
    insertion = 'insertion.txt'
    backbone = 'backbone.txt'
    output_file = 'pcrPrimer.txt'

    # Load template and backbone sequences
    with open(backbone, 'r', encoding='utf-8') as f:
        backbone_sequence = f.read().strip()

    with open(insertion, 'r', encoding='utf-8') as f:
        sense_strand = f.read().strip()

    # Generate and filter primer pairs
    #### Following parameters have been optimized and work well in practice.
    primer_pairs = get_prcPrimerPairs(
        sense_strand,   # top strand in snapgene
        backbone_sequence,
        overhang=10,
        length_range=(12, 50),
        GC_range=(20, 85),
        mode='KOD544',  # for Tm calculation
        Tm_range=(50, 80),
        GC_end=3,
        PTJ_range=(4, 8),
        one_GorC=False  # one_GorC at the end is not necessary
    )

    print(len(primer_pairs))
    pcrPrimer_writer(primer_pairs, output=output_file)
