# ================================
# PCR Primer Design Utility (v5)
# ================================
#
# Pipeline:
#   1. Generate forward / reverse primer candidates from the insertion sequence.
#   2. Filter each candidate:
#        a. Unique 3'-end binding site on the backbone
#        b. GC content within range
#        c. No all-GC 3' tail
#        d. Tm within range  (nearest-neighbor, KOD544 or IDT conditions)
#   3. Reject pairs that form 3'-end dimers (PTJ check).
#   4. Sort survivors by ΔTm (primary) then minimum GC (secondary).
#   5. Write results to a tab-delimited text file.

import os
import re
import random
from dataclasses import dataclass
from typing import List, Optional, Tuple

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


# ================================================================
# Data structure
# ================================================================

@dataclass
class Primer:
    """A primer sequence with its computed thermodynamic properties.

    Attributes:
        seq : sequence (5'→3', uppercase)
        gc  : GC content (%), rounded to 1 decimal place
        tm  : melting temperature (°C, nearest-neighbor), rounded to 1 decimal place
    """
    seq: str
    gc: float
    tm: float


# ================================================================
# Layer 1 — basic sequence calculations  (no side effects)
# ================================================================

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.  "ATGC" → "GCAT" """
    table = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(table[b] for b in reversed(seq.upper()))


def calculate_gc_content(sequence: str) -> float:
    """Return GC content (%) rounded to 1 decimal place."""
    s = sequence.upper()
    return round((s.count('G') + s.count('C')) / len(s) * 100, 1)


def calculate_tm(primer: str, mode: str = 'idt_default') -> float:
    """Return Tm (°C) using the nearest-neighbor method.

    Modes and reaction conditions:
        'idt_default' : [primer] 200 nM, [Na+] 44.2 mM, no Mg²⁺, no dNTPs
        'KOD544'      : [primer] 300 nM, [Na+] 50 mM, [Mg²⁺] 2 mM, [dNTPs] 0.16 mM
                        (KOD-Plus-Neo polymerase recommended conditions)

    Salt correction: SantaLucia 1998 (saltcorr=3).
    """
    seq = Seq(primer.upper())
    if mode == 'idt_default':
        tm = mt.Tm_NN(seq, dnac1=200, Na=44.2, Mg=0,   dNTPs=0,    saltcorr=3)
    elif mode == 'KOD544':
        tm = mt.Tm_NN(seq, dnac1=300, Na=50,   Mg=2,   dNTPs=0.16, saltcorr=3)
    else:
        raise ValueError(f"Unknown mode '{mode}'. Use 'idt_default' or 'KOD544'.")
    return round(tm, 1)


# ================================================================
# Layer 2 — single-primer filtering
# ================================================================

def _passes_params(primer: str,
                   gc_range: Tuple[int, int],
                   tm_range: Tuple[int, int],
                   tm_mode:  str,
                   gc_end:   int) -> Optional[Primer]:
    """Apply GC / Tm / 3'-end filters to one primer.

    Filters (early-exit order, cheapest first):
        1. All-GC 3' tail: rejects primers where the last `gc_end` bases are
           all G/C — such tails increase non-specific annealing risk.
        2. GC content must fall within gc_range (%).
        3. Tm must fall within tm_range (°C).
           Computed last because it is the most expensive step.

    Returns a Primer object on success, None on any filter failure.
    """
    if all(b.upper() in 'GC' for b in primer[-gc_end:]):
        return None
    gc = calculate_gc_content(primer)
    if not (gc_range[0] <= gc <= gc_range[1]):
        return None
    tm = calculate_tm(primer, mode=tm_mode)
    if not (tm_range[0] <= tm <= tm_range[1]):
        return None
    return Primer(seq=primer, gc=gc, tm=tm)


def _count_binding_sites(template: str, template_rc: str,
                         primer: str, overhang: int) -> int:
    """Count how many times the primer's 3' overhang appears across both strands."""
    pattern = re.compile(re.escape(primer[-overhang:]))
    return len(pattern.findall(template)) + len(pattern.findall(template_rc))


# ================================================================
# Layer 3 — pipeline functions  (generate → filter → pair → write)
# ================================================================

def generate_pcr_primers(template: str,
                         length_range: Tuple[int, int] = (17, 30),
                         one_gc_end: bool = False) -> Tuple[List[str], List[str]]:
    """Generate candidate forward and reverse primers from an insertion sequence.

    Forward primers  : 5' prefixes of the sense strand (template[:n]).
    Reverse primers  : 5' prefixes of the antisense strand (RC of template),
                       which correspond to the 3' end of the sense strand.

    Args:
        template     : sense strand of the insertion (5'→3', as shown in SnapGene)
        length_range : (min_nt, max_nt) primer length range to generate
        one_gc_end   : if True, keep only primers whose last base is G or C
    """
    antisense = reverse_complement(template)
    primers_f = [template[:i]   for i in range(length_range[0], length_range[1] + 1)]
    primers_r = [antisense[:i]  for i in range(length_range[0], length_range[1] + 1)]
    if one_gc_end:
        primers_f = [p for p in primers_f if p[-1].upper() in 'GC']
        primers_r = [p for p in primers_r if p[-1].upper() in 'GC']
    return primers_f, primers_r


def filter_primers(backbone: str,
                   primers:  List[str],
                   overhang: int = 6,
                   gc_range: Tuple[int, int] = (30, 60),
                   tm_range: Tuple[int, int] = (50, 65),
                   tm_mode:  str = 'KOD544',
                   gc_end:   int = 2) -> List[Primer]:
    """Filter primer candidates against backbone uniqueness and thermodynamic criteria.

    Each candidate is checked independently for:
        • Binding uniqueness : the 3' overhang must match exactly one site on the
          backbone (both strands).  backbone_rc is precomputed once for efficiency.
        • Parameter fitness  : GC content, Tm, and 3'-end composition (see _passes_params).

    Rejection counts are printed separately for each failure mode.

    Returns a list of Primer objects that pass all criteria.
    """
    backbone_rc = reverse_complement(backbone)   # computed once, reused for all primers

    passed        = []
    n_no_site     = 0   # 3' overhang absent from backbone
    n_multi_site  = 0   # 3' overhang found at >1 location
    n_bad_params  = 0   # GC / Tm / 3'-end check failed

    for primer in primers:
        n_sites = _count_binding_sites(backbone, backbone_rc, primer, overhang)
        result  = _passes_params(primer, gc_range, tm_range, tm_mode, gc_end)

        if n_sites == 1 and result is not None:
            passed.append(result)

        if   n_sites == 0: n_no_site    += 1
        elif n_sites  > 1:
            print(f"# of Binding sites: {n_sites}")
            n_multi_site += 1
        if result is None: n_bad_params += 1

    print(
        f"  {n_no_site} primers: no binding site on backbone.\n"
        f"  {n_multi_site} primers: multiple binding sites (rejected).\n"
        f"  {n_bad_params} primers: GC / Tm / 3'-end parameters not met."
    )
    return passed


def _forms_ptj_dimer(primer1: str, primer2: str,
                     ptj_range: Tuple[int, int]) -> bool:
    """Return True if the two primers can form a 3'-end dimer (Partial Terminal Junction).

    A PTJ dimer occurs when the 3' tail of one primer hybridises to the other,
    allowing mutual extension and generation of spurious products.

    Each primer's 3'-end fragments (lengths min_ptj … max_ptj) are tested against
    the reverse complement of the other primer, longest fragment first, so the most
    stable potential overlap is reported first.
    """
    min_ptj, max_ptj = ptj_range
    if min_ptj == max_ptj:
        return False

    def _3prime_fragments(seq: str) -> List[str]:
        """3'-end substrings from length max_ptj down to min_ptj."""
        cap = min(max_ptj, len(seq))
        # range(-min_ptj, -cap-1, -1) yields indices for suffixes of decreasing length
        return [seq[i:] for i in range(-min_ptj, -cap - 1, -1)]

    rc1, rc2 = reverse_complement(primer1), reverse_complement(primer2)

    for ext in reversed(_3prime_fragments(primer1)):
        if ext in rc2:
            print(f"3'end anneals with {len(ext)} base pairs.")
            print(f"Primer1:{primer1}-->{ext}\nPrimer2:{primer2}-->{reverse_complement(ext)}\n{ext}")
            return True

    for ext in reversed(_3prime_fragments(primer2)):
        if ext in rc1:
            print(f"3'end anneals with {len(ext)} base pairs.")
            print(f"Primer1:{primer1}-->{reverse_complement(ext)}\nPrimer2:{primer2}-->{ext}\n")
            return True

    return False


def get_pcr_primer_pairs(sense_strand: str,
                         backbone:     str,
                         length_range: Tuple[int, int] = (17, 30),
                         tm_mode:      str = 'KOD544',
                         overhang:     int = 6,
                         gc_range:     Tuple[int, int] = (30, 60),
                         tm_range:     Tuple[int, int] = (50, 65),
                         gc_end:       int = 2,
                         ptj_range:    Tuple[int, int] = (4, 8),
                         one_gc_end:   bool = False) -> List[Tuple[Primer, Primer]]:
    """Full pipeline: generate, filter, and rank PCR primer pairs.

    Steps:
        1. Generate forward / reverse candidates from sense_strand.
        2. Filter each set against backbone uniqueness + Tm/GC criteria.
        3. Reject pairs that form 3'-end dimers (PTJ check).
        4. Sort by ΔTm (ascending), then by the lower GC of the pair.
           Sort key: ΔTm × 1000 + min(GC_F, GC_R)  — the ×1000 weight ensures
           ΔTm dominates while GC acts as a tiebreaker.

    Args:
        sense_strand : top strand of the insertion (5'→3', as shown in SnapGene)
        backbone     : full vector sequence used for uniqueness checking
        length_range : (min_nt, max_nt) primer length range
        tm_mode      : Tm calculation conditions — 'idt_default' or 'KOD544'
        overhang     : number of 3'-end bases used for backbone uniqueness check
        gc_range     : acceptable GC content range (%)
        tm_range     : acceptable Tm range (°C)
        gc_end       : number of 3'-end bases checked for all-GC clamp
        ptj_range    : (min_bp, max_bp) 3'-end overlap range for dimer detection
        one_gc_end   : if True, require the last base of each primer to be G or C

    Returns:
        Sorted list of (forward_Primer, reverse_Primer) tuples.
    """
    fwd_candidates, rev_candidates = generate_pcr_primers(sense_strand, length_range, one_gc_end)
    print(f"{len(fwd_candidates)} forward primers generated.")
    print(f"{len(rev_candidates)} reverse primers generated.\n")

    print("Filtering forward primers:")
    filtered_f = filter_primers(backbone, fwd_candidates, overhang, gc_range, tm_range, tm_mode, gc_end)

    print("\nFiltering reverse primers:")
    filtered_r = filter_primers(backbone, rev_candidates, overhang, gc_range, tm_range, tm_mode, gc_end)

    pairs = [
        (fwd, rev)
        for fwd in filtered_f
        for rev in filtered_r
        if not _forms_ptj_dimer(fwd.seq, rev.seq, ptj_range)
    ]

    print(f"\n{len(pairs)} primer pairs passed all filters.")
    return sorted(pairs, key=lambda p: abs(p[0].tm - p[1].tm) * 1000 + min(p[0].gc, p[1].gc))


def write_primer_pairs(primer_pairs:   List[Tuple[Primer, Primer]],
                       product_length: int,
                       output:         str = 'pcrPrimers.txt',
                       mode:           str = 'w') -> Optional[bool]:
    """Write primer pairs to a tab-delimited text file.

    Columns:  Pair# | Primer | Length | GC content/% | Tm/C | Delta-Tm | Product/bp
    Delta-Tm and Product/bp appear only on the -R row (pair-level values).

    The PCR product length = len(sense_strand), because forward primers always
    start at position 0 and reverse primers always end at the last base of the
    sense strand (see generate_pcr_primers).

    Args:
        primer_pairs   : list of (forward Primer, reverse Primer) tuples
        product_length : expected PCR product length in bp (= len(sense_strand))
        output         : output filename
        mode           : file open mode — 'w' to overwrite, 'a' to append
    """
    if output in os.listdir('.'):
        if input("文件已存在，是否覆盖? [y/n]: ").strip().lower() not in ('y', 'yes'):
            print("取消")
            return None

    with open(output, mode, encoding='utf-8') as f:
        f.write("Pair#\tPrimer\tLength\tGC content/%\tTm/C\tDelta-Tm\tProduct/bp\n")
        for i, (pf, pr) in enumerate(primer_pairs):
            delta_tm = round(abs(pf.tm - pr.tm), 1)
            f.write(f"{i}-F\t{pf.seq}\t{len(pf.seq)} nt\t{pf.gc} %\t{pf.tm} C\t\t\n")
            f.write(f"{i}-R\t{pr.seq}\t{len(pr.seq)} nt\t{pr.gc} %\t{pr.tm} C\t{delta_tm}\t{product_length} bp\n")
            f.write('\n')

    print("文件已写入")
    return True


# ================================================================
# Utilities — random primer generation (testing / benchmarking)
# ================================================================

def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence of the given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def write_random_sequences(length: int = 25, count: int = 10) -> None:
    """Write `count` random sequences of `length` nt to 'random_primers.txt'.
    Output format (tab-delimited): <index>  <sequence>
    """
    with open('random_primers.txt', 'w', encoding='utf-8') as f:
        for i in range(count):
            f.write(f'{i}\t{generate_random_sequence(length)}\n')


def batch_tm(filepath: str = 'random_primers.txt') -> None:
    """Print Tm and GC content for each sequence in a tab-delimited primer file.
    Expected input format: <index>  <sequence>
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            seq = line.split('\t')[1].strip()
            print(f"{seq}\t{calculate_tm(seq)}\t{calculate_gc_content(seq)}")


# ================================================================
# Main
# ================================================================

if __name__ == '__main__':
    insertion_file = 'insertion.txt'
    backbone_file  = 'backbone.txt'
    output_file    = 'pcrPrimer.txt'

    with open(backbone_file,  'r', encoding='utf-8') as f:
        backbone_sequence = f.read().strip()
    with open(insertion_file, 'r', encoding='utf-8') as f:
        sense_strand = f.read().strip()

    # Parameters below have been experimentally validated and work well in practice.
    primer_pairs = get_pcr_primer_pairs(
        sense_strand,           # top strand of insertion as shown in SnapGene (5'→3')
        backbone_sequence,
        overhang     = 10,      # 3'-end bases used for backbone uniqueness check
        length_range = (12, 50),# primer length range (nt)
        gc_range     = (20, 85),# acceptable GC content (%)
        tm_mode      = 'KOD544',# Tm calculated under KOD-Plus-Neo conditions
        tm_range     = (50, 80),# acceptable Tm range (°C)
        gc_end       = 3,       # reject primers with all-GC in last 3 bases
        ptj_range    = (4, 8),  # 3'-end overlap range for dimer detection (bp)
        one_gc_end   = False,   # GC clamp at 3' end not required for KOD
    )

    print(len(primer_pairs))
    write_primer_pairs(primer_pairs, product_length=len(sense_strand), output=output_file)