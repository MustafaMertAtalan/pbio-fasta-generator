
# Album No : s26808
# Date     : 2026-05-06
# Desc     : Random DNA sequence generator in FASTA format.
#            Supports motif search, complement, transcription,
#            ORF finding, sliding-window GC analysis


import random
import csv



#1. CORE FUNCTIONS

def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length.
    Uses random.choices for uniform nucleotide distribution."""
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=length))


def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary with nucleotide percentages and GC content.

    Keys: 'A', 'C', 'G', 'T' (float, %), 'GC' (float, %).
    Only uppercase nucleotide characters are counted;
    embedded name letters (lowercase) are intentionally ignored.
    """
    # Filter out lowercase letters (inserted name) before counting
    clean = [ch for ch in sequence if ch.isupper()]
    total = len(clean)

    if total == 0:
        return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'GC': 0.0}

    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for ch in clean:
        if ch in counts:
            counts[ch] += 1

    stats = {nuc: round(counts[nuc] / total * 100, 2) for nuc in counts}
    stats['GC'] = round((counts['G'] + counts['C']) / total * 100, 2)
    return stats


def insert_name(sequence: str, name: str) -> str:
    """Inserts the user's name (lowercase) at a random position in the sequence.

    The name does NOT affect statistics because calculate_stats()
    ignores lowercase characters.
    """
    name_lower = name.lower()
    pos = random.randint(0, len(sequence))  # inclusive on both ends
    return sequence[:pos] + name_lower + sequence[pos:]


def format_fasta(seq_id: str, description: str,
                 sequence: str, line_width: int = 80) -> str:
    """Returns a properly formatted FASTA record as a string.

    Header: '>ID description' (description omitted if empty).
    Sequence is wrapped at line_width characters.
    """
    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    # Split sequence into lines of exactly line_width characters
    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]
    return header + '\n' + '\n'.join(lines) + '\n'


def validate_positive_int(prompt: str,
                          min_val: int = 1,
                          max_val: int = 100_000) -> int:
    """Repeatedly prompts the user until a valid integer in [min_val, max_val] is given."""
    while True:
        raw = input(prompt)
        try:
            value = int(raw)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def validate_id(prompt: str) -> str:
    """Prompts until a whitespace-free, non-empty sequence ID is entered."""
    while True:
        seq_id = input(prompt).strip()
        if seq_id and not any(ch.isspace() for ch in seq_id):
            return seq_id
        print("Error: ID cannot be empty or contain whitespace.")


# 2.ADDITIONAL FEATURES

# Feature A: Motif search
def find_motifs(sequence: str, motif: str) -> list[int]:
    """Finds all (1-based) positions of motif in the clean uppercase sequence.

    Biological convention: positions start at 1.
    Overlapping matches ARE reported.
    """
    clean = ''.join(ch for ch in sequence if ch.isupper())
    positions = []
    start = 0
    while True:
        pos = clean.find(motif.upper(), start)
        if pos == -1:
            break
        positions.append(pos + 1)   # convert to 1-based
        start = pos + 1             # allow overlapping hits
    return positions


# Feature B: Complement & reverse complement
def complement(sequence: str) -> str:
    """Returns the complementary DNA strand (5'→3' direction preserved).

    Only uppercase nucleotides are complemented;
    lowercase name letters are left unchanged.
    """
    pairs = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(pairs)


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of the sequence."""
    return complement(sequence)[::-1]


# Feature C: In silico transcription (DNA → mRNA)
def transcribe(sequence: str) -> str:
    """Returns the mRNA sequence by replacing T with U (uppercase and lowercase)."""
    return sequence.replace('T', 'U').replace('t', 'u')


# Feature D: ORF identification
def find_orfs(sequence: str, min_length: int = 100) -> list[dict]:
    """Finds all ORFs (ATG → stop codon) with length >= min_length nucleotides.

    Searches all three reading frames on the sense strand.
    Returns list of dicts: {frame, start (1-based), end, length}.
    Stop codons: TAA, TAG, TGA.
    """
    stop_codons = {'TAA', 'TAG', 'TGA'}
    clean = ''.join(ch for ch in sequence if ch.isupper())
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(clean) - 2:
            codon = clean[i:i + 3]
            if codon == 'ATG':
                # Scan forward for stop codon
                for j in range(i + 3, len(clean) - 2, 3):
                    stop = clean[j:j + 3]
                    if stop in stop_codons:
                        orf_len = (j + 3) - i
                        if orf_len >= min_length:
                            orfs.append({
                                'frame': frame + 1,
                                'start': i + 1,        # 1-based
                                'end': j + 3,          # 1-based inclusive
                                'length': orf_len,
                            })
                        i = j + 3   # resume after stop codon
                        break
                else:
                    i += 3
            else:
                i += 3

    return orfs


# Feature E: Sliding-window GC analysis → CSV
def sliding_window_gc(sequence: str, window: int, step: int = 1) -> list[dict]:
    """Calculates GC content in a sliding window along the clean sequence.

    Returns list of dicts: {start_position (1-based), gc_content (%)}.
    """
    clean = ''.join(ch for ch in sequence if ch.isupper())
    results = []
    for i in range(0, len(clean) - window + 1, step):
        window_seq = clean[i:i + window]
        gc = (window_seq.count('G') + window_seq.count('C')) / window * 100
        results.append({'start_position': i + 1, 'gc_content': round(gc, 2)})
    return results


def save_gc_csv(data: list[dict], filename: str) -> None:
    """Saves sliding-window GC data to a CSV file."""
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['start_position', 'gc_content'])
        writer.writeheader()
        writer.writerows(data)


#3. MAIN

def main():
    print("=== DNA Sequence Generator (FASTA) ===\n")

    # Required inputs
    length = validate_positive_int("Enter sequence length: ")
    seq_id = validate_id("Enter sequence ID: ")
    description = input("Enter a description of the sequence (optional): ").strip()
    name = input("Enter your name: ").strip()

    #  Generate & decorate
    sequence = generate_sequence(length)
    sequence_with_name = insert_name(sequence, name)

    # Write FASTA file
    fasta_content = format_fasta(seq_id, description, sequence_with_name)
    filename = f"{seq_id}.fasta"
    with open(filename, 'w') as f:
        f.write(fasta_content)

    print(f"\nSequence saved to file: {filename}")

    # Statistics (on clean sequence, without name)
    stats = calculate_stats(sequence_with_name)
    print(f"\nSequence statistics (n={length}):")
    for nuc in ('A', 'C', 'G', 'T'):
        print(f"  {nuc}: {stats[nuc]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")

    # Preview
    print(f"\nSample contents of {filename}:")
    with open(filename) as f:
        for i, line in enumerate(f):
            print(line, end='')
            if i >= 3:   # show first few lines only
                print("  ...")
                break

    # Additional features (interactive)

    # A: Motif search
    do_motif = input("\nSearch for a motif? (y/n): ").strip().lower()
    if do_motif == 'y':
        motif = input("Enter motif (e.g. ATG): ").strip().upper()
        positions = find_motifs(sequence_with_name, motif)
        if positions:
            print(f"Motif '{motif}' found at positions (1-based): {positions}")
        else:
            print(f"Motif '{motif}' not found.")

    # B: Complement & reverse complement
    do_comp = input("Add complement & reverse complement to FASTA? (y/n): ").strip().lower()
    if do_comp == 'y':
        comp_seq = complement(sequence)          # use clean sequence for biology
        rcomp_seq = reverse_complement(sequence)
        with open(filename, 'a') as f:
            f.write(format_fasta(f"{seq_id}_complement", "Complementary strand", comp_seq))
            f.write(format_fasta(f"{seq_id}_revcomp", "Reverse complement", rcomp_seq))
        print("Complement and reverse complement appended to FASTA file.")

    # C: Transcription
    do_trans = input("Add mRNA transcription to FASTA? (y/n): ").strip().lower()
    if do_trans == 'y':
        mrna = transcribe(sequence)
        with open(filename, 'a') as f:
            f.write(format_fasta(f"{seq_id}_mRNA", "In silico transcription", mrna))
        print("mRNA record appended to FASTA file.")

    # D: ORF finding
    do_orf = input("Search for ORFs? (y/n): ").strip().lower()
    if do_orf == 'y':
        min_orf = validate_positive_int("Minimum ORF length (nt): ", min_val=3, max_val=length)
        orfs = find_orfs(sequence_with_name, min_length=min_orf)
        if orfs:
            print(f"\nFound {len(orfs)} ORF(s):")
            for o in orfs:
                print(f"  Frame {o['frame']} | Start: {o['start']} | End: {o['end']} | Length: {o['length']} nt")
        else:
            print("No ORFs found with the given minimum length.")

    # E: Sliding-window GC
    do_sw = input("Sliding-window GC analysis? (y/n): ").strip().lower()
    if do_sw == 'y':
        window = validate_positive_int("Window size (nt): ", min_val=1, max_val=length)
        step = validate_positive_int("Step size (nt): ", min_val=1, max_val=length)
        gc_data = sliding_window_gc(sequence_with_name, window, step)
        csv_file = f"{seq_id}_gc_window.csv"
        save_gc_csv(gc_data, csv_file)
        print(f"GC sliding-window data saved to: {csv_file} ({len(gc_data)} windows)")

    print("\nDone.")


if __name__ == "__main__":
    main()