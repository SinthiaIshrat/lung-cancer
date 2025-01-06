import difflib
import os
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def load_reference_sequence(filepath):
    """
    Load the reference sequence from a file.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
    # Ignore the header and join the rest into a single sequence
    reference_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return reference_sequence

def calculate_similarity(seq1, seq2):
    """
    Calculate the similarity percentage between two DNA sequences using SequenceMatcher.
    """
    matcher = difflib.SequenceMatcher(None, seq1, seq2)
    return matcher.ratio() * 100

def perform_global_alignment(seq1, seq2):
    """
    Perform global alignment between two sequences using Needleman-Wunsch algorithm.
    """
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]  # Take the best alignment
    alignment_score = best_alignment[2]
    return best_alignment, alignment_score

def perform_local_alignment(seq1, seq2):
    """
    Perform local alignment between two sequences using Smith-Waterman algorithm.
    """
    alignments = pairwise2.align.localxx(seq1, seq2)
    best_alignment = alignments[0]  # Take the best alignment
    alignment_score = best_alignment[2]
    return best_alignment, alignment_score

def check_infection(reference, user_sequence):
    """
    Compare the user sequence with the reference to determine infection status.
    """
    similarity = calculate_similarity(reference, user_sequence)
    threshold = 90.0  # Infection threshold, set at 90% similarity
    is_infected = similarity < threshold
    return is_infected, similarity

def main():
    print("Human Lung Cancer-Associated Polyomavirus Detection")
    print("==================================================")
    
    # Load the reference genome sequence
    filepath = "cancer_lung.fna"
    if not os.path.exists(filepath):
        print(f"Error: The file '{filepath}' was not found in the current directory.")
        print("Please ensure the file is available and try again.")
        return

    try:
        reference_sequence = load_reference_sequence(filepath)
        print("Reference genome loaded successfully!")
    except Exception as e:
        print(f"Error: Could not load the reference genome. Details: {e}")
        return

    # Get the user's DNA sequence
    user_sequence = input("Enter the DNA sequence to analyze: ").strip().upper()

    # Validate the DNA sequence
    if not all(base in "ATGC" for base in user_sequence):
        print("Error: Invalid DNA sequence. Please only use A, T, G, and C bases.")
        return

    # Check for infection
    is_infected, similarity = check_infection(reference_sequence, user_sequence)

    # Perform global and local alignment
    global_alignment, global_score = perform_global_alignment(reference_sequence, user_sequence)
    local_alignment, local_score = perform_local_alignment(reference_sequence, user_sequence)

    # Display results
    print("\nAnalysis Results:")
    print(f"Similarity with reference genome: {similarity:.2f}%")
    if is_infected:
        print("Status: Infected (Significant variation detected)")
    else:
        print("Status: Not Infected (Sequence matches reference genome closely)")

    # Display global alignment results
    print("\nGlobal Alignment (Needleman-Wunsch):")
    print(format_alignment(*global_alignment))
    print(f"Global Alignment Score: {global_score:.2f}")

    # Display local alignment results
    print("\nLocal Alignment (Smith-Waterman):")
    print(format_alignment(*local_alignment))
    print(f"Local Alignment Score: {local_score:.2f}")

if __name__ == "__main__":
    main()
