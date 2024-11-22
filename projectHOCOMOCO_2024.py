import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import motifs


def plot_binding_sites(sequence, binding_sites):
    """
    Plots the DNA sequence highlighting the binding sites.

    Args:
        sequence (str): The DNA sequence to analyze.
        binding_sites (list[tuple[int, float]]): List of binding site positions and scores.
    """
    sequence_length = len(sequence)
    positions = [position for position, _ in binding_sites]
    scores = [score for _, score in binding_sites]
    plt.figure(figsize=(12, 4))
    plt.plot(range(sequence_length), [0] * sequence_length, color="lightgrey", lw=0.5,
             label="DNA sequence")  # Sequenza come linea piatta
    plt.scatter(positions, [0] * len(positions), c=scores, cmap="viridis", edgecolor="black", s=100,
                label="Binding sites")
    for pos, score in zip(positions, scores):
        plt.text(pos, 0.1, f"Similarity score: {score:.2f}", fontsize=8, ha="center", rotation=90, color="red")

    plt.xlabel("Position in DNA sequence")
    plt.ylabel("Binding Sites")
    plt.title("Highlighted Binding Sites on the DNA Sequence")
    plt.colorbar(label="PWM score", orientation="vertical")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


def load_pwm(file_path):
    with open(file_path) as f:
        motifs_list = motifs.parse(f, "transfac")  # Reads the motifs from the file in Transfac format
        pwm = motifs_list[0].pssm  # Take the first one and convert it to PWM
    return pwm


def find_binding_sites(sequence, pwm, threshold=0.8):
    """
    Finds the binding sites of a given PWM in a DNA sequence.

    Args:
        sequence (str): The DNA sequence to search for binding sites.
        pwm (motifs.PSSM): The PWM to use for searching.
        threshold (float, optional): The minimum score required for a match. Defaults to 0.8.

    Returns:
        list[tuple[int, float]]: A list of tuples, where the first element of each tuple is the position of the binding site, and the second element is the score of the match.
    """
    matches = []
    for position, score in pwm.search(sequence, threshold):
        matches.append((position, score))
    return matches


def load_sequence(dna_file_path):
    with open(dna_file_path) as f:
        fasta_sequences = SeqIO.parse(f, "fasta")
        for fasta in fasta_sequences:
            return str(fasta.seq)


def main():
    dna_file_path = "BCL6-FASTA-ENSEMBL-Promotore.fasta"
    pwm_file_path = "BCL6.H12CORE.0.PSM.A_transfac_format.txt"
    dna_sequence = load_sequence(dna_file_path)
    pwm = load_pwm(pwm_file_path)
    threshold = 0.9
    binding_sites = find_binding_sites(dna_sequence, pwm, threshold)

    # Visualizzazione dei siti di legame
    if binding_sites:
        plot_binding_sites(dna_sequence, binding_sites)
        for position, score in binding_sites:
            print(f"Motif founded into the DNA sequence at position {position} with a similarity score of {score:.2f}")
    else:
        print("No binding sites found above the threshold, No binding sites found.")


if __name__ == "__main__":
    main()
