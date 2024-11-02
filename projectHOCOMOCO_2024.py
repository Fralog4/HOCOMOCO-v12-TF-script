import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import motifs


def plot_score_distribution(binding_sites):
    scores = [score for _, score in binding_sites]
    plt.figure(figsize=(8, 5))
    plt.hist(scores, bins=10, color="grey", edgecolor="black")
    plt.xlabel("Punteggio PWM")
    plt.ylabel("Frequenza")
    plt.title("Distribuzione dei punteggi dei motivi di legame trovati")
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
    dna_file_path = "Homo_Sapiens_chromosome_3_BCL6_protein_gene_exon1_sequence.fasta"
    pwm_file_path = "BCL6.H12CORE.0.PSM.A_transfac_format.txt"
    dna_sequence = load_sequence(dna_file_path)
    pwm = load_pwm(pwm_file_path)
    threshold = 0.8
    binding_sites = find_binding_sites(dna_sequence, pwm, threshold)
    plot_score_distribution(binding_sites)
    if binding_sites:
        for position, score in binding_sites:
            print(f"Motif founded into the DNA sequence at position {position} with a similarity score of {score:.2f}")
    else:
        print("No binding sites found above the threshold, No binding sites found.")


if __name__ == "__main__":
    main()
