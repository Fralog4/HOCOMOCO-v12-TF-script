from Bio import motifs
from Bio import SeqIO
import matplotlib.pyplot as plt

def plot_binding_sites(sequence, pwm, threshold):
    scores = [pwm.calculate(sequence[i:i+len(pwm)]) for i in range(len(sequence) - len(pwm) + 1)]
    positions = list(range(len(scores)))

    plt.figure(figsize=(10, 6))
    plt.plot(positions, scores, label="Punteggio PWM", color="blue")
    plt.axhline(y=threshold, color="red", linestyle="--", label="Threshold")
    plt.xlabel("Posizione nella Sequenza")
    plt.ylabel("Punteggio PWM")
    plt.legend()
    plt.title("Punteggi PWM lungo la sequenza")
    plt.show()

def plot_score_distribution(binding_sites):
    scores = [score for _, score in binding_sites]
    plt.figure(figsize=(8, 5))
    plt.hist(scores, bins=10, color="grey", edgecolor="black")
    plt.xlabel("Punteggio PWM")
    plt.ylabel("Frequenza")
    plt.title("Distribuzione dei punteggi dei motivi di legame trovati")
    plt.show()

# Funzione per caricare PWM da un file in formato Trasfac
def load_pwm(file_path):
    with open(file_path) as f:
        motifs_list = motifs.parse(f, "transfac")  # Leggi tutti i motivi dal file Transfac
        pwm = motifs_list[0].pssm  # Prendi il primo motivo, convertito in Position-Specific Scoring Matrix (PSSM)
    return pwm

# Funzione per cercare siti di legame in una sequenza
def find_binding_sites(sequence, pwm, threshold=0.8):
    matches = []
    for position, score in pwm.search(sequence, threshold):
        matches.append((position, score))
    return matches

# Definisci la sequenza di DNA da analizzare
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
    #plot_binding_sites(dna_sequence, pwm, threshold)
    plot_score_distribution(binding_sites)
    if binding_sites:
        print("Siti di legame trovati:")
        for position, score in binding_sites:
            print(f"Motivo trovato nella sequenza di DNA del gene alla posizione {position} con un punteggio di similarità di {score:.2f}")
    else:
        print("Nessun sito di legame trovato sopra la soglia, Nessun motivo di legame trovato.")

if __name__ == "__main__":
    main()