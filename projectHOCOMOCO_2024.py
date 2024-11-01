from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO

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

# Definisci il percorso del file Transfac e carica il modello PWM
transfac_file = "BCL6.H12CORE.0.PSM.A_transfac_format.txt"  # File Transfac nella directory del progetto
pwm = load_pwm(transfac_file)

# Definisci la sequenza di DNA da analizzare
#dna_sequence = Seq("AGCTTAGCTTTCAGGAATTAGGCTTAGGCTT")# Sostituisci con la sequenza da analizzare
def load_sequence(dna_file_path):
    with open(dna_file_path) as f:
        fasta_sequences = SeqIO.parse(f, "fasta")
        for fasta in fasta_sequences:
            return str(fasta.seq)

# Carica la sequenza da analizzare
dna_file_path = "Homo_Sapiens_chromosome_3_BCL6_protein_gene_exon1_sequence.fasta"  # File FASTA nella directory del progetto
dna_sequence = load_sequence(dna_file_path)

# Cerca i siti di legame nella sequenza
threshold = 0.8  # Regola il threshold per il punteggio
binding_sites = find_binding_sites(dna_sequence, pwm, threshold)

# Stampa i risultati
if binding_sites:
    print("Siti di legame trovati:")
    for position, score in binding_sites:
        print(f"Posizione: {position}, Punteggio: {score:.2f}")
else:
    print("Nessun sito di legame trovato sopra la soglia.")
