import os
import requests
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

# Standard Genetic Code Dictionary
CODON_TABLE = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
}

AMINO_ACID_NAMES = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic Acid', 'E': 'Glutamic Acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine (Start)', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine',
    '_': 'Stop Codon'
}

def detect_sequence_type(seq):
    seq = seq.upper()
    valid_chars = set('ACGTU')
    if not set(seq).issubset(valid_chars):
        return "Invalid", "Sequence contains invalid characters. Only A, C, G, T, and U are allowed."
    if 'T' in seq and 'U' in seq:
        return "Invalid", "Sequence contains both T and U, which is biologically invalid for a single strand."
    if 'U' in seq:
        return "RNA", "RNA detected because the sequence contains Uracil (U)."
    if 'T' in seq:
        return "DNA", "DNA detected because the sequence contains Thymine (T)."
    return "DNA", "Sequence contains only A, C, G. Defaulting to DNA."

def transcribe_dna(seq, strand_type):
    seq = seq.upper()
    if strand_type == 'non-template':
        # Coding strand: just replace T with U
        return seq.replace('T', 'U')
    else:
        # Template strand: reverse complement then T->U 
        # (A->U, T->A, C->G, G->C)
        complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in seq)

def translate_mrna(mrna):
    codons = [mrna[i:i+3] for i in range(0, len(mrna) - len(mrna)%3, 3)]
    amino_acids = []
    polypeptide = ""
    for codon in codons:
        if len(codon) == 3:
            aa_abbr = CODON_TABLE.get(codon, '?')
            aa_name = AMINO_ACID_NAMES.get(aa_abbr, 'Unknown')
            amino_acids.append({'codon': codon, 'abbr': aa_abbr, 'name': aa_name})
            if aa_abbr == '_':
                break # Stop codon halts translation
            polypeptide += aa_abbr
    return amino_acids, polypeptide

def fetch_uniprot_data(protein_sequence):
    if not protein_sequence or len(protein_sequence) < 5:
        return {"error": "Sequence too short for a reliable database match."}
    
    url = f"https://rest.uniprot.org/uniprotkb/search?query=sequence:{protein_sequence}&size=1&format=json"
    try:
        response = requests.get(url, timeout=5)
        data = response.json()
        if data.get('results'):
            result = data['results'][0]
            name = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown Protein')
            organism = result.get('organism', {}).get('scientificName', 'Unknown Organism')
            return {"name": name, "organism": organism, "match": "Found"}
        return {"error": "No exact match found in UniProt for this specific polypeptide sequence."}
    except Exception as e:
        return {"error": "Could not connect to the UniProt database."}

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    data = request.json
    sequence = data.get('sequence', '').strip().upper().replace(" ", "").replace("\n", "")
    strand_type = data.get('strand_type', 'non-template')

    # 1. Detection
    seq_type, det_explanation = detect_sequence_type(sequence)
    if seq_type == "Invalid":
        return jsonify({"error": det_explanation})

    # 2. Transcription
    mrna = sequence if seq_type == 'RNA' else transcribe_dna(sequence, strand_type)
    
    # 3 & 4. Translation & Amino Acids
    amino_acids, polypeptide = translate_mrna(mrna)

    # 5. Protein & Database Lookup
    db_result = fetch_uniprot_data(polypeptide)

    return jsonify({
        "input_sequence": sequence,
        "type": seq_type,
        "detection_exp": "Biological sequences are made of specific nucleotides. DNA uses A, C, G, and Thymine (T), while RNA uses Uracil (U) instead of Thymine. " + det_explanation,
        "mrna": mrna,
        "transcription_exp": "Transcription is the cellular process where DNA is copied into messenger RNA (mRNA). If a template strand was provided, we created its exact chemical complement. If a non-template strand was provided, we simply replaced Thymine (T) with Uracil (U).",
        "amino_acids": amino_acids,
        "translation_exp": "Translation is the process where cellular ribosomes read mRNA in groups of three letters, called 'codons'. Each codon acts like a specific word that translates into a single amino acid.",
        "polypeptide": polypeptide,
        "polypeptide_exp": "Amino acids are the core building blocks of life. When linked together in a specific order, they form a 'polypeptide chain'. This chain is the raw material that will eventually fold into a functional protein.",
        "db_result": db_result,
        "protein_exp": "A protein is a complex 3D structure made from folded polypeptide chains that does the physical work inside cells. We queried the UniProt global biological database to see if your newly generated sequence matches any known proteins in nature."
    })

if __name__ == '__main__':
    app.run(debug=True)