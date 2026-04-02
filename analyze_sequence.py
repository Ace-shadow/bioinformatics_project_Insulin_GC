from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
import os

def analyze_sequence(filepath, species_name):
    """Analyze a DNA sequence and return metrics"""
    record = SeqIO.read(filepath, "fasta")
    sequence = record.seq
    
    analysis = {
        "Species": species_name,
        "Length": len(sequence),
        "GC_Content": gc_fraction(sequence) * 100,
        "AT_Content": (1 - gc_fraction(sequence)) * 100,
        "G_count": sequence.count("G"),
        "C_Count": sequence.count("C"),
        "A_Count": sequence.count("A"),
        "T_Count": sequence.count("T"),
        "Description": record.description[:50] + "..."
    }
    return analysis

def main():
    files = [
    ("data/human.sequence.fasta", "Human"),
    ("data/mouse.sequence.fasta", "Mouse"),
    ("data/zebrafish.sequence.fasta", "Zebrafish")
]
    
    results = []
    for filepath, species in files:
        if os.path.exists(filepath):
            print(f"Analyzing {species}...")
            results.append(analyze_sequence(filepath, species))
        else:
            print(f"Warning: {filepath} not found. Run fetch_sequences.py first!")
    
    if results:
        df = pd.DataFrame(results)
        print("\n=== SEQUENCE ANALYSIS RESULTS ===")
        print(df.to_string(index=False))
        
        # Save results
        os.makedirs("results", exist_ok=True)
        output_file = "results/sequence_analysis_results.csv"
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        return df
    else:
        print("No data to analyze!")
        return None

if __name__ == "__main__":
    main()