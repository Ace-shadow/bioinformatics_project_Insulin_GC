import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

def main():
    # Load CSV safely
    csv_path = "results/sequence_analysis_results.csv"
    if not os.path.exists(csv_path):
        print("Error: Run analyze_sequence.py first!")
        return
    
    df = pd.read_csv(csv_path)
    
    # Clean column names and convert numeric columns
    df.columns = df.columns.str.strip()  # remove extra spaces
    numeric_cols = ["A_Count", "T_Count", "G_Count", "C_Count", "Length", "GC_Content", "AT_Content"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        else:
            # If a column is missing, fill with zeros
            df[col] = 0
    
    print("Columns after cleaning:", df.columns)
    print(df.dtypes)
    
    # Start plotting
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: GC Content
    ax1 = axes[0, 0]
    colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D"]
    bars = ax1.bar(df["Species"], df["GC_Content"], color=colors[:len(df)], edgecolor='black')
    ax1.set_ylabel("GC Content (%)", fontsize=12)
    ax1.set_title("GC Content Across Species\n(Insulin Gene)", fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 100)
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # Plot 2: Sequence Length
    ax2 = axes[0, 1]
    bars2 = ax2.bar(df["Species"], df["Length"], color=colors[:len(df)], edgecolor='black')
    ax2.set_ylabel("Sequence Length (bp)", fontsize=12)
    ax2.set_title("Sequence Length Comparison", fontsize=14, fontweight='bold')
    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 3: Nucleotide Composition (stacked bar)
    ax3 = axes[1, 0]
    nucleotides = ["A_Count", "T_Count", "G_Count", "C_Count"]
    colors_nt = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"]
    labels = ["A", "T", "G", "C"]
    bottom = [0] * len(df)
    for nt, color, label in zip(nucleotides, colors_nt, labels):
        values = df[nt].values
        ax3.bar(df["Species"], values, bottom=bottom, label=label, color=color, edgecolor='black')
        bottom = [sum(x) for x in zip(bottom, values)]
    ax3.set_ylabel("Count", fontsize=12)
    ax3.set_title("Nucleotide Composition", fontsize=14, fontweight='bold')
    ax3.legend(title="Nucleotide")
    
    # Plot 4: GC vs AT Content Pie chart (Human)
    ax4 = axes[1, 1]
    human_data = df[df["Species"].str.lower() == "human"]
    if not human_data.empty:
        human_gc = human_data["GC_Content"].values[0]
        sizes = [human_gc, 100 - human_gc]
        colors_pie = ["#FF6B6B", "#4ECDC4"]
        explode = (0.05, 0)
        ax4.pie(
            sizes,
            explode=explode,
            labels=["GC Content", "AT Content"],
            colors=colors_pie,
            autopct='%1.1f%%',
            startangle=90,
            textprops={'fontsize': 11, 'fontweight': 'bold'}
        )
        ax4.set_title("Human Insulin Gene\nNucleotide Composition", fontsize=14, fontweight='bold')
    else:
        ax4.text(0.5, 0.5, "Human data not available", ha='center', va='center', fontsize=12, fontweight='bold')
        ax4.axis('off')
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs("results", exist_ok=True)
    output_path = "results/bioinformatics_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_path}")
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    main()