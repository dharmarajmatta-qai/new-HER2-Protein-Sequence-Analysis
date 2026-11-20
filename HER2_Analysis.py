
import requests
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from io import StringIO
import os

# Create project folder
project_folder = "Project01_HER2"
os.makedirs(project_folder, exist_ok=True)

# Download HER2 sequence
url = "https://rest.uniprot.org/uniprotkb/P04626.fasta"
response = requests.get(url, timeout=30)
response.raise_for_status()

# Read protein sequence
record = SeqIO.read(StringIO(response.text), "fasta")
sequence = str(record.seq)

# Calculate molecular weight
analysis = ProteinAnalysis(sequence)
molecular_weight = analysis.molecular_weight()

# Count amino acids
counts = analysis.count_amino_acids()

# Create results table
amino_df = pd.DataFrame(
    list(counts.items()),
    columns=["Amino Acid", "Count"]
)

amino_df["Percentage"] = (
    amino_df["Count"] / len(sequence)
) * 100

amino_df = amino_df.sort_values(
    by="Count", ascending=False
)

# Save CSV
amino_df.to_csv(
    os.path.join(project_folder, "HER2_Amino_Acid_Composition.csv"),
    index=False
)

# Save FASTA
with open(
    os.path.join(project_folder, "HER2_Protein_Sequence.fasta"), "w"
) as file:
    file.write(f">{record.id}\n")
    file.write(sequence + "\n")

# Generate and save graph
plt.figure(figsize=(12, 5))
plt.bar(amino_df["Amino Acid"], amino_df["Percentage"])
plt.xlabel("Amino Acid")
plt.ylabel("Percentage (%)")
plt.title("HER2 Protein Amino Acid Composition")
plt.tight_layout()

plt.savefig(
    os.path.join(project_folder, "HER2_Amino_Acid_Graph.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Protein ID:", record.id)
print("Sequence length:", len(sequence))
print("Molecular weight (Da):", molecular_weight)
print("Molecular weight (kDa):", molecular_weight / 1000)
print("Project completed successfully.")
