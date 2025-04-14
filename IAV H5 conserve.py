from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===== STEP 1: Load Aligned FASTA =====
fasta_path = "/Users/chenchi/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/influenza/np alignment/np_H7N9_H5N6_H9N2_H3N8_H6N1sequences_20250414_4010632.fasta"  # <-- Change if needed
fasta_lines = Path(fasta_path).read_text().splitlines()

# ===== STEP 2: Parse Aligned Sequences =====
sequences = []
current_seq = []

for line in fasta_lines:
    if line.startswith(">"):
        if current_seq:
            sequences.append(''.join(current_seq))
            current_seq = []
    else:
        current_seq.append(line.strip())
if current_seq:
    sequences.append(''.join(current_seq))

# ===== STEP 3: Convert to Matrix =====
max_len = max(len(seq) for seq in sequences)
aligned_array = np.array([list(seq.ljust(max_len, '-')) for seq in sequences])

# ===== STEP 4: Calculate Percent Identity Per Position =====
conservation_scores = []
for i in range(aligned_array.shape[1]):
    column = aligned_array[:, i]
    column = column[column != '-']  # exclude gaps
    if len(column) == 0:
        conservation_scores.append(0)
    else:
        most_common = max(set(column), key=list(column).count)
        score = list(column).count(most_common) / len(column)
        conservation_scores.append(score)

# ===== STEP 5: Save Table =====
df = pd.DataFrame({
    "Position": np.arange(1, len(conservation_scores) + 1),
    "Conservation_Score": conservation_scores
})
df.to_csv("NCBI H7N9_H5N6_H9N2_H3N8_H6N1 NP 10 years H5 sequences.csv", index=False)

# ===== STEP 6: Plot Conservation =====
plt.figure(figsize=(12, 4))
plt.plot(conservation_scores, color='red', linewidth=1)
plt.title("NCBI H7N9_H5N6_H9N2_H3N8_H6N1 NP 10 years H5 sequences")
plt.xlabel("Amino Acid Position")
plt.ylabel("Conservation Score")
plt.ylim(0, 1.05)
plt.grid(True)
plt.tight_layout()
plt.show()