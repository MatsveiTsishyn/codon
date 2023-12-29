
# Codon

A small and simple module to **manage** and **mutate** **amino acids** and **codons** with `python3`. \
Manage amino acids one-letter code vs. three-letters code expressions. \
Manage codon translation to amino acid by standard genetic code and synonymous codons.

## Usage examples

Code to list all amino acids that are accessible by **one nucleotide mutation from alanine**.

```python
# Imports
from Codon import AminoAcid, Codon

# Define alanine and all its possible codons
ala = AminoAcid("A")
ala_codons = ala.get_codons() # [Codon("GCA"), Codon("GCC"), ...], len=4

# List all possible codons accessible by one nucleotide mutation (from any alanine codon)
mutant_codons = []
for ala_codon in ala_codons:
    mutant_codons += ala_codon.get_mutant_codons(n_mutations=1)
mutant_codons = list(set(mutant_codons)) # Clear Codons repetitions, len=28

# List all possible amino acids represented by those codons
mutant_aminoacids = [codon.aminoacid for codon in mutant_codons] # Codon -> AminoAcid, len=28
mutant_aminoacids = list(set(mutant_aminoacids)) # Clear AminoAcid repetitions, len=8
```

More usage examples can be found in `usage_examples.py`.