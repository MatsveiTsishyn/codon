
# Imports ----------------------------------------------------------------------
from Codon import AminoAcid, Codon



# Use AminoAcid class ----------------------------------------------------------

# Initialize an Amino Acid
ala = AminoAcid("A") # Use upper-case one- or three-letters code for standard proteinogenic amino acids
stop_aa = AminoAcid("*") # Stop codon translation product can also be an instance of AminoAcid

# Access base AminoAcid properties
ala.one # "A"
ala.three # "ALA"
ala.type # "hydrophobic"
ala.name # "alanine"
ala.id # 1
ala.is_aminoacid # True
ala.is_stop # False

# List all AminoAcids
# Set only_aminoacids=True to exclude stop codon AminoAcid
all_aa_list = AminoAcid.get_all(only_aminoacids=True) # [AminoAcid("A"), AminoAcid("C"), ...], len=20

# Get all possible codons of an AminoAcid
gly_codons_list = AminoAcid("G").get_codons() # [Codon("GGA"), Codon("GGC"), ...], len=4

# Use AminoAcid class as mapping: one-letter-code <-> three-letters-code
one_str = "H"
three_str = AminoAcid.one2three(one_str) # "HIS"
one_copy_str = AminoAcid.three2one(three_str) # "H"

# AminoAcid instances are hashable
aa_set = set([AminoAcid("A"), AminoAcid("A"), AminoAcid("A"), AminoAcid("C")]) # {AminoAcid("A"), AminoAcid("C")}, len=2



# Use Codon class --------------------------------------------------------------

# Initialize a Codon
codon = Codon("AAT") # Use 3-nucleotides code

# Access base Codon properties
codon.seq # "AAT"
codon.aminoacid # Codon('ATT') (I, ILE)
codon.is_aminoacid # True
codon.is_stop # False

# List all Codons
# Set only_aminoacids=True to exclude stop codons
all_codons_list = Codon.get_all(only_aminoacids=False) # [Codon("AAA"), Codon("AAC"), ...], len=64

# Codons instances are hashable
codons_set = set([Codon("AAT"), Codon("AAT"), Codon("AGT")]) # {Codon("AAT"), Codon("AGT")}, len=2

# List synonimous codons
synonymous_codons_list = Codon("TTA").get_synonymous_codons() # [Codon("CTA"), Codon("CTC"), ...], len=6



# Explore mutations relations between Codons -----------------------------------

# Get mutant codon by specifying a nucleotide mutation
# Nucleotide positions are indexed by 0, 1 and 2
codon_aat = Codon("AAT")
codon_aag = codon_aat.get_codon_with_mutant_nucleotide("G", 2) # Codon("AAG")

# List mutant codons by exactly 1 mutation (at any nucleotide position)
mutant_codons_list = codon_aat.get_mutant_codons(n_mutations=1) # [Codon("CAT"), Codon("TAT"), ...], len=9

# List mutant codons by exactly 1 mutation at a given nucleotide position
mutant_codons_list = codon_aat.get_mutant_codons(n_mutations=1, positions=[0]) # [Codon("CAT"), Codon("TAT"), Codon("GAT")], len=3

# List mutant codons by exactly 2 mutations (at any pair of nucleotide positions)
mutant_codons_list = codon_aat.get_mutant_codons(n_mutations=2) # [Codon("CCT"), Codon("CTT"), ...], len=27

# List mutant codons by exactly 2 mutations at a given pair of nucleotide positions
mutant_codons_list = codon_aat.get_mutant_codons(n_mutations=2, positions=[0, 2]) # [Codon("CAA"), Codon("CAC"), ...], len=9

# List mutant codons by exactly 3 mutations (at all nucleotide positions)
mutant_codons_list = codon_aat.get_mutant_codons(n_mutations=3) # [Codon("CCA"), Codon("CCC"), ...], len=27



# Example of application -------------------------------------------------------
# Code to list all amino acids that are accessible by exactly 1 nucleotide mutation from alanine

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
