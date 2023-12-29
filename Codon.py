
# Imports ----------------------------------------------------------------------
from typing import Union, List
from typing_extensions import Self

# AminoAcid --------------------------------------------------------------------

class _AminoAcidPropties():
    """Container for AA properties (private class)."""
    def __init__(self, id: int, one: str, three: str, name: str, aa_type: str):
        self.id = id
        self.one = one
        self.three = three
        self.name = name
        self.type = aa_type

class AminoAcid():
    """Amino Acid Class (stop-codon is also considered as an AminoAcid here)"""

    # AminoAcid Static Properties ----------------------------------------------

    # Amino Acids information
    _aa_list = [
        _AminoAcidPropties( 0, "*", "***", "stopcodon",     "stopcodon"),
        _AminoAcidPropties( 1, "A", "ALA", "alanine",       "hydrophobic"),
        _AminoAcidPropties( 2, "C", "CYS", "cysteine",      "special"),
        _AminoAcidPropties( 3, "D", "ASP", "aspartate",     "negative"),
        _AminoAcidPropties( 4, "E", "GLU", "glutamate",     "negative"),
        _AminoAcidPropties( 5, "F", "PHE", "phenylalanine", "hydrophobic"),
        _AminoAcidPropties( 6, "G", "GLY", "glycine",       "special"),
        _AminoAcidPropties( 7, "H", "HIS", "histidine",     "positive"),
        _AminoAcidPropties( 8, "I", "ILE", "isoleucine",    "hydrophobic"),
        _AminoAcidPropties( 9, "K", "LYS", "lysine",        "positive"),
        _AminoAcidPropties(10, "L", "LEU", "leucine",       "hydrophobic"),
        _AminoAcidPropties(11, "M", "MET", "methionine",    "hydrophobic"),
        _AminoAcidPropties(12, "N", "ASN", "asparagine",    "polar"),
        _AminoAcidPropties(13, "P", "PRO", "proline",       "special"),
        _AminoAcidPropties(14, "Q", "GLN", "glutamine",     "polar"),
        _AminoAcidPropties(15, "R", "ARG", "arginine",      "positive"),
        _AminoAcidPropties(16, "S", "SER", "serine",        "polar"),
        _AminoAcidPropties(17, "T", "THR", "thréonine",     "polar"),
        _AminoAcidPropties(18, "V", "VAL", "valine",        "hydrophobic"),
        _AminoAcidPropties(19, "W", "TRP", "tryptophane",   "hydrophobic"),
        _AminoAcidPropties(20, "Y", "TYR", "tyrosine",      "hydrophobic"),
    ]
    # Mappings
    _one_map = {aa.one: aa for aa in _aa_list}
    _three_map = {aa.three: aa for aa in _aa_list}
    _aa_map = _one_map | _three_map
    _three2one = {aa.three: aa.one for aa in _aa_list}
    _one2three = {aa.one: aa.three for aa in _aa_list}

    # Genetic code information
    _genetic_code = {
        "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S", "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
        "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C", "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
    }
    # Inverse genetic code: AA -> list of codons
    _inverse_genetic_code = {}
    for condon, aa_one in _genetic_code.items():
        if aa_one not in _inverse_genetic_code:
            _inverse_genetic_code[aa_one] = []
        _inverse_genetic_code[aa_one].append(condon)

    # Class methods ------------------------------------------------------------

    # Check existance of one/three/codon
    @classmethod
    def one_exists(cls, aa_one: str) -> bool:
        return aa_one in cls._one_map
    
    @classmethod
    def three_exists(cls, aa_three: str) -> bool:
        return aa_three in cls._three_map
    
    @classmethod
    def aa_exists(cls, aa: str) -> bool:
        return aa in cls._aa_map
    
    @classmethod
    def codon_exists(cls, codon_seq: str) -> bool:
        return codon_seq in cls._genetic_code
    
    # Mappings functions
    @classmethod
    def three2one(cls, aa_three: str) -> str:
        return cls._three2one[aa_three]
    
    @classmethod
    def one2three(cls, aa_one: str) -> str:
        return cls._one2three[aa_one]
    
    @classmethod
    def condon2one(cls, codon_seq: str) -> str:
        return cls._genetic_code[codon_seq]
    
    @classmethod
    def one2codons(cls, aa_one: str) -> List[str]:
        return cls._inverse_genetic_code[aa_one]
    
    @classmethod
    def get_all(cls, only_aminoacids: bool=False) -> List[Self]:
        """Get the list of all AminoAcids"""
        aa_list = [AminoAcid(aa_info.one) for aa_info in cls._aa_list]
        if only_aminoacids:
            aa_list = aa_list[1:]
        return aa_list

    # Construcor ---------------------------------------------------------------
    def __init__(self, aa_str: str):
        aa_info = self._aa_map.get(aa_str, None)
        assert aa_info is not None, f"ERROR in AminoAcid('{aa_str}'): invalid amino acid one/three letter-code."
        self._id = aa_info.id
        self._one = aa_info.one
        self._three = aa_info.three
        self._name = aa_info.name
        self._type = aa_info.type

    # Properties and Methods ---------------------------------------------------
    @property
    def id(self) -> int:
        return self._id
    
    @property
    def one(self) -> str:
        return self._one
    
    @property
    def three(self) -> str:
        return self._three
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def type(self) -> str:
        return self._type
    
    @property
    def is_aminoacid(self) -> bool:
        return self.id != 0
    
    @property
    def is_stop(self) -> bool:
        return self.id == 0

    def __str__(self) -> str:
        return f"Amino Acid: {self.one} ({self.three}, {self.name}) [{self.type}]"
    
    def __int__(self) -> int:
        return self.id
    
    def __eq__(self, other: Self) -> bool:
        return self.id == other.id
    
    def __hash__(self) -> int:
        return self.id
    
    def get_codons(self) -> list:
        """Returns the list of Codons that represents this AminoAcid"""
        return [Codon(condon_str) for condon_str in self.one2codons(self.one)]
    

# Condon -----------------------------------------------------------------------
class Codon():
    """Codon Class"""

    # Codon Static Properties --------------------------------------------------
    _nucleotides = ["A", "C", "T", "G"]
    _positions = [0, 1, 2]

    # Class Methods ------------------------------------------------------------
    @classmethod
    def nucleotide_exists(cls, nucleotide_str: str) -> bool:
        return nucleotide_str in cls._nucleotides
    
    @classmethod
    def position_exists(cls, position_int: int) -> bool:
        return position_int in cls._positions
    
    @staticmethod
    def codon_exists(codon_seq: str) -> bool:
        return AminoAcid.codon_exists(codon_seq)
    
    @classmethod
    def get_all(cls, only_aminoacids: bool=False) -> List[Self]:
        """Get the list of all Codons"""
        condons_list = [Codon(condon_seq) for condon_seq in AminoAcid._genetic_code.keys()]
        if only_aminoacids:
            condons_list = [codon for codon in condons_list if codon.is_aminoacid]
        return condons_list

    # Constructor --------------------------------------------------------------
    def __init__(self, codon_seq: str):
        aa_one = AminoAcid._genetic_code.get(codon_seq, None)
        assert aa_one is not None, f"ERROR in Codon('{codon_seq}'): invalid condon input."
        self._seq = codon_seq
        self._aminoacid = AminoAcid(aa_one)

    # Properties and Methods ---------------------------------------------------
    @property
    def seq(self) -> str:
        return self._seq
    
    @property
    def aminoacid(self) -> AminoAcid:
        return self._aminoacid
    
    @property
    def is_aminoacid(self) -> bool:
        return self.aminoacid.is_aminoacid
    
    @property
    def is_stop(self) -> bool:
        return self.aminoacid.is_stop

    def __str__(self) -> str:
        return f"Codon('{self.seq}') ({self.aminoacid.one}, {self.aminoacid.three})"
    
    def __int__(self) -> int:
        return int(self.aminoacid)
    
    def __eq__(self, other: Self) -> bool:
        return self.seq == other.seq
    
    def __hash__(self) -> int:
        return hash(self.seq)
    
    @property
    def _errorstr(self) -> str:
        return f"ERROR in Codon('{self.seq}')"
    
    # Mutations Methods --------------------------------------------------------
    def get_mutant_nucleotides(self, position: int) -> List[str]:
        """Returns the list of nucleotides that are different from wild-type nucleotide in a given position."""
        assert self.position_exists(position), f"{self._errorstr}.get_mutant_nucleotides(), invalid position={position}."
        wildtype_nucleotide = self.seq[position]
        return [nucleotide for nucleotide in self._nucleotides if nucleotide != wildtype_nucleotide]

    def get_codon_with_mutant_nucleotide(self, mutant_nucleotide: str, position: int) -> Self:
        """Returns the Codon in which we have mutated one nucleotide"""
        assert self.nucleotide_exists(mutant_nucleotide), f"{self._errorstr}.get_codon_with_mutant_nucleotide(): invalid mutant_nucleotide='{mutant_nucleotide}'."
        assert self.position_exists(position), f"{self._errorstr}.get_codon_with_mutant_nucleotide(): invalid position='{position}'."
        assert mutant_nucleotide != self.seq[position], f"{self._errorstr}.get_codon_with_mutant_nucleotide(): mutant_nucleotide='{mutant_nucleotide}' is same as wild-type at position {position}."
        mutant_codon_seq = self.seq[:position] + mutant_nucleotide + self.seq[position+1:]
        return Codon(mutant_codon_seq)

    def get_synonymous_codons(self) -> List[Self]:
        """Returns the list of all synonymous Codons"""
        return self.aminoacid.get_codons()
    
    def get_mutant_codons(self, n_mutations: int=1, positions: Union[None, List[int]]=None) -> List[Self]:
        """Returns the list of all Codons that are separated by n_mutations at a given position (if None, we consider all possible positions)."""

        # Guardians
        assert n_mutations in [1, 2, 3], f"{self._errorstr}.get_mutant_codons(): n_mutations={n_mutations} should be 1, 2 or 3."
        if positions is not None:
            assert len(positions) == n_mutations, f"{self._errorstr}.get_mutant_codons(): positions={positions} length should match with n_mutations={n_mutations}."
            for position in positions:
                assert self.position_exists(position), f"{self._errorstr}.get_mutant_codons(): invalid position={position}."
            assert len(positions) == len(set(positions)), f"{self._errorstr}.get_mutant_codons(): positions={positions} contain some repetitions."
        
        # Case of 1 mutation
        if n_mutations == 1:
            if positions is None:
                return self.get_mutant_codons(1, [0]) + self.get_mutant_codons(1, [1]) + self.get_mutant_codons(1, [2])
            position = positions[0]
            mutant_codons = []
            for mutant_nucleotide in self.get_mutant_nucleotides(position):
                codon = self.get_codon_with_mutant_nucleotide(mutant_nucleotide, position)
                mutant_codons.append(codon)
            return mutant_codons
        
        # Case of 2 mutations
        elif n_mutations == 2:
            if positions is None:
                return self.get_mutant_codons(2, [0, 1]) + self.get_mutant_codons(2, [0, 2]) + self.get_mutant_codons(2, [1, 2])
            position1, position2 = positions
            mutant_codons = []
            for mutant_nucleotide1 in self.get_mutant_nucleotides(position1):
                for mutant_nucleotide2 in self.get_mutant_nucleotides(position2):
                    codon = self.get_codon_with_mutant_nucleotide(mutant_nucleotide1, position1) \
                        .get_codon_with_mutant_nucleotide(mutant_nucleotide2, position2)
                    mutant_codons.append(codon)
            return mutant_codons
        
        # Case of 3 mutations
        else:
            position1, position2, position3 = 0, 1, 2
            mutant_codons = []
            for mutant_nucleotide1 in self.get_mutant_nucleotides(position1):
                for mutant_nucleotide2 in self.get_mutant_nucleotides(position2):
                    for mutant_nucleotide3 in self.get_mutant_nucleotides(position3):
                        codon = self.get_codon_with_mutant_nucleotide(mutant_nucleotide1, position1) \
                            .get_codon_with_mutant_nucleotide(mutant_nucleotide2, position2) \
                            .get_codon_with_mutant_nucleotide(mutant_nucleotide3, position3)
                        mutant_codons.append(codon)
            return mutant_codons
