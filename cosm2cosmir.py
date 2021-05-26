"""Simulate radiation damage to a COSM structure.

This module simulates ionizing radiation-induced DNA damage in the context of DNA
origami represented by COSM model. It does so by converting a COSM pdb file
into a COSM_IR pdb file which represents the damaged structure.

Attributes:
    END_ATOMS: Names of the atoms that mark start of a new chain.
    TERMINAL_ATOMS: Names of the atoms that mark start/end of a DS region.
    ATOMS: All valid names of COSM atoms.
    DAMAGED_ATOM_NAMES: Templates for special atoms' names (needed to apply damage).
    DAMAGED_ATOMS: Names of the special atoms.
    CHAIN_NAMES: Chain names to use when writing output pdb file.
"""

import numpy as np
import argparse
import dataclasses
import string
import damage_probabilities
import itertools

from collections import namedtuple, defaultdict


# TODO: parser
parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', required=True,
#                     action='store', dest='pdb',
#                     help='input cosm pdb file')
# parser.add_argument('-o', '--output', required=True,
#                     action='store', dest='output',
#                     help='output pdb file')
# parser.add_argument('-r', '--restr', required=True,
#                     action='store', dest='restr',
#                     help='input cosm distance/angle restraints')
# parser.add_argument('-s', '--seq', required=True,
#                     action='store', dest='seq',
#                     help='scaffold sequence')
# parser.add_argument('-m', '--map', required=True,
#                     action='store', dest='map',
#                     help='2d coords map file')
# parser.add_argument('-d', '--dose',
#                     action='store', dest='dose',
#                     help='radiation dose (in Gy)')
args = parser.parse_args()

args.pdb = "example/out/example.pdb"
args.output = "example/example_ir.pdb"
args.restr = "example/out/example.r"
# args.restr_out = "example/example_ir.r"
args.sequence = "example/inp/example_scaffold.txt"
args.map = "example/out/example.map"
args.dose = 100


END_ATOMS = {
    "ST",
    "TT",
    "T1T",
    "T2T",
    "T3T",
    "T4T",
    "T5T",
    "T6T",
    "T7T",
    "B1T",
    "B2T",
    "B3T",
    "B4T",
    "B5T",
    "B6T",
    "B7T",
    "BT",
}

TERMINAL_ATOMS = {
    "T",
    "T1",
    "T2",
    "T3",
    "T4",
    "T5",
    "T6",
    "T7",
    "TT",
    "T1T",
    "T2T",
    "T3T",
    "T4T",
    "T5T",
    "T6T",
    "T7T",
    "B",
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "BT",
    "B1T",
    "B2T",
    "B3T",
    "B4T",
    "B5T",
    "B6T",
    "B7T",
}

ATOMS = TERMINAL_ATOMS | {"S", "ST", "N", "H", "PT"}
DAMAGED_ATOM_NAMES = {
    "ssb": "SSB",
    "rev_ssb": "SSB",
    "abasic": "AP",
    "rev_abasic": "AP",
    "base": {"A": ["DAO", "DAF"],
             "G": ["DGO", "DGF"],
             "T": ["DT"],
             "C": [],
             },
    "rev_base": {"T": ["DAO", "DAF"],
                 "C": ["DGO", "DGF"],
                 "A": ["DT"],
                 "G": [],
                 },
}
DAMAGED_ATOMS = {"SSB", "AP", "DT", "DAO", "DAF", "DGO", "DGF"}

CHAIN_NAMES = itertools.chain(
    string.ascii_uppercase, string.ascii_lowercase, string.digits, itertools.repeat("_")
)


Restraint = namedtuple("Restraint", ["type", "atom", "parameters"])
MapCoord = namedtuple("MapCoord", ["row", "column"])


@dataclasses.dataclass
class NucleotideAtoms:
    prev: "Atom"
    current: "Atom"
    next: "Atom"


@dataclasses.dataclass
class Atom:
    resn: str
    coords: np.array
    nucleotide: "Nucleotide" = dataclasses.field(default=None)
    restrs: str = dataclasses.field(default_factory=list)
    mapping: "MapCoord" = dataclasses.field(default=None)

    atomid: int = dataclasses.field(init=False)
    chain: str = dataclasses.field(init=False)

    def get_mode(self):
        """Determine the orientation of a terminal atom."""
        assert self.resn in TERMINAL_ATOMS
        assert self.nucleotide is not None

        if (
            prev_atom := self.nucleotide.atoms.prev,
            next_atom := self.nucleotide.atoms.next,
        ) == (None, None):
            return "0"
        if prev_atom is None:
            assert next_atom not in ["S", "ST"]
            return "0 +1"
        if next_atom is None:
            assert prev_atom is not None
            assert prev_atom not in ["S", "ST"]
            return "-1 0"
        if prev_atom.resn in ["S", "ST"]:
            assert next_atom.resn not in ["S", "ST"]
            return "0 +1"
        if next_atom.resn in ["S", "ST"]:
            assert prev_atom.resn not in ["S", "ST"]
            return "-1 0"
        if self.mapping.row == prev_atom.mapping.row:
            assert self.mapping.row != next_atom.mapping.row
            return "-1 0"
        assert self.mapping.row == next_atom.mapping.row
        return "0 +1"

    def _generate_pdbline(self):
        name = "D"
        altloc = ""
        resi = self.atomid
        insert = ""
        occupancy = 1
        temperature = 0
        segi = ""
        symbol = "D"
        charge = ""

        pdbline = (
            f"ATOM  {self.atomid:>5} {name:^4}{altloc:>1}{self.resn:>3} {self.chain:>1}"
            f"{resi:>4}{insert:>1}   {self.coords[0]:>8.3f}"
            f"{self.coords[1]:>8.3f}{self.coords[2]:>8.3f}{occupancy:>6.2f}"
            f"{temperature:>6.2f}      {segi:<4}{symbol:>2}{charge:>2}\n"
        )
        return pdbline

    @classmethod
    def parse_pdbline(cls, line):
        """Create an Atom object from the pdb line."""
        atomid = line[6:11].strip()
        resn = line[17:20].strip()
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()

        atom = cls(resn, coords=np.array([x, y, z], dtype=float))
        atom.atomid = int(atomid)
        return atom

    @staticmethod
    def write_pdb(atoms, outfile):
        """Write a pdb file from a list of Atom objects."""
        for atom in atoms:
            outfile.write(atom._generate_pdbline())
        for atom, next_atom in zip(atoms[:-1], atoms[1:]):
            if next_atom.resn not in END_ATOMS:
                outfile.write(f"CONECT{atom.atomid:>5d}{next_atom.atomid:>5d}\n")
            for restr in atom.restrs:
                if restr.type == "staple":
                    if atom.atomid < restr.atom.atomid:
                        outfile.write(
                            f"CONECT{atom.atomid:>5d}{restr.atom.atomid:>5d}\n"
                        )


@dataclasses.dataclass
class Nucleotide:
    letter: str
    index: int
    atoms: NucleotideAtoms

    damage_prob: "DNADamage" = dataclasses.field(init=False)
    damage: str = dataclasses.field(default=None, init=False)

    def conjure_atom(self, resn):
        """Create an atom which corresponds to the nucleotide.
        Note: this method does not assign next and prev atoms to any nucleotides."""
        assert self.atoms.current is None
        assert (prev_atom := self.atoms.prev) is not None
        assert (next_atom := self.atoms.next) is not None
        assert (row := prev_atom.mapping.row) == next_atom.mapping.row
        prev_index = prev_atom.nucleotide.index
        index = self.index
        next_index = next_atom.nucleotide.index
        chunk_length = next_index - prev_index
        coord = np.linspace(prev_atom.coords, next_atom.coords, chunk_length + 1)[
            index - prev_index
        ]
        map_col = np.linspace(
            prev_atom.mapping.column, next_atom.mapping.column, chunk_length + 1
        )[index - prev_index]
        self.atoms.current = Atom(
            resn, coord, nucleotide=self, mapping=MapCoord(row, map_col)
        )

    @classmethod
    def parse_sequence(cls, seq, cosm_atoms, deletions):
        """Create a list of Nucleotide objects.
        Args:
            seq: A nucleotide sequence.
            cosm_atoms: A list of Atom objects that correspond to the sequnce.
            deletions: A list of MapCoord objects with deletion coordinates.
        """
        nucleotides = []
        prev_atom = None
        seq_index = 0
        for atom, next_atom in zip(cosm_atoms[:-1], cosm_atoms[1:]):
            nucleotides.append(
                cls(
                    letter=seq[seq_index],
                    index=seq_index,
                    atoms=NucleotideAtoms(prev_atom, atom, next_atom),
                )
            )
            prev_atom = atom
            seq_index += 1
            if atom.mapping.row != next_atom.mapping.row:
                assert atom.mapping.column == next_atom.mapping.column
                continue
            if atom.mapping.column == next_atom.mapping.column:
                continue
            step = np.sign(next_atom.mapping.column - atom.mapping.column)
            column = atom.mapping.column + step
            while column != next_atom.mapping.column:
                if (atom.mapping.row, column) not in deletions:
                    nucleotides.append(
                        cls(
                            letter=seq[seq_index],
                            index=seq_index,
                            atoms=NucleotideAtoms(atom, None, next_atom),
                        )
                    )
                    seq_index += 1
                column += step
        nucleotides.append(
            cls(
                letter=seq[seq_index],
                index=seq_index,
                atoms=NucleotideAtoms(atom, next_atom, None),
            )
        )
        for nucleotide in nucleotides:
            if nucleotide.atoms.current is not None:
                nucleotide.atoms.current.nucleotide = nucleotide
        return nucleotides

    @classmethod
    def apply_damage(cls, nucleotides):
        """Resolve structural changes according to the simulated damage."""
        # Strand breaks (DSB; SSB at terminal atoms).
        for nucl in nucleotides:
            if (atom := nucl.atoms.current) is not None and atom.resn in ["S", "ST"]:
                continue
            if nucl.damage == "dsb":
                if (
                    atom is not None
                    and atom.resn in TERMINAL_ATOMS
                    and atom.get_mode() == "-1 0"
                ):
                    for restr in nucl.atoms.current.restrs:
                        if restr.type == "scaffold":
                            restr.atom.restrs.remove(
                                Restraint(
                                    restr.type, nucl.atoms.current, restr.parameters
                                )
                            )
                            nucl.atoms.current.restrs.remove(restr)
                            if restr.atom.nucleotide.damage == "rev_ssb":
                                restr.atom.nucleotide.damage = None
                            if (next_atom := nucl.atoms.next) is not None:
                                next_atom.resn += "T"
                                next_atom.nucleotide.atoms.prev = None
                                nucl.atoms.next = None
                            break
                    else:
                        raise RuntimeError("DSB at the end of a DS region.")
                else:
                    cls._convert_T_to_B(atom or nucl.atoms.prev, mode="prev")
                    cls._convert_T_to_B(nucl.atoms.next, mode="next")

                    cls._break_left(nucl, nucleotides)
                    cls._break_right(nucl, nucleotides)
                    pass
            elif (
                nucl.damage == "ssb"
                and atom is not None
                and atom.resn in TERMINAL_ATOMS
                and atom.get_mode() == "-1 0"
            ):
                if nucl.atoms.next is None:
                    continue
                nucl.atoms.next = None
                nucleotides[nucl.index + 1].atoms.prev = None
                nucleotides[nucl.index + 1].atoms.current.resn += "T"
            elif (
                nucl.damage == "rev_ssb"
                and atom is not None
                and atom.resn in TERMINAL_ATOMS
                and atom.get_mode() == "0 +1"
                and nucl.index != 0
            ):
                for restr in atom.restrs:
                    if restr.type == "scaffold":
                        restr.atom.restrs.remove(
                            Restraint(restr.type, nucl.atoms.current, restr.parameters)
                        )
                        nucl.atoms.current.restrs.remove(restr)
                        if restr.atom.nucleotide.damage == "dsb":
                            restr.atom.nucleotide.damage = "ssb"
                        break
        # Introducing special atoms (SSB in DS regions; AP sites; base damage).
        for nucl in nucleotides:
            if nucl.damage is None or (
                (atom := nucl.atoms.current) is not None
                and atom.resn in TERMINAL_ATOMS | {"S", "ST"}
            ):
                continue
            cls._fill_with_N(
                nucleotides,
                nucl.atoms.prev.nucleotide.index,
                nucl.atoms.next.nucleotide.index,
            )
            new_resn = DAMAGED_ATOM_NAMES[nucl.damage]
            if type(new_resn) == dict:
                possible_names = new_resn[nucl.letter]
                if len(possible_names) == 0:
                    continue
                new_resn = random.choice(possible_names)
            nucl.atoms.current.resn = new_resn
            
        # Fix potential N-PT-N.
        for nucl in nucleotides:
            if (nucl.atoms.current is not None
                and nucl.atoms.current.resn == "PT"
                and nucl.atoms.prev is not None
                and nucl.atoms.prev.resn in {"N"} | DAMAGED_ATOMS
                and nucl.atoms.next is not None
                and nucl.atoms.next.resn in {"N"} | DAMAGED_ATOMS
            ):
                nucl.atoms.current.resn = "N"

    @classmethod
    def _fill_with_N(cls, nucleotides, left_index, right_index):
        """Fill all nucleotides between left and right indexes with N atoms.
        Note: this method also fixes resn at the ends."""
        assert left_index <= right_index
        assert (left_atom := nucleotides[left_index].atoms.current) is not None
        assert (right_atom := nucleotides[right_index].atoms.current) is not None
        assert left_atom.resn in TERMINAL_ATOMS | {"H", "PT", "N"} | DAMAGED_ATOMS
        if left_atom.resn in ["H", "PT"]:
            left_atom.resn = "PT"
        elif left_atom.resn in TERMINAL_ATOMS:
            new_name = left_atom.resn[0] + "1"
            if left_atom.resn in END_ATOMS:
                new_name += "T"
            left_atom.resn = new_name
        assert right_atom.resn in TERMINAL_ATOMS | {"H", "PT", "N"} | DAMAGED_ATOMS
        if right_atom.resn in ["H", "PT"]:
            right_atom.resn = "PT"
        elif right_atom.resn in TERMINAL_ATOMS:
            assert right_atom.resn not in END_ATOMS
            right_atom.resn = right_atom.resn[0] + "1"

        for index in range(left_index + 1, right_index):
            if nucleotides[index].atoms.current is None:
                nucleotides[index].conjure_atom(resn="N")
            else:
                nucleotides[index].atoms.current.resn = "N"
        cls._fix_near_atoms(nucleotides, left_index, right_index)

    @classmethod
    def _fix_near_atoms(cls, nucleotides, left_index, right_index):
        """Assign prev and next atoms for nucleotides between indexes (ends included).
        Note: the chunk should be connected (one chain).
        """
        assert left_index < right_index
        assert (left_atom := nucleotides[left_index].atoms.current) is not None
        for index in range(left_index + 1, right_index + 1):
            nucl = nucleotides[index]
            nucl.atoms.prev = left_atom
            left_atom = nucl.atoms.current or left_atom
        assert (right_atom := nucleotides[right_index].atoms.current) is not None
        for index in range(right_index - 1, left_index, -1):
            nucl = nucleotides[index]
            nucl.atoms.next = right_atom
            right_atom = nucl.atoms.current or right_atom
        nucleotides[left_index].atoms.next = right_atom

    @staticmethod
    def _convert_T_to_B(atom, mode):
        """Resolve terminal atoms' names (T -> B if no staples)."""
        assert mode in ["prev", "next"], f"{mode=}"
        if atom is None:
            return
        staple_restrs = set()  # only row numbers
        row = atom.mapping.row
        next_atom = atom
        while next_atom.mapping.row == row:
            atom = next_atom
            for restr in atom.restrs:
                if restr.type == "staple":
                    staple_restrs.add(restr.atom.mapping.row)
            if (next_atom := getattr(atom.nucleotide.atoms, mode)) is None:
                return
        if atom.resn not in TERMINAL_ATOMS or next_atom.resn not in TERMINAL_ATOMS:
            return
        assert atom.resn[0] == next_atom.resn[0]
        for restr in atom.restrs:
            if restr.type == "scaffold":
                return
        if atom.resn[0] == "T" and next_atom.mapping.row not in staple_restrs:
            atom.resn = f"B{atom.resn[1:]}"
            next_atom.resn = f"B{next_atom.resn[1:]}"

    @classmethod
    def _break_left(cls, nucl, nucleotides):
        """Resolve DSB for nucleotides to the left of nucl and fix resn of nucl."""
        if (atom := nucl.atoms.current) is None:
            nucl.conjure_atom(resn="_")
            cls._fix_near_atoms(
                nucleotides, nucl.atoms.prev.nucleotide.index, nucl.index
            )
            atom = nucl.atoms.current
        if atom.resn in TERMINAL_ATOMS:
            if (mode := atom.get_mode()) == "-1 0":
                raise RuntimeError("DSB at the end of a DS region.")
            if mode == "0 +1":
                new_name = atom.resn[0] + "1"  # B1 or T1
                if atom.resn in END_ATOMS:
                    new_name += "T"
                atom.resn = new_name
                return
            raise RuntimeError(f"Unexpected mode: {mode}.")
        if (prev_atom := nucl.atoms.prev).resn in ["H", "PT"]:
            if atom.resn in ["H", "PT"]:
                atom.resn = "T"
                return
            prev_atom.resn = "PT"
            atom.resn = f"T{nucl.index - prev_atom.nucleotide.index}"
            return
        if prev_atom.resn in TERMINAL_ATOMS:
            atom.resn = "T1"
            cls._fill_with_N(nucleotides, prev_atom.nucleotide.index, nucl.index)
            if nucl.index == prev_atom.nucleotide.index + 1:
                atom.resn = "T1" if prev_atom.resn.startswith("B") else "B1"
                return
        if nucl.atoms.prev.resn in ["N", "S", "ST"]:
            assert atom is not None
            assert nucl.index == nucl.atoms.prev.nucleotide.index + 1
            atom.resn = "T1"
            return
        raise RuntimeError("Could not resolve DSB to the left.")

    @classmethod
    def _break_right(cls, nucl, nucleotides):
        """Resolve DSB for nucleotides to the right of nucl and fix nucl.atoms.next."""
        assert nucl.atoms.next is not None, "DSB at the end of a chain."
        next_nucl = nucleotides[nucl.index + 1]
        if (atom := next_nucl.atoms.current) is None:
            next_nucl.conjure_atom(resn="_")
            cls._fix_near_atoms(
                nucleotides, next_nucl.index, next_nucl.atoms.next.nucleotide.index
            )
            atom = next_nucl.atoms.current
        assert atom.resn not in END_ATOMS
        if atom.resn in TERMINAL_ATOMS:
            mode = atom.get_mode()
        nucl.atoms.next = next_nucl.atoms.prev = None
        if atom.resn in TERMINAL_ATOMS:
            if mode == "0 +1":
                atom.resn += "T"
                return
            if mode == "-1 0":
                atom.resn = atom.resn[0] + "1T"
                return
            raise RuntimeError(f"Unexpected mode: {mode}.")
        if (next_atom := next_nucl.atoms.next).resn in ["H", "PT"]:
            if atom.resn in ["H", "PT"]:
                atom.resn = "TT"
                return
            next_atom.resn = "PT"
            atom.resn = f"T{next_atom.nucleotide.index - next_nucl.index}T"
            return
        if next_atom.resn in TERMINAL_ATOMS:
            atom.resn = "T1T"
            cls._fill_with_N(nucleotides, next_nucl.index, next_atom.nucleotide.index)
            if next_nucl.index + 1 == next_atom.nucleotide.index:
                atom.resn = "T1T" if next_atom.resn.startswith("B") else "B1T"
                return
        if next_nucl.atoms.next.resn in ["N", "S"]:
            assert atom is not None
            assert next_nucl.index + 1 == next_nucl.atoms.next.nucleotide.index
            atom.resn = "T1T"
            return
        raise RuntimeError("Could not resolve DSB to the right.")


@dataclasses.dataclass
class DNADamage:
    dsb_prob: float
    ssb_prob: float
    abasic_prob: float
    base_prob: float
    rev_ssb_prob: float
    rev_abasic_prob: float
    rev_base_prob: float

    _DAMAGE_TYPE_PRIORITY = [
        "dsb",
        "ssb",
        "rev_ssb",
        "abasic",
        "rev_abasic",
        "base",
        "rev_base",
    ]
    _doublet_probabilities = {"dsb": {}, "ssb": {}}
    _triplet_probabilities = {}

    def simulate(self):
        """Return random damage type based on assigned probabilities."""
        for damage_type in self._DAMAGE_TYPE_PRIORITY:
            damage_prob = getattr(self, f"{damage_type}_prob")
            if np.random.random() < damage_prob:
                return damage_type
        return None

    @classmethod
    def adjust_probabilities(cls, dose):
        """Adjust damage probabilities based on the dose."""
        prior_triplet_probabilities = damage_probabilities.get_triplet_probabilities()
        for damage_type, prior_dose in damage_probabilities.PRIOR_DOSES.items():
            adjusted_probabilities = {
                seq: 1 - (1 - prob) ** (dose / prior_dose)
                for seq, prob in prior_triplet_probabilities[damage_type].items()
            }
            cls._triplet_probabilities[damage_type] = defaultdict(
                lambda: 0, adjusted_probabilities
            )

        # Set DSB and SSB doublet probabilities.
        prior_doublet_probabilities = damage_probabilities.get_doublet_probabilities()
        for damage_type in ["dsb", "ssb"]:
            prior_dose = damage_probabilities.PRIOR_DOSES[damage_type]
            for mode in ["-1 0", "0 +1"]:
                adjusted_probabilities = {
                    seq: 1 - (1 - prob) ** (dose / prior_dose)
                    for seq, prob in prior_doublet_probabilities[damage_type][
                        mode
                    ].items()
                }
                cls._doublet_probabilities[damage_type][mode] = defaultdict(
                    lambda: 0, adjusted_probabilities
                )

    @classmethod
    def generate_prob(cls, triplet, mode):
        """Calculate damage probabilities for the triplet."""
        assert len(triplet) == 3, f"Invalid sequence length: {len(triplet)}."
        assert cls._triplet_probabilities, (
            "DNADamage.adjust_probabilities has to be invoked before generating"
            " probabilities."
        )
        if mode == "-1 +1":
            return cls(
                cls._triplet_probabilities["dsb"][triplet],
                cls._triplet_probabilities["ssb"][triplet],
                cls._triplet_probabilities["abasic"][triplet],
                cls._triplet_probabilities["base"][triplet],
                cls._triplet_probabilities["ssb"][rev_compl(triplet)],
                cls._triplet_probabilities["abasic"][rev_compl(triplet)],
                cls._triplet_probabilities["base"][rev_compl(triplet)],
            )
        if mode == "-1 0":
            doublet = triplet[:-1]
            return cls(
                dsb_prob=0,
                ssb_prob=cls._doublet_probabilities["ssb"][mode][doublet],
                abasic_prob=0,
                base_prob=0,
                rev_ssb_prob=0,
                rev_abasic_prob=0,
                rev_base_prob=0,
            )
        if mode == "0 +1":
            doublet = triplet[1:]
            return cls(
                dsb_prob=cls._doublet_probabilities["dsb"][mode][doublet],
                ssb_prob=0,
                abasic_prob=0,
                base_prob=0,
                rev_ssb_prob=0,
                rev_abasic_prob=0,
                rev_base_prob=0,
            )
        raise ValueError(f"Invalid mode: {mode}.")

    @classmethod
    def assign_probs(cls, nucleotides):
        """Calculate and assign damage probabilities for nucleotides."""
        for prev_nucl, nucl, next_nucl in zip(
            [None] + nucleotides[:-1], nucleotides, nucleotides[1:] + [None]
        ):
            mode = "-1 +1"
            if (atom := nucl.atoms.current) is not None and atom.resn in TERMINAL_ATOMS:
                mode = atom.get_mode()
                for restraint in atom.restrs:
                    if restraint.type == "scaffold":
                        if mode == "-1 0":
                            next_nucl = restraint.atom.nucleotide
                        elif mode == "0 +1":
                            prev_nucl = restraint.atom.nucleotide
                        else:
                            raise ValueError(f"Invalid mode: {mode}.")
                        mode = "-1 +1"
            if prev_nucl is None:
                triplet = "_" + nucl.letter + next_nucl.letter
            elif next_nucl is None:
                triplet = prev_nucl.letter + nucl.letter + "_"
            else:
                triplet = prev_nucl.letter + nucl.letter + next_nucl.letter
            nucl.damage_prob = cls.generate_prob(triplet, mode)


def rev_compl(seq):
    """Return reverse complement of the sequence."""
    compl = str.maketrans("ATGC", "TACG")
    return seq[::-1].translate(compl)


# Read input pdb.
initial_atoms = []
with open(args.pdb) as infile:
    for line in infile.readlines():
        if line.startswith("ATOM"):
            initial_atoms.append(Atom.parse_pdbline(line))

# Check that atom ids are consecutive and start with 1.
for i, atom in enumerate(initial_atoms):
    if atom.atomid != i + 1:
        raise ValueError(f"Not consecutive atom ids in pdb at atom id {atom.atomid}.")

# Check that there is only one chain.
assert initial_atoms[0].resn in TERMINAL_ATOMS
for atom in initial_atoms[1:]:
    if atom.resn in END_ATOMS:
        raise NotImplementedError("More than one chain in input pdb.")

# Check atom types
for atom in initial_atoms:
    assert atom.resn in ATOMS

# Read input restraints.
with open(args.restr) as infile:
    restr_type = None
    for line in infile.readlines():
        if not line.strip() or line.startswith("["):
            continue
        elif line.startswith(";"):
            restr_type = line.split()[1]
            if restr_type not in ["staple", "scaffold", "lattice"]:
                raise ValueError(f"Got an unexpected restraint type: {restr_type}")
        else:
            atomid1, atomid2, _, _, _, *parameters = line.split()
            atom1 = initial_atoms[int(atomid1) - 1]
            atom2 = initial_atoms[int(atomid2) - 1]
            atom1.restrs.append(Restraint(restr_type, atom2, parameters))
            atom2.restrs.append(Restraint(restr_type, atom1, parameters))

# Read input mapping.
with open(args.map) as infile:
    deletions = set()
    map_index = 0
    for line in infile.readlines():
        if line.startswith("D"):
            column, row = line.split()[1].split(":")
            deletions.add(MapCoord(int(row), int(column)))
        elif line.startswith(("T", "I")):
            row, column = line.split(":")[-2:]
            initial_atoms[map_index].mapping = MapCoord(int(row), int(column))
            map_index += 1

# Read input sequence.
with open(args.sequence) as infile:
    seq = "".join(infile.read().split())

DNADamage.adjust_probabilities(args.dose)
nucleotides = Nucleotide.parse_sequence(seq, initial_atoms, deletions)
DNADamage.assign_probs(nucleotides)
for nucleotide in nucleotides:
    nucleotide.damage = nucleotide.damage_prob.simulate()
Nucleotide.apply_damage(nucleotides)

final_atoms = [
    nucl.atoms.current for nucl in nucleotides if nucl.atoms.current is not None
]

# Fix atom_ids and chains.
chain_index = -1
for i, atom in enumerate(final_atoms):
    atom.atomid = i + 1
    if atom.resn in END_ATOMS:
        try:
            chain = next(CHAIN_NAMES)
        except StopIteration:
            raise OverflowError("Too many chains in the output.")
    atom.chain = chain


# Write output pdb.
with open(args.output, "w") as outfile:
    Atom.write_pdb(final_atoms, outfile)

Write output restraints.
final_restrs = {"staple": [], "scaffold": [], "lattice": []}
for atom in final_atoms:
    for restr_type, other_atom, parameters in atom.restrs:
        if atom.atomid < other_atom.atomid:
            final_restrs[restr_type].append(
                (atom.atomid, other_atom.atomid, parameters)
            )

# Write output restraints.
restr_index = 1
with open(args.restr_out, "w") as outfile:
    outfile.write("[ distance_restraints ]\n")

    outfile.write("; staple crossovers\n")
    for atomid1, atomid2, parameters in final_restrs["staple"]:
        parameters = "\t".join(parameters)
        outfile.write(f"{atomid1}\t{atomid2}\t1\t{restr_index}\t1\t{parameters}\n")
        restr_index += 1

    outfile.write("; scaffold crossovers\n")
    for atomid1, atomid2, parameters in final_restrs["scaffold"]:
        parameters = "\t".join(parameters)
        outfile.write(f"{atomid1}\t{atomid2}\t1\t{restr_index}\t1\t{parameters}\n")
        restr_index += 1

    outfile.write("; lattice restraints\n")
    for atomid1, atomid2, parameters in final_restrs["lattice"]:
        parameters = "\t".join(parameters)
        outfile.write(f"{atomid1}\t{atomid2}\t1\t{restr_index}\t1\t{parameters}\n")
        restr_index += 1
