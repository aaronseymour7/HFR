#!/usr/bin/env python3

import os
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from pulp import LpProblem, LpVariable, lpSum, LpMinimize, LpInteger, PULP_CBC_CMD, LpStatus
from rdkit.Chem import Descriptors, rdMolDescriptors
from itertools import combinations
from AaronTools.geometry import Geometry
from AaronTools.atoms import Atom
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from AaronTools.fileIO import FileWriter
from core import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic


def isogyric_count(mol):
    atom_counts = {}
    if isinstance(mol, Chem.Atom):
        key = (mol.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    for atom in mol.GetAtoms():

        key = (atom.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
    return atom_counts

def isodesmic_count(mol):
    bond_counts = {}
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        atom_pair = tuple(sorted([atom1.GetSymbol(), atom2.GetSymbol()]))

        bond_type = (atom_pair, bond.GetBondType())
        bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1


    atom_counts = {}
    if isinstance(mol, Chem.Atom):
        key = (mol.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    for atom in mol.GetAtoms():
        if atom.GetSymbol() !='H':
            key = (atom.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
    return bond_counts,atom_counts

def hypohomodesmotic_count(mol):
    mol = Chem.AddHs(mol)
    hydrogen_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetNumExplicitHs())
        hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
    mol = Chem.RemoveHs(mol)
    atom_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetHybridization())
        atom_counts[key] = atom_counts.get(key, 0) + 1
    return hydrogen_counts, atom_counts

def homodesmotic_count(mol):
    bond_counts = {}
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol()== 'H' or bond.GetEndAtom().GetSymbol()== 'H':
            continue 
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        atom_pair = tuple(sorted(
            [(atom1.GetSymbol(), atom1.GetHybridization()), 
             (atom2.GetSymbol(), atom2.GetHybridization())]
        ))

        bond_type = (atom_pair, bond.GetBondType())
        bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1


    hydrogen_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetHybridization(), atom.GetNumExplicitHs())
        hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
    return bond_counts, hydrogen_counts

def display_reaction_counts(input_mol, reaction_fn):
    def format_key(k):
        if isinstance(k, tuple):
            return "(" + ", ".join(format_key(i) for i in k) + ")"
        if hasattr(k, 'name'):  # for RDKit enums like HybridizationType
            return k.name
        return str(k)

    def print_dict(title, d):
        if not d:
            return
        print(f"\n{title}:")
        print("-" * len(title))
        for key, val in sorted(d.items()):
            print(f"{format_key(key):<60} : {val}")    

    if reaction_fn is Isogyric:
        atom_counts = isogyric_count(input_mol)
        print_dict("Isogyric Atom Counts", atom_counts)

    elif reaction_fn is Isodesmic:
        bond_counts, atom_counts = isodesmic_count(input_mol)
        print_dict("Isodesmic Bond Counts", bond_counts)
        print_dict("Isodesmic Atom Counts", atom_counts)

    elif reaction_fn is Hypohomodesmotic:
        hydrogen_counts, atom_counts = hypohomodesmotic_count(input_mol)
        print_dict("Hypohomodesmotic Hydrogen Counts", hydrogen_counts)
        print_dict("Hypohomodesmotic Atom Counts", atom_counts)

    elif reaction_fn is Homodesmotic:
        bond_counts, hydrogen_counts = homodesmotic_count(input_mol)
        print_dict("Homodesmotic Bond Counts", bond_counts)
        print_dict("Homodesmotic Hydrogen Counts", hydrogen_counts)


def geom_from_rdkit(rdkitmol):
    """
    Takes an RDKit molecule (already embedded) and returns an AaronTools Geometry along with
    the weighted adjacency matrix
    """
    result = AllChem.EmbedMolecule(rdkitmol)
    if result != 0:
        raise ValueError("Embedding failed")
        
    atom_list = []
    for i, atom in enumerate(rdkitmol.GetAtoms()):
        positions = rdkitmol.GetConformer().GetAtomPosition(i)
        atom_list.append(Atom(element=atom.GetSymbol(), coords=positions))
    return Geometry(atom_list)

method = Theory(
    method =  "M062X",
    basis= "6-31G",
    #method="B3LYP",
    #basis="6-31G(d)",
    job_type=[OptimizationJob(), FrequencyJob()]
)

def parse_smiles_list(arg):
    if isinstance(arg, list):
        smiles_list = []
        for item in arg:
            if ',' in item:
                smiles_list.extend(item.split(','))
            else:
                smiles_list.append(item)
    else:
        smiles_list = arg.split(',') if ',' in arg else [arg]

    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise argparse.ArgumentTypeError(f"Invalid SMILES: {smi}")
        mols.append(mol)
    return mols

reaction_map = {
    "homodesmotic": Homodesmotic,
    "isodesmic": Isodesmic,
    "isogyric": Isogyric,
    "hypohomodesmotic": Hypohomodesmotic,
}

def parse_smiles_list(arg):
    if isinstance(arg, list):
        smiles_list = []
        for item in arg:
            if ',' in item:
                smiles_list.extend(item.split(','))
            else:
                smiles_list.append(item)
    else:
        smiles_list = arg.split(',') if ',' in arg else [arg]

    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise argparse.ArgumentTypeError(f"Invalid SMILES: {smi}")
        mols.append(mol)
    return mols

parser = argparse.ArgumentParser(description="Run reaction balancing (e.g. Homodesmotic, Isodesmic, etc.)")
parser.add_argument("action_type", choices=["write", "view", "count"], help="Type of action")
parser.add_argument("reaction_type", choices=reaction_map.keys(), help="Type of reaction")
parser.add_argument("input", type=str, help="Input molecule SMILES")

parser.add_argument("--lhs", nargs="+", type=str, help="LHS required molecules (SMILES or comma-separated list)")
parser.add_argument("--rhs", nargs="+", type=str, help="RHS required molecules (SMILES or comma-separated list)")
parser.add_argument("--substruct", nargs="+", type=str, help="Substructures to replace (SMILES or list)")
parser.add_argument("--replacement", nargs="+", type=str, help="Replacement structures (SMILES or list)")
parser.add_argument("--outfolder", type=str, default=None, help="Output folder to write files for the 'write' action.")
args = parser.parse_args()

input_smiles = str(args.input)  
input_mol = Chem.AddHs(Chem.MolFromSmiles(input_smiles))
if input_mol is None:
    raise argparse.ArgumentTypeError(f"Invalid SMILES: {args.input}")

lhs_required = parse_smiles_list(args.lhs) if args.lhs else None
rhs_required = parse_smiles_list(args.rhs) if args.rhs else None
Substruct = parse_smiles_list(args.substruct) if args.substruct else None
Replacement = parse_smiles_list(args.replacement) if args.replacement else None
Outfolder = args.outfolder

reaction_fn = reaction_map[args.reaction_type.lower()]


rhs, lhs, status = reaction_fn(input_mol, lhs_required, rhs_required, Substruct, Replacement)
if status == 'Optimal':
    Li = 1
    Ri = 1
    if args.action_type == 'write':
        if args.outfolder is None:
            parser.error("The --outfolder argument is required when action_type is 'write'.")
        else:
            outfolder = args.outfolder
            os.makedirs(outfolder, exist_ok=True)
            index_file = os.path.join(outfolder, "index.txt")
            with open(index_file, "w") as idx:
                idx.write(f"Level:\t {reaction_fn.__name__}\n")
                idx.write("Filename\tInChI\tSMILES\n")
                print("-----------Reactants-----------")
                for mol, coeff in lhs:
                    geom = geom_from_rdkit(mol)
                    smiles = Chem.MolToSmiles(mol)
                    inchi = Chem.MolToInchi(mol)
                    filename = f"R{Li}_{coeff}.com"
                    name = f"R{Li}_{coeff}"
                    outfile = os.path.join(outfolder, filename)
                    print(f"({coeff})*{smiles} = {filename}")
                    FileWriter.write_file(geom=geom, style="com", outfile=outfile, theory=method)
                    idx.write(f"{name}\t{inchi}\t{smiles}\n")
                    Li += 1

                print("-----------Products-----------")
                for mol, coeff in rhs:
                    geom = geom_from_rdkit(mol)
                    smiles = Chem.MolToSmiles(mol)
                    inchi = Chem.MolToInchi(mol)
                    filename = f"P{Ri}_{coeff}.com"
                    name = f"P{Ri}_{coeff}"
                    outfile = os.path.join(outfolder, filename)
                    print(f"({coeff})*{smiles} = {filename}")
                    FileWriter.write_file(geom=geom, style="com", outfile=outfile, theory=method)
                    idx.write(f"{name}\t{inchi}\t{smiles}\n")
                    Ri += 1
    if args.action_type  == 'view':
        print("-----------Reactants-----------")
        for mol, coeff in lhs:
            name = Chem.MolToSmiles(mol)
            print(f"({coeff})*{name}")
            Li+=1
        print("-----------Products-----------")
        for mol, coeff in rhs:
            name = Chem.MolToSmiles(mol)
            print(f"({coeff})*{name}")
            Ri+=1
        
    elif args.action_type  == 'count':
        display_reaction_counts(input_mol, reaction_fn)
    print("Complete")
else:
    print("Infeasible")
(base) ads09449@c4-20 bin$ cat hfr_driver.py 
#!/usr/bin/env python3

import os
import subprocess

reaction_map = {
    "1": "isogyric",
    "2": "isodesmic",
    "3": "hypohomodesmotic",
    "4": "homodesmotic"
}

def main(input_file):
    with open(input_file, "r") as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"[SKIP] Line {i} malformed: {line}")
                continue

            level, smiles = parts
            reaction_type = reaction_map.get(level)
            if not reaction_type:
                print(f"[SKIP] Invalid level {level} on line {i}")
                continue

            outfolder = f"{i}.mhfr"
            print(f"[RUN] Line {i}: level={level}, type={reaction_type}, SMILES={smiles}, folder={outfolder}")

            cmd = [
                "python3", "HFR.py",
                "write",
                reaction_type,
                smiles,
                "--outfolder", outfolder
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[ERROR] Failed on line {i}:\n{result.stderr}")
            else:
                print(f"[DONE] Folder written: {outfolder}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python3 batch_hfr_driver.py input_file.txt")
        sys.exit(1)
    main(sys.argv[1])
