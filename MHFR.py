#!/usr/bin/env python3

import sys
import argparse
from core import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def prepare_mol(mol):
    try:
        mol = Chem.Mol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in mol.GetBonds():
            bond.SetIsAromatic(False)
        mol.UpdatePropertyCache(strict=False)
        

        mol = Chem.AddHs(mol)
        return mol 
    except Exception as e:
        print(f"Kekulization failed: {Chem.MolToSmiles(mol)} - {e}")
        return None  

function_map = {
    1: Isogyric,
    2: Isodesmic,
    3: Hypohomodesmotic,
    4: Homodesmotic
}
name_map ={
    1: "Isogyric",
    2: "Isodesmic",
    3: "Hypohomodesmotic",
    4: "Homodesmotic"
}
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_reactions.py infile outfile")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]


    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue

            try:
                u, v = line.split(maxsplit=1)
                level = int(u)
                
                mol = Chem.MolFromSmiles(v)

                if mol is None:
                    fout.write(f"INPUT: {v} LEVEL: {level} STATUS: ERROR - Invalid SMILES\n\n")
                    continue
                if prepare_mol(mol) is None:
                    fout.write(f"INPUT: {v} LEVEL: {level} STATUS: ERROR - Unkekulizable molecule\n\n")
                    continue
                if level not in function_map:
                    fout.write(f"INPUT: {v} LEVEL: {level} STATUS: ERROR - Invalid level\n\n")
                    continue

                fn = function_map[level]
                p, q, r = fn(mol)
                name = name_map[level]
                fout.write(f"INPUT: {v} LEVEL: {name} STATUS: {r}\n")
                fout.write("REACTANTS:\n")
                for g, f in q:
                    fout.write(f"-{f}*{Chem.MolToSmiles(g)}\n")

                fout.write("\nPRODUCTS:\n")
                for g, f in p:
                    fout.write(f"-{f}*{Chem.MolToSmiles(g)}\n")

                fout.write("\n")

            except Exception as e:
                fout.write(f"INPUT: {line} STATUS: ERROR - {str(e)}\n\n")
