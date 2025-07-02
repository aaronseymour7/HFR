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
from AaronTools.fileIO import FileReader
from core import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic




logfile = sys.argv[1]
reader = FileReader(logfile, just_geom=False)
energy = reader["energy"]
enthalpy = reader["enthalpy"]
free_energy = reader["free_energy"]
zpe = reader["E_ZPVE"]

print(f"SCF Energy: {energy} Eh")
print(f"Enthalpy: {enthalpy} Eh")
print(f"Free Energy: {free_energy} Eh")
print(f"ZPVE: {zpe} Eh")
