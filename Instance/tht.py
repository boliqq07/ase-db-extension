# -*- coding: utf-8 -*-

# @Time     : 2021/8/30 15:26
# @Software : PyCharm
# @License  : GNU General Public License v3.0
# @Author   : xxx
import os
import json
from multiprocessing import Pool
import numpy as np
import pandas as pd
from pymatgen.core import Structure
from pymatgen.io import ase
from tqdm import tqdm

from ase_db_extension import db_from_structure_json, aaa, db_from_structure_dict

os.chdir(r"F:\data_new")
os.chdir("Thermodynamic-High-throuput")
f = open("perovskites-HT-ComputedStructureEntry.json")
entries = json.load(f)["entries"]
f.close()


def fmt(dct):
    si = Structure.from_dict(dct)
    atoms = aaa.get_atoms(si)
    return atoms


data = db_from_structure_dict(entries, new_database_name="Thermodynamic-High-throughput.db",
                              structure_index_name="structure", fmt=fmt)
