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
os.chdir("JARVIS-3D")

files = os.listdir()

file = "d3-5-16-2021.json"
new_file = file[:-5]
new_file = new_file + ".db"


def fmt(dct):
    si = Structure(lattice=np.array(dct["lattice_mat"]), species=dct['elements'], coords=dct['coords'], )
    atoms = aaa.get_atoms(si)
    return atoms


data = db_from_structure_json(file, new_database_name=new_file, structure_index_name="atoms", fmt=fmt)

# def fmt(dct):
#     si = Structure.from_dict(dct)
#     atoms = aaa.get_atoms(si)
#     return atoms
#
# data = db_from_structure_json(file, new_database_name=new_file, structure_index_name= "final_str", fmt = fmt,pop_name=["initial_str"])
