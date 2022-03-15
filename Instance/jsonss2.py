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
from ase.db import connect
from pymatgen.core import Structure
from pymatgen.io import ase
from tqdm import tqdm

from ase_db_extension import db_from_structure_json, aaa, db_from_structure_dict

os.chdir(r"F:\data_new")
os.chdir("JARVIS-2D")

files = os.listdir()

file = "d3-5-16-2021.json"
new_file = file[:-5]
new_file = new_file + ".db"

db = connect("jdft_2d-7-7-2018.db")
