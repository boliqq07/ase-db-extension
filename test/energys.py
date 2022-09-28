# -*- coding: utf-8 -*-

# @Time  : 2022/9/14 20:06
# @Author : boliqq07
# @Software: PyCharm
# @License: MIT License
import os

import numpy as np
import pandas as pd
from pymatgen.core import Structure

from ase_db_extension import db_from_structure_dict
from featurebox.utils.general import AAA

os.remove("st_and_energy.db")

def fmt(st):
    st.lattice._pbc=np.array([True,True,True])
    try:
        atoms = AAA.get_atoms(st)
    except ValueError:
        st.remove_site_property("selective_dynamics")
        atoms = AAA.get_atoms(st)
    return atoms

data = pd.read_pickle("st_and_energy.pkl_pd")

db = db_from_structure_dict(data, new_database_name="st_and_energy.db",
                           structure_index_name="structure", fmt=fmt,
                           pop_name =None)