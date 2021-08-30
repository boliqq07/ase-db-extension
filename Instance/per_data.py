# -*- coding: utf-8 -*-

# @Time     : 2021/8/30 10:48
# @Software : PyCharm
# @License  : GNU General Public License v3.0
# @Author   : xxx
import os

from ase.db import connect

os.chdir(r"F:\data_new")
# database1 = connect("CMR-absorption-spetra-perovskites/absorption_perovskites.db")
# database2 = connect("CMR-Functional-Perovskites/funct_perovskites.db")
# database3 = connect("CMR-low-sym-perovskites/low_symmetry_perovskites.db")
# database4 = connect("CMR-Organometal-Halide-Perovskites/organometal.db")
# database5 = connect("CMR-Perovskite-Water-Splitting/cubic_perovskites.db")

os.chdir("Discovery of Pb-free perovskite")
database6 = connect()