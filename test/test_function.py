import os
import unittest
import numpy as np
import pandas as pd
from ase.db import connect

from ase_db_extension import db_row2dct, db_dct2row, db_to_csv, db_from_ase_csv, db_from_structure_csv, db_cat, \
    db_rename
from ase_db_extension.ext_row import atoms_row_matcher, HashAtomsRow, atoms_row_matcher_log

os.chdir("../Instance")


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.db = connect("organometal.db")

    def test_matcher(self):
        with self.db:
            self.assertTrue(all([not isinstance(i,np.ndarray) for i in db_row2dct(self.db[1]).keys()]))

    def test_reback(self):
        with self.db:
            ar1 = self.db[1]
            dct = db_row2dct(self.db[1])
            atomsrow, new_kv, new_data = db_dct2row(dct)
            print(atoms_row_matcher_log(atomsrow,self.db[1]))
            self.assertTrue(atoms_row_matcher(atomsrow,self.db[1])>0.9)

    def test_csv(self):
        with self.db:

            db_to_csv(self.db, csv_file_name = "")
            file = db_to_csv(self.db, csv_file_name = None)
            self.assertIsInstance(file,pd.DataFrame)

    def test_from_ase_csv(self):

        if os.path.isfile("temp.db"):
            os.remove("temp.db")

        db = db_from_ase_csv("organometal.csv",new_database_name="temp.db")
        try:
            assert isinstance(db[1].get("space_group"),str)
            print(db[1].get("space_group"))
        except FileExistsError:
            pass

    def test_from_structure_csv(self):

        if os.path.isfile("temp.db"):
            os.remove("temp.db")

        db = db_from_structure_csv("kim_raw_data.csv",new_database_name="temp.db")
        try:
            assert isinstance(db[1].data["Label"],str)
            print(db[1].data["Label"])
        except FileExistsError:
            pass

    def test_data_cat(self):

        db2 = self.db
        db1 = connect("temp.db")
        db3 = db_cat(db1,db2,new_database_name="new2.db")

    def test_rename(self):

        ar = db_rename(self.db, name_pair = (("space_group2", "space_group"),), check = True)
        print(ar[1].key_value_pairs)
        print(ar[1].data)


        # ar = db_rename(self.db, name_pair = (("space_group2", "space_group1"),), check = True)


if __name__ == '__main__':
    unittest.main()
