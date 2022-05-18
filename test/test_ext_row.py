import os
import unittest
import numpy as np
from ase.db import connect

from ase_db_extension.ext_row import atoms_row_matcher, HashAtomsRow, atoms_row_rename

os.chdir("../Instance")


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.db = connect("organometal.db")

    def test_matcher(self):
        with self.db:
            self.assertEqual(atoms_row_matcher(self.db[1], self.db[1]), 1.0)
            self.assertLess(atoms_row_matcher(self.db[1], self.db[2]), 1.0)

    def test_hash(self):
        list1 = []
        with self.db:
            for i in self.db.select():
                list1.append(HashAtomsRow.from_atomsrow(i))
            a = set(list1)
            list1.append(HashAtomsRow.from_atomsrow(self.db[1]))
            b = set(list1)
            self.assertEqual(len(a), len(b))

    def test_rename(self):
        with self.db:
            ar = atoms_row_rename(self.db[1], name_pair=(("space_group", "space_group2"),), check=True)
            print(ar.key_value_pairs)
            print(ar.data)


if __name__ == '__main__':
    unittest.main()
