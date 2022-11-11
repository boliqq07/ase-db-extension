# -*- coding: utf-8 -*-

# @Time     : 2021/8/29 18:55
# @Software : PyCharm
# @License  : GNU General Public License v3.0
# @Author   : xxx
"""Interface for ase.db and csv"""
import json
import numbers
import os
import re
import shutil
import warnings
from collections import Callable
from pathlib import PurePath
from typing import Dict, Union, Tuple, Any, List

import numpy as np
import pandas as pd
from ase.db import connect
from ase.db.core import Database, reserved_keys
from ase.db.row import AtomsRow, atoms2dict
from pymatgen.core import Structure

from ase_db_extension.ext_row import check_name_tup, atoms_row_rename
from ase_db_extension.general import _aaa as aaa

reserved_keys_extend = reserved_keys

special_kv_keys = {'project', 'name', 'symmetry', 'space_group'}

word = re.compile('[_a-zA-Z][_0-9a-zA-Z]*$')


# data just accept float, ndarray, str.

# #############store########################

def _code_dct(dct: Dict) -> Dict:
    """Change the numpy as list for store."""
    for k, v in dct.items():
        if isinstance(v, np.ndarray):
            dct[k] = v.tolist()
        elif isinstance(v, dict):
            dct[k] = _code_dct(v)
        elif isinstance(v, (list, tuple)):
            dct[k] = [_code_dct(i) if isinstance(i, dict) else v for i in v]
        else:
            dct[k] = v
    return dct


def db_row2dct(row: AtomsRow
               ) -> Dict[str, Any]:
    """Convert row to dict,and change the numpy as list for store."""
    temp = {}
    temp.update(row.data)
    temp.update(row.key_value_pairs)
    dct = {i: row.__dict__[i] for i in reserved_keys if i in row.__dict__}
    temp.update(dct)

    temp = _code_dct(temp)

    return temp


# #############load########################


def _decode_dct(dct: Dict) -> Dict:
    """Change the list as np.ndarray for store."""
    for k, v in dct.items():
        if k == "unique_id":
            dct[k] = v

        elif isinstance(v, float):
            if np.isnan(v):
                dct[k] = None
            else:
                dct[k] = v
        else:
            try:
                if isinstance(v, str):
                    v = eval(v)
                if isinstance(v, list):
                    if v == []:
                        dct[k] = None
                    else:
                        nda = np.array(v)
                        if nda.dtype in (np.float64, np.int32) and nda.ndim >= 1:
                            dct[k] = nda
                        else:
                            dct[k] = v
                elif isinstance(v, dict):
                    dct[k] = _decode_dct(v)
                else:
                    dct[k] = v
            except (ValueError, NameError, TypeError, SyntaxError):
                pass
    return dct


def collect_dct3(dct: Dict, kv_names=None) -> Tuple:
    """Take part the dict to 3 part."""
    if "data" in dct:
        data = dct["data"]
        del dct["data"]
    else:
        data = {}

    if "key_value_pairs" in dct:
        key_value_pairs = dct["key_value_pairs"]
        del dct["key_value_pairs"]
    else:
        key_value_pairs = {}

    temp = {}

    temp.update(data)
    temp.update(key_value_pairs)
    temp.update(dct)

    new_dct, new_data, new_kv = {}, {}, {}

    kv_keys = special_kv_keys if kv_names is None else special_kv_keys & set(kv_names)

    for k, v in temp.items():
        if k in reserved_keys:
            new_dct.update({k: v})
        elif k in kv_keys:
            new_kv.update({k: v})
        elif not word.match(k) or not isinstance(v, (numbers.Real, np.bool_)):
            new_data.update({k: v})
        else:
            new_kv.update({k: v})

    new_dct["key_value_pairs"] = new_kv
    new_dct["data"] = new_data

    return new_dct, new_kv, new_data


def db_dct2row(dct: Dict, kv_names=None) -> [AtomsRow, Dict, Dict]:
    """Convert dict to row."""
    dct = _decode_dct(dct)
    new_dct, new_kv, new_data = collect_dct3(dct, kv_names=kv_names)
    atomsrow = AtomsRow(new_dct)
    return atomsrow, new_kv, new_data


def collect_dct1(dct: Dict) -> Dict:
    new_dct, new_kv, new_data = collect_dct3(dct)
    new_data.update(new_kv)
    new_data.update(new_dct)
    return new_data


def db_to_csv(database: Union[str, Database], csv_file_name=""):
    """Store the db to csv."""
    if isinstance(database, str):
        database = connect(database)

    data = []

    for row in database.select(selection=None):
        try:
            data.append(db_row2dct(row))
        except (KeyError, StopIteration, AttributeError) as e:
            print(e, row)

    data = pd.DataFrame.from_dict(data)

    if csv_file_name:
        data.to_csv("{}.csv".format(os.path.splitext(csv_file_name)[0]))
    elif csv_file_name == "":
        csv_file_name = os.path.split(database.filename)[-1]
        csv_file_name = os.path.splitext(csv_file_name)[0]
        data.to_csv("{}.csv".format(csv_file_name))
    else:
        return data


def db_from_ase_csv(csv_file: Union[str, pd.DataFrame], new_database_name="new_database.db", index_col=0,
                    pop_name: List[str] = None
                    ) -> Database:
    """Convert csv(with strand Atoms table) to ase.db.
    see also: db_to_csv"""
    if pop_name is not None:
        assert isinstance(pop_name, list)
    if new_database_name == "" and isinstance(csv_file, str):
        new_database_name = "".join((os.path.splitext(csv_file)[0], ".db"))
    assert ".db" in new_database_name or ".json" in new_database_name
    if os.path.isfile(new_database_name):
        raise FileExistsError(new_database_name, "is exist, please change the `new_database_name`.")
    database = connect(new_database_name)
    if isinstance(csv_file, str):
        file = pd.read_csv(csv_file, index_col=index_col)
    else:
        file = csv_file

    data = file.to_dict(orient="index")

    with database:
        for k in data.keys():
            atomsrow, new_kv, new_data = db_dct2row(data[k])
            if pop_name is not None:
                for i in pop_name:
                    if i in new_kv:
                        del new_kv[i]
                    elif i in new_data:
                        del new_data[i]
            database.write(atomsrow, key_value_pairs=new_kv, data=new_data)
    return database


def db_from_structure_csv(csv_file: Union[str, pd.DataFrame], new_database_name="new_database.db", index_col=0,
                          structure_index_name: str = None, fmt: Union[str, Callable] = "json",
                          pop_name: List[str] = None
                          ):
    """Convert csv to ase.db. The structure must be offered."""

    if isinstance(csv_file, str):
        file = pd.read_csv(csv_file, index_col=index_col)
    else:
        file = csv_file

    data = file.to_dict(orient="index")

    if new_database_name == "" and isinstance(csv_file, str):
        new_database_name = "".join((os.path.splitext(csv_file)[0], ".db"))
    assert ".db" in new_database_name or ".json" in new_database_name

    if structure_index_name and structure_index_name in file.columns:
        name = structure_index_name
    elif "structure" in file.columns:
        name = "structure"
    else:
        raise NameError("There must be structure index name.")

    database = db_from_structure_dict(data, new_database_name=new_database_name,
                                      structure_index_name=name, fmt=fmt, pop_name=pop_name)

    return database


def db_from_structure_json(json_file: Union[str, dict], new_database_name="new_database.db",
                           structure_index_name: str = "structure", fmt: Union[str, Callable] = "json",
                           pop_name: List[str] = None
                           ):
    """Convert csv to ase.db. The structure must be offered."""

    def read_json(i):
        f = open(i)
        entries = json.load(f)
        f.close()
        return entries

    if isinstance(json_file, str):
        data = read_json(json_file)
    else:
        data = json_file

    if new_database_name == "" and isinstance(json_file, str):
        new_database_name = "".join((os.path.splitext(json_file)[0], ".db"))
    assert ".db" in new_database_name or ".json" in new_database_name

    if structure_index_name is None:
        structure_index_name = "structure"
    name = structure_index_name

    database = db_from_structure_dict(data, new_database_name=new_database_name,
                                      structure_index_name=name, fmt=fmt, pop_name=pop_name)

    return database


def db_from_structure_dict(data, new_database_name="new_database.db",
                           structure_index_name: str = None, fmt: Union[str, Callable] = "json",
                           pop_name: List[str] = None):
    """Convert csv to ase.db. The structure must be offered."""
    if pop_name is not None:
        assert isinstance(pop_name, list)

    assert ".db" in new_database_name or ".json" in new_database_name
    if os.path.isfile(new_database_name):
        raise FileExistsError(new_database_name, "is exist, please change the `new_database_name`.")

    if isinstance(data, list):
        number = range(len(data))
    else:
        number = data.keys()

    database = connect(new_database_name)

    with database:
        for k in number:
            try:
                datak = data[k]
                structure = datak.pop(structure_index_name)
                if isinstance(fmt, str):
                    st = Structure.from_str(structure, fmt=fmt)
                    atoms = aaa.get_atoms(st)
                elif isinstance(fmt, Callable):
                    atoms = fmt(structure)
                dct = _decode_dct(datak)
                new_dct, new_kv, new_data = collect_dct3(dct)
                dct = atoms2dict(atoms)
                dct.update(new_dct)
                if pop_name is not None:
                    for i in pop_name:
                        if i in new_kv:
                            del new_kv[i]
                        elif i in new_data:
                            del new_data[i]

                database.write(dct, key_value_pairs=new_kv, data=new_data)
            except NameError as e:
                # except BaseException as e:
                warnings.warn("The {} sample is can't be analysis.".format(k))
                print(e)

    return database


# ####function#######


def check_postfix(name: Union[Database, str, dict]):
    if isinstance(name, PurePath):
        name = str(name)

    if isinstance(name, dict):
        postfix = '.json'
    elif isinstance(name, str):
        postfix = os.path.splitext(name)[-1]
    else:
        from ase.db.jsondb import JSONDatabase
        from ase.db.sqlite import SQLite3Database

        if isinstance(name, JSONDatabase):
            postfix = '.json'
        elif isinstance(name, SQLite3Database):
            postfix = '.db'
        else:
            postfix = ''

    return postfix


def db_append(database1: Database, database2: Database, parallel=False) -> None:
    """Replace the in situ."""
    if parallel:
        with database1:
            with database2:
                for ar in database2.select():
                    id = database1.reserve(name=str(ar))
                    if id is None:
                        continue
                    database1.write(ar, ar.key_value_pairs, ar.data, id=id)
    else:
        with database1:
            with database2:
                for ar in database2.select():
                    database1.write(ar, ar.key_value_pairs, ar.data)


def db_remove(database: Database, sli: [slice, List, Tuple]) -> None:
    """Delete the in situ!!!!!,Keep the id in sli!!!."""
    if isinstance(sli, Tuple):
        assert len(sli) <= 3, 'just accept (start,stop,[step])'
        sli = slice(*sli)
    if isinstance(sli, slice):
        assert sli.start >= 1, "slice for database start from 1"
        if sli.step is None:
            sli = slice(sli.start, sli.stop, 1)
        sli_list = list(range(sli.start, sli.stop, sli.step))
    else:
        sli_list = sli
    with database:
        database.delete(sli_list)
    print(sli_list, "is deleted")


def db_slice(database: Database, sli: [slice, Tuple], new_database_name: str = None, paraller=False) -> Union[
    Database, List[AtomsRow]]:
    """Get slice of database, in new database, Keep the id in sli!!!"""
    if isinstance(database, str):
        database = connect(database)

    if isinstance(sli, Tuple):
        assert len(sli) <= 3, 'just accept (start,stop,[step])'
        sli = slice(*sli)

    assert isinstance(sli, slice) and sli.start >= 1, "slice for database start from 1"
    if sli.step is None:
        sli = slice(sli.start, sli.stop, 1)
    if isinstance(new_database_name, str):
        new_database_name = os.path.splitext(new_database_name)[0]
        new_database_name = "".join((new_database_name, "{}"))
        database1_ = database.__class__(new_database_name.format(check_postfix(database)))
        if paraller:
            with database:
                with database1_:
                    for i in range(sli.start, sli.stop, sli.step):
                        id = database1_.reserve(name=i)
                        if id is None:
                            continue
                        ar = database[i]
                        database1_.write(ar, ar.key_value_pairs, ar.data, id=id)
        else:
            with database:
                with database1_:
                    for i in range(sli.start, sli.stop, sli.step):
                        ar = database[i]
                        database1_.write(ar, ar.key_value_pairs, ar.data)
    else:
        database1_ = []
        with database:
            for i in range(sli.start, sli.stop, sli.step):
                database1_.append(database[i])
    return database1_


def db_cat(database1: Union[Database, str], database2: Union[Database, str], new_database_name: str) -> Database:
    """new database, in new database."""
    new_database_name = os.path.splitext(new_database_name)[0]
    new_database_name = "".join((new_database_name, "{}"))

    if isinstance(database1, Database):
        database1 = database1.filename

    if isinstance(database1, str):
        if isinstance(database2, str):
            database2 = connect(database2)

        assert os.path.isfile(database1)
        new = new_database_name.format(check_postfix(database1))
        shutil.copy(database1, new)
        new = connect(new)
        db_append(new, database2)
        return database1

    else:
        # old deprecated
        if isinstance(database1, str):
            database1 = connect(database1)
        if isinstance(database2, str):
            database2 = connect(database2)
        database1_ = database1.__class__(new_database_name.format(check_postfix(database1)))

        with database1_:
            with database1:
                for ar in database1.select(selection=None):
                    database1_.write(ar, ar.key_value_pairs, ar.data)
            with database2:
                for ar in database2.select(selection=None):
                    database1_.write(ar, ar.key_value_pairs, ar.data)
        return database1_


def db_ids(database: Database, selection=None, **kwargs) -> Tuple:
    """Get ids.
    """
    ids = []
    with database:
        for row in database.select(selection, **kwargs):
            ids.append(row.id)

    t = ids[0]
    i = ids[0]
    pstr = "({}-".format(t)
    for i in ids:
        if i - t <= 1:
            pass
        else:
            pstr = pstr + "{},{}-".format(t, i)
        t = i
    pstr = pstr + "{})".format(i)
    print(pstr)
    return tuple(ids)


def db_rename(database: Database, name_pair=(("old_name1", "new_name1"),), check=True, selection=None,
              **kwargs) -> Database:
    if check:
        check_name_tup(name_pair)

    if os.path.splitext(database.filename)[1] == '.json':
        js = True
    else:
        js = False
    with database:
        for row in database.select(selection, **kwargs):
            new_row = atoms_row_rename(row, name_pair=name_pair, check=False)
            kvp = new_row.key_value_pairs
            data = new_row.get('data', {})

            if js:
                database._write(new_row, kvp, data, row.id)
            else:
                database._update(row.id, kvp, data)

    return database


def db_transform(database, new_database_name="new_database.json"):
    """Convert ase.json to ase.db, or inverse."""

    assert ".db" in new_database_name or ".json" in new_database_name
    if os.path.isfile(new_database_name):
        raise FileExistsError(new_database_name, "is exist, please change the `new_database_name`.")

    to_type = os.path.splitext(new_database_name)[1]

    assert to_type in [".json", ".db"]

    from ase.db.jsondb import JSONDatabase
    from ase.db.sqlite import SQLite3Database

    if isinstance(database, str):
        database = connect(database)
    database1_ = connect(new_database_name)

    if isinstance(database, JSONDatabase) and to_type == ".json":
        print("Initial type is json, do nothing.")
        return
    elif isinstance(database, SQLite3Database) and to_type == ".db":
        print("Initial type is db (SQLite3), do nothing.")
        return

    else:
        with database:
            with database1_:
                for ar in database.select():
                    database1_.write(ar, ar.key_value_pairs, ar.data)
        return database1_
