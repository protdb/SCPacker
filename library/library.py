import glob
import os.path
import pickle
from pathlib import Path

import numpy as np
from biotite.structure import stack, AtomArray, from_template
from biotite.structure.io.pdb import PDBFile

from config_.config import Config
from data.templates import NCAA_TEMPLATES_
from utils.pdb_utils import read_pdb
import chilife as pack


def build_ptm_lib(ptm_type):
    config = Config()
    pdb_path = Path(config.pdb_folder)

    with open(config.candidates_file, 'rb') as fh:
        data = pickle.load(fh)

    data = data[ptm_type]
    struct_arr = []
    out_pdb = PDBFile()
    output_file = config.lib_folder / f'{ptm_type}.pdb'
    template_len = len(NCAA_TEMPLATES_[ptm_type]['atoms'])
    template = AtomArray(template_len)
    template.chain_id = ['A' for _ in range(template_len)]
    template.element = [t[0] for t in NCAA_TEMPLATES_[ptm_type]['atoms']]
    template.res_id = [1 for _ in range(template_len)]
    template.atom_name = NCAA_TEMPLATES_[ptm_type]['atoms']

    for i, item in enumerate(data):
        file_id, ptm_list = item
        file = pdb_path / f'{file_id}.pdb'

        if not os.path.exists(file):
            continue

        filedata = read_pdb(file)
        try:
            for ptm_data in ptm_list:
                chain, res_id = ptm_data.split('_')
                chain_mask = filedata.chain_id == chain
                res_mask = filedata.res_id == int(res_id)
                ptm_mask = chain_mask * res_mask
                residues_name = filedata.res_name[ptm_mask]

                test = set(residues_name)
                assert len(test) == 1
                test = test.pop()
                assert test == ptm_type
                model = filedata[ptm_mask]
                assert template_len == len(model.atom_name)

                is_success = True
                coo_arr = {a: [] for a in NCAA_TEMPLATES_[ptm_type]['atoms']}

                for idx, atom in enumerate(model.atom_name):
                    if atom == 'H':
                        continue
                    if atom not in coo_arr:
                        is_success = False
                        break

                    coo_arr[atom].append(model.coord[idx])

                assert is_success
                coords = []

                for a in coo_arr:
                    if not coo_arr[a]:
                        is_success = False
                        break
                    if len(coo_arr[a]) != 1:
                        is_success = False
                        break
                    coords.append(*coo_arr[a])
                assert is_success
                coords = np.array(coords)
                template_model = template.copy()
                template_model.coord = coords.copy()
                struct_arr.append(template_model)

        except AssertionError:
            continue

    stack_ = stack(struct_arr)
    out_pdb.set_structure(stack_)
    out_pdb.write(str(output_file))


def build_rot_lib(ptm_type, site=1):
    config = Config()
    pdb_file = config.lib_folder / f'{ptm_type}.pdb'
    assert os.path.exists(pdb_file)

    print(f"Building rotamer for {pdb_file}")

    pack.create_library(libname=ptm_type,
                        resname=ptm_type,
                        site=site,
                        pdb=str(pdb_file),
                        dihedral_atoms=NCAA_TEMPLATES_[ptm_type]['dihedral_atoms'],
                        spin_atoms=NCAA_TEMPLATES_[ptm_type]['spin_atoms'])


def make_library(force=False):
    config = Config()
    lib_files = glob.glob(f'{config.custom_rotlib_folder}/*.npz')
    exists_libs = [Path(f).stem.split('_')[0] for f in lib_files]
    to_build = set(NCAA_TEMPLATES_.keys()).difference(exists_libs) if not force else NCAA_TEMPLATES_

    for ptm in to_build:
        build_ptm_lib(ptm_type=ptm)
        build_rot_lib(ptm)


if __name__ == '__main__':
    make_library()
