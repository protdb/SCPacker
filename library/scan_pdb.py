import glob
import os
import pickle
from pathlib import Path

import numpy as np
from tqdm import tqdm
from config_.config import Config
from utils.pdb_utils import read_pdb

SCAN_PTM_SET = {'SAH', 'MLY', 'ALY', 'SAM', 'CSO', 'HYP'}


class PDBDatabaseScan(object):
    def __init__(self):
        config = Config()
        self.candidates_file = config.candidates_file
        pdb_folder = config.pdb_folder
        assert os.path.exists(pdb_folder), f"Can't load PDB folder: {pdb_folder}"

        self.pdb_files = glob.glob(f'{pdb_folder}/*')
        self.candidates = {code: [] for code in SCAN_PTM_SET}

    def scan_database(self):
        for file in tqdm(self.pdb_files, total=(len(self.pdb_files))):
            try:
                structure = read_pdb(file)
                residues_names = structure.res_name
            except (AssertionError, Exception) as e:
                continue

            intersection = set(residues_names).intersection(SCAN_PTM_SET)

            if not intersection:
                continue

            for res_ in intersection:
                mask = residues_names == res_
                try:
                    residues_ids = structure.res_id
                    chains = structure.chain_id
                    candidate_residues_ids = residues_ids[mask]
                    candidate_chains = chains[mask]
                    assert len(candidate_chains) == len(candidate_residues_ids)
                    positions = map(lambda x: f'{x[0]}_{x[1]}', zip(candidate_chains, candidate_residues_ids))
                    positions = sorted(set(positions))
                    pdb_id = Path(file).stem
                    self.candidates[res_].append((pdb_id, positions))
                except (AssertionError, Exception):
                    continue

    def save_results(self):
        with open(self.candidates_file, 'wb') as fh:
            pickle.dump(self.candidates, fh)


def scan_pdb_database():
    scanner = PDBDatabaseScan()
    scanner.scan_database()
    scanner.save_results()


if __name__ == '__main__':
    scan_pdb_database()
