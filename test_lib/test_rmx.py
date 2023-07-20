import os.path
import shutil
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile

from config_.config import Config
from data.templates import NCAA_TEMPLATES_
from utils.pdb_utils import read_pdb

BACKBONE_MASK = ['N', 'CA', 'C', 'O']
backbone_len = 4
HYDROGEN = 'H'

SAMPLE_NAMES = ['source', 'subj']
SAMPLE_PARAMS = [('MCMC', False), ('MCMC', True), ('GA', True), ('GA', False)]

#SAMPLE_PARAMS = [('GA', False)]


class RMXTest(object):
    def __init__(self,
                 source_file,
                 chain,
                 output_folder=None,
                 ):

        self.source_file = source_file
        self.chain = chain
        config = Config()
        output_folder = config.test_folder if not output_folder else output_folder
        self.output_folder = output_folder / Path(source_file).stem
        self.output_folder.mkdir(exist_ok=True)
        self.test_files = {}

    def remove_results(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(str(self.output_folder))

    def create_sample(self):
        assert os.path.exists(self.source_file)
        structure, backbone_mask = self.get_file_records(self.source_file)
        full_mask = backbone_mask | ~backbone_mask

        source_structure = structure[full_mask]
        subj_structure = structure[backbone_mask]

        for i, structure in enumerate([source_structure, subj_structure]):
            outfile = self.output_folder / f'{Path(self.source_file).stem}_{SAMPLE_NAMES[i]}.pdb'
            outData = PDBFile()
            outData.set_structure(structure)
            outData.write(outfile)
            self.test_files.update({SAMPLE_NAMES[i]: str(outfile)})

    def get_sample_params(self, mode):
        assert mode in SAMPLE_PARAMS

        sample_type, with_library = mode
        file_postfix = sample_type if not with_library else f'{sample_type}_lib'

        input_file = self.test_files['subj']
        output_file = self.output_folder / f'{Path(input_file).stem}_{file_postfix}.pdb'

        self.test_files.update({file_postfix: output_file})

        return input_file, output_file

    def get_file_records(self, file):
        assert os.path.exists(file)
        structure = read_pdb(file)
        hyd_mask = structure.element == HYDROGEN
        structure = structure[~hyd_mask]
        chains = structure.chain_id
        chain_mask = chains == self.chain

        het_mask = structure.hetero

        ptm_mask = np.zeros_like(het_mask)
        res_names = structure.res_name

        for i, res_name in enumerate(res_names):
            if res_name in NCAA_TEMPLATES_:
                ptm_mask[i] = True

        het_mask *= ~ptm_mask
        mask = ~het_mask * chain_mask
        structure = structure[mask]
        atom_names = structure.atom_name
        backbone_mask = np.zeros(len(atom_names), dtype=bool)

        for i in range(len(atom_names)):
            atom = atom_names[i]

            if atom in BACKBONE_MASK:
                backbone_mask[i] = True

        return structure, backbone_mask

    def compare(self, repacked_residues=None):
        files = self.test_files.copy()
        source_file = files['source']

        del files['subj']
        del files['source']

        results = {}

        for key in files:
            file_x = files[key]
            rmsd, aa_rms = self.pairwise_distance(source_file, file_x, repacked_residues)
            results.update({key: {'rmsd': rmsd, 'aa_rmsd': aa_rms}})
        return results

    def pairwise_distance(self, source_file, file_x, repacked_residues):
        source_structure, source_backbone_mask = self.get_file_records(source_file)
        x_structure, x_backbone_mask = self.get_file_records(file_x)

        source_res_ids = source_structure.res_id
        x_res_ids = x_structure.res_id

        distances = []
        aa_rms = []

        for res_id in repacked_residues:
            source_res_mask = source_res_ids == res_id
            x_res_mask = x_res_ids == res_id
            res_name = source_structure.res_name[source_res_mask]
            source_sch_mask = ~source_backbone_mask[source_res_mask]
            x_sch_mask = ~x_backbone_mask[x_res_mask]
            source_xyz = source_structure.coord[source_res_mask][source_sch_mask]
            x_xyz = x_structure.coord[x_res_mask][x_sch_mask]

            if len(source_xyz) == len(x_xyz) and len(x_xyz) > 0:
                dist = np.sum(np.linalg.norm(source_xyz - x_xyz, axis=1)) / len(source_xyz)
                distances.append(dist)
                aa_rms.append((set(res_name).pop(), res_id, dist))

        assert distances

        return np.mean(distances), aa_rms
