import os.path
from pathlib import Path

import chilife as pack
import numpy as np
from chilife import get_lj_energy, get_lj_rep

from config_.config import Config
from packer.mutation import MutationPacker
from test_lib.test_rmx import RMXTest, SAMPLE_PARAMS

N_SAMPLES = 200


def my_sfxn(rotlib):
    weight = 0.078
    f = 0.94
    lj_E = get_lj_energy(rotlib, forgive=f)
    SASA = rotlib.get_sasa()
    return lj_E - weight * SASA


def test_repack(source_file, ptm_code, site, chain):
    assert os.path.exists(source_file)

    rmx = RMXTest(source_file=source_file, chain=chain)
    rmx.create_sample()
    repacked_residues = None

    for sample_type in SAMPLE_PARAMS:
        input_file, output_file = rmx.get_sample_params(mode=sample_type)
        repack_mode, with_library = sample_type

        print(f'Starting {repack_mode} only library: {with_library} ...')

        packer = MutationPacker(energy_fnx=get_lj_energy)
        packer.load_from_pdb(input_file)

        packer.mutate(ptm_code,
                      site=site,
                      chain=chain,
                      n_samples=N_SAMPLES,
                      with_library=with_library,
                      repack_mode=repack_mode
                      )

        packer.save(filepath=output_file)
        repacked_residues = packer.repacked_residues()

    results = rmx.compare(repacked_residues=repacked_residues)

    print('*' * 50)

    for key in results:
        print(f"{key}: rmsd: {results[key]['rmsd']} per_aa: {results[key]['aa_rmsd']}")


test_file = "/home/dp/Data/PDB/4bag.pdb"

if __name__ == '__main__':
    test_repack(source_file=test_file,
                ptm_code='SEP',
                site=954,
                chain='B'
                )
