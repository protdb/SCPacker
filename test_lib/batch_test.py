import glob
import os.path
import pickle
import shutil
import time
from pathlib import Path

import numpy as np
from chilife import get_lj_energy, get_sasa, get_lj_rep, get_lj_scwrl
from numba import njit

from config_.config import Config
from data.templates import NCAA_TEMPLATES_
from packer.mutation import MutationPacker
from test_lib.foldx_repack import FoldXRepack
from test_lib.rosetta_repack import RosettaRepack
from test_lib.statistic import StatLogger
from test_lib.test_rmx import RMXTest, SAMPLE_PARAMS


def custom_sfx(rotlib):
    weight = 0.0178
    f = 0.94
    lj_E = get_lj_energy(rotlib, forgive=f)
    sasa = rotlib.get_sasa()
    return lj_E - weight * sasa


class BatchTest(object):
    REPETITION = 200

    def __init__(self, n_process):
        self.n_process = n_process

        self.config = Config()
        test_files = glob.glob(f'{self.config.dataset_folder}/*')
        self.test_data = {}

        for filepath in test_files:
            with open(filepath, 'rb') as fh:
                data = pickle.load(fh)
                self.test_data.update(data)

        for key in self.test_data.copy():
            if key not in NCAA_TEMPLATES_:
                del self.test_data[key]

        self.logger = StatLogger()

        self.rose_repack = RosettaRepack()
        self.fx_repack = FoldXRepack()

    def test_ptm(self, ptm):
        assert ptm in self.test_data

        output_folder = self.config.test_folder / ptm
        output_folder.mkdir(exist_ok=True)
        test_data = self.test_data[ptm]
        rnd_idx = np.random.randint(0, len(test_data))
        test_record = test_data[rnd_idx]
        test_file_id, ptm_points = test_record
        test_file = Path(self.config.pdb_folder) / f'{test_file_id}.pdb'
        assert os.path.exists(test_file)
        ptm_point = np.random.choice(ptm_points, 1)[0]
        chain, site = ptm_point.split('_')

        assert not self.logger.is_exist(pdb_id=test_file_id, chain=chain, site=site, ptm=ptm)

        results = self.repack(
            source_file=test_file,
            output_folder=output_folder,
            ptm=ptm,
            chain=chain,
            site=site,
            n_repetition=self.REPETITION,
            verbose=True,

        )

        assert results

        self.logger.add_results(
            pdb_id=test_file_id,
            chain=chain,
            site=site,
            ptm=ptm,
            results=results
        )

    def process_batch(self):
        n_processed = 0

        while n_processed <= self.n_process:
            ptm = np.random.choice(['TPO'], 1)[0]#np.random.choice(list(self.test_data.keys()), 1)[0]
            try:
                self.test_ptm(ptm)
                n_processed += 1
            except (AssertionError, Exception, RuntimeError) as e:
                continue

    def repack(self,
               source_file,
               output_folder,
               ptm,
               chain,
               site,
               n_repetition,
               verbose=True,
               ):

        rmx = None
        results = {}

        try:
            assert os.path.exists(source_file)

            rmx = RMXTest(source_file=source_file,
                          output_folder=output_folder,
                          chain=chain
                          )
            rmx.create_sample()

            if verbose:
                print(f'Starting {Path(source_file).stem} PTM: {ptm} chain: {chain} site: {site}')

            repacked_residues = []

            times = {}

            for ic, sample_type in enumerate(SAMPLE_PARAMS):

                input_file, output_file = rmx.get_sample_params(mode=sample_type)
                repack_mode, with_library = sample_type
                start_time = time.time()

                if verbose:
                    print(f'Starting {repack_mode} only library: {with_library} ...')

                packer = MutationPacker(energy_fnx=get_lj_scwrl)
                packer.load_from_pdb(input_file)
                packer.mutate(ptm,
                              site=site,
                              chain=chain,
                              n_samples=n_repetition,
                              with_library=with_library,
                              repack_mode=repack_mode
                              )
                repacked_residues = packer.repacked_residues()
                packer.save(filepath=output_file)
                end_time = time.time()
                times.update({sample_type: end_time - start_time})

                if ic == 0:
                    # rosetta repack
                    print('Rosetta repack')
                    start_time = time.time()
                    subj_file = rmx.test_files['subj']
                    rose_file = self.rose_repack.repack(
                        subj_file=subj_file,
                        site=site,
                        repack_residues=repacked_residues
                    )
                    rmx.test_files.update({'rosetta': rose_file})
                    end_time = time.time()
                    times.update({'rosetta': (end_time - start_time)})

            # start_time = time.time()
            # fx_file = self.fx_repack.repack(subj_file)
            # end_time = time.time()
            # times.update({'foldx': end_time - start_time})

            # rmx.test_files.update({'fx': fx_file})

            results = rmx.compare(repacked_residues=repacked_residues)

            assert results
            self.logger.log_msg({f'LEN:{len(repacked_residues)}'})
            self.logger.log_times(times)

            if verbose:
                print('*' * 50)
                for key in results:
                    print(f"{key}: rmsd: {results[key]['rmsd']} per_aa: {results[key]['aa_rmsd']}")

        except (AssertionError, Exception, RuntimeError) as e:
            self.logger.log_msg(e)

            if rmx is not None:
                rmx.remove_results()

            raise e

        return results

    def print_samples(self):
        st = {}

        for ptm in self.test_data:
            data = self.test_data[ptm]
            samples_ids = np.random.choice(len(data), 20)
            st.update({ptm: []})
            for i in samples_ids:
                test_ = data[i]
                pdb_id, sites = test_
                site = np.random.choice(np.array(sites), 1)[0]
                st[ptm].append((pdb_id, site))

        ic = 1
        for k in st:
            st_data = st[k]

            for item in st_data:
                pdb_id, site = item
                site = site.split('_')

                print(f'{ic} {k} {pdb_id} {site[0]} {site[1]}')

                ic += 1


def batch_test():
    test_ = BatchTest(n_process=20)
    test_.print_samples()

    #test_.process_batch()


if __name__ == '__main__':
    batch_test()
