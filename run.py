import os.path

from packer.mutation import MutationPacker
# from test_lib.batch_test import batch_test

TEST_INPUT_FILE = 'test_data/1edi.pdb'
PTM = 'PTR'
CHAIN = 'A'
SITE = 12
OUTPUT_FILE = 'test_data/1edi_ptm.pdb'

if __name__ == '__main__':
    assert os.path.exists(TEST_INPUT_FILE)

    print(f'Processing {TEST_INPUT_FILE}...')

    packer = MutationPacker()
    packer.load_from_pdb(TEST_INPUT_FILE)
    print('loaded!')

    packer.mutate(aa_code=PTM,
                  site=SITE,
                  chain=CHAIN,
                  n_samples=200,
                  with_library=True,
                  repack_mode='MCMC' # MCMC, GA
                  )
    packer.save(filepath=OUTPUT_FILE)
