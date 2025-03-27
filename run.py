import os.path

from packer.mutation import MutationPacker

TEST_INPUT_FILE = '/home/dp/Data/3b_tmp/hd_ib_md2.pdb'
PTM = 'SEP'
CHAIN = 'H'
SITE = 15
OUTPUT_FILE = '/home/dp/Data/3b_tmp/3bwj_reconstruct.pdb'

if __name__ == '__main__':
    assert os.path.exists(TEST_INPUT_FILE)

    print(f'Processing {TEST_INPUT_FILE}...')

    packer = MutationPacker()
    packer.load_from_pdb(TEST_INPUT_FILE)

    '''
        args: position info: array of tuples (AA_CODE, CHAIN', SITE)
        n_samples: MCMC iterations (recommended > 1000, 1000 or 2000 - optimal) 
    '''
    packer.mutate(
                  position_info=[('SEP', 'H', 15), ('TPO', 'L', 150)],
                  n_samples=200,
                  with_library=True,
                  repack_mode='MCMC' # MCMC, GA
                  )
    packer.save(filepath=OUTPUT_FILE)
