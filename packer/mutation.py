import os.path

import chilife as pack

from packer.ga_repack import ga_repack_neighbours
from packer.mcmc_repack import mcmc_repack_neighbours
from packer.ensemble import RTEnsemble


class MutationPacker(object):
    def __init__(self, energy_fnx=None):
        self.protein = None
        self.energy_fxn = energy_fnx if energy_fnx is not None else pack.get_lj_scwrl
        self.repacked_ensembles = None
        self.repacked_neighbors = None

    def load_from_pdb(self, pdb_path):
        assert os.path.exists(pdb_path), f"PDB file not found {pdb_path}"
        self.protein = pack.Protein.from_pdb(str(pdb_path))

    def mutate(self,
               aa_code,
               chain,
               site,
               repack=True,
               n_samples=500,
               with_library=True,
               repack_radius=10.0,
               temp=298,
               repack_mode='MCMC'
               ):
        assert self.protein is not None, "Protein data not loaded"

        mutation_ensemble = RTEnsemble(res=aa_code,
                                       chain=chain,
                                       site=site,
                                       protein=self.protein)

        if repack:
            if repack_mode == 'MCMC':
                self.protein, _ = mcmc_repack_neighbours(self.protein,
                                                         mutation_ensemble,
                                                         energy_func=self.energy_fxn,
                                                         repetitions=n_samples,
                                                         off_rotamer=not with_library,
                                                         repack_radius=repack_radius,
                                                         temp=temp,
                                                         add_missing_atoms=True
                                                         )

            elif repack_mode == 'GA':
                self.protein = ga_repack_neighbours(self.protein,
                                                    mutation_ensemble,
                                                    energy_func=self.energy_fxn,
                                                    repack_radius=repack_radius,
                                                    add_missing_atoms=True,
                                                    init_mode='rotlib' if with_library else 'random'
                                                    )


            else:
                raise NotImplementedError(f'Invalid repack mode: {repack_mode}')

        else:
            self.protein = pack.mutate(self.protein, mutation_ensemble, add_missing_atoms=True)

    def repacked_residues(self):
        try:
            repack_residues = self.protein.repack_residues_ids
        except:
            repack_residues = []

        return repack_residues

    def save(self, filepath):
        pack.save(filepath, self.protein)
