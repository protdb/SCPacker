import os.path
from pathlib import Path

import chilife as pack
import MDAnalysis as mda
from chilife import fmt_str

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
               position_info:list,
               repack=True,
               n_samples=500,
               with_library=True,
               repack_radius=10.0,
               temp=298,
               repack_mode='MCMC'
               ):
        assert self.protein is not None, "Protein data not loaded"

        for record in position_info:
            aa_code, chain, site = record


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
        protein = self.protein

        with open(filepath, 'a+') as pdb_file:

            if isinstance(protein, (mda.AtomGroup, mda.Universe)):
                traj = protein.universe.trajectory
                name = Path(protein.universe.filename) if protein.universe.filename is not None else pdb_file.name
                name = name.name
            else:
                traj = protein.trajectory
                name = protein.fname

            if name is None:
                name = Path(pdb_file.name).name

            name = name[:-4] if name.endswith(".pdb") else name

            pdb_file.write(f'HEADER {name}\n')
            for mdl, ts in enumerate(traj):
                pdb_file.write(f"MODEL {mdl}\n")
                [
                    pdb_file.write(
                        fmt_str.format(
                            atom.index,
                            atom.name,
                            atom.resname[:3],
                            atom.segid,
                            atom.resnum,
                            *atom.position,
                            1.00,
                            1.0,
                            atom.type,
                        )
                    )
                    for atom in protein.atoms
                ]
                pdb_file.write("TER\n")
                pdb_file.write("ENDMDL\n")
