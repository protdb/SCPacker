import os.path

import numpy as np
from chilife import RotamerEnsemble, read_library, batch_ic2cart
import chilife as packer

from config_.config import Config
from data.templates import NCAA_TEMPLATES_


class RTEnsemble(RotamerEnsemble):
    def __init__(self, res, site, protein, chain, rotlib=None, **kwargs):

        self.Psi = None
        self.Phi = None
        self.config_ = Config()
        self.custom_library = True if res in NCAA_TEMPLATES_ else False
        self.residues_type = res

        super().__init__(res=res,
                         site=site,
                         protein=protein,
                         chain=chain,
                         rotlib=rotlib,
                         **kwargs)

    def get_lib(self, rotlib):
        if not self.custom_library:
            return super().get_lib(rotlib)

        rotlib_path = self.config_.custom_rotlib_folder / f'{self.residues_type}_rotlib.npz'

        assert os.path.exists(rotlib_path), f"Can't find library for {self.residues_type}"

        phi_sel, psi_sel = None, None

        if self.protein is not None:
            phi_sel = self.protein.select_atoms(f"resnum {self.site} and segid {self.chain}").residues[
                0].phi_selection()
            psi_sel = self.protein.select_atoms(f"resnum {self.site} and segid {self.chain}").residues[
                0].psi_selection()

        phi = None if phi_sel is None else np.rad2deg(packer.get_dihedral(phi_sel.positions))
        psi = None if psi_sel is None else np.rad2deg(packer.get_dihedral(psi_sel.positions))

        self.Phi, self.Psi = phi, psi

        library = read_library(rotlib_path, phi, psi)

        library = {key: value.copy() if hasattr(value, 'copy') else value for key, value in library.items()}
        library['internal_coords'] = [a.copy() for a in library['internal_coords']]
        library['_coords'] = library.pop('coords')
        library['_dihedrals'] = library.pop('dihedrals')

        if 'skews' not in library:
            library['skews'] = None

        return library

    def _base_copy(self, rotlib=None):
        return RTEnsemble(
            res=self.res,
            site=self.site,
            protein=self.protein,
            chain=self.chain,
            rotlib=rotlib)

    def sample(self, n=1, off_rotamer=False, **kwargs):
        return super().sample(n, off_rotamer=off_rotamer, **kwargs)

    def _off_rotamer_sample(self, idx, off_rotamer, **kwargs):

        if len(self._weights) == 1 or np.all(np.isinf(self.sigmas)):
            new_dihedrals = np.random.random((len(idx), len(off_rotamer))) * 2 * np.pi
            new_weights = np.ones(len(idx))

        else:
            new_dihedrals = np.random.vonmises(self._rdihedrals[idx][:, off_rotamer],
                                               self._rkappas[idx][:, off_rotamer])
            new_weights = np.ones(len(idx))

        new_weights = self._weights[idx] * new_weights

        internal_coords = [self._lib_IC[i_idx].copy() for i_idx in idx]
        internal_coords = [IC.set_dihedral(new_dihedrals[i], 1,
                                           self.dihedral_atoms[off_rotamer]) for i, IC in enumerate(internal_coords)]

        Zmat_idxs = np.array([IC.zmat_idxs[1] for IC in internal_coords])
        Zmats = np.array([IC.zmats[1] for IC in internal_coords])
        coords = batch_ic2cart(Zmat_idxs, Zmats)
        mx, ori = internal_coords[0].chain_operators[1]["mx"], internal_coords[0].chain_operators[1]["ori"]
        coords = np.einsum("ijk,kl->ijl", coords, mx) + ori

        for i, IC in enumerate(internal_coords):
            IC._coords = coords[i]

        coords = coords[:, self.ic_mask]

        if kwargs.setdefault("return_dihedrals", False):
            return coords, new_weights, internal_coords, new_dihedrals, self._rkappas[idx][:, off_rotamer]
        else:
            return coords, new_weights
