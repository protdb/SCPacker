import numpy as np
from chilife import GAS_CONST, mutate, chilife, RotamerEnsemble
from tqdm import tqdm

from data.templates import NCAA_TEMPLATES_, NON_REPACKED_RESIDUES
from packer.ensemble import RTEnsemble


def mcmc_repack_neighbours(protein,
                           *spin_labels,
                           repetitions=200,
                           temp=50,
                           energy_func,
                           off_rotamer=False,
                           verbose=True,
                           **kwargs):
    temp = np.atleast_1d(temp)
    KT = {t: GAS_CONST * t for t in temp}

    repack_radius = kwargs.pop("repack_radius") if "repack_radius" in kwargs else None  # Angstroms
    if repack_radius is None:
        repack_radius = max([SL.clash_radius for SL in spin_labels])


    spin_label_str = " or ".join(
        f"( {spin_label.selstr} )" for spin_label in spin_labels
    )
    protein = mutate(protein, *spin_labels, **kwargs).atoms

    repack_residues = protein.select_atoms(
        f"(around {repack_radius} {spin_label_str} ) or {spin_label_str}"
    ).residues

    repack_res_kwargs = spin_labels[0].input_kwargs
    repack_res_kwargs['eval_clash'] = False

    supported_residues = chilife.SUPPORTED_RESIDUES.union(set(NCAA_TEMPLATES_))

    repack_residue_libraries = [
        RTEnsemble.from_mda(res, **repack_res_kwargs)
        for res in repack_residues
        if res.resname not in NON_REPACKED_RESIDUES and res.resname in supported_residues
    ]

    protein = mutate(protein, *repack_residue_libraries, **kwargs).atoms

    repack_residues = protein.select_atoms(
        f"(around {repack_radius} {spin_label_str} ) " f"or {spin_label_str}"
    ).residues

    repack_residue_libraries = [
        RTEnsemble.from_mda(res, **repack_res_kwargs)
        for res in repack_residues if res.resname not in ["GLY", "ALA"] and res.resname in supported_residues
    ]

    traj = np.empty((repetitions, *protein.positions.shape))
    deltaEs = []

    sample_freq = np.array([len(res._weights) for res in repack_residue_libraries], dtype=np.float64)
    sample_freq /= sample_freq.sum()

    dummy_ = []

    for res in repack_residue_libraries:
        sl = res.copy()
        sl.protein = protein
        sl.mask = np.isin(protein.ix, sl.clash_ignore_idx)
        sl._coords = np.atleast_3d([protein.atoms[sl.mask].positions])
        dummy_.append(sl)

    count = 0
    a_count = 0
    b_count = 0
    bidx = 0
    schedule = repetitions / (len(temp) + 1)

    repack_residues_ids = [res.site for res in repack_residue_libraries]

    while count < repetitions:
        select_idx = np.random.choice(len(repack_residue_libraries), p=sample_freq)
        SiteLibrary = repack_residue_libraries[select_idx]

        if not hasattr(SiteLibrary, "dummy_label"):
            SiteLibrary.dummy_label = SiteLibrary.copy()
            SiteLibrary.dummy_label.protein = protein
            SiteLibrary.dummy_label.mask = np.isin(protein.ix, SiteLibrary.clash_ignore_idx)
            SiteLibrary.dummy_label._coords = np.atleast_3d([protein.atoms[SiteLibrary.dummy_label.mask].positions])

            with np.errstate(divide="ignore"):
                SiteLibrary.dummy_label.E0 = energy_func(SiteLibrary.dummy_label) - \
                                             KT[temp[bidx]] * np.log(SiteLibrary.current_weight)

        DummyLabel = SiteLibrary.dummy_label
        coords, weight = SiteLibrary.sample(off_rotamer=off_rotamer)

        with np.errstate(divide="ignore"):
            DummyLabel._coords = np.atleast_3d([coords])
            E1 = energy_func(DummyLabel) - KT[temp[bidx]] * np.log(weight)

        deltaE = E1 - DummyLabel.E0
        deltaE = np.maximum(deltaE, -10.0)

        a_count += 1
        # Metropolis-Hastings criteria

        if E1 < DummyLabel.E0 or np.exp(-deltaE / KT[temp[bidx]]) > np.random.rand():
            deltaEs.append(deltaE)

           # dummy_[select_idx] = DummyLabel.copy()
            energies = []

            # for sl in dummy_:
            #     with np.errstate(divide="ignore"):
            #         energy_sl = energy_func(sl)
            #     energies.append(energy_sl)
            #
            # total_energies = np.sum(energies)
            #
            # print(f'Step: {count} LJ: {total_energies}')

            try:
                protein.atoms[DummyLabel.mask].positions = coords
                DummyLabel.E0 = E1

            except ValueError as e:
                raise e

            traj[count] = protein.atoms.positions
            SiteLibrary.update_weight(weight)
            count += 1
            b_count += 1

            if b_count > schedule:
                b_count = 0
                bidx = np.minimum(bidx + 1, len(temp) - 1)
        else:
            continue
    # protein.protein.trajectory = chilife.Trajectory(traj, protein)
    protein.repack_residues_ids = repack_residues_ids.copy()

    return protein, np.squeeze(deltaEs)
