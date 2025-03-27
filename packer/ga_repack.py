from copy import deepcopy

import numpy as np
from chilife import mutate, SUPPORTED_RESIDUES, batch_ic2cart, GAS_CONST, get_sasa
from deap.tools import mutPolynomialBounded, cxSimulatedBinaryBounded, cxBlend


from data.templates import NCAA_TEMPLATES_, NON_REPACKED_RESIDUES
from packer.ensemble import RTEnsemble
from deap import base
from deap import creator
from deap import tools
from utils.ga_utils import eaSimpleWithElitism, GAIndividual

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", GAIndividual, fitness=creator.FitnessMin)


class GARoutine(object):
    def __init__(self,
                 protein,
                 repack_residues_libraries,
                 energy_func,
                 init_population_mode='rotlib'
                 ):

        assert init_population_mode in ('rotlib', 'random')

        self.repack_residues_libraries = []
        self.init_population_mode = init_population_mode
        self.E0 = []
        self.coo = []

        for x in repack_residues_libraries:
            dummy_ = x
            dummy_.protein = protein
            dummy_.mask = np.isin(protein.ix, x.clash_ignore_idx)
            dummy_.coords = np.atleast_3d([protein.atoms[dummy_.mask].positions])
            self.repack_residues_libraries.append(dummy_)

            with np.errstate(divide="ignore"):
                xE0 = energy_func(dummy_)[0]
            self.E0.append(xE0)
            self.coo.append(dummy_.coords)

        self.toolbox = base.Toolbox()
        self.POPULATION_SIZE = 300
        self.CROWDING_FACTOR = 12.0
        self.HALL_OF_FAME_SIZE = int(0.15 * self.POPULATION_SIZE)
        self.P_CROSSOVER = 0.9
        self.P_MUTATION = 0.3
        self.MAX_GENERATIONS = 100
        self.IND_PB = 0.8
        self.KAPPA_INIT = (10, 60)

        self.dihedral_mask = []
        self.protein = protein
        self.energy_func = energy_func
        self.first_individual = True

    def __len__(self):
        return len(self.repack_residues_libraries)

    def population_creator(self):
        population = []

        from_library = self.init_population_mode == 'rotlib'

        for lib_idx, sl in enumerate(self.repack_residues_libraries):
            dihedral_mask = ([True] * len(sl.dihedral_atoms) if len(sl.dihedral_atoms) > 0 else [True])
            self.dihedral_mask.append(dihedral_mask)

            if from_library:
                coords, _, _, dihedrals, kappas = sl.sample(off_rotamer=True,
                                                            return_dihedrals=True)

                try:
                    _ = len(dihedrals)
                except TypeError:
                    dihedrals = [dihedrals]
                    kappas = [kappas]
            else:
                dh_size = len(dihedral_mask)

                if self.first_individual:
                    dihedrals = np.radians(sl.dihedrals)
                    dihedrals = dihedrals[0]
                    try:
                        _ = len(dihedrals)
                    except TypeError:
                        dihedrals = [dihedrals]

                else:
                    dihedrals = np.random.uniform(-np.pi, np.pi, dh_size)

                kappas = np.random.uniform(*self.KAPPA_INIT, dh_size)

            individuals_pairs = list(zip(dihedrals, kappas))

            for item in individuals_pairs:
                dh, kappas = item
                population.append(np.float32(dh))
                population.append(np.float32(kappas))

        self.first_individual = False
        return population

    def initialize(self):
        self.toolbox.register("PCreator", self.population_creator)
        self.toolbox.register("individualCreator", tools.initIterate, creator.Individual, self.toolbox.PCreator)
        self.toolbox.register("populationCreator", tools.initRepeat, list, self.toolbox.individualCreator)
        self.toolbox.register("evaluate", self.evaluator)

        # genetic operators:
        self.toolbox.register("select", tools.selTournament, tournsize=self.POPULATION_SIZE // 20)
        self.toolbox.register("mate", self.crossover)
        self.toolbox.register("mutate", self.mutate)

    def run(self):
        population = self.toolbox.populationCreator(n=self.POPULATION_SIZE)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", np.min)
        stats.register("avg", np.mean)
        stats.register("max", np.max)

        hof = tools.HallOfFame(self.HALL_OF_FAME_SIZE)

        population, logbook = eaSimpleWithElitism(population,
                                                  self.toolbox,
                                                  cxpb=self.P_CROSSOVER,
                                                  mutpb=self.P_MUTATION,
                                                  ngen=self.MAX_GENERATIONS,
                                                  stats=stats,
                                                  halloffame=hof,
                                                  verbose=True)

        best = hof.items[0]
        print("-- Best Individual = ", best)
        print("-- Best Fitness = ", best.fitness.values[0])

        for i, sl in enumerate(self.repack_residues_libraries):
            dh = best[i]
            dihedral_mask = self.dihedral_mask[i]
            ic = sl._lib_IC[0]
            ic.set_dihedral(dh, 1, sl.dihedral_atoms[dihedral_mask])

            z_mat_ids = np.array([ic.zmat_idxs[1]])
            z_mat = np.array([ic.zmats[1]])
            coo = batch_ic2cart(z_mat_ids, z_mat)
            mx, ori = ic.chain_operators[1]["mx"], ic.chain_operators[1]["ori"]
            coo = np.einsum("ijk,kl->ijl", coo, mx) + ori
            coo = coo[:, sl.ic_mask]
            coo = np.atleast_3d(coo)

            self.protein.atoms[sl.mask].positions = coo

        return self.protein

    def evaluator(self, individuals):
        energies = []
        dh_mask = self._mask_individuals(individuals)
        dh_items = np.array(individuals)[dh_mask].tolist()

        for i, sl in enumerate(self.repack_residues_libraries):
            individual_item = dh_items[i]
            dh = individual_item
            dihedral_mask = self.dihedral_mask[i]
            ic = sl._lib_IC[0]
            back_coo = sl.coords.copy()
            back_dh = ic.get_dihedral(1, sl.dihedral_atoms[dihedral_mask])
            ic.set_dihedral(dh, 1, sl.dihedral_atoms[dihedral_mask])
            z_mat_ids = np.array([ic.zmat_idxs[1]])
            z_mat = np.array([ic.zmats[1]])
            coo = batch_ic2cart(z_mat_ids, z_mat)
            mx, ori = ic.chain_operators[1]["mx"], ic.chain_operators[1]["ori"]
            coo = np.einsum("ijk,kl->ijl", coo, mx) + ori
            coo = coo[:, sl.ic_mask]
            sl.coords = np.atleast_3d(coo)

            with np.errstate(divide="ignore"):
                sl_energy = self.energy_func(sl)[0]

            if self.E0[i] > sl_energy:
                self.E0[i] = sl_energy
                sl.protein.atoms[sl.mask].positions = coo

            sl.coords = back_coo
            ic.set_dihedral(back_dh, 1, sl.dihedral_atoms[dihedral_mask])

            energies.append(sl_energy)
        energy = np.sum(energies)

        return energy,

    def mutate(self, individual):
        ind_values = np.array(individual)
        dh_mask = self._mask_individuals(ind_values)
        individual_dh = ind_values[dh_mask]
        individual_kappa = ind_values[~dh_mask]
        is_mutate_kappa = self.init_population_mode == 'random'

        if np.random.rand() < 0.5:
            individual_dh = mutPolynomialBounded(individual_dh,
                                                 eta=self.CROWDING_FACTOR,
                                                 low=-np.pi,
                                                 up=np.pi,
                                                 indpb=self.IND_PB)

            if is_mutate_kappa:
                individual_kappa = mutPolynomialBounded(individual_kappa,
                                                        eta=self.CROWDING_FACTOR,
                                                        low=self.KAPPA_INIT[0],
                                                        up=self.KAPPA_INIT[1],
                                                        indpb=self.IND_PB
                                                        )
                individual_kappa = individual_kappa[0]

            individual_dh = individual_dh[0]

        else:
            mutation_mask = [np.random.rand() < 0.3 for _ in range(len(individual_dh))]
            mutation_dh = individual_dh[mutation_mask]
            mutation_kappa = individual_kappa[mutation_mask]
            mutation_dh = self.vonMisses_sample(mutation_dh, mutation_kappa)
            individual_dh[mutation_mask] = mutation_dh

        individual_dh = individual_dh
        ind_values[dh_mask] = individual_dh
        ind_values[~dh_mask] = individual_kappa
        individual.reassign(ind_values)

        return individual,

    def crossover(self, ind1, ind2):
        ind_values1 = np.array(ind1)
        ind_values2 = np.array(ind2)

        dh_mask1 = self._mask_individuals(ind_values1)
        dh_mask2 = self._mask_individuals(ind_values2)
        individual_dh1 = ind_values1[dh_mask1]
        individual_dh2 = ind_values2[dh_mask2]
        individual_kappa1 = ind_values1[~dh_mask1]
        individual_kappa2 = ind_values2[~dh_mask2]
        is_cross_over_kappa = self.init_population_mode == 'random'

        if np.random.rand() < 0.5:
            individual_dh1, individual_dh2 = cxSimulatedBinaryBounded(individual_dh1,
                                                                      individual_dh2,
                                                                      eta=self.CROWDING_FACTOR,
                                                                      low=-np.pi,
                                                                      up=np.pi)
            if is_cross_over_kappa:
                individual_kappa1, individual_kappa2 = cxSimulatedBinaryBounded(individual_kappa1,
                                                                                individual_kappa2,
                                                                                eta=self.CROWDING_FACTOR,
                                                                                low=self.KAPPA_INIT[0],
                                                                                up=self.KAPPA_INIT[1])
        else:
            individual_dh1, individual_dh2 = cxBlend(individual_dh1, individual_dh2, alpha=0.5)
            individual_dh1 = self.assert_dihedrals(individual_dh1)
            individual_dh2 = self.assert_dihedrals(individual_dh2)

        ind_values1[dh_mask1] = individual_dh1
        ind_values2[dh_mask2] = individual_dh2

        ind_values1[~dh_mask1] = individual_kappa1
        ind_values2[~dh_mask2] = individual_kappa2

        ind1.reassign(ind_values1)
        ind2.reassign(ind_values2)

        return ind1, ind2

    @staticmethod
    def assert_dihedrals(x):
        low_mask = x < -np.pi
        x[low_mask] = -np.pi
        high_mask = x > np.pi
        x[high_mask] = np.pi
        return x

    @staticmethod
    def vonMisses_sample(dh, kappa):
        new_dh = np.random.vonmises(dh, kappa)
        return new_dh

    @staticmethod
    def _mask_individuals(individual):
        dh_mask = [i % 2 == 0 for i in range(len(individual))]
        return np.array(dh_mask)


def ga_repack_neighbours(protein,
                         *spin_labels,
                         repack_radius,
                         init_mode,
                         energy_func,
                         **kwargs):
    protein = mutate(protein,
                     *spin_labels,
                     add_missing_atoms=True,
                     random_rotamers=False
                     ).atoms

    spin_label_str = " or ".join(
        f"( {spin_label.selstr} )" for spin_label in spin_labels
    )

    repack_residues = protein.select_atoms(
        f"(around {repack_radius} {spin_label_str} ) or {spin_label_str}"
    ).residues

    repack_res_kwargs = {'eval_clash': False}
    supported_residues = SUPPORTED_RESIDUES.union(set(NCAA_TEMPLATES_))

    repack_residue_libraries = [RTEnsemble.from_mda(res, **repack_res_kwargs)
                                for res in repack_residues
                                if res.resname not in NON_REPACKED_RESIDUES and res.resname in supported_residues
                                ]

    repack_residues_ids = [res.site for res in repack_residue_libraries]

    ga_routine = GARoutine(protein=protein,
                           repack_residues_libraries=repack_residue_libraries,
                           energy_func=energy_func,
                           init_population_mode=init_mode
                           )
    ga_routine.initialize()
    updated_protein = ga_routine.run()
    updated_protein.repack_residues_ids = repack_residues_ids

    return updated_protein
