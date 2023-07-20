import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover

subj_mask = '_subj'


class RosettaRepack(object):
    def __init__(self):
        pyrosetta.init(extra_options="-mute core -mute basic")
        self.mm_std_sf = get_fa_scorefxn()

    def repack(self, subj_file, site, repack_residues):
        assert os.path.exists(subj_file)

        input_file = str(subj_file)
        target_file = input_file.replace(subj_mask, f'{subj_mask}_rose')
        initial_pose = pose_from_pdb(input_file)
        info = initial_pose.pdb_info()

        repacked_residue_map = []
        site_map = -1
        repack_residues = [str(r) for r in repack_residues]

        for i in range(info.nres()):
            num = str(info.number(i))

            if num == site:
                site_map = i
                continue

            if num in repack_residues:
                repacked_residue_map.append(i)

        assert site_map != -1
        assert repacked_residue_map

        task = pyrosetta.standard_packer_task(initial_pose)
        task.restrict_to_repacking()

        for i in range(1, len(initial_pose.residues) + 1):
            if i == site:
                continue

            if i in repacked_residue_map:
                task.nonconst_residue_task(i).restrict_to_repacking()
            else:
                task.nonconst_residue_task(i).prevent_repacking()

        pack = RotamerTrialsMover(self.mm_std_sf, task)

        pack.apply(initial_pose)

        pyrosetta.dump_pdb(initial_pose, target_file)

        return target_file
