import logging
import os.path

from worker_framework import Runner, TaskTransfer
from worker_framework.tools.workdir import CreateDirMode
from worker_framework.message_managers import HttpPollMsgManager
from pymol import cmd
from packer.mutation import MutationPacker
import requests

from worker_config import ModificationConfig


class ModificationRunner(Runner):
    cfg = ModificationConfig()
    create_dir_mode = CreateDirMode.always
    auto_store_file = False
    msg_manager_class = HttpPollMsgManager
    worker_name = 'scpacker'

    @staticmethod
    def run_mutation(in_path, chain, code, no, out_path):
        packer = MutationPacker()
        packer.load_from_pdb(in_path)
        packer.mutate(
            aa_code=code,
            chain=chain,
            site=no,
            n_samples=200,
            with_library=True,
            repack_mode='MCMC'  # MCMC, GA
        )
        packer.save(out_path)

    def handler(self, task: TaskTransfer) -> TaskTransfer:
        logging.debug(task.json())
        with self.stage_writer('download'):
            if task.model_id > 0:
                task.model_id -= 1
            raw_path = os.path.join(task.workdir, 'raw.pdb')
            in_path = os.path.join(task.workdir, 'backbone.pdb')
            with open(raw_path, 'w') as out_file:
                out_file.write(requests.get(task.url).text)
            cmd.reinitialize()
            cmd.load(raw_path)
            cmd.remove('not polymer')
            cmd.remove('elem H')
            cmd.remove('not backbone and resn HIS')
            resns = []
            cmd.iterate(f'resi {task.amino_replacements[0].no} and n. CA', lambda atom: resns.append(atom.resn))
            source_residue = resns[0]
            cmd.save(in_path)

        logging.debug('Initiating packer')
        with self.stage_writer('refine_intact', True):
            self.run_mutation(
                in_path,
                task.apfid[5],
                source_residue,
                task.amino_replacements[0].no,
                os.path.join(task.workdir, 'reference.pdb')
            )
        with self.stage_writer('apply_modification', True):
            self.run_mutation(
                in_path,
                task.apfid[5],
                task.amino_replacements[0].replace,
                task.amino_replacements[0].no,
                os.path.join(task.workdir, 'modified.pdb')
            )
        return task


if __name__ == "__main__":
    ModificationRunner().run()