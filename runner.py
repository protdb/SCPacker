import logging
import os.path

from worker_framework import Runner, TaskTransfer
from worker_framework.tools.workdir import CreateDirMode
from worker_framework.message_managers import HttpPollMsgManager
from worker_framework.tools.recorder import Recorder
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

    def handler(self, task: TaskTransfer) -> TaskTransfer:
        logging.debug(task.json())
        self.stage_writer.write_stage('download', 'start')
        if task.model_id > 0:
            task.model_id -= 1
        in_path = os.path.join(task.workdir, 'reference.pdb')
        with open(in_path, 'w') as out_file:
            out_file.write(requests.get(task.url).text)
        cmd.reinitialize()
        cmd.load(in_path)
        cmd.remove('not polymer')
        cmd.save(in_path)
        self.stage_writer.write_stage('download', 'finish', 'OK')
        logging.debug('Initiating packer')
        self.stage_writer.write_stage('mutation', 'start')
        rec = Recorder()
        rec.start()
        packer = MutationPacker()
        packer.load_from_pdb(in_path)
        packer.mutate(
            aa_code=task.amino_replacements[0].replace,
            chain=task.apfid[5],
            site=task.amino_replacements[0].no,
            n_samples = 200,
            with_library = True,
            repack_mode = 'MCMC'  # MCMC, GA
        )
        packer.save(os.path.join(task.workdir, 'modified.pdb'))
        rec.stop()
        self.stage_writer.write_stage('mutation', 'finish', 'OK', rec.get_result())
        return task


if __name__ == "__main__":
    ModificationRunner().run()