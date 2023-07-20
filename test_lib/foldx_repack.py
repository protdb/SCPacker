import glob
import os.path
import shutil
import subprocess
import sys
import pwd
from pathlib import Path

subj_mask = '_subj'
result_mask = 'Rebuilt_'


class FoldXRepack(object):
    FOLD_X_WORKING_FOLDER = Path("/home/dp/Data/Foldx/")
    FOLD_X_CFG = "fx.cfg"
    USER = 'dp'

    def __init__(self):
        pass

    def repack(self, filepath):
        pdb_id = Path(filepath).stem.replace(subj_mask, '')
        try:
            assert os.path.exists(filepath)
            folder = Path(filepath).parent
            input_filename = f'{pdb_id}.pdb'
            fx_input_filepath = self.FOLD_X_WORKING_FOLDER / input_filename
            shutil.copy(filepath, fx_input_filepath)
            self.make_config(input_filename)
            self.run_command()
            target_file = self.get_results(pdb_id=pdb_id, output_folder=folder)
        except (RuntimeError, Exception, AssertionError) as e:
            target_file = None
        finally:
            self.clean_ws(pdb_id)
        return target_file

    def run_command(self):
        os.setuid(pwd.getpwnam(self.USER).pw_uid)
        os.chdir(str(self.FOLD_X_WORKING_FOLDER))
        fx_exec = self.FOLD_X_WORKING_FOLDER / 'foldx'

        process = subprocess.Popen(f'{fx_exec} -f {self.FOLD_X_WORKING_FOLDER / self.FOLD_X_CFG}',
                                   user='dp',
                                   shell=True)
        process.communicate()

    def make_config(self, filename):
        original_stdout = sys.stdout
        cfg_path = self.FOLD_X_WORKING_FOLDER / self.FOLD_X_CFG

        with open(cfg_path, 'w') as fh:
            sys.stdout = fh
            print('command=ReconstructSideChains')
            print(f'pdb={filename}')

        sys.stdout = original_stdout

    def get_results(self, pdb_id, output_folder):
        result_file = self.FOLD_X_WORKING_FOLDER / f'{result_mask}{pdb_id}.pdb'
        assert os.path.exists(result_file)
        target_file = output_folder / f'{pdb_id}_fx.pdb'
        shutil.move(result_file, target_file)
        return target_file

    def clean_ws(self, pdb_id):
        files = glob.glob(f'{self.FOLD_X_WORKING_FOLDER}/*{pdb_id}*.*')
        cfg_path = self.FOLD_X_WORKING_FOLDER / self.FOLD_X_CFG
        if os.path.exists(cfg_path):
            os.remove(cfg_path)

        for file in files:
            os.remove(file)




# TEST

if __name__ == '__main__':
    testfile = "/home/dp/Data/sch_lib/test/MLY/2zpm/2zpm_subj.pdb"
    test_ = FoldXRepack()
    test_.repack(testfile)
