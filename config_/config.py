import os.path
from pathlib import Path


class Config(object):
    def __init__(self):
        self.pdb_folder = "/home/dp/Data/PDB"
        self.working_folder = "/home/dp/Data/sch_lib"
        self.log_folder = "log"
        self.lib_folder = 'library'
        self.test_folder = 'test'
        self.dataset_folder = 'dataset'
        self.stat_file = 'stats.pkl'
        self.test_log_file = 'processed.log'
        self.log_every = 2
        self.candidates_file = 'candidates.pkl'
        self.custom_rotlib_folder = Path(__file__).parent.parent / 'data' / 'custom_lib'
        self.create_workspace()

    def create_workspace(self):
        self.working_folder = Path(self.working_folder)
        assert os.path.exists(self.working_folder)

        self.log_folder = self.working_folder / self.log_folder
        self.log_folder.mkdir(exist_ok=True)

        # Dataset folders
        self.dataset_folder = self.working_folder / self.dataset_folder
        self.dataset_folder.mkdir(exist_ok=True)
        self.candidates_file = self.dataset_folder / self.candidates_file

        self.lib_folder = self.working_folder / self.lib_folder
        self.lib_folder.mkdir(exist_ok=True)

        self.test_folder = self.working_folder / self.test_folder
        self.test_folder.mkdir(exist_ok=True)

        self.stat_file = self.log_folder / self.stat_file
        self.test_log_file = self.log_folder / self.test_log_file

