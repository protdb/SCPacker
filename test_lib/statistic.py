import os.path
import pickle
import sys

from config_.config import Config
from data.templates import NCAA_TEMPLATES_


class StatLogger(object):
    def __init__(self):
        self.config = Config()
        self.records = {k: [] for k in NCAA_TEMPLATES_}
        self.processed = []
        self.n_total = 0

        self.load()

    def load(self):
        if not os.path.exists(self.config.stat_file):
            return

        with open(self.config.stat_file, 'rb') as fh:
            self.records = pickle.load(fh)

        print('Loaded...')

        for ptm in self.records:
            records = self.records[ptm]
            for item in records:
                key, _ = item
                self.processed.append(key)

            print(f'{ptm}: {len(records)}')

    def add_results(self,
                    pdb_id,
                    chain,
                    site,
                    ptm,
                    results,
                    ):
        key = f'{pdb_id}_{chain}_{site}_{ptm}'
        self.processed.append(key)

        self.records[ptm].append((key, results))
        self.n_total += 1

        self.dump_results()

        if self.n_total % self.config.log_every == 0:
            self.dump_log()

    def dump_results(self):
        with open(self.config.stat_file, 'wb') as fh:
            pickle.dump(self.records, fh)

    def log_times(self, times_):
        original_stdout = sys.stdout
        with open(self.config.test_log_file, 'a') as fh:
            sys.stdout = fh

            for key in times_:
                print(f'{key}: {times_[key] / 60:2f}')

            print()

        sys.stdout = original_stdout

    def log_msg(self, msg):
        original_stdout = sys.stdout

        with open(self.config.test_log_file, 'a') as fh:
            sys.stdout = fh
            print(msg)

            print()

        sys.stdout = original_stdout

    def is_exist(self,
                 pdb_id,
                 chain,
                 site,
                 ptm):
        key = f'{pdb_id}_{chain}_{site}_{ptm}'

        return key in self.processed

    def dump_log(self):
        original_stdout = sys.stdout

        with open(self.config.test_log_file, 'a') as fh:
            sys.stdout = fh

            print(f'Total processed: {self.n_total}')

            for key in self.records:
                print(f'{key}: {len(self.records[key])}')

            print()

        sys.stdout = original_stdout
