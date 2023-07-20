import pickle

import numpy as np

from config_.config import Config
from Bio.PDB.Polypeptide import standard_aa_names

from data.templates import NCAA_TEMPLATES_


class TestMetrics(object):
    def __init__(self):
        config = Config()
        stat_file = config.stat_file

        with open(stat_file, 'rb') as fh:
            self.stat_data = pickle.load(fh)

        aa_keys = list(self.stat_data.keys()) + list(standard_aa_names) + ['TOTAL']
        samples = self.stat_data[aa_keys[0]]
        samples_keys = [list(s[1].keys()) for s in samples]
        samples_type = set([item for sublist in samples_keys for item in sublist])

        self.rms_metrics = {s: {aa: [] for aa in aa_keys} for s in samples_type}

    def group_metrics(self):
        for ptm in self.stat_data:
            items = self.stat_data[ptm]

            for rec in items:
                pdb_id, stat = rec
                for sample_type in stat:
                    total_rms = stat[sample_type]['rmsd']
                    self.rms_metrics[sample_type]['TOTAL'].append(total_rms)
                    aa_rmsd = stat[sample_type]['aa_rmsd']
                    for aa_rec in aa_rmsd:
                        aa_name, _, aa_rm = aa_rec
                        self.rms_metrics[sample_type][aa_name].append(aa_rm)
                    print(pdb_id, sample_type, total_rms)

    def calc_rms(self):
        self.group_metrics()

        for sample_type in self.rms_metrics:
            print(sample_type)
            for aa in self.rms_metrics[sample_type]:

                data = self.rms_metrics[sample_type][aa]
                if not data:
                    continue
                mean = np.round(np.mean(data), 2)
                std = np.round(np.std(data), 2)
                count = len(data)
                print(f'{aa} {count} {mean}')

    def calc_err(self):
        tres_ = {'SEP': 1.5, 'TPO': 1.5, 'PTR': 1.5, 'MLY': 3.0, 'CSO': 2.5}

        self.group_metrics()

        for sample_type in self.rms_metrics:
            print(sample_type)
            for aa in self.rms_metrics[sample_type]:
                if aa not in NCAA_TEMPLATES_:
                    continue
                data = self.rms_metrics[sample_type][aa]
                if not data:
                    continue
                count = len(data)
                true = np.array(data) <= tres_[aa]
                correct = np.sum(true) / count
                print(f'{aa}, {(1.0 - correct) * 100}')


def calc_metrics():
    metrics_ = TestMetrics()
    metrics_.calc_err()


if __name__ == '__main__':
    calc_metrics()
