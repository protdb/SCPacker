

NCAA_TEMPLATES_ = {'SEP': {'atoms': ['N', 'CA', 'C', 'O', 'CB', 'OG', 'P', 'O1P', 'O2P', 'O3P'],
                           'dihedral_atoms': [['N', 'CA', 'CB', 'OG'],
                                              ['CA', 'CB', 'OG', 'P'],
                                              ['CB', 'OG', 'P', 'O1P'],
                                              ['CB', 'OG', 'P', 'O2P'],
                                              ['CB', 'OG', 'P', 'O3P']],
                           'spin_atoms': ['OG', 'P', 'O1P', 'O2P', 'O3P']},
                   'TPO': {'atoms': ['N', 'CA', 'C', 'O', 'CB', 'CG2', 'OG1', 'P', 'O1P', 'O2P', 'O3P'],
                           'dihedral_atoms': [['N', 'CA', 'CB', 'OG1'],
                                              ['CA', 'CB', 'OG1', 'P'],
                                              ['CB', 'OG1', 'P', 'O1P'],
                                              ['CB', 'OG1', 'P', 'O2P'],
                                              ['CB', 'OG1', 'P', 'O3P']],
                           'spin_atoms': ['OG', 'P', 'O1P', 'O2P', 'O3P']
                           },

                   'PTR': {'atoms': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
                                     'OH', 'P', 'O1P', 'O2P', 'O3P'],
                           'dihedral_atoms': [['N', 'CA', 'CB', 'CG'],
                                              ['CA', 'CB', 'CG', 'CD1'],
                                              ['CE2', 'CZ', 'OH', 'P'],
                                              ['CZ', 'OH', 'P', 'O1P']],
                           'spin_atoms': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
                                          'OH', 'P', 'O1P', 'O2P', 'O3P']},
                   'MLY': {'atoms': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'CH1', 'CH2'],
                           'dihedral_atoms': [['N', 'CA', 'CB', 'CG'],
                                              ['CA', 'CB', 'CG', 'CD'],
                                              ['CB', 'CG',  'CD', 'CE'],
                                              ['CG', 'CD', 'CE', 'NZ'],
                                              ['CD', 'CE', 'NZ', 'CH1']],
                           'spin_atoms': ['CG', 'CD', 'CE', 'NZ', 'CH1', 'CH2']},
                   'CSO': {
                       'atoms': ['N', 'CA', 'CB', 'SG', 'C', 'O', 'OD'],
                       'dihedral_atoms': [['N', 'CA', 'CB', 'SG'],
                                          ['CA', 'CB', 'SG', 'OD']],
                       'spin_atoms': ['SG', 'OD']
                   }
                   }

NON_REPACKED_RESIDUES = ["GLY", "ALA"]
