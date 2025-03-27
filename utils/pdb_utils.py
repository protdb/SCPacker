from biotite.structure.io.pdb import PDBFile


def read_pdb(filepath):
    data = PDBFile.read(filepath)
    structure = data.get_structure(model=1)
    # chains = structure.chain_id
    # residues_ids = structure.res_id
    # residues_names = structure.res_name
    # assert len(chains) == len(residues_ids) == len(residues_names)
    #
    # atom_names = structure.atom_name
    # atom_coo = structure.coord
    # elements = structure.element
    #
    # assert len(atom_names) == len(atom_coo) == len(residues_names) == len(elements)

    return structure
