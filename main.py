"""
Shridhar Athinarayanan
Metropolis-Monte Carlo using Spring Potential Energy Function
(CS279 - requires CS279 PyRosetta Built Environment)

$ python protein_folder.py <input PDB file> <atom stripped XYZ file> <output XYZ file>
Example: $ python protein_folder.py hras.pdb hras_init.xyz hras_out.xyz
"""

import utils
import monte_carlo



def main():
    import sys, os

    pdbfile = sys.argv[1]
    name = pdbfile[0: pdbfile.index('.')]

    atoms, x, y, z = utils.get_coordinates_from_pdb(pdbfile, name)
    atoms, x, y, z = utils.filter_by_name(atoms, x, y, z, 'CA')  # filter to just get the alpha carbons

    utils.write_xyz_file(atoms, x, y, z, name,
                         sys.argv[2])  # write coordinates to an XYZ file which can be converted to PDB


    iters = int(input("Enter amount of Monte Carlo iterations: "))
    k = int(input("Enter spring constant (how rigid is the polymer?): "))
    kBt = int(input("Enter kbT value (Boltzmann constant times temperature of system): "))
    rad = int(input("Enter how far to move atoms radially: "))

    utils.graph_protein(x, y, z, 'Initial Protein')  # graph initial protein
    monte_carlo.monte_carlo(iters, k, kBt, rad, x, y, z)

    utils.write_xyz_file(atoms, x, y, z, name, sys.argv[3])  # write XYZ file of new structure's coordinates
    utils.graph_protein(x, y, z, 'Final Protein')  # graph final protein


# Boiler plate invokes the main function
if __name__ == "__main__":
    main()



