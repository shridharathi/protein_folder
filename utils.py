"""
Shridhar Athinarayanan

This file includes the utility functions for the program (extracting info from PDB coordinates, writing XYZ files,
calculating spring constants and average distances between atoms, etc.)
"""


from Bio import PDB
import matplotlib.pyplot as plt

"""
FUNCTION: get_coordinates_from_pdb
This function takes in any PDB file and writes the coordinates of each atom into returned arrays

PARAMETERS - pdbfile: PDB file to get atom coordinates of
             struct_name: the name of the protein structure
             
RETURN - atoms - list of atom names
         x_arr - list of atom x-coordinates 
         y_arr - list of atom y-coordinates 
         z_arr - list of atom z-coordinates 
"""
def get_coordinates_from_pdb(pdbfile, struct_name):
  parser = PDB.PDBParser()
  io = PDB.PDBIO()
  struct = parser.get_structure(struct_name, pdbfile)
  #f = open(textfile, "r+")

  atoms = []
  x_arr = []
  y_arr = []
  z_arr = []

  for model in struct:
      for chain in model:
          for residue in chain:
              for atom in residue:
                  atoms.append(atom.get_id())
                  x,y,z = atom.get_coord()
                  x_arr.append(x)
                  y_arr.append(y)
                  z_arr.append(z)
                  #f.write(line + '\n')
  return atoms, x_arr, y_arr, z_arr



"""
FUNCTION: filter_by_name
This function will filter out the lists of atoms and coordinates such that only atoms with a certain ID's
coordinates are returned as lists.

PARAMETERS - atoms: list of atom IDs
             x: list of x-coordinates of the protein
             y: list of y-coordinates of the protein
             z: list of z-coordinates of the protein
             name: ID of desired atom to filter by
             
RETURN - new_atoms: list with filtered atom names
         new_x: list with filtered atom x-coordinates 
         new_y: list with filtered atom y-coordinates 
         new_z: list with filtered atom z-coordinates 
"""
def filter_by_name(atoms, x, y, z, name):
    new_atoms = []
    new_x = []
    new_y = []
    new_z = []
    for i in range(len(x)):
        if atoms[i] == name:
            new_atoms.append(atoms[i])
            new_x.append(x[i])
            new_y.append(y[i])
            new_z.append(z[i])
    return new_atoms, new_x, new_y, new_z


"""
FUNCTION: write_xyz_file
This function will write the coordinates of each atom in the protein to an XYZ file format.

PARAMETERS - textfile: path to the textfile in which to save coordinates

RETURN - lists of x, y, and z-coordinates
"""
def write_xyz_file(atoms, x, y, z, name, file):
    f = open(file, "r+")
    f.write(str(len(x)) + '\n')
    #f.write('\n')
    f.write(name + '\n')

    for i in range(len(x)):
        line = str(atoms[i]) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i])
        f.write(line + '\n')
    return x, y, z


"""
FUNCTION: graph_protein
This function will display the graph of any protein given its coordinates.

PARAMETERS - x: list of x-coordinates of the protein
             y: list of y-coordinates of the protein
             z: list of z-coordinates of the protein
"""
def graph_protein(x, y, z, title):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='r', marker='o')

    plt.title(title)
    plt.show()



