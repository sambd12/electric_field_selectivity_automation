import numpy as np
import argparse

def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs=1, action='store')

    args=parser.parse_args()
    return args
    #args is the entire namespace
    

def read_xyz(filename):
    """ Read in an xyz file.  Returns a list of the atom column and 
    numpy array containing the atomic coordinates
    """
    f = open(filename,'r')

    n_atoms = int(f.readline())
    f.readline()
    atoms_col = []
    coords = [[],[],[]]

    for line in f:
        words = line.split()
        atoms_col.append(words[0])
        coords[0].append(float(words[1]))
        coords[1].append(float(words[2]))
        coords[2].append(float(words[3]))
    return atoms_col, np.array(coords).T

def write_xyz(filename,atom_col,coords,write_file=True):
    """ Write an xyz file, given a list containing the atom names and a
    numpy array containing the atomic coordinates
    """
    coords_to_write = coords.T
    natoms = len(atom_col)
    file_str = ""
    file_str += str(natoms) + '\n'
    file_str += 'generated by a python script \n'
    if write_file:
        f = open(filename,'w')
        f.write(file_str)

    for i in range(natoms):
        atom = atom_col[i]
        x = str(coords_to_write[0][i])
        y = str(coords_to_write[1][i])
        z = str(coords_to_write[2][i])
        sep = '     '
        line = atom + sep + x + sep + y + sep + z + '\n'
        file_str += line
        if write_file:
            f.write(line)
    return file_str


def normalize(v):
    """Normalize a vector."""
    norm = np.linalg.norm(v)
    if norm == 0: 
        return v
    return v / norm

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation 
    about the given axis by theta radians.
    """
    axis = normalize(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def align_bond_to_z_axis(coords):
    bond_vector = coords[1] - coords[0]     # Calculate bond vector
    z_axis = np.array([0, 0, 1])            # Z-axis vector
    
    axis = np.cross(bond_vector, z_axis)    # Axis of rotation 
    angle = np.arccos(np.dot(normalize(bond_vector), z_axis))
    
    R = rotation_matrix(axis, angle)        # Rotation matrix
    
    # Apply rotation to all coordinates
    # transpose R for right-hand multiplication
    rotated_coords = np.dot(coords, R.T)  
    
    return rotated_coords

def rotation_matrix_z(theta):
    """Return the rotation matrix for a counterclockwise rotation around 
    the z-axis by theta radians."""
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta),  np.cos(theta), 0],
                     [0,              0,             1]])

def align_bond_to_x_axis(coords):
    # Calculate the bond vector between the first and third atoms
    bond_vector = coords[2] - coords[0]
    
    # Angle between the projected bond vector and the x-axis
    angle = np.arctan2(bond_vector[1], bond_vector[0])
    
    #  Since we want the bond vector to align with the x-axis,
    #  we rotate by -angle
    R = rotation_matrix_z(-angle)
    
    # Apply rotation to all coordinates
    # transpose R for right-hand multiplication
    rotated_coords = np.dot(coords, R.T)  
    
    return rotated_coords

def round_and_clean_up(coords, tol=1e-15, decimals = 6):
    result = np.zeros_like(coords)
    for i, atom_vector in enumerate(coords):
        for j, component in enumerate(atom_vector):
            if np.abs(component) <= tol:
                result[i,j] = 0. 
            else:
                result[i,j] = np.round(component, decimals)
    return result

def align_xyz(filename):
    atom_col, coords = read_xyz(filename)
    coords -= coords[0]
    rotated_coords_z = align_bond_to_z_axis(coords)
    rotated_coords = align_bond_to_x_axis(rotated_coords_z)
    rotated_coords = round_and_clean_up(rotated_coords)
    
    name=filename.split('.')[0] 
    aligned_filename = name + '_aligned.xyz'
    write_xyz(aligned_filename, atom_col, rotated_coords)
    return aligned_filename

    

def main():
    args=get_arguments()
    filename = args.filename[0]
    align_xyz(filename) 

if __name__ == "__main__":
    main()
