import os
import sys
import subprocess
import tempfile
import re
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import *
import shutil

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
# Working directory
os.chdir(script_dir)

# Declare Variables
protein = '3B67'
nativeligand = "B67"
vina_path = "/home/hwcopeland/Sandbox/autodock_vina_1_1_2_linux_x86/bin/vina"
num_modes = 10


##############################################################################
######################### Dont Edit Below this line ##########################
##############################################################################

# Create a PDBList object & Download the PDB file

pdbl = PDBList()
pdbl.retrieve_pdb_file(protein, pdir='./Protein/', file_format='pdb')
os.chdir('./Protein/')
os.rename(f'pdb{protein.lower()}.ent', f'{protein}.pdb')

# Parse the structure
parser = PDBParser()
structure = parser.get_structure(protein, f"{protein}.pdb")

water_residues = []
ligand_residues = []
for model in structure:
    for chain in model:
        for res in chain:
            if res.get_resname() == "HOH":
                water_residues.append((chain, res.get_id()))
            elif res.get_resname() == nativeligand:
                ligand_residues.append((chain, res.get_id()))

# Create a new structure for ligand
ligand_structure = Structure.Structure('Ligand')
# Add a model to the ligand structure
ligand_model = Model.Model(0)
ligand_structure.add(ligand_model)
# Create a new chain
ligand_chain = Chain.Chain('A')
ligand_model.add(ligand_chain)
ligand_files = os.listdir(f"../Ligands/PDB/")

# Cleaning the protein structure
with open(f"{script_dir}/Ligands/PDBQT/all_ligands.pdbqt", "w") as out_file:
    # Write each ligand to the output file
    for i, ligand_file in enumerate(ligand_files):
        ligand_structure = parser.get_structure(ligand_file, f"{script_dir}/Ligands/PDB/{ligand_file}") # Parse the ligand structure
        out_file.write(f"MODEL {i+1}\n") # Write the MODEL line
        io = PDBIO() # Write the ligand structure
        io.set_structure(ligand_structure)
        io.save(out_file)
        out_file.write("ENDMDL\n") # Write the ENDMDL line
for chain, res_id in ligand_residues: # Add ligand residues to the ligand structure
    res = chain[res_id]
    ligand_chain.add(res.copy())
for chain, res_id in ligand_residues: # Remove ligand residues from the protein structure
    chain.detach_child(res_id)
for chain, res_id in water_residues: # Remove water molecules from the protein structure
    chain.detach_child(res_id)
io = PDBIO() # Save the cleaned protein structure
io.set_structure(structure)
io.save(f"{protein}_clean.pdb")
io.set_structure(ligand_structure) # Save the ligand structure with a different filename
ligand_filename = os.path.join(script_dir, 'Ligands','PDB', 'ligand_native.pdb')
io.save(ligand_filename)

def convert_pdb_to_pdbqt(protein):
    # Define the command
    command = f'obabel {protein}_clean.pdb -O {protein}_clean.pdbqt'

    # Run the command
    subprocess.run(command, shell=True)

def make_rigid(protein):
    with open(f"{protein}_clean.pdbqt", "r") as input_file, open(f"{protein}_rigid.pdbqt", "w") as output_file:
        for line in input_file:
            if not line.startswith(("ROOT", "BRANCH", "ENDBRANCH","ENDROOT", "TORSDOF")):
                output_file.write(line)

# Call the function with the protein name
convert_pdb_to_pdbqt(protein)
# Call the function with the protein name
make_rigid(protein)

os.chdir('../Ligands/')

def convert_ligand_to_pdbqt(ligand_path):
    # Define the command
    command = f'obabel {ligand_path} -O {script_dir}/Ligands/PDBQT/{ligand}.pdbqt -h'

    # Run the command and capture the output
    result = subprocess.run(command, shell=True)
    print(result)

# Call the function with the ligand path
for ligand in ligand_files:
    ligand_path = f"{script_dir}/Ligands/PDB/{ligand}"
    convert_ligand_to_pdbqt(ligand_path)

# Update the receptor path
receptor_path = f"{script_dir}/Protein/{protein}_rigid.pdbqt"
ligand_path = f"{script_dir}/Ligands/PDBQT/"  # Use the correct ligand structure filename
output_path = f"{script_dir}/Results/output.pdbqt"

os.chdir('../Protein/')

def run_fpocket(protein):
    # Run fpocket
    command = f"fpocket -f {protein}.pdb -o {protein}_out"
    subprocess.run(command, shell=True)


def get_pocket_center(protein):
    run_fpocket(protein)
    time.sleep(5)
    file_path = f"{protein}_out/{protein}_pockets.pqr"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist. Check if fpocket ran correctly and the protein file is correctly formatted.")
    with open(file_path, "r") as file:
        lines = file.readlines()
    coords = [re.findall(r"[\d\.-]+", line) for line in lines if "ATOM" in line]
    for coord in coords:
        if len(coord) < 7:
            print(f"Unexpected line format: {coord}")
            continue
    x_coords = [float(coord[2]) for coord in coords if len(coord) >= 7]
    y_coords = [float(coord[3]) for coord in coords if len(coord) >= 7]
    z_coords = [float(coord[4]) for coord in coords if len(coord) >= 7]  
    center_x = sum(x_coords) / len(x_coords) if x_coords else None
    center_y = sum(y_coords) / len(y_coords) if y_coords else None
    center_z = sum(z_coords) / len(z_coords) if z_coords else None
    return center_x, center_y, center_z

center_x, center_y, center_z = get_pocket_center(protein)

# Define the size of the search space
size_x = 20
size_y = 20
size_z = 20
# Construct the command
for ligand in ligand_files:
    ligand_path = f"{script_dir}/Ligands/PDBQT/{ligand}.pdbqt"
    output_path = f"{script_dir}/Results/{ligand[:-4]}_output.pdbqt"
    command = f"{vina_path} --receptor {receptor_path} --ligand {ligand_path} --out {output_path} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --num_modes {num_modes}"

    # Print the command for debugging
    print(f"Running command: {command}")

    # Run the command
    subprocess.call(command, shell=True)


