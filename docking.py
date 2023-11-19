import os
import sys
import subprocess
import re
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import *
import shutil
from Bio.PDB import Structure, Model, Chain

sys.path.append('/usr/lib/python3/dist-packages')

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
# Working directory
os.chdir(script_dir)

# User input for protein
protein = input("Please enter the protein name: ")
# Define directories
dirs = [f"{script_dir}/Ligands/PDB", f"{script_dir}/Ligands/PDBQT", "Protein", "Results"]

# Ensure directories exist and are clean
for dir_name in dirs:
    # Create directory if it doesn't exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    # If it's Protein or Results directory, clean it
    elif dir_name in ["Protein", "Results"]:
        for file_name in os.listdir(dir_name):
            file_path = os.path.join(dir_name, file_name)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")

# Declare Variables
vina_path = "/home/hwcopeland/Sandbox/autodock_vina_1_1_2_linux_x86/bin/vina"
num_modes = "20"

# Create a PDBList object & Download the PDB file

pdbl = PDBList()
pdbl.retrieve_pdb_file(protein, pdir='./Protein/', file_format='pdb')
os.chdir('./Protein/')
os.rename(f'pdb{protein.lower()}.ent', f'{protein}.pdb')

# Define a list of common metals and ions
metals_and_ions = ['NA', 'MG', 'K', 'CA', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'CD', 'W', 'AU', 'HG', 'CL', 'BR', 'F', 'I', 'SO4', 'PO4', 'NO3', 'CO3']

# Load the structure
parser = PDBParser()
structure = parser.get_structure(protein, f"{protein}.pdb")

# Iterate over all atoms in the structure
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() in metals_and_ions:
                # Remove the residue if it's a metal or ion
                chain.detach_child(residue.get_id())

def get_ligand_residues(structure):
    standard_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 
                   'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
                   'TYR', 'VAL']
    ligand_residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname() not in standard_aa and res.get_resname() != "HOH":
                    ligand_residues.append((chain, res))
    return ligand_residues

def write_ligands_to_files(ligand_residues, output_dir, metals_and_ions_dir, metals_and_ions):
    for chain, ligand in ligand_residues:
        io = PDBIO()
        s = Structure.Structure('Ligand')
        m = Model.Model(0)
        s.add(m)
        c = Chain.Chain(chain.id)
        m.add(c)
        c.add(ligand.copy())
        io.set_structure(s)
        # Include the residue name in the filename
        filename = f'ligand_{ligand.get_resname()}.pdb'
        if ligand.get_resname() in metals_and_ions:
            io.save(os.path.join(metals_and_ions_dir, filename))
        else:
            io.save(os.path.join(output_dir, filename))


def clean_protein_structure(structure, ligand_residues, water_residues):
    for chain, ligand in ligand_residues:  # Remove ligand residues from the protein structure
        chain.detach_child(ligand.id)
    for chain, res_id in water_residues:  # Remove water molecules from the protein structure
        chain.detach_child(res_id)
    io = PDBIO()  # Save the cleaned protein structure
    io.set_structure(structure)
    io.save(f"{protein}_clean.pdb")

def convert_pdb_to_pdbqt(protein):
    # Define the command
    command = f'obabel {protein}_clean.pdb -O {protein}_clean.pdbqt'
    # Run the command
    subprocess.run(command, shell=True)

# Define the variables
water_residues = []

# Get ligand residues
ligand_residues = get_ligand_residues(structure)

# Write each ligand to a separate PDB file
write_ligands_to_files(ligand_residues, f'{script_dir}/Ligands/PDB', f'{script_dir}/Ligands/Ions', metals_and_ions)

def pdb_to_smiles(script_dir):
    ligand_dir = os.path.join(script_dir, 'Ligands', 'PDB')
    ligand_files = [os.path.join(ligand_dir, f) for f in os.listdir(ligand_dir) if f.endswith('.pdb')]

    smiles_list = []

    for ligand_file in ligand_files:
        # Define the command
        command = f'obabel -ipdb "{ligand_file}" -osmi'
        # Run the command
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        # Extract the SMILES string from the output
        output_parts = result.stdout.strip().split()
        if output_parts:  # Check that the output contains at least one part
            smiles = output_parts[0]  # SMILES string is the first part of the output
            smiles_list.append(smiles)

    return smiles_list

smiles_list = pdb_to_smiles(script_dir)
print(smiles_list)

## Generate Pharmacophore Model with RDKit
def generate_pharmacophore_models(smiles_list, output_dir):
    pharmacophore_models = []
    for i, smiles in enumerate(smiles_list):
        if smiles:  # Check that the SMILES string is not empty
            m = Chem.MolFromSmiles(smiles)
            if m is not None:  # Check that the molecule is valid
                AllChem.EmbedMolecule(m)
                pharmacophore_models.append(m)
                # Save the molecule as an SDF file
                Chem.rdmolfiles.MolToMolFile(m, os.path.join(output_dir, f'molecule_{i}.sdf'))
    return pharmacophore_models

pharmacophore_models = generate_pharmacophore_models(smiles_list, f'{script_dir}/Ligands/Pharmacophores/')

print(pharmacophore_models)
## Scan the molport database for similar compounds

## Generate 3D conformations of ligands with RDKit

## save top 100 results in /ligands/PDB

# Clean the protein structure
clean_protein_structure(structure, ligand_residues, water_residues)

# Convert the cleaned protein structure to PDBQT format
convert_pdb_to_pdbqt(protein)

def make_rigid(protein):
    with open(f"{protein}_clean.pdbqt", "r") as input_file, open(f"{protein}_rigid.pdbqt", "w") as output_file:
        for line in input_file:
            if not line.startswith(("ROOT", "BRANCH", "ENDBRANCH","ENDROOT", "TORSDOF", "REMARK", "MODEL", "ENDMDL", "HETATM")):
                output_file.write(line)

# Call the function with the protein name
convert_pdb_to_pdbqt(protein)
# Call the function with the protein name
make_rigid(protein)

os.chdir('../Ligands/')

def convert_ligand_to_pdbqt(ligand_path):
    # Define the command
    command = f'obabel {ligand_path} -O {output_path} -h'

    # Run the command and capture the output
    result = subprocess.run(command, shell=True)
    print(result)

# Define the variables
ligand_dir = f'{script_dir}/Ligands/PDB/'
ligand_files = [os.path.join(ligand_dir, file) for file in os.listdir(ligand_dir) if file.endswith('.pdb')]

# Call the function with the ligand path
for ligand in ligand_files:
    ligand_p = f"{ligand}"
    output_path = f"{script_dir}/Ligands/PDBQT/{os.path.basename(ligand)[:-4]}.pdbqt"
    convert_ligand_to_pdbqt(ligand_p)

os.chdir('../Protein/')

def run_fpocket(protein):
    # Run fpocket
    command = f"fpocket -f {protein}.pdb -o {protein}_out"
    subprocess.run(command, shell=True)

# Call the function with the protein name

def get_pocket_center(protein, pocket_number):
    run_fpocket(protein)
    time.sleep(10)

    file_path = f"{script_dir}/Protein/{protein}_out/{protein}_pockets.pqr"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist. Check if fpocket ran correctly and the protein file is correctly formatted.")
    with open(file_path, "r") as file:
        lines = file.readlines()
    # Only consider lines that belong to the specified pocket
    lines = [line for line in lines if line.startswith("ATOM") and int(line[22:26].strip()) == int(pocket_number)]
    coords = [list(map(float, re.findall(r"[\d\.-]+", line[30:54]))) for line in lines]
    for coord in coords:
        if len(coord) < 3:
            print(f"Unexpected line format: {coord}")
            continue
    x_coords = [coord[0] for coord in coords if len(coord) >= 3]
    y_coords = [coord[1] for coord in coords if len(coord) >= 3]
    z_coords = [coord[2] for coord in coords if len(coord) >= 3]

    if not x_coords or not y_coords or not z_coords:
        print("Error: No coordinates were extracted for the specified pocket.")
        sys.exit(1)

    center_x = sum(x_coords) / len(x_coords) if x_coords else None
    center_y = sum(y_coords) / len(y_coords) if y_coords else None
    center_z = sum(z_coords) / len(z_coords) if z_coords else None
    return center_x, center_y, center_z
def get_ligand_center(ligand):
    file_path = f"{ligand}"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist.")
    with open(file_path, "r") as file:
        lines = file.readlines()
    lines = [line for line in lines if line.startswith("HETATM")]
    coords = [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in lines]
    for coord in coords:
        if len(coord) < 3:
            print(f"Unexpected line format: {coord}")
            continue
    x_coords = [coord[0] for coord in coords if len(coord) >= 3]
    y_coords = [coord[1] for coord in coords if len(coord) >= 3]
    z_coords = [coord[2] for coord in coords if len(coord) >= 3]

    if not x_coords or not y_coords or not z_coords:
        print("Error: No coordinates were extracted for the ligand.")
        sys.exit(1)

    center_x = sum(x_coords) / len(x_coords) if x_coords else None
    center_y = sum(y_coords) / len(y_coords) if y_coords else None
    center_z = sum(z_coords) / len(z_coords) if z_coords else None
    return center_x, center_y, center_z
  

# Request user input for the pocket
if input("Would you like to view the pockets? (y/n): ") == "y":
    # Run the pymol script
    run_fpocket(protein)
    time.sleep(10)
    subprocess.run(f"{script_dir}/Protein/{protein}_out/{protein}_PYMOL.sh", shell=True)
    pocket_number = input("Please identify the pocket you would like to target: ")
    center_x, center_y, center_z = get_pocket_center(protein, pocket_number)
else:
    if input("Would you like to select a pocket based on ligand center? (y/n): ") == "y":
        ligpocket = input("What Ligand would you like to use?: ")
        center_x, center_y, center_z = get_ligand_center(f"{script_dir}/Ligands/PDB/ligand_{ligpocket}.pdb")
    else:
        print("No pocket selection method chosen. Exiting...")
        sys.exit(1)

# Use the values of center_x, center_y, and center_z in the subsequent code

# Define the size of the search space
size_x = 20
size_y = 20
size_z = 20
# Construct the command

if input("Would you like to use a custom box size? (y/n): ") == "y":
    size_x = input("Please enter the size of the box in the x direction: ")
    size_y = input("Please enter the size of the box in the y direction: ")
    size_z = input("Please enter the size of the box in the z direction: ") 
if input("Would you like to use a custom number of modes? (y/n): ") == "y":
    num_modes = input("Please enter the number of modes: ")

sdf_file = f"{script_dir}/Ligands/Pharmacophores/query_results.sdf"

# Load the SDF file
supplier = Chem.SDMolSupplier(sdf_file)

for i, mol in enumerate(supplier):
    if mol is not None:
        output_dir = f"{script_dir}/Ligands/PDBQT"
        # Get the MolPort ID from the molecule's title
        molport_id = mol.GetProp('_Name').split()[0]

        # Write the molecule to a separate SDF file
        output_sdf = f"{output_dir}/{molport_id}.sdf"
        Chem.MolToMolFile(mol, output_sdf)

        # Convert the SDF file to PDBQT format using Open Babel
        output_pdbqt = f"{output_dir}/{molport_id}.pdbqt"
        command = f"obabel {output_sdf} -O {output_pdbqt}"
        subprocess.run(command, shell=True)
        
# Update the receptor path
receptor_path = f"{script_dir}/Protein/{protein}_rigid.pdbqt"


ligand_dir = f'{script_dir}/Ligands/PDBQT/'
ligand_files = [os.path.join(ligand_dir, file) for file in os.listdir(ligand_dir) if file.endswith('.pdbqt')]

for ligand in ligand_files:
    ligandqt_path = f"{ligand}"
    output_path = f"{script_dir}/Results/{os.path.basename(ligand)[:-6]}_output.pdbqt"
    command = f"{vina_path} --receptor {receptor_path} --ligand {ligandqt_path} --out {output_path} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --num_modes {num_modes}"

    # Print the command for debugging
    print(f"Running command: {command}")

    # Run the command
    subprocess.call(command, shell=True)


