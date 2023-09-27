import subprocess
import numpy as np
from copy import Error

import molecular_model.atom

class Molecule_Analizer():

    def __init__(self, protein_name, aminoacids, input_filename, output_filename, basis, energy_method, protein_id=-1):

        self.protein_name = protein_name
        self.aminoacids = aminoacids
        self.protein_id = protein_id

        self.input_filename = input_filename
        self.output_filename = output_filename
        self.basis = basis
        self.energy_method = energy_method

    #Get the atoms (and the properties) of a protein
    def extract_atoms(self):

        print('    ⬤ Extracting atoms from proteins')
        #call psi4 to get the atoms of the protein
        atoms = self.get_atoms_from_protein()
        
        # if atoms length is 0 means that the proteins was not find in the database
        if len(atoms) == 0:
            raise Exception("Protein name not found. There is no atoms for that protein")

        print('    ⬤ Calculating connections between atoms')
        #Calculate the connection between atoms
        atoms, backbone = self.calculate_atom_connection(atoms)

        return atoms, backbone

    def get_atoms_from_protein(self):

        #create input file
        self.create_input_file()

        #execute psi4
        self.execute_psi_command()

        #read/parse outputfile
        [atoms, protein_id] = self.parse_psi_output_file()

        # if protein_id is not -1 means that psi4 was not able to find the protein but multiples ids for the protein
        # the solution is to create an input file with the name and the id
        if protein_id != -1:
            self.create_input_file(self.protein_name, protein_id)
            self.execute_psi_command()
            [atoms, protein_id] = self.parse_psi_output_file(self.protein_name)

        if atoms == []:
            raise Error('No atoms have been found!')

        return atoms

    def create_input_file(self):

        inputFile = open(self.input_filename+'.dat', 'w')

        inputFile.write('molecule ' + self.protein_name + '{\n')

        inputFile.write(' pubchem: '+ self.protein_name+'\n') if self.protein_id == -1 else inputFile.write(' pubchem: '+ self.protein_id+'\n')

        inputFile.write('}\n\n')

        inputFile.write('set basis ' +  self.basis + '\n')
        inputFile.write('set reference rhf\n')
        inputFile.write("energy('" + self.energy_method + "')\n")

        inputFile.close()

    def execute_psi_command(self):

        # execute psi4 by command line (it generates the file output.dat with the information)
        subprocess.run(['./psi4/psi4conda/bin/psi4','-n', str(8), self.input_filename+".dat", self.output_filename+".dat"], stdout=subprocess.DEVNULL)

    def main_chain_builder(self, nitrogen_starts, aminoacids):
        '''Takes all the nitrogens that are only connected to a single C and returns the backbone of the protein'''
        best_chains = []
        len_best_chain = 0
        for nitro in nitrogen_starts:
            candidate_chain = []
            nit = nitro

            for amino_index in range(len(aminoacids)):

                aminolist = []
                aminolist.append(nit)

                # Searching for C-alpha
                carbons = nit.linked_to_dict['C']
                carbons_not_in_chain = [carbon for carbon in carbons if (carbon not in candidate_chain and carbon not in aminolist)]
                if (len(carbons_not_in_chain)==1 and aminoacids[amino_index] != 'P'):
                    car_alpha = carbons_not_in_chain[0]
                    aminolist.append(car_alpha)
                elif (len(carbons_not_in_chain)==2 and aminoacids[amino_index] == 'P'):
                    car_alpha = (carbons_not_in_chain[0] if (len(carbons_not_in_chain[0].linked_to_dict['N']) == 1 and len(carbons_not_in_chain[0].linked_to_dict['C']) == 2 and len(carbons_not_in_chain[0].linked_to_dict['H']) == 1) else carbons_not_in_chain[1])
                    aminolist.append(car_alpha)
                else:
                    break

                # Searching for Carboxy
                carbons = car_alpha.linked_to_dict['C']
                carboxys_not_in_chain = [carbon for carbon in carbons if (carbon not in candidate_chain and carbon not in aminolist and len(carbon.linked_to_dict['O']) > 0)]
                if amino_index+1 < len(aminoacids):
                    carboxys_not_in_chain = [carbox for carbox in carboxys_not_in_chain if len(carbox.linked_to_dict['N']) > 0]
                if len(carboxys_not_in_chain)==1:
                    carbox = carboxys_not_in_chain[0]
                    aminolist.append(carbox)
                else:
                    break

                #We have a full aminoacid, so we save it to the candidate list
                candidate_chain += aminolist

                # Searching for next aminoacid Nitrogen
                nitrogens = carbox.linked_to_dict['N']
                nitrogens_not_in_chain = [n for n in nitrogens if (n not in candidate_chain and n not in aminolist)]
                if len(nitrogens_not_in_chain)==1:
                    nit = nitrogens_not_in_chain[0]
                else:
                    break

            # Is the found chain longer than the one we already had?
            if len(candidate_chain) > len_best_chain:
                len_best_chain = len(candidate_chain)
                best_chains = [candidate_chain]
            elif len(candidate_chain) == len_best_chain:
                best_chains.append(candidate_chain)
            else: 
                pass

        if len(best_chains) != 1 or len(best_chains[0])//3 != len(aminoacids):
            raise ValueError('There should be a single lengthy chain!', best_chains, nitrogen_starts)
        else:
            return best_chains[0]

    def parse_psi_output_file(self):

        atomId = 0
        protein_id = -1
        with open(self.output_filename+'.dat', 'r') as filehandle:

            isDataLine = False
            isInfoLine = False
            atoms = []
            for line in filehandle:

                #if line is an empty string after reading data
                if isDataLine and line.isspace():
                    break
                
                # Data has ------ and it is necessary to avoid it
                if isDataLine and not '--' in line:

                    lineChunks = line.split()
                    atoms += [molecular_model.atom.Atom(atomId, lineChunks[0], float(lineChunks[1]), float(lineChunks[2]), float(lineChunks[3]), float(lineChunks[4]))]
                    atomId += 1

                if isInfoLine and not 'Chemical ID' in line:
                    protein_id = line.split()[0]
                    break

                if 'Center' in line:
                    isDataLine = True

                if 'Chemical ID' in line:
                    isInfoLine = True
                        
        return [atoms, protein_id]
    
    def distance(self, atom, atom2):
        return np.sqrt((atom.x-atom2.x)**2+(atom.y-atom2.y)**2+(atom.z-atom2.z)**2)
    
    def calculate_atom_connection(self, atoms):

        #Let us first map the topology. Currently cost is O(N^2). Some other algorithm could be desirable
        for at1 in atoms:
            for at2 in atoms:
                if at1 != at2:
                    if at1.element != 'H' and at2.element != 'H' and self.distance(at1,at2)<2  and (at1 not in at2.linked_to): 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to

                    elif at1.element != 'H' and at2.element == 'H' and self.distance(at1,at2)<1.3  and (at1 not in at2.linked_to):
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to

        # Next we give an structure to each linked_to list
        for at in atoms:
            at.linked_to_dict = {'N': [], 'O': [], 'C': [], 'H': [], 'Other': []}
            for at1 in at.linked_to:
                if at1.element == 'N':
                    at.linked_to_dict['N'].append(at1)
                elif at1.element == 'O':
                    at.linked_to_dict['O'].append(at1)
                elif at1.element == 'C':
                    at.linked_to_dict['C'].append(at1)
                elif at1.element == 'H':
                    at.linked_to_dict['H'].append(at1)
                else:
                    at.linked_to_dict['Other'].append(at1)

        #self.plotting(list_of_atoms = atoms, title = 'Peptide_plot')

        #make a list of nitrogen atoms where one could start the main chain
        nitrogen_starts = []

        # For any aminoacid except proline
        if self.aminoacids[0] != 'P':
            for at in atoms:
                # This allows to identify any initial N to start except for Proline which has a weird structure
                if at.element == 'N' and len(at.linked_to_dict['C']) == 1 and len(at.linked_to_dict['H'])==2:
                    nitrogen_starts.append(at)

        # For the protein starting at proline
        elif self.aminoacids[0] == 'P':
            for at in atoms:
                # This allows to identify any initial N to start except for Proline which has a weird structure
                if at.element == 'N' and self.is_proline_N(at):
                    nitrogen_starts.append(at)

        # Find main_chain
        backbone = self.main_chain_builder(nitrogen_starts, self.aminoacids)

        # Name the atoms
        for (atom,i) in zip(backbone, range(len(backbone))):
            if atom.element == 'N' and (i % 3 == 0):
                atom.c_type = 'N_backbone'
            elif atom.element == 'C' and (i % 3 == 1) and (atom.linked_to_dict['O'] == []):
                atom.c_type = 'C_alpha'
            elif atom.element == 'C' and (i % 3 == 2) and (atom.linked_to_dict['O'] != []):
                atom.c_type = 'Carboxy'
            else:
                raise TypeError('The atom', atom.element, 'does not fulfill the requirements to be part of the backbone')

        return atoms, backbone