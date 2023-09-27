import copy
import json
import math
import subprocess
import numpy as np
import tqdm as tqdm

class Energy_Calculator():

    def __init__(self, input_filename, output_filename, precalculated_energies_path, basis, energy_method):

        self.input_filename = input_filename
        self.output_filename = output_filename
        self.precalculated_energies_path = precalculated_energies_path
        self.basis = basis
        self.energy_method = energy_method

    def calculate_energy_of_rotations(self, copied_atoms):

        #Write the file with the actual rotations
        self.write_file_energies(copied_atoms)

        #Calculate the energy of the actual rotations using PSI4
        self.execute_psi_command()

        #Read the PSI4 output file and get the energy
        energy = self.read_energy_from_psi4_file()

        return energy
    
    def execute_psi_command(self):

        # execute psi4 by command line (it generates the file output.dat with the information)
        subprocess.run(['./psi4/psi4conda/bin/psi4','-n', str(8), self.input_filename+".dat", self.output_filename+".dat"], stdout=subprocess.DEVNULL)
    
    def write_file_energies(self, atoms):

        #Write file with all atoms rotated
        rotationHandle = open(self.input_filename+'.dat', 'w')

        rotationHandle.write('molecule glycylglycine{\n')
        #write input.dat with all rotated atoms
        for at in atoms:
            rotationHandle.write(" " + at.element + " " + str(at.x) + " " + str(at.y) + " " + str(at.z)+'\n')
        
        rotationHandle.write('}\n\n')
        rotationHandle.write('set basis ' +  self.basis + '\n')
        rotationHandle.write("set reference rhf\n")
        rotationHandle.write("energy('" + self.energy_method + "')\n")

        rotationHandle.close()

    def read_energy_from_psi4_file(self):

        energy = 0
        with open(self.output_filename+'.dat', 'r') as filehandle:
            for line in filehandle:

                #If the PSI4 algorithm converges
                if 'Final Energy' in line:
                    energy = float(line.split(':')[1])
                
                #If the PSI4 algorithm does not converge, the energy used is the calculated in the last iteration (iteration 100)
                if 'iter 100:' in line:
                    energy = float(line.split()[3])

        return energy
    
    #This method returns the json with all rotations and energies associated to these rotations
    def calculate_all_deltas_of_rotations(self, atoms, aminoacids, min_energy_psi4, proteinName, numberBitsRotation, method_rotations_generation, backbone):

        rotationSteps = pow(2, int(numberBitsRotation))
        
        # it calculates the number of necessary bits to represent the number of angles
        # example, 4 aminoacids: 3 phis/psis => 2 bits
        bits_number_angles = math.ceil(np.log2(len(aminoacids)-1))

        print('    ⬤ Calculating energies for all posible rotations')
        number_angles = 2*(len(aminoacids)-1)
        energies = self.calculate_all_energies(atoms, rotationSteps, number_angles, number_angles, aminoacids)

        #Write the headers of the energies json that is going to be returned
        deltas_json = {}
        deltas_json['protein'] = proteinName
        deltas_json['numberBitsRotation'] = numberBitsRotation
        deltas_json['psi4_min_energy'] = min_energy_psi4
        deltas_json['deltas'] = {}

        print('    ⬤ Calculating deltas for all possible combinations of rotations')

        min_energy = 99999
        index_min_energy = -1  
        # iterates over all calculated energies using the keys (contains the values of the phi/psi angles)      
        for e_key in energies.keys():

            old_energy = energies[e_key]

            # check if the energy is lower than the minimum
            if old_energy < min_energy:
                min_energy = old_energy
                index_min_energy = e_key

            angle_keys = e_key.split(' ')
            # iterate over all angles keys
            for index_a_key in range(len(angle_keys)):

                # calculate the plus/minus 1 rotation delta
                for plusminus in [0,1]:
                    pm = (-2)*plusminus + 1

                    new_value = (int(angle_keys[index_a_key]) + pm) % (2**numberBitsRotation)

                    # create a key with the values of the angles (if the index is equal to the index of the modified angle, insert the modified one)
                    angle_key = ''
                    binary_key = ''
                    for index_key in range(len(angle_keys)):

                        binary_key += np.binary_repr(int(angle_keys[index_key]), width = numberBitsRotation)

                        if index_key == index_a_key:
                            angle_key += str(new_value)+ ' '
                        
                        else:
                            angle_key += angle_keys[index_key] + ' '

                    new_energy = energies[angle_key.strip()]

                    # add 0/1 for phi/psi
                    if index_a_key % 2 == 0:
                        # if it is even number is phi (add 0)
                        binary_key += str(0)
                    else:
                        # if it is odd number is psi (add 1)
                        binary_key += str(1)

                    # add the index of the phi/psi angle (with more than 2 aminoacids, there are more than one phi/psi)
                    binary_key += np.binary_repr(int(index_a_key/2), width = bits_number_angles)

                    # add 0/1 for plus/minus
                    binary_key += str(plusminus)

                    delta = new_energy - old_energy

                    #Add the values to the file with the precalculated energies
                    deltas_json['deltas'][binary_key] = delta

        deltas_json['initial_min_energy'] = min_energy
        deltas_json['index_min_energy'] = index_min_energy.replace(' ', '-')

        return deltas_json
    
    # RECURSIVE function to calculate all energies of each possible rotation 
    def calculate_all_energies(self, atoms, rotation_steps, protein_sequence_length, max_lenght, aminoacids, index_sequence='', pbar=None, energies={}):

        if pbar is None: pbar = tqdm.tqdm(total=rotation_steps**protein_sequence_length, desc=f"Calculating energy of each position")

        # iterate to calculate all possible rotations
        # for example if there are 4 rotation steps, it executes the loop 4 times, but in each iteration, it calls recursively to all rotations starting with 0 (first iteration)  
        for index in range(rotation_steps):

            if protein_sequence_length > 0:
                # returned energy is added to a data structure (this structure is multi-dimensional)
                # index_sequence contains the accumulated index (it helps to know the general index_sequence)
                energies = self.calculate_all_energies(atoms, rotation_steps, protein_sequence_length-1, max_lenght, aminoacids, index_sequence+str(index)+' ', pbar, energies)
            
            else:
                
                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                for at in copied_atoms:
                    if at.c_type == 'N_backbone' and ((len(at.linked_to_dict['C']) == 1 and len(at.linked_to_dict['H']) == 2) or self.is_proline_N(at)):
                        nitro_start = at
                        break

                copied_backbone = self.main_chain_builder([nitro_start], aminoacids)

                x_values = []
                y_values = []

                # remove last whitespace
                index_sequence = index_sequence.strip()
                index_values = index_sequence.split(' ')
                for index in range(len(index_values)):

                    if index%2 == 0:
                        # rotation sequence even (0, 2, 4, ...)
                        x_values.append(int(index_values[index])) 

                    if index%2 != 0:
                        # rotation sequence odd (1, 3, 5, ...)
                        y_values.append(int(index_values[index]))

                for index in range(len(x_values)):

                    #Always rotate from state (0,0) angle_type, angle, starting_atom, backbone
                    self.rotate(angle_type='psi', angle=(y_values[index]/rotation_steps) * 2*math.pi, starting_atom = copied_backbone[3*index + 2], backbone = copied_backbone)
                    self.rotate(angle_type='phi', angle=(x_values[index]/rotation_steps) * 2*math.pi, starting_atom = copied_backbone[3*index + 4], backbone = copied_backbone) 
                    
                
                #Calculate the energy of the protein structure after the previous rotations
                energies[index_sequence] = self.calculate_energy_of_rotations(copied_atoms)
                pbar.update(1)

                # We eliminate previous copies
                del copied_atoms
                del copied_backbone

                break

        return energies

    def is_proline_N(self, atom):

        carbon_ring = []

        if atom.element != 'N' or len(atom.linked_to_dict['C']) != 2 or len(atom.linked_to_dict['H'])!=1:
            return False
        else:
            carbons = atom.linked_to_dict['C']
            if len(carbons[0].linked_to_dict['C']) == 1  and len(carbons[0].linked_to_dict['N']) == 1 and len(carbons[1].linked_to_dict['C']) == 2 and len(carbons[1].linked_to_dict['N']) == 1:
                current_carbon = carbons[0]
                ending_carbon = carbons[1]
            elif len(carbons[1].linked_to_dict['C']) == 1  and len(carbons[1].linked_to_dict['N']) == 1 and len(carbons[0].linked_to_dict['C']) == 2 and len(carbons[0].linked_to_dict['N']) == 1:
                current_carbon = carbons[1]
                ending_carbon = carbons[0]
            else:
                return False
            
            for _ in range(2):
                carbon_ring.append(current_carbon)
                current_carbon = (current_carbon.linked_to_dict['C'][0] if current_carbon.linked_to_dict['C'][0] not in carbon_ring else current_carbon.linked_to_dict['C'][1])
                if len(current_carbon.linked_to_dict['C']) != 2 or len(current_carbon.linked_to_dict['N']) != 0 or len(current_carbon.linked_to_dict['O']) != 0 or len(current_carbon.linked_to_dict['H']) != 2:
                    return False
                
            return (True if current_carbon in ending_carbon.linked_to else False)

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
    
    def rotate(self, angle_type, angle, starting_atom, backbone):

        previous_atom = backbone[backbone.index(starting_atom)-1]

        if angle_type == 'phi':
            if previous_atom.c_type != 'N_backbone' or starting_atom.c_type != 'C_alpha':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type,'or wrong previous atom',previous_atom.c_type )
                    
        elif angle_type == 'psi':
            if previous_atom.c_type != 'C_alpha' or starting_atom.c_type != 'Carboxy':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )

        else:
            raise Exception('Angle not recognised!:',angle_type)
        
        # Define the list of atoms to rotate and then rotate them
        
        backbone2rotate = backbone[backbone.index(starting_atom)+1:]
        ##self.backbone_to_rotate(angle_type,starting_atom, backbone)
        list_of_atoms_to_rotate = self.decorations_to_rotate(backbone2rotate,backbone)

        for atom in list_of_atoms_to_rotate:
            # The axis is defined by the starting atom and the atom prior to the starting atom in the backbone
            atom.rotate(previous_atom, starting_atom, angle, angle_type)

    def decorations_to_rotate(self, backbone2rotate, backbone):
        
        atoms2rotate = backbone2rotate

        newly_added = backbone2rotate

        while newly_added != []:
            previously_added = newly_added
            newly_added = []
            for atom in previously_added:
                for at2 in atom.linked_to:
                    if at2 not in atoms2rotate and at2 not in newly_added and at2 not in backbone:
                        newly_added.append(at2)

            atoms2rotate += newly_added 
        
        return atoms2rotate

    def write_json(self, json_data, file_name, proteinName, numberBitsRotation, method_rotations_generation):

        #Create json with calculated energies
        #TODO: extract the path to a config file
        with open(self.precalculated_energies_path+file_name+'_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json', 'w') as outfile:
            json.dump(json_data, outfile)

    def read_energy_json(self, proteinName, numberBitsRotation, method_rotations_generation):

        with open(self.precalculated_energies_path + 'delta_energies_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json') as json_file:
            data = json.load(json_file)

            return [data['deltas'], data['psi4_min_energy'], data['initial_min_energy'], data['index_min_energy'], data['initialization_stats']]

        