from sympy.combinatorics.graycode import GrayCode, gray_to_bin, bin_to_gray
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit, Qubit
from collections import OrderedDict
from qiskit import Aer, transpile

import math
import numpy as np

class Oracle():
    '''Outputs the binary coord of rotation to get the correct probability. Tested ok'''
    def __init__(self, deltas_dictionary, n_angles, angle_precision_bits, out_bits, precision_coords, optimization=False, mct_mode='noancilla'):

        self.out_bits = out_bits
        self.deltas_dictionary = OrderedDict(sorted(deltas_dictionary.items()))
        self.n_angles = n_angles
        self.angle_precision_bits = angle_precision_bits
        self.precision_coords = precision_coords

    def generate_angles_codification(self, beta):

        angles = {}

        for key in self.deltas_dictionary.keys():

            if self.deltas_dictionary[key] >= 0:
                probability = math.exp(-beta * self.deltas_dictionary[key])
            else: 
                probability = 1
            # Instead of encoding the angle corresponding to the probability, we will encode the angle theta such that sin^2(pi/2 - theta) = probability.
            # That way 1 -> 000, but if probability is 0 there is some small probability of acceptance
            
            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(pi/2-theta/2).
            # The theta/2 is because if you input theta, qiskits rotates theta/2. Also normalised (divided between pi the result)
            angle = 1 - 2/math.pi * math.asin(math.sqrt(probability))

            # Ensure that the angle stays minimally away from 1
            angle = np.minimum(angle, 1-2**(-self.out_bits-1))
            # Convert it into an integer and a string
            if angle == 1.:
                print('probability = ',probability)
                print('angle',angle)
                raise ValueError('Warning: angle seems to be pi/2, and that should not be possible')
            
            # angle will be between 0 and 1, so we move it to between 0 and 2^out_bits. Then calculate the integer and the binary representation
            angles[key] = np.binary_repr(int(angle*2**self.out_bits), width= self.out_bits)

            #if key[:10] == '1101000101':
            #    print('<DEBUG> For key:', key, 'angle is:', angles[key])

        self.angles = angles

        return angles

    def generate_oracle(self, beta):

        angles = self.generate_angles_codification(beta)

        # create a quantum circuit with the same length than the key of the deltas energies

        oracle_key = QuantumRegister(len(list(self.deltas_dictionary.keys())[0]))
        oracle_value = QuantumRegister(len(list(angles.values())[0]))
        oracle_circuit = QuantumCircuit(oracle_key, oracle_value)

        len_key = len(list(self.angles.keys())[0])

        for key in self.angles.keys():
            
            # apply x gates to the 0s in the key
            for key_bit_index in range(len(key)):
                if key[(len_key-1) - key_bit_index] == '0':
                    oracle_circuit.x(oracle_key[key_bit_index])

            # apply mcx gates with the 1s in the angles (control the whole key)
            angle = angles[key]
            for angle_bit_index in range(len(angle)):
                if angle[(len(angle)-1) - angle_bit_index] == '1':
                    oracle_circuit.mcx(oracle_key, oracle_value[angle_bit_index])

            # apply x gates to the 0s in the key
            for key_bit_index in range(len(key)):
                if key[(len_key-1) - key_bit_index] == '0':
                    oracle_circuit.x(oracle_key[key_bit_index])

        return oracle_circuit.to_instruction()

    # convert the interger key of the deltas dict in binary key
    # format of deltas key coord1-coord2-coord3-...|id(integer)|plusminus
    def key_to_binary(self, key):

        binary_key = ''

        coords = key.split('|')[0].split('-')
        id_integer = int(key.split('|')[1])
        plus_minus = key.split('|')[2]

        for coord_index in range(len(coords)):
            binary_key += np.binary_repr(int(coords[coord_index]), width = int(self.precision_coords[coord_index]))

        if len(self.precision_coords) == 1:
            binary_key += '0'
        else:
            number_bits_to_represent_id = int(math.ceil(np.log2(len(self.precision_coords))))
            binary_key += np.binary_repr(id_integer, width=number_bits_to_represent_id)

        binary_key += plus_minus

        return binary_key

    def execute_x_multiple_register(self, circuit, target, target_value, index, direction):

        if direction == 'fordward':

            for i in range(len(index)):
                if index[i] == target_value:
                    circuit.x(target[i])

        elif direction == 'backward':

            for i in range(len(index)):
                if index[len(index)-i-1] == target_value:
                    circuit.x(target[i])
        else:
            raise Exception("Wrong direction value")