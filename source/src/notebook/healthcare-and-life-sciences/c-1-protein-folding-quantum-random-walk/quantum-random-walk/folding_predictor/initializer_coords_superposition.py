###############################################################################################################################################
#                                                                                                                                             #
#                                                           Code from Qiskit Terra                                                            #
# (https://github.com/Qiskit/qiskit-terra/blob/b42142eeaf0a72fc83b1bacb8be8c4e36fda9b5f/qiskit/extensions/quantum_initializer/initializer.py) #
#                        Modification necessary to remove the resets in the code (an invertible circuit was necessary)                        #
#                                                                                                                                             #
###############################################################################################################################################


from qiskit import QuantumCircuit, execute
from qiskit.circuit import QuantumRegister
from qiskit.quantum_info import Statevector
from qiskit import Aer
import math
import numpy as np
from qiskit.circuit.library.standard_gates.x import CXGate
from qiskit.circuit.library.standard_gates.rz import RZGate
from qiskit.circuit.library.standard_gates.ry import RYGate

EPS = 1e-10  # global variable used to chop very small numbers to zero
def _bloch_angles(pair_of_complex):
    """
    Static internal method to work out rotation to create the passed-in
    qubit from the zero vector.
    """
    [a_complex, b_complex] = pair_of_complex
    # Force a and b to be complex, as otherwise numpy.angle might fail.
    a_complex = complex(a_complex)
    b_complex = complex(b_complex)
    mag_a = np.absolute(a_complex)
    final_r = float(np.sqrt(mag_a ** 2 + np.absolute(b_complex) ** 2))
    if final_r < EPS:
        theta = 0
        phi = 0
        final_r = 0
        final_t = 0
    else:
        theta = float(2 * np.arccos(mag_a / final_r))
        a_arg = np.angle(a_complex)
        b_arg = np.angle(b_complex)
        final_t = a_arg + b_arg
        phi = b_arg - a_arg

    return final_r * np.exp(1.0j * final_t / 2), theta, phi


def _multiplex(target_gate, list_of_angles, last_cnot=True):
    """
    Return a recursive implementation of a multiplexor circuit,
    where each instruction itself has a decomposition based on
    smaller multiplexors.

    The LSB is the multiplexor "data" and the other bits are multiplexor "select".

    Args:
        target_gate (Gate): Ry or Rz gate to apply to target qubit, multiplexed
            over all other "select" qubits
        list_of_angles (list[float]): list of rotation angles to apply Ry and Rz
        last_cnot (bool): add the last cnot if last_cnot = True

    Returns:
        DAGCircuit: the circuit implementing the multiplexor's action
    """
    list_len = len(list_of_angles)
    local_num_qubits = int(math.log2(list_len)) + 1

    q = QuantumRegister(local_num_qubits)
    circuit = QuantumCircuit(q, name="multiplex" + local_num_qubits.__str__())

    lsb = q[0]
    msb = q[local_num_qubits - 1]

    # case of no multiplexing: base case for recursion
    if local_num_qubits == 1:
        circuit.append(target_gate(list_of_angles[0]), [q[0]])
        return circuit

    # calc angle weights, assuming recursion (that is the lower-level
    # requested angles have been correctly implemented by recursion
    angle_weight = np.kron([[0.5, 0.5], [0.5, -0.5]], np.identity(2 ** (local_num_qubits - 2)))

    # calc the combo angles
    list_of_angles = angle_weight.dot(np.array(list_of_angles)).tolist()

    # recursive step on half the angles fulfilling the above assumption
    multiplex_1 = _multiplex(target_gate, list_of_angles[0 : (list_len // 2)], False)
    circuit.append(multiplex_1.to_instruction(), q[0:-1])

    # attach CNOT as follows, thereby flipping the LSB qubit
    circuit.append(CXGate(), [msb, lsb])

    # implement extra efficiency from the paper of cancelling adjacent
    # CNOTs (by leaving out last CNOT and reversing (NOT inverting) the
    # second lower-level multiplex)
    multiplex_2 = _multiplex(target_gate, list_of_angles[(list_len // 2) :], False)
    if list_len > 1:
        circuit.append(multiplex_2.to_instruction().reverse_ops(), q[0:-1])
    else:
        circuit.append(multiplex_2.to_instruction(), q[0:-1])

    # attach a final CNOT
    if last_cnot:
        circuit.append(CXGate(), [msb, lsb])

    return circuit

def _rotations_to_disentangle(local_param):
    """
    Static internal method to work out Ry and Rz rotation angles used
    to disentangle the LSB qubit.
    These rotations make up the block diagonal matrix U (i.e. multiplexor)
    that disentangles the LSB.

    [[Ry(theta_1).Rz(phi_1)  0   .   .   0],
        [0         Ry(theta_2).Rz(phi_2) .  0],
                                .
                                    .
        0         0           Ry(theta_2^n).Rz(phi_2^n)]]
    """
    remaining_vector = []
    thetas = []
    phis = []

    param_len = len(local_param)

    for i in range(param_len // 2):
        # Ry and Rz rotations to move bloch vector from 0 to "imaginary"
        # qubit
        # (imagine a qubit state signified by the amplitudes at index 2*i
        # and 2*(i+1), corresponding to the select qubits of the
        # multiplexor being in state |i>)
        (remains, add_theta, add_phi) = _bloch_angles(
            local_param[2 * i : 2 * (i + 1)]
        )

        remaining_vector.append(remains)

        # rotations for all imaginary qubits of the full vector
        # to move from where it is to zero, hence the negative sign
        thetas.append(-add_theta)
        phis.append(-add_phi)

    return remaining_vector, thetas, phis
    
def _gates_to_uncompute(amplitudes):
    """Call to create a circuit with gates that take the desired vector to zero.

    Returns:
        QuantumCircuit: circuit to take self.params vector to :math:`|{00\\ldots0}\\rangle`
    """

    num_qubits = int(math.log2(len(amplitudes)))

    q = QuantumRegister(num_qubits)
    circuit = QuantumCircuit(q, name="disentangler")

    # kick start the peeling loop, and disentangle one-by-one from LSB to MSB
    remaining_param = amplitudes

    for i in range(num_qubits):
        # work out which rotations must be done to disentangle the LSB
        # qubit (we peel away one qubit at a time)
        (remaining_param, thetas, phis) = _rotations_to_disentangle(remaining_param)

        # perform the required rotations to decouple the LSB qubit (so that
        # it can be "factored" out, leaving a shorter amplitude vector to peel away)

        add_last_cnot = True
        if np.linalg.norm(phis) != 0 and np.linalg.norm(thetas) != 0:
            add_last_cnot = False

        if np.linalg.norm(phis) != 0:
            rz_mult = _multiplex(RZGate, phis, last_cnot=add_last_cnot)
            circuit.append(rz_mult.to_instruction(), q[i : num_qubits])

        if np.linalg.norm(thetas) != 0:
            ry_mult = _multiplex(RYGate, thetas, last_cnot=add_last_cnot)
            circuit.append(ry_mult.to_instruction().reverse_ops(), q[i : num_qubits])
    circuit.global_phase -= np.angle(sum(remaining_param))
    return circuit

# this method generate a list of normalized amplitudes with the parameter number. If the number is not a power of 2, it fills with zeros until a power of 2
def _generate_amplitude_values_for_number(number, number_qubits):
    return [math.sqrt(1/number) if index_amplitudes < number else 0 for index_amplitudes in range(2**number_qubits)]

def generate_uniform_superposition(number_coordinates, number_qubits):

    # call to generate the circuit that takes the desired vector to zero
    disentangling_circuit = _gates_to_uncompute(_generate_amplitude_values_for_number(number_coordinates, number_qubits))

    # invert the circuit to create the desired vector from zero (assuming
    # the qubits are in the zero state)
    initialize_instr = disentangling_circuit.to_instruction().inverse()

    return initialize_instr