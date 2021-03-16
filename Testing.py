import numpy as np
import qiskit
from qiskit import Aer

import quantum_circuit


class Testing:
    def __init__(self, part):
        self.part = str(part)

    def test(self):
        return exec('self.' + self.part + '()')

    def tensor_product(self):
        n = np.random.randint(1, 10)
        W = []
        V = []
        for i in range(n):
            V_row = []
            W_row = []
            for j in range(n):
                V_row.append([np.random.randint(0, 50)])
                W_row.append([np.random.randint(0, 50)])
            V.append(V_row)
            W.append(W_row)
        tensor_product_test = qiskit.aqua.utils.tensorproduct(np.array(V), np.array(W))
        tensor_product_func = quantum_circuit.QuantumCircuit(2).tensor_product(np.array(V), np.array(W))
        all_zeros = not np.any(
            tensor_product_func.reshape(n ** 2, n ** 2) - tensor_product_test.reshape(n ** 2, n ** 2))
        if (all_zeros == True):
            print('Success')
        else:
            print('Fail')

    def get_initial_state(self):
        testing = []
        for qubit_no in range(1, self.test_no):
            circ = qiskit.QuantumCircuit(qubit_no)
            job = qiskit.execute(circ, Aer.get_backend('statevector_simulator'))
            initial_state_correct = np.abs(job.result().get_statevector(circ, decimals=3))
            initial_state_to_test = quantum_circuit.QuantumCircuit(qubit_no).get_initial_state()
            all_zeros = not np.any(initial_state_correct - np.reshape(initial_state_to_test, (1, 2 ** qubit_no)))
            if all_zeros == True:
                testing.append(0)
            else:
                testing.append(1)
        if sum(testing) == 0:
            print('Success')
        else:
            print('Fail')

    def apply_hadamard(self):
        testing = []
        for qubit_no in range(1, self.test_no):
            for wire_index in range(0, qubit_no):
                circ = qiskit.QuantumCircuit(qubit_no)
                circ.h(wire_index)
                job = qiskit.execute(circ, Aer.get_backend('statevector_simulator'))
                hadamard_correct = np.abs(job.result().get_statevector(circ, decimals=3))
                hadamard_to_test = quantum_circuit.QuantumCircuit(qubit_no)
                hadamard_to_test.apply_hardmard(wire_index)
                print(qubit_no, wire_index, hadamard_correct, np.reshape(hadamard_to_test.state, (1, 2 ** qubit_no)))
                all_zeros = not np.any(
                    np.round(hadamard_correct - np.reshape(hadamard_to_test.state, (1, 2 ** qubit_no))))
                if all_zeros == True:
                    testing.append(0)
                else:
                    testing.append(1)
        if sum(testing) == 0:
            print('Success')
        else:
            print('Fail')
