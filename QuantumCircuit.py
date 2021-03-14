"""
This module holds the Quantum Circuit class
"""

import numpy as np
import QuantumGate as QG
from matplotlib import pyplot as plt
from sparse_matrix import SparseMatrix


class QuantumCircuit:
    """
    This is a class for implimenting grover's algorithm
    as well as various quantum gates. \n
    It hold a quantum state as a sparse matrix and the number of qubits. \n
    The sparse matrices are implimented via a seperate class SparseMatrix. \n
    """

    def __init__(self, qubit_number):
        """
        This is the constructor method for Quantum Circuit \n
        @param qubit_number: Number of qubits in state \n
        """
        assert (qubit_number > 0), 'Qubit number should be more than 0'
        self.qn = qubit_number
        self.state = self.get_initial_state()

    def show_state(self):
        '''
        Method to show current state in terminal \n
        '''
        print(SparseMatrix.numpy(self.state))

    def get_initial_state(self):
        '''
        Initialize the qubit state given by the number of qubits \n
        @return: Quantum state as a sparse matrix \n
        '''
        state = np.zeros(2 ** self.qn)
        state[0] = 1
        return SparseMatrix.sparsify(state.reshape(len(state), 1))

    def apply_hardmard(self, wire_index):
        '''
        Method applying the hardman gate to the quantum state \n
        @param wire_index: Integer indicating location of hardman gate \n
        @return: State changed via going through hardman gate \n
        '''
        assert -1 < wire_index < self.qn, (
            'Input argument should be between wire 0 to ' + str(self.qn - 1))

        if self.qn == 1:
            self.state = SparseMatrix.dot(QG.H, self.state)
        else:
            gate_list = []
            for i in range(self.qn):
                if i == wire_index:
                    gate_list.append(QG.H)
                else:
                    gate_list.append(QG.eye)

            gate_M = gate_list[0]
            for i in range(1, self.qn):
                gate_M = SparseMatrix.tensordot(gate_M, gate_list[i])
            self.state = SparseMatrix.dot(gate_M, self.state)

    def apply_pauliX(self, wire_index):
        '''
        Method applying the Pauli X gate to the quantum state \n
        @param wire_index: Integer indicating location of gate \n
        @return: State changed via going through gate \n
        '''
        assert -1 < wire_index < self.qn, (
            'Input argument should be between wire 0 to ' + str(self.qn - 1))

        if self.qn == 1:
            self.state = SparseMatrix.dot(QG.PX, self.state)
        else:
            gate_list = []
            for i in range(self.qn):
                if i == wire_index:
                    gate_list.append(QG.PX)
                else:
                    gate_list.append(QG.eye)

            gate_M = gate_list[0]
            for i in range(1, self.qn):
                gate_M = SparseMatrix.tensordot(gate_M, gate_list[i])
            self.state = SparseMatrix.dot(gate_M, self.state)

    def apply_pauliZ(self, wire_index):
        '''
        Method applying the Pauli Z gate to the quantum state \n
        @param wire_index: Integer indicating location of gate \n
        @return: State changed via going through gate \n
        '''
        assert -1 < wire_index < self.qn, (
            'Input argument should be between wire 0 to ' + str(self.qn - 1))

        if self.qn == 1:
            self.state = SparseMatrix.dot(QG.PZ, self.state)
        else:
            gate_list = []
            for i in range(self.qn):
                if i == wire_index:
                    gate_list.append(QG.PZ)
                else:
                    gate_list.append(QG.eye)
            gate_M = gate_list[0]
            for i in range(1, self.qn):
                gate_M = SparseMatrix.tensordot(gate_M, gate_list[i])
            self.state = SparseMatrix.dot(gate_M, self.state)

    def apply_swap(self, wire_index1, wire_index2):
        '''
        Method applying the swap gate to the quantum state \n
        @param wire_index1: Integer indicating location of gate \n
        @param wire_index2: Integer indicating location of gate \n
        @return: State changed via going through gate \n
        '''
        assert wire_index1 < self.qn or wire_index2 < self.qn, (
            'Input argument should be between wire 0 to ' + str(self.qn - 1))

        if self.qn == 2:
            self.state = SparseMatrix.dot(QG.SWAP, self.state)
        else:
            if wire_index1 < wire_index2:
                a = wire_index1
            else:
                a = wire_index2
            gate_list = []
            for i in range(self.qn - 1):
                if i == a:
                    gate_list.append(QG.SWAP)
                else:
                    gate_list.append(QG.eye)
            gate_M = gate_list[0]
            for i in range(1, self.qn - 1):
                gate_M = SparseMatrix.tensordot(gate_M, gate_list[i])
            self.state = SparseMatrix.dot(gate_M, self.state)

    def apply_grover_oracle(self, marks):
        '''
        Method to apply grover oracle \n
        @param marks: Integer or list of location for oracle \n
        @return: Changed qubit state with grover oracle \n
        '''
        eye = np.eye(2 ** self.qn)
        oracle = eye
        if isinstance(marks, int):
            oracle[marks][marks] = -1
        else:
            for mark in marks:
                oracle[mark][mark] = -1
        oracle = SparseMatrix.sparsify(oracle)
        self.state = SparseMatrix.dot(oracle, self.state)

    def apply_amplification(self):
        '''
        Method applying amplitude amplification of marked item \n
        @return: Changed qubit state via amplification of marked item \n
        '''
        s = np.zeros(2 ** self.qn)
        s[0] = 1
        s = SparseMatrix.sparsify(s.reshape(len(s), 1))
        gate_list = []
        for i in range(self.qn):
            gate_list.append(QG.H)

        H = gate_list[0]
        for i in range(1, self.qn):
            H = SparseMatrix.tensordot(H, gate_list[i])
        s = SparseMatrix.dot(H, s)
        s_T = SparseMatrix.transpose(s)
        s_s = SparseMatrix.tensordot(s, s_T)
        s_s_2 = SparseMatrix.multiply(s_s, 2)
        eye = SparseMatrix.sparsify(np.eye(2 ** self.qn))
        diffuser = SparseMatrix.minus(s_s_2, eye)
        self.state = SparseMatrix.dot(diffuser, self.state)

    def plot_pr(self):
        '''
        Plotting probabilities of qubit states
        @result: Bar graph of probabilities
        '''
        temp_x = range(1, 2 ** self.qn + 1)
        x = []
        for elem in temp_x:
            x.append(str(elem - 1))
        y = []

        ss = SparseMatrix.numpy(self.state)
        for i in range(ss.shape[0]):
            y.append((ss[i][0]) ** 2)
        plt.style.use('seaborn')
        plt.bar(x, y, width=0.5)
        plt.tick_params(axis='both', labelsize=15)
        if self.qn > 4:
            plt.tick_params(axis='x', labelsize=10)
        plt.ylim(0, 1)
        plt.ylabel('Probability', fontsize=15)
        plt.xlabel('State', fontsize=15)
        plt.show()
