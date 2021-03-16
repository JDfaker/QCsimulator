import cmath
import numpy as np
from matplotlib import pyplot as plt


class QFT:
    '''

    '''
    def __init__(self, no_qubits):
        self.n = no_qubits

    def generate_states(self):
        '''

        Returns:

        '''
        state_list = []
        for i in range(self.n):
            state_list.append([1, 0])
        return np.array(state_list)

    def qft_circuit(self):
        '''

        Returns:

        '''
        states = self.generate_states()
        had = (1 / np.sqrt(2)) * np.array([[1, 1],
                                           [1, -1]])
        output = []
        for i in range(1, len(states) + 1):
            initial = np.matmul(had, states[i - 1])
            for j in range(2, len(states) - i):
                initial = np.matmul(self.r_n(j), initial)
            output.append(initial)
        return np.array(output)

    def r_n(self, n):
        '''

        Args:
            n:

        Returns:

        '''
        val = cmath.exp((2 * np.pi * 1j) / (2 ** n))
        return np.array([[1, 0], [0, val]])

    def generate_superposition(self):
        '''

        Returns:

        '''
        states = self.qft_circuit()
        for i in range(len(states)):
            initial = np.kron(states[0], states[1])
            for j in range(2, len(states)):
                initial = np.kron(states[j], initial)
        return initial

    def plot_qft(self):
        '''

        Returns:

        '''
        plt.title(' Probaility Amplitudes for a ' + str(self.n) + ' qubit register')
        plt.ylabel('$|\psi|^{2}$')
        plt.xlabel('Number')
        return plt.bar(np.arange(0, 2 ** self.n), (2 * (np.real(self.generate_superposition()) ** 2)))
