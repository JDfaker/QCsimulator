"""
This module contains the class QFT (Quantum Fourier Transform)
"""
import cmath
import numpy as np
from matplotlib import pyplot as plt
import quantum_gate as QG


class QFT:
    """
    This class impliments the quantum fourier transform \n
    within the Shor Algorith \n
    """
    def __init__(self, no_qubits):
        """
        This is the constructor method for Shors \n
        @param no_qubits: Integer number of qubits \n
        """
        self.n = no_qubits

    def generate_states(self):
        """
        Method for generating states for fourier transform \n
        @return: A list of numpy arrays of states \n
        """
        state_list = []
        for i in range(self.n):
            state_list.append([1, 0])
        return state_list

    def qft_circuit(self):
        """
        Method for implimenting the quantum fourier transform in circuit \n
        @return: A list of numpy arrays of states gone through gate \n
        """
        states = self.generate_states()
        output = []
        for i in range(1, len(states) + 1):
            initial = np.dot(QG.H_Shors, states[i - 1])
            for j in range(2, len(states) - i):
                initial = np.dot(self.r_n(j), initial)
            output.append(initial)
        return output

    def r_n(self, n):
        """
        Method for the r gate \n
        @param n: Integer of nth qubit in circuit \n
        @return: Complex numpy array of r gate \n
        """
        val = cmath.exp((2 * np.pi * 1j) / (2 ** n))
        return np.array([[1, 0], [0, val]])

    def generate_superposition(self):
        """
        Method for implimenting a superposition of qubit states \n
        @return: Numpy array of superposition \n
        """
        states = self.qft_circuit()
        initial = np.kron(states[0], states[1])
        for j in range(2, len(states)):
            initial = np.kron(states[j], initial)
        return initial

    def plot_qft(self):
        """
        Method for plotting the quantum fourier transform \n
        @return: Plot of states as bar graph \n
        """
        plt.title(' Probaility Amplitudes for a '
                  + str(self.n) + ' qubit register')
        plt.ylabel('$|\psi|^{2}$')
        plt.xlabel('Number')
        return plt.bar(np.arange(0, 2 ** self.n),
                       (2 * (np.real(self.generate_superposition()) ** 2)))
