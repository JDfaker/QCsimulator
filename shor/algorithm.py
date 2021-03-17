"""
This modeule holds the Shors class
"""
import numpy as np
from second_register import ShorsSecondRegister
from quantum_fourier_transform import QFT


class Shors:
    """
    This is a class to impliment Shor"s algorithm \n
    It is linked with the classes ShorsSecondRegister
    and QFT (quantum fourier transform) \n
    """
    def __init__(self, a, N):
        """
        This is the constructor method for Shors \n
        @param a: A number coprime of N \n
        @param N: A number we wish to factor \n
        """
        self.a = int(a)
        self.N = int(N)
        self.t = int(np.ceil(np.log2(N)))
        # minimum_number of classical bits
        self.state_no = int(int(2 ** (np.ceil(np.log2(N)))))
        # min no quantum states
        self.second_register = ShorsSecondRegister(int(a), int(N))
        self.qft = QFT(self.t)

    def compute_r(self):
        """
        Method for computing probability amplitudes of the QFT \n
        @return: Numerical indexes of the states \n
        """
        indexes = ShorsSecondRegister.pick_out_states(self)
        probs = np.abs((np.real(self.qft.generate_superposition()))) ** 2
        states = np.empty([1, len(indexes[0]), 2])
        for i in range(len(indexes[0])):
            states[0, i, 0] = probs[indexes[0, i]]
            states[0, i, 1] = indexes[0, i]
        point = len(states[0, :, 0])
        p = states[0, :, 0] / sum(states[0, :, 0])
        return states[0, :, 1][np.random.choice(point)].astype(int)

    def cf(self, n):
        """
        Method for computing partial expansion of n/2^N \n
        @param n: Number of qubits in the system \n
        """
        d = self.state_no
        res = []
        q, r = divmod(n, d)
        while r != 0:
            res = res + [q]
            prev_r = r
            q, r = divmod(d, r)
            d = prev_r
        return sum(res + [q])

    def compute_shors(self):
        """
        Putting the states through shor's\n
        @return: Output of Shors \n
        """
        thing = Shors.compute_r(self)
        r = Shors.cf(self, thing)
        while r % 2 != 0:
            r = Shors.cf(self, Shors.compute_r(self))
            a1 = int((self.a ** (r / 2)) - 1)
            a2 = int((self.a ** (r / 2)) + 1)
        a1 = int((self.a ** (r / 2)) - 1)
        a2 = int((self.a ** (r / 2)) + 1)
        return np.gcd(a1, self.N), np.gcd(a2, self.N)
