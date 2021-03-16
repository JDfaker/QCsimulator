import numpy as np

from shor.second_register import ShorsSecondRegister
from shor.quantum_fourier_transform import QFT


class Shors:
    '''

    '''
    def __init__(self, a, N):
        '''

        Args:
            a:
            N:
        '''
        self.a = int(a)
        self.N = int(N)
        self.t = int(np.ceil(np.log2(N)))  # minimum_number of classical bits
        self.state_no = int(int(2 ** (np.ceil(np.log2(N)))))  # minimum number of quantum states
        self.second_register = ShorsSecondRegister(int(a), int(N))
        self.qft = QFT(self.t)

    def compute_r(self):
        '''

        Returns:

        '''
        indexes = self.second_register.work_out_original_indexes()
        probs = np.abs((np.real(self.qft.generate_superposition()))) ** 2
        states = np.empty([1, len(indexes[0]), 2])
        for i in range(len(indexes[0])):
            states[0, i, 0] = probs[indexes[0, i]]
            states[0, i, 1] = indexes[0, i]
        return states[0, :, 1][np.random.choice(len(states[0, :, 0]), p=states[0, :, 0] / sum(states[0, :, 0]))].astype(
            int)

    def cf(self, n):
        '''

        Args:
            n:

        Returns:

        '''
        d = self.state_no
        res = []
        q, r = divmod(n, d)
        while r != 0:
            res = res + [q]
            prev_r = r
            q, r = divmod(d, r)
            d = prev_r
        return sum(res + [q])

