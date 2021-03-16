import numpy as np

from shor.algorithm import Shors


class ClassicalPart:
    '''

    '''

    def __init__(self, test_range):
        '''

        Args:
            test_range:
        '''
        self.range = test_range

    def main(self, N):
        '''

        Args:
            N:

        Returns:

        '''
        gcd = []
        for i in range(N):
            if np.gcd(i, N) == 1:
                gcd.append(Shors(i, N).main())
            else:
                pass
        return gcd

    def testing(self):
        '''

        Returns:

        '''
        factors = []
        for i in range(6, self.range):
            factors.append(self.main(i))
        return factors

    def continued_fractions(self, n):
        '''

        Args:
            self:
            n:

        Returns:

        '''
        remainders = []
        d = self.t
        for i in range(10):
            remainder = np.floor(n / d)
            if remainder > 0:
                diff = n - remainder * d
                if diff == 0:
                    remainders.append(remainder)
                    break
                new_frac = n / diff
                n, d = d, diff
            else:
                n, d = d, n
                if d == 0:
                    break
            remainders.append(remainder)
        return np.array(remainders)

    def compute_fractions(self, n):
        '''

        Args:
            self:
            n:

        Returns:

        '''
        remainders = self.continued_fractions(n)
        # new_value = 0
        # fractions = []
        new_denominator = 0
        new_numerator = 0
        for i in range(len(remainders[0:]) - 1, 1, -1):
            if i == (len(remainders[0:])) - 1:
                new_denominator = int(remainders[i])
                new_numerator = 1 + new_denominator
                new_denominator, new_numerator = new_numerator, new_denominator
            else:
                new_value = int(remainders[i])
                new_numerator = new_numerator + (new_value * new_denominator)
                new_denominator, new_numerator = new_numerator, new_denominator
        return new_denominator
