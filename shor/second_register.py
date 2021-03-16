import numpy as np


class ShorsSecondRegister:
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

    def generate_states(self):
        '''

        Returns:

        '''
        binaries = []
        registers = []
        my_dict = {0: np.array([1, 0]), 1: np.array([0, 1])}
        for i in range(0, self.state_no):
            binaries.append(format(i, ''.join('0' + str(self.t) + 'b')))
        for j in range(0, self.state_no):
            temp = 0
            initial = np.kron(my_dict[int(binaries[j][0])], my_dict[int(binaries[j][1])])
            if self.t == 2:
                registers.append(initial)
            else:
                for k in range(2, self.t):
                    temp = np.kron(my_dict[int(binaries[j][k])], initial)
                    initial = temp
            registers.append(temp)
        return np.array(registers)[::-1]

    def toffoli_gate(self):
        '''

        Returns:

        '''
        input_value = np.empty([self.t, self.t])
        output_value = np.empty([self.state_no, self.t])
        function = np.empty([self.state_no, 2, self.t])
        for i in range(0, self.state_no):
            value_input = map(int, str(format(i, ''.join('0' + str(self.t) + 'b'))))
            value_output = map(int, str(format((self.a ** i) % self.N, ''.join('0' + str(self.t) + 'b'))))
            function[i, 0, 0:self.t] = np.array(list(value_input))
            function[i, 1, 0:self.t] = np.array(list(value_output))
        return function.astype(int)

    def binary_to_quantum(self, function):
        '''

        Args:
            function:

        Returns:

        '''
        function_quantum_state = np.empty([self.state_no, 2, self.state_no])
        my_dict = {0: np.array([1, 0]), 1: np.array([0, 1])}
        for i in range(0, self.state_no):
            temp_input = 0
            temp_output = 0
            initial_input = np.kron(my_dict[int(function[i, 0, 0])], my_dict[int(function[i, 0, 1])])
            initial_output = np.kron(my_dict[int(function[i, 1, 0])], my_dict[int(function[i, 1, 1])])
            if self.t == 2:
                function_quantum_state[i, 0, 0:self.state_no] = initial_input
                function_quantum_state[i, 1, 0:self.state_no] = initial_output
            else:
                for k in range(2, self.t):
                    temp_input = np.kron(my_dict[int(function[i, 0, k])], initial_input)
                    temp_output = np.kron(my_dict[int(function[i, 1, k])], initial_output)
                    initial_input = temp_input
                    initial_output = temp_output
            function_quantum_state[i, 0, 0:self.state_no] = initial_input
            function_quantum_state[i, 1, 0:self.state_no] = initial_output
        return function_quantum_state[::-1].astype(int)

    def perform_mod_fn(self):
        '''

        Returns:

        '''
        function = self.binary_to_quantum(self.toffoli_gate())
        initial_states = self.generate_states()
        combined_registers = np.empty([len(initial_states), 2, self.state_no])
        for i in range(0, len(initial_states)):
            find = initial_states[i]
            index = np.int(np.where((function[:, 0, :]).dot(find) == sum(find))[0])
            combined_registers[i, 0, 0:self.state_no] = find
            combined_registers[i, 1, 0:self.state_no] = function[index, 1]
        return combined_registers

    def pick_out_states(self):
        '''

        Returns:

        '''
        register = self.perform_mod_fn()
        unique_states = np.unique(register[:, 1, :], axis=0)
        prob = []
        for i in range(len(unique_states)):
            value = unique_states[i]
            occurances = np.sum((register[:, 1]).dot(value) == sum(value))
            prob.append(occurances / self.state_no)
        prob = np.array(prob)
        value = np.random.choice(len(unique_states), p=prob)
        find = unique_states[value, :].astype(int)
        positions = np.where((register[:, 1, :]).dot(find) == 1)
        return np.array(positions)

    def work_out_original_indexes(self):
        '''

        Returns:

        '''
        return self.pick_out_states()