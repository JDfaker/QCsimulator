#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import cmath
from matplotlib import pyplot as plt 
import random


# In[3]:


class QFT:
    def __init__(self,no_qubits):
        self.n = no_qubits
    def generate_states(self):
        state_list = []
        for i in range(self.n):
            state_list.append([1,0])
        return np.array(state_list)
    def qft_circuit(self):
        states = self.generate_states()
        had = (1/np.sqrt(2))*np.array([[1,1],
                        [1,-1]])
        output = []
        for i in range(1,len(states)+1):
            initial = np.matmul(had,states[i-1])
            for j in range(2,len(states) - i):
                initial = np.matmul(self.r_n(j), initial)
            output.append(initial)
        return np.array(output)
    def r_n(self,n):
        val = cmath.exp((2*np.pi*1j)/(2**n))
        return np.array([[1,0],[0,val]])
    def generate_superposition(self):
        states = self.qft_circuit()
        for i in range(len(states)):
            initial = np.kron(states[0],states[1])
            for j in range(2,len(states)):
                initial = np.kron(states[j],initial)
        return initial
    def plot_qft(self):
        plt.title(' Probaility Amplitudes for a '+str(self.n)+ ' qubit register')
        plt.ylabel('$|\psi|^{2}$')
        plt.xlabel('Number')
        return plt.bar(np.arange(0,2**self.n),(2*(np.real(self.generate_superposition())**2)))


# In[4]:


QFT(10).plot_qft()


# In[3]:


class Shors_Second_Register:
    def __init__(self, a, N):
        self.a = int(a)
        self.N = int(N)
        self.t = int(np.ceil(np.log2(N)))#minimum_number of classical bits
        self.state_no = int(int(2**(np.ceil(np.log2(N)))))#minimum number of quantum states 
    def generate_states(self):
        binaries = []
        registers = []
        my_dict = {0: np.array([1,0]), 1: np.array([0,1])}
        for i in range(0,self.state_no):
            binaries.append(format(i, ''.join('0'+str(self.t)+'b')))
        for j in range(0,self.state_no):
            temp = 0
            initial = np.kron(my_dict[int(binaries[j][0])],my_dict[int(binaries[j][1])])
            if self.t == 2:
                registers.append(initial)
            else:
                for k in range(2,self.t):
                    temp = np.kron(my_dict[int(binaries[j][k])],initial)
                    initial = temp
            registers.append(temp)
        return np.array(registers)[::-1]
    def toffoli_gate(self):
        input_value = np.empty([self.t,self.t])
        output_value  = np.empty([self.state_no,self.t])
        function = np.empty([self.state_no,2,self.t])
        for i in range(0,self.state_no):
            value_input = map(int,str(format(i, ''.join('0'+str(self.t)+'b'))))
            value_output = map(int,str(format((self.a**i)%self.N, ''.join('0'+str(self.t)+'b'))))
            function[i,0,0:self.t] = np.array(list(value_input))
            function[i,1,0:self.t] = np.array(list(value_output))
        return function.astype(int)
    def binary_to_quantum(self,function):
        function_quantum_state = np.empty([self.state_no,2,self.state_no])
        my_dict = {0: np.array([1,0]), 1: np.array([0,1])}
        for i in range(0,self.state_no):
            temp_input = 0
            temp_output = 0
            initial_input = np.kron(my_dict[int(function[i,0,0])],my_dict[int(function[i,0,1])])
            initial_output = np.kron(my_dict[int(function[i,1,0])],my_dict[int(function[i,1,1])])
            if self.t == 2:
                function_quantum_state[i,0,0:self.state_no] = initial_input
                function_quantum_state[i,1,0:self.state_no] = initial_output
            else:
                for k in range(2,self.t):
                    temp_input = np.kron(my_dict[int(function[i,0,k])],initial_input)
                    temp_output = np.kron(my_dict[int(function[i,1,k])],initial_output)
                    initial_input = temp_input
                    initial_output = temp_output
            function_quantum_state[i,0,0:self.state_no] = initial_input
            function_quantum_state[i,1,0:self.state_no] = initial_output
        return function_quantum_state[::-1].astype(int)
    def perform_mod_fn(self):
        function = self.binary_to_quantum(self.toffoli_gate())
        initial_states = self.generate_states()
        combined_registers = np.empty([len(initial_states),2,self.state_no])
        for i in range(0,len(initial_states)):
            find = initial_states[i]
            index = np.int(np.where((function[:,0,:]).dot(find) == sum(find))[0])
            combined_registers[i,0,0:self.state_no] = find
            combined_registers[i,1,0:self.state_no] = function[index,1]
        return combined_registers
    def pick_out_states(self):
        register = self.perform_mod_fn()
        unique_states = np.unique(register[:,1,:],axis = 0 )
        prob = []
        for i in range(len(unique_states)):
            value = unique_states[i]
            occurances = np.sum((register[:,1]).dot(value) == sum(value))
            prob.append(occurances/self.state_no)
        prob = np.array(prob)
        value = np.random.choice(len(unique_states), p = prob)
        find = unique_states[value,:].astype(int)
        positions = np.where((register[:,1,:]).dot(find) == 1)
        return np.array(positions)

    def work_out_original_indexes(self):
        return self.pick_out_states()


# In[4]:


Shors_Second_Register(4,9).work_out_original_indexes()


# In[5]:


class Shors:
    def __init__(self, a, N):
        self.a = int(a)
        self.N = int(N)
        self.t = int(np.ceil(np.log2(N)))#minimum_number of classical bits
        self.state_no = int(int(2**(np.ceil(np.log2(N)))))#minimum number of quantum states 
        self.second_register = Shors_Second_Register(int(a),int(N))
        self.qft = QFT(self.t)
    def compute_r(self):
        indexes = self.second_register.work_out_original_indexes()
        probs = np.abs((np.real(self.qft.generate_superposition())))**2
        states = np.empty([1,len(indexes[0]),2])
        for i in range(len(indexes[0])):
            states[0,i,0] = probs[indexes[0,i]]
            states[0,i,1] = indexes[0,i]
        return states[0,:,1][np.random.choice(len(states[0,:,0]),p = states[0,:,0]/sum(states[0,:,0]))].astype(int)
    def cf(self,n):
        d = self.state_no
        res = []
        q, r = divmod(n, d)
        while r != 0:
            res = res + [q]
            prev_r = r
            q, r = divmod(d, r)
            d = prev_r
        return sum(res + [q])
    def main(self):
        thing = self.compute_r()
        r = self.cf(thing)
        while r%2!=0:
            r = self.cf(self.compute_r())
            a1 = int((self.a**(r/2))-1)
            a2 = int((self.a**(r/2))+1)
        a1 = int((self.a**(r/2))-1)
        a2 = int((self.a**(r/2))+1)
        return np.gcd(a1,self.N),np.gcd(a2,self.N)


# In[95]:


factors = []
n = 78
for i in range(0,500):
    random = np.random.randint(0,n)
    random
    while np.gcd(random,n)!=1:
        random = np.random.randint(0,n)
    factors.append(Shors(random,n).main())
factors = np.array(factors).flatten()


# In[96]:


factors.shape


# In[100]:


plt.hist(factors, bins = range(0,n+1,1))
plt.xlabel('Factors')
plt.ylabel('Frequencies')
plt.title('Histogram of Shors Algorithm Output for  N = '+str(n))


# In[43]:


for i in range(100):
    Shors(5,9).main()


# In[688]:


Shors(7,9).main()


# In[486]:


def cf(n, d):
    res = []
    q, r = divmod(n, d)
    while r != 0:
        res = res + [q]
        prev_r = r
        q, r = divmod(d, r)
        d = prev_r
    return res + [q]


# In[488]:


cf(3,8)


# In[530]:


class Classical_Part():
    def __init__(self,test_range):
        self.range = test_range
    def main(self,N):
        gcd = []
        for i in range(N):
            if np.gcd(i,N)==1:
                gcd.append(Shors(i,N).main())
            else:
                pass
        return gcd
    def testing(self):
        factors = []
        for i in range(6,self.range):
            factors.append(self.main(i))
        return factors


# In[532]:


Classical_Part(21).testing()


# In[ ]:


def continued_fractions(self,n):
    remainders = []
    d = self.t
    for i in range(10):
        remainder = np.floor(n/d)
        if remainder > 0:
            diff = n - remainder*d
            if diff ==0:
                remainders.append(remainder)
                break
            new_frac = n/diff
            n, d = d, diff
        else:
            n,d = d,n
            if d==0:
                break
        remainders.append(remainder)
    return np.array(remainders)
def compute_fractions(self,n):
    remainders = self.continued_fractions(n)
    new_value = 0
    fractions = [] 
    new_denominator = 0
    new_numerator = 0
    for i in range(len(remainders[0:])-1,1,-1):
        if i == (len(remainders[0:]))-1 :
            new_denominator = int(remainders[i])
            new_numerator = 1 +  new_denominator
            new_denominator, new_numerator = new_numerator, new_denominator
        else:
            new_value = int(remainders[i])
            new_numerator = new_numerator + (new_value * new_denominator)
            new_denominator, new_numerator = new_numerator, new_denominator
    return new_denominator

