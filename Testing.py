#!/usr/bin/env python
# coding: utf-8

# In[1]:


import qiskit
import numpy as np 
import QuantumCircuit
import sparse_matrix
import scipy.sparse


# In[10]:


class Tensor_Testing:
    def __init__(self, seed=None):
        self.seed = np.random.seed(seed)
        
    def get_random_matrices(self):
        W = np.random.randint(1,50,(np.random.randint(2,10),np.random.randint(2,10)))
        V = np.random.randint(1,50,np.shape(np.transpose(W)))
        return V, W
    
    def tensor_product_test(self, V=None, W=None):
        if V==None or W==None:
            V, W = self.get_random_matrices()
        tensor_product_test = qiskit.aqua.utils.tensorproduct(np.array(V),np.array(W))
        tensor_product_func = QuantumCircuit.QuantumCircuit(2).tensor_product(np.array(V),np.array(W))
        assert not np.any(tensor_product_func-tensor_product_test) == True, "The tensor products do not match. {} != {}".format(tensor_product_func, tensor_product_test)
    
    #def dot_product_test(self):
    #    V, W = self.get_random_matrices()
    #    dot_product_test = np.matmul(np.array(V), np.array(W))
    #    dot_product_func = QuantumCircuit.QuantumCircuit(2).dot_product(np.array(V), np.array(W))
    #    all_zeros = not np.any(np.array([dot_product_func-dot_product_test]))
    #    assert not np.any(dot_product_func-dot_product_test) == True, "The dot products do not match"


# In[11]:


class Gate_For_Test:
    '''
    Defines gates to be tested.
    
    qiskit_name - defines the name of the method for the gate in the qiskit library. Must be str.
    our_name - defines the corresponding method name in our simulator. Must be str.
    num_qubits - defines how many qubits the gate operates on and therefore how many must be provided to not raise an error.
                Must be int.
    '''
    def __init__(self, qiskit_name, our_name, num_qubits):
        self.qiskit_name = qiskit_name
        self.our_name = our_name
        self.num_qubits = num_qubits


# In[21]:


class Gate_Testing:
    '''
    Tests the given gate with both qiskit and our simulator and compares the results
    
    Quick gate name reference: 
    Hadamard - h
    Pauli-X - x
    Pauli-Y - y
    Pauli-Z - z
    SWAP - swap
    Controlled-Z - cz
    Controlled-NOT - cx
    
    Parameters:
    qubits - size of test circuits. Type should be int.
    
    gate_input - gate to be tested. Should be a str corresponding to the qiskit name of the gate (for brevity; this can be 
                changed easily if required).
    test_qubit_1 - first qubit to apply gate to. Will raise an error if outside the registry. Should be int.
    test_qubit_2 - second target for use with 2- and 3-qubit gates ONLY. Type should be int. For gates with a control and a 
                target, the first test qubit is the control and the second is the target.
    test_qubit_3 - third target for use with 3-qubit gates ONLY. Type should be int.
    test_angle - angle for use with Rotation gates. Type can be int or float.
    '''
    def __init__(self, qubits):
        self.num_qubits = qubits
        self.qiskit_circ = qiskit.QuantumCircuit(qubits)
        self.our_circ = QuantumCircuit.QuantumCircuit(qubits)
        
        self.gate_database = [Gate_For_Test("h", "apply_hardmard", 1), Gate_For_Test("x", "apply_pauliX", 1), 
                              Gate_For_Test("y", "apply_pauliY", 1), Gate_For_Test("z", "apply_pauliZ", 1), 
                              Gate_For_Test("swap", "apply_swap", 2), Gate_For_Test("cz", "apply_controlZ", 2),
                              Gate_For_Test("cx", "apply_cnot", 2)]
    
    def run_qiskit_circuit(self,circ):
        backend = qiskit.Aer.get_backend('statevector_simulator')
        job = qiskit.execute(circ, backend)
        result = job.result()
        outputstate = result.get_statevector(circ, decimals=3)
        return outputstate
    
    def qiskit_gate_test(self, qiskit_gate, test_qubits, test_angle):
        #if qiskit_gate == "crz":
        #    exec("self.qiskit_circ." + str(qiskit_gate) + "(" str(test_angle) + ", " + str(test_qubits[0]) + ", " + str(test_qubits[1]) + ")")
        if len(test_qubits) == 3:
            exec("self.qiskit_circ." + str(qiskit_gate) + "(" + str(test_qubits[0]) + ", " + str(test_qubits[1]) + ", " + str(test_qubits[2]) + ")")
        elif len(test_qubits) == 2:
            exec("self.qiskit_circ." + str(qiskit_gate) + "(" + str(test_qubits[0]) + ", " + str(test_qubits[1]) + ")")
        elif len(test_qubits) == 1:
            exec("self.qiskit_circ." + str(qiskit_gate) + "(" + str(test_qubits[0]) + ")")
        else:
            raise ValueError("Invalid number of target qubits provided to qiskit_gate_test ({})".format(len(test_qubits)))
        qiskit_output = self.run_qiskit_circuit(self.qiskit_circ)
        return qiskit_output
   
    def our_gate_test(self, our_gate, test_qubits):
        if len(test_qubits) == 3:
            exec("self.our_circ." + str(our_gate) + "(" + str(test_qubits[0]) + ", " + str(test_qubits[1]) + ", " + str(test_qubits[2]) + ")")
        elif len(test_qubits) == 2:
            exec("self.our_circ." + str(our_gate) + "(" + str(test_qubits[0]) + ", " + str(test_qubits[1]) + ")")
        elif len(test_qubits) == 1:
            exec("self.our_circ." + str(our_gate) + "(" + str(test_qubits[0]) + ")")
        else:
            raise ValueError("Invalid number of target qubits provided to our_gate_test ({})".format(len(test_qubits)))
        our_output = np.transpose(self.our_circ.state)[0].astype(complex)
        return our_output
    
    def gate_test(self, gate_input, test_qubit_1, test_qubit_2=None, test_qubit_3=None, test_angle=None):
        assert np.size(np.unique(np.array([x for x in [test_qubit_1, test_qubit_2, test_qubit_3] if x != None]))) == np.size(np.array([x for x in [test_qubit_1, test_qubit_2, test_qubit_3] if x != None])), "Test qubits may not contain any duplicates. {} != {}".format(np.unique(np.array([x for x in [test_qubit_1, test_qubit_2, test_qubit_3] if x != None])), np.array([x for x in [test_qubit_1, test_qubit_2, test_qubit_3] if x != None]))
        
        test_qubits_qiskit = [x for x in [test_qubit_1, test_qubit_2, test_qubit_3] if x != None]
        test_qubits_ours = np.subtract(np.dot(np.ones_like(test_qubits_qiskit), self.num_qubits-1), test_qubits_qiskit)
        
        gate_this_test = next((x for x in self.gate_database if x.qiskit_name == gate_input), None)
        
        assert len(test_qubits_qiskit)<=gate_this_test.num_qubits, "Number of test qubits provided is greater than number given gate operates on. {} > {}".format(len(test_qubits_qiskit), gate_this_test.num_qubits)
        assert np.max(np.array(test_qubits_qiskit))<self.num_qubits, "Test qubit index greater than largest register index. {} >= {}".format(np.max(np.array(test_qubits_qiskit)), self.num_qubits)
        
        qiskit_output = self.qiskit_gate_test(gate_this_test.qiskit_name, test_qubits_qiskit, test_angle)
        our_output = self.our_gate_test(gate_this_test.our_name, test_qubits_ours)
        
        #tolerance implemented to deal with miniscule errors that were throwing the test off
        assert (np.abs(np.subtract(qiskit_output, our_output)) <= 0.00005).all(), "The states after the gate's application do not match. {} != {}".format(qiskit_output, our_output)


# In[43]:


class Sparse_Testing:
    def __init__(self, seed=None):
        self.seed = np.random.seed(seed)
        
    def get_random_matrices(self):
        test_matrix_1 = np.random.randint(1,50,(np.random.randint(1,10),np.random.randint(1,10)))
        test_matrix_2 = np.random.randint(1,50,np.shape(np.transpose(test_matrix_1)))
        return test_matrix_1, test_matrix_2
    
    def basic_sparsify_test(self, test_matrix=None):
        '''
        Ensures that the sparsify and numpy methods are consistent for the given matrix by converting to sparse and back again
        and comparing the result with the original matrix passed to the method.
        
        Parameters:
        test_matrix - matrix to convert. Must be a numpy array.  Will be generated randomly if not provided.
        '''
        if test_matrix==None:
            test_matrix = self.get_random_matrices()[0]
        original_matrix = test_matrix
        test_matrix = sparse_matrix.SparseMatrix.sparsify(test_matrix)
        test_matrix = sparse_matrix.SparseMatrix.numpy(test_matrix)
        assert np.equal(test_matrix, original_matrix).all(), "The matrix does not match its original form. {} != {}".format(test_matrix, original_matrix)
    
    def sparse_dot_test(self, test_matrix_1=None, test_matrix_2=None):
        '''
        Performs a dot product with both scipy.sparse and our class and confirm that the results are identical.
        
        Parameters:
        test_matrix_1, test_matrix_2 - matrices to multiply. Must be numpy arrays. Will be generated randomly if not provided.
        '''
        if test_matrix_1==None or test_matrix_2==None:
            test_matrix_1, test_matrix_2 = self.get_random_matrices()
            
        our_dot = sparse_matrix.SparseMatrix.numpy(sparse_matrix.SparseMatrix.sparsify(test_matrix_1).dot(sparse_matrix.SparseMatrix.sparsify(test_matrix_2)))
        scipy_dot = scipy.sparse.csc_matrix(test_matrix_1).dot(scipy.sparse.csc_matrix(test_matrix_2)).toarray()
        assert np.equal(our_dot, scipy_dot).all(), "The output matrices do not match. {} != {}".format(our_dot, scipy_dot)
        
    def sparse_tensor_dot_test(self, test_matrix_1=None, test_matrix_2=None):
        '''
        Performs a tensor dot product with both scipy.sparse and our class and confirm that the results are identical.
        
        Parameters:
        test_matrix_1, test_matrix_2 - matrices to multiply. Must be numpy arrays. Will be generated randomly if not provided.
        '''
        if test_matrix_1==None or test_matrix_2==None:
            test_matrix_1, test_matrix_2 = self.get_random_matrices()
            
        our_tensor_dot = sparse_matrix.SparseMatrix.numpy(sparse_matrix.SparseMatrix.sparsify(test_matrix_1).tensordot(sparse_matrix.SparseMatrix.sparsify(test_matrix_2)))
        scipy_tensor_dot = scipy.sparse.kron(scipy.sparse.csc_matrix(test_matrix_1), scipy.sparse.csc_matrix(test_matrix_2)).toarray()
        assert np.equal(our_tensor_dot, scipy_tensor_dot).all(), "The output matrices do not match.{} != {}".format(our_tensor_dot, scipy_tensor_dot)
        
    def get_attribute_test(self, operation, *args, test_matrix=None):
        '''
        Performs the 'get' method corresponding to the operation given on the sparse form of test_matrix, passing the given 
        arguments, then compares the output to the original dense matrix to ensure correctness.
        
        Operation quick reference:
        row - Gets all nonzero entries of stated row as a dictionary with their positions in the row as keys.
        col - Gets all nonzero entries of stated column as a dictionary with their positions in the column as keys.
        value - Gets value at the given row and column
        nonzero_rows - Gets all rows with nonzero elements
        nonzero_cols - Gets all columns with nonzero elements
        
        Parameters:
        test_matrix - matrix to operate on. Must be a numpy array. Will be generated randomly if one is not provided.
        operation - operation to perform. Must be one of those stated above or the test will return an error.
        args - arguments of the test function. Should be a single int for row and col, a pair of ints (row, col) for value, and 
        nothing for nonzero_rows and nonzero_cols.
        '''
        if test_matrix==None:
            test_matrix = self.get_random_matrices()[0]
        
        output = []
        exec("output.append(sparse_matrix.SparseMatrix.sparsify(test_matrix).get_" + str(operation) + "(*args))")
        output = output[0]
        if operation == "col":
            out_keys = output.keys()
            for key in out_keys:
                assert output[key] == test_matrix[int(key)][args[0]], "An incorrect nonzero value has been retrieved by get_col. {} != {}".format(output[key], test_matrix[int(key)][args[0]])
            assert not np.any(np.array([test_matrix[x][args[0]] for x in range(len(test_matrix)) if x not in out_keys])), "Nonzero values in the column have not been retrived by get_col: {}.".format(np.array([test_matrix[x][args[0]] for x in range(len(test_matrix)) if x not in out_keys]))
        elif operation == "row":
            out_keys = output.keys()
            for key in out_keys:
                assert output[key] == test_matrix[args[0]][int(key)], "An incorrect nonzero value has been retrieved by get_row. {} != {}".format(output[key], test_matrix[args[0]][int(key)])
            assert not np.any(np.array([test_matrix[args[0]][x] for x in range(len(test_matrix[args[0]])) if x not in out_keys])), "Nonzero values in the row have not been retrived by get_row: {}.".format(np.array([test_matrix[args[0]][x] for x in range(len(test_matrix[args[0]])) if x not in out_keys]))
        elif operation == "value":
            assert output == test_matrix[args], "The found value does not match the corresponding value in the dense matrix. {} != {}".format(output, test_matrix[args])
        elif operation == "nonzero_rows":
            assert output == [x for x in range(len(test_matrix)) if test_matrix[x].any() == True], "Rows with nonzero elements exist that have not been retrived by get_nonzero_rows: {}.".format([x for x in range(len(test_matrix)) if test_matrix[x].any() == True and x not in output])
        elif operation == "nonzero_cols":
            assert output == [x for x in range(len(test_matrix[0])) if np.array([test_matrix[i][x] for i in range(len(test_matrix))]).any() == True],"Columns with nonzero elements exist that have not been retrived by get_nonzero_cols: {}.".format([x for x in range(len(test_matrix[0])) if np.array([test_matrix[i][x] for i in range(len(test_matrix))]).any() == True and x not in output])
        else:
            raise ValueError("Invalid operation provided to get_attribute_test ({})".format(operation))


# In[57]:


class Grover_Testing:
    def __init__(self, qubits):
        '''
        Initialises Grover test by creating circuits and then putting them into a state of superposition by applying a hadamard
        gate to each qubit.
        
        Parameters:
        qubits - number of qubits for test circuits. Must be int.
        '''
        self.qiskit_circ = qiskit.QuantumCircuit(qubits)
        self.our_circ = QuantumCircuit.QuantumCircuit(qubits)
        self.num_qubits = qubits
        
        #initialise superposition
        for i in range(self.num_qubits):
            self.qiskit_circ.h(i)
            self.our_circ.apply_hardmard(i)
    
    def our_grover_test(self):
        '''
        Performs a Grover test on our simulator, stopping iterations when the target state has been located to a suitable
        accuracy.
        '''
        while np.transpose(self.our_circ.state)[0][self.target] < 0.999:
            self.our_circ.apply_grover_oracle(self.target)
            self.our_circ.apply_amplification()
        
        return np.transpose(self.our_circ.state)[0]
        
    def run_qiskit_circuit(self):
        '''
        Same qiskit circuit process as used in Gate_Testing.
        '''
        backend = qiskit.Aer.get_backend('statevector_simulator')
        job = qiskit.execute(self.qiskit_circ, backend)
        result = job.result()
        outputstate = result.get_statevector(self.qiskit_circ, decimals=3)
        return outputstate
        
    def qiskit_oracle(self):
        '''
        Uses the same method as our simulator to construct an oracle matrix, then converts it to a qiskit gate.
        '''
        #construct a gate matrix the same way as in our simulator
        I = np.eye(2 ** self.num_qubits)
        oracle = I
        if isinstance(self.target, int):
            oracle[self.target][self.target] = -1
        else:
            for mark in self.target:
                oracle[mark][mark] = -1

        #convert oracle to a qiskit gate
        return qiskit.extensions.UnitaryGate(oracle)
    
    def qiskit_diffuser(self):
        '''
        General Grover diffuser converted from the version given at https://qiskit.org/textbook/ch-algorithms/grover.html.
        Takes the number of qubits and constructs a diffuser of that size.
        '''
        for qubit in range(self.num_qubits):
            self.qiskit_circ.h(qubit)
        for qubit in range(self.num_qubits):
            self.qiskit_circ.x(qubit)
        self.qiskit_circ.h(self.num_qubits-1)
        self.qiskit_circ.mct(list(range(self.num_qubits-1)), self.num_qubits-1)  # multi-controlled-toffoli
        self.qiskit_circ.h(self.num_qubits-1)
        for qubit in range(self.num_qubits):
            self.qiskit_circ.x(qubit)
        for qubit in range(self.num_qubits):
            self.qiskit_circ.h(qubit)
        return self.qiskit_circ  
    
    def qiskit_grover_test(self):
        '''
        Performs a Grover test on our simulator, stopping iterations when the target state has been located to a suitable
        accuracy. Qiskit builds circuits and then runs them, rather than running automatically at every step like our simulator,
        meaning that the circuit must be manually run at each stage.
        '''
        out = self.run_qiskit_circuit()
        while out[self.target]*np.conj(out[self.target]) < complex(0.999)*np.conj(complex(0.999)):
            self.qiskit_circ.append(self.qiskit_oracle(), range(self.num_qubits))
            self.qiskit_diffuser()
            out = self.run_qiskit_circuit()
        return out
    
    def grover_test(self,target):
        '''
        Performs a Grover test with both simulators searching for the target state. Compares the index of the found state, and
        the value of the state itself within a given tolerance (fairly large, as the discrepency between the two simulators
        increases significantly with the number of qubits in the circuit). Because of differences in output format, the squares
        of the results are compared instead of the raw results.
        
        Parameters:
        target - target state expressed as a decimal int.
        '''
        self.target = target
        qiskit_result = self.qiskit_grover_test()
        our_result = self.our_grover_test()
        assert np.where(qiskit_result**2 >= complex(0.999)**2) == np.where(our_result >= 0.999), "The simulators did not find the same state. {} != {}".format(np.where(qiskit_result**2 >= complex(0.999)**2), np.where(our_result >= 0.999))
        assert np.abs(np.real(np.amax(qiskit_result*np.conj(qiskit_result))) - np.amax(our_result)**2) <= 0.005, "The converted values of the found states do not match to within +/- 0.005. {} != {}".format(np.real(np.amax(qiskit_result*np.conj(qiskit_result))), np.amax(our_result)**2)


# In[58]:


Tensor_Testing().tensor_product_test()

Gate_Testing(5).gate_test("h", 0)
Gate_Testing(5).gate_test("x", 0)
#Gate_Testing(5).gate_test("y", 0)
Gate_Testing(5).gate_test("z", 0)
Gate_Testing(5).gate_test("swap", 0, 1)
Gate_Testing(5).gate_test("cz", 0, 1)
Gate_Testing(5).gate_test("cx", 0, 1)

Sparse_Testing().basic_sparsify_test()
Sparse_Testing().sparse_dot_test()
Sparse_Testing().sparse_tensor_dot_test()
Sparse_Testing().get_attribute_test("col", 0)
Sparse_Testing().get_attribute_test("row", 0)
Sparse_Testing().get_attribute_test("value", 0, 0)
Sparse_Testing().get_attribute_test("nonzero_rows")
Sparse_Testing().get_attribute_test("nonzero_cols")

Grover_Testing(5).grover_test(2)


# In[ ]:




