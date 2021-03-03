import qiskit
import numpy as np 
import QuantumCircuit

class Tensor_Dot_Testing:
    def __init__(self, part):
        self.part = str(part)
        #self.seed = np.random.seed()
    
    def test(self):
        return exec('self.'+self.part+'()')
    
    def get_random_matrices(self):
        W = np.random.randint(1,50,(np.random.randint(1,10),np.random.randint(1,10)))
        V = np.random.randint(1,50,np.shape(W))
        return V, W
    
    def tensor_product(self):
        V, W = self.get_random_matrices()
        tensor_product_test = qiskit.aqua.utils.tensorproduct(np.array(V),np.array(W))
        tensor_product_func = QuantumCircuit.QuantumCircuit(2).tensor_product(np.array(V),np.array(W))
        assert not np.any(tensor_product_func-tensor_product_test) == True, "The tensor products do not match"
    
    def dot_product(self):
        V, W = self.get_random_matrices()
        dot_product_test = np.matmul(np.array(V), np.array(W))
        dot_product_func = QuantumCircuit.QuantumCircuit(2).dot_product(np.array(V), np.array(W))
        all_zeros = not np.any(np.array([dot_product_func-dot_product_test]))
        assert not np.any(dot_product_func-dot_product_test) == True, "The dot products do not match"

class Gate_For_Test:
    '''
    Defines gates to be tested.
    
    qiskit_name defines the name of the method for the gate in the qiskit library
    our_name defines the corresponding method name in our simulator
    num_qubits defines how many qubits the gate operates on and therefore how many must be provided to not raise an error
    '''
    def __init__(self, qiskit_name, our_name, num_qubits):
        self.qiskit_name = qiskit_name
        self.our_name = our_name
        self.num_qubits = num_qubits
        
class Gate_Testing:
    '''
    Tests the given gate with both qiskit and our simulator and compares the results
    
    Quick gate name reference: 
    Hadamard - h
    Pauli-X - x
    Pauli-Y - y
    Pauli-Z - z
    SWAP - swap
    
    Parameters:
    gate_input - gate to be tested. Should be the qiskit name of the gate (for brevity; this can be changed easily if required).
    qubits - size of test circuits
    test_qubits - qubits to apply gate to. Will throw an error if outside the registry or an incorrect number of targets.
    '''
    def __init__(self, gate_input, qubits, *test_qubits):
        self.gate_input = str(gate_input)
        self.qiskit_circ = qiskit.QuantumCircuit(qubits)
        self.our_circ = QuantumCircuit.QuantumCircuit(qubits)
        self.num_qubits = qubits
        self.test_qubits = test_qubits
        assert len(self.test_qubits)<=self.num_qubits, "More test qubits provided than number of qubits present."
        
        self.gate_database = [Gate_For_Test("h", "apply_hardmard", 1), Gate_For_Test("x", "apply_pauliX", 1), 
                              Gate_For_Test("y", "apply_pauliY", 1), Gate_For_Test("z", "apply_pauliZ", 1), 
                              Gate_For_Test("swap", "apply_swap", 2)]
    
    def run_qiskit_circuit(self,circ):
        backend = qiskit.Aer.get_backend('statevector_simulator')
        job = qiskit.execute(circ, backend)
        result = job.result()
        outputstate = result.get_statevector(circ, decimals=3)
        return outputstate
    
    def qiskit_gate_test(self, qiskit_gate, *test_qubits):
        exec("self.qiskit_circ." + str(qiskit_gate) + "(" + str(*test_qubits) + ")")
        qiskit_output = self.run_qiskit_circuit(self.qiskit_circ)
        return qiskit_output
   
    def our_gate_test(self, our_gate, *test_qubits):
        exec("self.our_circ." + str(our_gate) + "(" + str(*test_qubits) + ")")
        our_output = self.our_circ.state 
        return our_output
    
    def test_gate(self):
        gate_this_test = next((x for x in self.gate_database if x.qiskit_name == self.gate_input), None)
        
        assert len(self.test_qubits)<=gate_this_test.num_qubits, "Number of test qubits provided is greater than number given gate operates on"
        assert np.max(np.array(self.test_qubits))<self.num_qubits, "Test qubit index greater than largest register index"
        
        qiskit_output = self.qiskit_gate_test(gate_this_test.qiskit_name, *self.test_qubits)
        our_output = self.our_gate_test(gate_this_test.our_name, *self.test_qubits)
        
        assert qiskit_output.all() == our_output.all(), "The states after the gate's application do not match."
        print("Success!"
