#!/usr/bin/env python
# coding: utf-8

# In[ ]:



def run_time(qubit_no):
    times = []
    for i in range(3,qubit_no):
        value = 2**i
        random = 0
        datapoints = 0 
        start = time.process_time()
        random = np.random.randint(0,value)
        while np.gcd(random,n)!=1:
            random = np.random.randint(0,value)
        datapoints = np.array(Shors(random,value).main())
        times.append(time.process_time() - start)
        print(i)
    plt.figure(figsize = (10,5))
    plt.plot(np.arange(3,12),np.array(times))
    plt.ylabel('Time(s)')
    plt.xlabel('Qubit Number')
    plt.title('Time for a single itteration as a function of # qubits')


# In[ ]:


def test_outputs(n):
success = []
factors_true = []
accuracy = []
inaccuracy = []
for i in range (0, n):
        value = np.fromiter(sympy.ntheory.factorint(i).keys(), dtype=int)
        factors_true.append(value.tolist())
for i in range (5,n): 
    factors_new = []
    success_temp = 0 
    for j in range(0,50):
        random = np.random.randint(0,i)
        while np.gcd(random,n)!=1:
            random = np.random.randint(0,i)
        datapoints = np.array(Shors(random,i).main())
        factors_new.append(datapoints)
    factors_new = np.unique(np.array(factors_new).flatten())
    for k in range(len(factors_true[i])):
        if np.sum(np.where(factors_new % factors_true[i][k] == 0))>=1:
            success_temp += 1
        else:
            pass
    accuracy.append(success_temp/len(factors_new)*100)
    inaccuracy.append((1-(success_temp/len(factors_new)))*100)
    fig, ax = plt.subplots(figsize = (15,5))
    labels = np.arange(5,n)
    ax.bar( labels, accuracy, width = 0.5, label='Correct')
    ax.bar( labels, inaccuracy, width = 0.5, bottom = accuracy, label='Incorrect')
    ax.legend()
    ax.set_title('Number of Correct Factors determined for $ 5 < N \leq 50$')
    ax.set_ylabel('Percentage of factors')
    ax.set_xlabel('N')

