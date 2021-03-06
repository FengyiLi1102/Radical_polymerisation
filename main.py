# Radical polymerisation simulation
# Author: Joanna Jiang, Angela Pan, Rachael Wang, Hanjia Zhang 
#--------------------------------
import numpy as np
from timeit import default_timer as timer
import math
import scipy.constants as spc
import copy as cp
import matplotlib.pyplot as plt


class myDict(dict):
    """
    Two functions for this simulation to run more effectively.
    """
    def check_existence(self, target):
        """
        Check if a chain is existed in the dictionary or not.
        Input parameters
            self: dictionary required to be checked
            target: key for checking
        Output parameters
            if not existing, dictionary with added {key:0}, 
            and otherwise return the original dictionary
        """
        keys = list(self.keys())
        if target in keys:
            pass
        else:
            self.update({target:0})

    def delItem(self, target):
        """
        Delete a pair if the value = 0
        Input parameters
            self: dictionary required to be checked
            target: key to check
        """
        if self[target] == 0:
            del self[target]
        else:
            pass


def chain_prob_0(chainDict):
    """
    Calculate the probability for which chain will be chosen to react
    Applicable for reaction 2.1, 3.1, 3.2, 3.3, 4.2, 4.3 and 5.3 (one polymer with n monomers)
    reaction number = 1, 2, 3, 4, 6, 7, 10
    Input parameters
        chainDict: dictionary storing n and chain number
    Output parameters
        n of the chain reacting
    """
    totalNum = sum(chainDict.values()) #calculate the total number of chains
    probChain = np.asarray(list(chainDict.values())) / totalNum #calculate the probability
    chainTypes = np.asarray(list(chainDict.keys())) #make all of the keys into a numpy array
    chainKey = np.random.choice(chainTypes, size=1, replace=True, p=probChain)[0] #random choose a chain to react according to the probabilty
    
    return chainKey


def chain_prob_1(chainDict_m, chainDict_n, flag=False):
    """
    Calculate the probability for which chains combination will react
    Applicable for reaction 4.1, 5.1, 5.2 (two polymers with n and m monomers respectively)
    reaction number = 5, 8, 9, 11
    Input parameters
        chainDict_m: dictionary storing chain type and chain number
        chainDict_n: dictionary storing chain type and chain number
        flag: True for reaction containing two variables and False for only one
    Output parameters
        [one key for chainDict_m, one key for chainDict_n]
    """
    chainNums_m = np.asarray(list(chainDict_m.values()))
    chainTypes_m = np.asarray(list(chainDict_m.keys()))
    chainNums_n = np.asarray(list(chainDict_n.values()))
    chainTypes_n = np.asarray(list(chainDict_n.keys()))

    # Combination of different chains with specific n, m values
    combinations = np.array([[i, j] for i in chainTypes_m for j in chainTypes_n])
    num_product = np.array([[i, j] for i in chainNums_m for j in chainNums_n])
        
    # reaction 5.1 and 5.2
    # Find the combinations with the same length of free radicals
    if flag:
        mask = np.where(combinations[:,0] == combinations[:,1], True, False)
        mask = np.vstack([mask, mask]).T
        num_product = np.where(mask, np.vstack([num_product[:,0], 
                               num_product[:,1]-1]).T, num_product)
    else: # reaction 4.1 and 5.4
        pass

    #calculate the probability
    prob = num_product[:,0] * num_product[:,1]
    probabilities = prob / sum(num_product[:,0] * num_product[:,1])
    combination_index = np.random.choice(np.arange(len(combinations)), size=1, 
                                         replace=True, p=probabilities)[0]

    return combinations[combination_index] # [n, m]


# reaction 1.1
def reaction11():
    global I, Rn_d
    I -= 1
    #check the existence of R0 and update the value of it
    Rn_d.check_existence(0)
    Rn_d.update({0: Rn_d[0]+2})

#reaction 2.1, 3.1, 3.2, 3.3, 4.2, 4.3
def reaction_prob_0(num):
    # Combine reaction 2.1, 3.1, 3.2 and 3.3 together
    # Reaction 2.1 and 3.1
    global M, T, Rn_d, d_TRn, TRn, RnTRm
    if num in [1, 2]:
        reacts = chain_prob_0(Rn_d) # returns which n to react
        Rn_d.update({reacts: Rn_d[reacts]-1})
        Rn_d.delItem(reacts)
        if num == 1: # reaction 2.1
            M -= 1
            newChain = reacts + 1
            Rn_d.check_existence(newChain)
            Rn_d.update({newChain: Rn_d[newChain]+1})
        else: # reaction 3.1
            T -= 1
            d_TRn.check_existence(reacts)
            d_TRn.update({reacts: d_TRn[reacts]+1})
    # reaction 3.2 and 3.3
    elif num in [3, 4]:
        reacts = chain_prob_0(d_TRn)
        d_TRn.update({reacts: d_TRn[reacts]-1})
        d_TRn.delItem(reacts)
        if num == 3: # reaction 3.2
            Rn_d.check_existence(reacts)
            Rn_d.update({reacts: Rn_d[reacts]+1})
            T += 1
        else: # reaction 3.3
            Rn_d.check_existence(0)
            Rn_d.update({0: Rn_d[0]+1})
            TRn.check_existence(reacts)
            TRn.update({reacts: TRn[reacts]+1})
    else: # reaction 4.2 and 4.3
        reacts = chain_prob_0(RnTRm) # string 'n,m'
        RnTRm.update({reacts: RnTRm[reacts]-1})
        RnTRm.delItem(reacts)
        n, m = int(list(reacts)[0]), int(list(reacts[2]))
        if num == 6: # reaction 4.2
            Rn_d.check_existence(n)
            Rn_d.update({n: Rn_d[n]+1})
            TRn.check_existence(m)
            TRn.update({m: TRn[m]+1})
        else: # reaction 4.3
            Rn_d.check_existence(m)
            Rn_d.update({m: Rn_d[m]+1})
            TRn.check_existence(n)
            TRn.update({n: TRn[n]+1})

#reaction 4.1, 5.1, 5.2
def reaction_prob_1(num):
    global Rn_d, TRn, RnTRm
    if num == 5: # reaction 4.1
        reacts = chain_prob_1(Rn_d, TRn)
        n, m = reacts[0], reacts[1]
        Rn_d.update({n: Rn_d[n]-1})
        Rn_d.delItem(n)
        TRn.update({m: TRn[m]-1})
        TRn.delItem(m)
        newChain = str(n)+','+str(m)
        RnTRm.check_existence(newChain)
        RnTRm.update({newChain: RnTRm[newChain]+1})
    else: # reaction 5.1 or 5.2
        reacts = chain_prob_1(Rn_d, Rn_d, True)
        n, m = reacts[0], reacts[1]
        Rn_d.update({n: Rn_d[n]-1})
        Rn_d.delItem(n)
        Rn_d.update({m: Rn_d[m]-1})
        Rn_d.delItem(m)

        if num == 8: # reaction 5.1
            """
            Places for updating the concentration of the P(n+m)
            """
        else: # reaction 5.2
            """
            For updating the concentration of the Pn and Pm
            """


global M, Rn_d, d_TRn, TRn, RnTRm
# Global constants required to be tuned before starting the simulation
# Kinetic constant (in the unit of mM and h)
k11 = 0.36/3600 * 10E3 # 0.36
k21 = 3.6*10**15*3600 # E7
k31 = k41 = 3.6*10**18*3600 # E10
k32 = k33 = k42 = k43 = 18*10**25/3600 # E8
k51 = k52 = k53 = k54 = 3.6*10**20*3600 # E11
t_max = 36000           # Total time for the simulation
interval = 100          # time interval for n and P(n)
M = 5                   # Monomer concentration (mM/L)
V = 0.1                 # Reaction volume       (L)
NA = spc.N_A            # Avogadro's constant   (per mole)
alpha = 1/(V*NA)        # Conversion factor from N/L to n/L
I_arr, T_arr = np.array([10E-1]), np.array([10E-1])
initial_concentrations = np.array([[I_arr[i], T_arr[j]] for i in range(len(I_arr)) 
                                    for j in range(len(T_arr))])

log = open('settings.txt', 'w')
log.write('K11: {0}\n'.format(k11))
log.write('K21: {0}\n'.format(k21))
log.write('K31, K41: {0}\n'.format(k31))
log.write('K32, K33, K42, K43: {0}\n'.format(k32))
log.write('K51, K52, K53, K54: {0}\n'.format(k51))
log.write('Total time: {0}s\n'.format(t_max))
log.write('[M]: {0}mM\n'.format(M))
log.write('V: {0}L\n'.format(V))
log.write('[I]: {0}mM\n'.format(I_arr[0]))
log.write('[T]: {0}mM\n'.format(T_arr[0]))
log.write('Data recorded interval: {0}s\n'.format(interval))
log.close()

kinetic_constant = {"k11": k11, 
                    "k21": k21, 
                    "k31": k31, 
                    "k32": k32, 
                    "k33": k33, 
                    "k41": k41, 
                    "k42": k42, 
                    "k43": k43, 
                    "k51": k51, 
                    "k52": k52, 
                    "k53": k53, 
                    "k54": k54}
settingNum = 0

# Starts the simulation for each initial setting
for IandT in initial_concentrations:
    I0 = int(IandT[0] * 1E-3 * V * NA)
    T0 = int(IandT[1] * 1E-3 * V * NA)
    M0 = int(M * V * NA)
    # print(I0)
    t = 0                       # System time  (s)
    I = I0                      # Initial initiator concentration  (#)
    T = T0                      # Initial agent concentration      (#)
    M = M0                      # Initial monomer concentration    (#)
    Rn_d = myDict({0:0})        # Rn dot : number of chains
    d_TRn = myDict()            # dot TRn : number of chains
    TRn = myDict()    
    RnTRm = myDict()  
    times = 0                   # variable to update the name of the txt file
    wDict = myDict()            # store the data required for distribution plots
    t_and_R0 = myDict({0:0})    # store the data required for concentration plots

    # A specific initial setting
    while(t < t_max):
        ratesList = [] # a list of reaction rate

        # Step 1: calculate reaction rate
        R11 = 2 * kinetic_constant["k11"] * I * alpha
        ratesList.append(R11)

        R21 = kinetic_constant["k21"] * sum(Rn_d.values()) * alpha**2 * M
        ratesList.append(R21)

        R31 = kinetic_constant["k31"] * alpha**2 * sum(Rn_d.values()) * T
        ratesList.append(R31)

        R32 = kinetic_constant["k32"] * alpha * sum(d_TRn.values())
        ratesList.append(R32)

        R33 = kinetic_constant["k33"] * alpha * sum(d_TRn.values())
        ratesList.append(R33)

        Rnd_values = np.array(list(Rn_d.values()))
        TRn_values = np.array(list(TRn.values()))
        Rnd_TRn = np.array([i*j for i in Rnd_values for j in TRn_values], dtype=np.uint64).sum()
        R41 = kinetic_constant["k41"] * alpha**2 * Rnd_TRn
        ratesList.append(R41)

        R42 = kinetic_constant["k42"] * alpha * sum(RnTRm.values())
        ratesList.append(R42)

        R43 = kinetic_constant["k43"] * alpha * sum(RnTRm.values())
        ratesList.append(R43)

        Rnd_Rmd = np.array([i*j for i in Rnd_values for j in Rnd_values], dtype=np.uint64).sum()
        R51 = kinetic_constant["k51"] * alpha**2 * Rnd_Rmd
        ratesList.append(R51)

        R52 = kinetic_constant["k52"] * alpha**2 * Rnd_Rmd
        ratesList.append(R52)

        if 0 not in list(Rn_d.keys()):
            R53 = 0
        else:
            R53 = kinetic_constant["k53"] * alpha**2 * sum(RnTRm.values()) * Rn_d[0]
        ratesList.append(R53)

        RnTRm_values = np.array(list(RnTRm.values()))
        RnTRm_Rnd = np.array([i*j for i in RnTRm_values for j in Rnd_values], dtype=np.uint64).sum()
        R54 = kinetic_constant["k54"] * alpha**2 * RnTRm_Rnd
        ratesList.append(R54)
    
        ratesArr = np.array(ratesList)
        total_R = ratesArr.sum() # sum all the reaction rate
   

        #Step 2: calculate the probability
        probArr = ratesArr / total_R
        
        #Step 3: choose a type of reaction to happen
        reaction_index = np.arange(12)
        reaction_chosen = np.random.choice(reaction_index, size = 1, 
                                           replace = True, p = probArr)[0]
        r = np.random.uniform(0, 1, 1)
        tau = 1/total_R * math.log(1/r)
        t += tau # Update the system time
        print('t', t)
        print('tau', tau)
        

        # Step 4: choose the exact chain to react
        if reaction_chosen == 0:
            reaction11()
        elif reaction_chosen in [5, 8, 9, 11]:
            reaction_prob_1(reaction_chosen)
        else:
            reaction_prob_0(reaction_chosen)
        
        # collect concentration of the R0
        if 0 not in Rn_d.keys():
            t_and_R0.update({t: 0})
        else:
            t_and_R0.update({t:Rn_d[0]})

        if round(t) % interval == 0:
            times += 1
            name = times * interval
            chainNum_Rnd = list(Rn_d.values())
            chainType_Rnd = list(Rn_d.keys())
            total_chainRnd = sum(chainNum_Rnd)

            file=open(r'Rnd_fraction\system_{1:05d}.txt'.format(settingNum, name),'w+')
            file.write('R_n, n fraction(%)\n')
            ii = 0
            for i in chainNum_Rnd:
                fraction = i/total_chainRnd
                file.write('{0},{1:5.3f}\n'.format(chainType_Rnd[ii], fraction))
                ii += 1
            file.close()

            file=open(r'Rnd_fraction\system_{0:05d}.txt'.format(name),'r')
            dataList=[]
            for line in file.readlines()[1:]:
                dataList.append(line.strip())
                dataArr=np.zeros((len(dataList),2))
            for i in range(len(dataList)):
                dataArr[i,:]=[int(dataList[i].split(',')[0]),float(dataList[i].split(',')[1])]

            TotalPx=probArr[1] + probArr[2] + probArr[5] + probArr[8] + probArr[9]
            fig1 = plt.figure(dpi=300)
            plt.bar(dataArr[:,0], TotalPx*dataArr[:,1], width=0.5, fc='b')
            plt.xlabel('n')
            plt.ylabel('P(n)')
            plt.show()
            fig1.savefig(r'P(n)_n\P(n)_n_{0:05d}.png'.format(name))
            plt.close(fig1)

            avgn = sum(dataArr[:,0]) / len(dataArr[:,0])
            W=0
            for i in range(len(dataArr[:,0])):
                data = dataArr[:,0][i]
                Pn = dataArr[:,1][i]*TotalPx
                W += Pn*(data-avgn)**2
                wDict.update({t:W})
        
        print('I\n', I)
        print('Rn_d\n', Rn_d)
        print('TRn\n', TRn)
        print('d_TRn\n', d_TRn)
        print('RnTRm\n', RnTRm)
        print('Prob\n', probArr)
    
    fig2 = plt.figure(dpi=300)
    plt.plot(list(wDict.keys()), list(wDict.values()), color='green')
    plt.xlabel('t')
    plt.ylabel('W')
    fig2.savefig(r'W_t\W_t_{0:05d}.png'.format(name))
    plt.close(fig2)

    fig3 = plt.figure(dpi=300)
    ax = fig3.add_subplot(111)
    # fig.add_axes(ax)
    ax.plot([i for i in t_and_R0.keys()], [i*alpha for i in t_and_R0.values()], '-', color="r", lw=1)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.xaxis.major.formatter._useMathText = True
    ax.yaxis.major.formatter._useMathText = True
    ax.set_ylabel("concentration (n/L)")
    ax.set_xlabel("t(s)")
    fig3.savefig(r'concentration\Rnd_t_{0:05d}.png'.format(name))
    plt.close(fig3)

    settingNum +=1
    print('---------------------------------------------------------------')