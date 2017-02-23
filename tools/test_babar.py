import matplotlib.pylab as plt
import numpy as np

import babar_tools as bbr 

import sys

################################################################################
def inv_mass(p4s):

    tot_p4 = np.array([0.0, 0.0, 0.0, 0.0])
    for p4 in p4s:
        #print(tot_p4)
        tot_p4 += p4
        #print(tot_p4)

    m2 = tot_p4[0]*tot_p4[0] - \
           (tot_p4[1]*tot_p4[1] + \
            tot_p4[2]*tot_p4[2] + \
            tot_p4[3]*tot_p4[3]) 

    m = None
    if m2>=0:
        m = np.sqrt(m2)
    else:
        m = -np.sqrt(np.abs(m2))

    #print(m)
    return m

################################################################################

collisions = bbr.get_collisions(open(sys.argv[1]))

print(len(collisions))

masses = []
for count,collision in enumerate(collisions):

    if count%10000==0:
        print(count)

    pions,kaons,protons,muons,electrons,photons = collision

    '''
    npions = len(pions)
    for i in range(0,npions-1):
        p4i = np.array(pions[i][0:4])
        for j in range(i+1,npions):
            p4j = np.array(pions[j][0:4])
            m = inv_mass([p4i,p4j])
            masses.append(m)
    '''

    '''
    nmuons = len(muons)
    for i in range(0,nmuons-1):
        p4i = np.array(muons[i][0:4])
        for j in range(i+1,nmuons):
            p4j = np.array(muons[j][0:4])
            m = inv_mass([p4i,p4j])
            masses.append(m)
    '''
    '''
    nkaons = len(kaons)
    for i in range(0,nkaons-1):
        p4i = np.array(kaons[i][0:4])
        for j in range(i+1,nkaons):
            p4j = np.array(kaons[j][0:4])
            m = inv_mass([p4i,p4j])
            masses.append(m)
    '''

    '''
    # Look for D-->Kpipi
    nkaons = len(kaons)
    npions = len(pions)
    for k in range(0,nkaons):
        p4k = np.array(kaons[k][0:4])
        qk = kaons[k][4]
        for i in range(0,npions-1):
            p4i = np.array(pions[i][0:4])
            qi = pions[i][4]
            for j in range(i+1,npions):
                p4j = np.array(pions[j][0:4])
                qj = pions[j][4]
                if (qk==-qi and qi==qj):
                    m = inv_mass([p4i,p4j,p4k])
                    masses.append(m)
    '''
    nphotons = len(photons)
    #print(nphotons)
    for i in range(0,nphotons-1):
        p4i = np.array(photons[i][0:4])
        for j in range(i+1,nphotons):
            p4j = np.array(photons[j][0:4])
            m = inv_mass([p4i,p4j])
            masses.append(m)
    

print(masses)
plt.figure()
plt.hist(masses,bins=100,range=(0.0,3.0))
plt.show()
