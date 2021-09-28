# branch and cut applied in EVRP
import math
import random
from gurobipy import *
from scipy.io import loadmat
import numpy as np
import sys
import csv
import time

# tlimit = 200  # 3600*2
# map = 100
# vehicle = 3
# "Export into csv file "
# index = 2
# run = 12
# node =13


###### Main function
def BVRP(map, tlimit, node, vehicle,run):
    data = loadmat('DataMap/Vehicles_{}/RealMap_{}_test.mat'.format(run, node))

    Emax = data['Emax'][0][0]
    rmax = Emax
    e = data['e']
    T = 1e0 * data['T']

    M = int(data['M'])
    p = data['p'][0]
    g = data['g'][0]
    t = 1e0 * data['t'][0]
    c = data['c'][0]
    K = int(vehicle)
    N = int(data['N'])

    # Create variables
    "Gurobipy model formulation"
    m = Model()
    m.setParam('TimeLimit', tlimit)
    # m.setParam('Heuristics', 0.00)
    # m.setParam('Threads', 1)

    ini_time = time.time()

    x = m.addVars(K, N, N, vtype=GRB.BINARY, name="x")  # ****************
    Z = m.addVars(K, vtype=GRB.CONTINUOUS, name="Z")
    omega1 = m.addVars(N, vtype=GRB.BINARY, name="omega1")
    omega2 = m.addVars(N, vtype=GRB.BINARY, name="omega2")
    eta = m.addVars(K, N, N, vtype=GRB.CONTINUOUS, name="eta")  # ****************
    epsilon = m.addVars(K, N, N, vtype=GRB.CONTINUOUS, name="epsilon")  # ****************
    # biterm = m.addVars(N, vtype=GRB.CONTINUOUS, name="biterm")
    r = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="r")
    E = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="E")
    E0 = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="E0")
    tau = m.addVars(N, vtype=GRB.CONTINUOUS, name="tau")

    m.addConstrs(x.sum('*', i, '*') <= 1 for i in range(N) if i != 0)
    m.addConstrs(E[0, k] == 0.8 * Emax for k in range(K))

    m.addConstrs(x[k, i, i] == 0 for i in range(N) for k in range(K))
    m.addConstrs(x[k, N - 1, i] == 0 for i in range(N) for k in range(K))
    m.addConstrs(x[k, i, 0] == 0 for i in range(N) for k in range(K))

    # Flow conservation
    m.addConstrs(x.sum(k, 0, '*') - x.sum(k, '*', 0) == 1 for k in range(K))
    m.addConstrs(x.sum(k, N - 1, '*') - x.sum(k, '*', N - 1) == -1 for k in range(K))
    m.addConstrs(x.sum(k, i, '*') - x.sum(k, '*', i) == 0 for i in range(N - 1) if i != 0 for k in range(K))

    # SOC bound
    m.addConstrs(E[j, k] <= Emax for j in range(N) for k in range(K))
    m.addConstrs(0 <= E[j, k] for j in range(N) for k in range(K))
    m.addConstrs(r[j, k] <= rmax for j in range(N) for k in range(K))
    m.addConstrs(0 <= r[j, k] for j in range(N) for k in range(K))
    m.addConstrs(E[j, k] == E0[j, k] for j in range(N) for k in range(K))

    #  SOC cons
    m.addConstrs(
        E[i, k] - e[i, j] * x[k, i, j] + r[i, k] - E0[j, k] <= M * (1 - x[k, i, j]) for i in range(N) for j in range(N)
        for k in range(K))
    m.addConstrs(
        -M * (1 - x[k, i, j]) <= E[i, k] - e[i, j] * x[k, i, j] + r[i, k] - E0[j, k] for i in range(N) for j in range(N)
        for k in range(K))
    # Time cons
    m.addConstrs(t[j] - tau[j] <= 0 for j in range(N))
    m.addConstrs(tau[j] - t[j] <= 0 for j in range(N))
    m.addConstrs(
        tau[j] >= tau[i] + T[i, j] - M * (1 - x[k, i, j]) for j in range(N) for i in range(N) for k in range(K))

    # self.constraints.c8a = m.addConstrs(t[j] - tau[j] <= 0 for j in range(N))
    # self.constraints.c8b = m.addConstrs(tau[j] - t[j] <= 0 for j in range(N))
    #
    # self.constraints.c9 = m.addConstrs(
    #     tau[j] >= tau[i] + T[i, j] - M * (1 - x[k, i, j]) for j in range(N) for i in range(N) for k in range(K))

    # Linearized obj. terms
    m.addConstrs(
        eta[k, i, j] >= g[i] * r[i, k] - M * (1 - x[k, i, j]) for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(
        epsilon[k, i, j] >= p[i] * r[i, k] - M * (1 - x[k, i, j]) for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(eta[k, i, j] >= 0 for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(epsilon[k, i, j] >= 0 for i in range(N) for j in range(N) for k in range(K))

    end_time = time.time()
    Time = end_time - ini_time
    print("Running time %s" % Time)
    #  Z distance + cost
    for k in range(K):
        zz = 0
        for i in range(N):
            for j in range(N):
                zz += T[i, j] * x[k, i, j] + c[i] * x[k, i, j]
        m.addConstr(Z[k] >= zz)

    qq = 0
    for k in range(K):
        for i in range(N):
            for j in range(N):
                qq += eta[k, i, j] + epsilon[k, i, j]

    m.setObjective(qq + quicksum(Z))  # quicksum(Z)
    m.update()
    m._vars = m.getVars()
    m.Params.PreCrush = 1

    print('\n' + '#' * 20)
    print('Optimize: sol_map{}_N{}_K{}.txt'.format(map, N, K))
    print('\n' + '#' * 20)
    m.optimize()

    " Results analysis"

    print(m.getVars())
    varInfo = [(v.varName, v.X) for v in m.getVars() if v.X > 0]
    delta_val = m.getAttr('x', x)  # [(v.varName, v.X) for v in [x]]
    print(delta_val)
    qq = []
    for i in range(N):
        for j in range(N):
            for k in range(vehicle):
                qq.append(delta_val[k, i, j])  # delta_val[i, j, k]
    print(qq)

    charging_cost = []
    for k in range(K):
        for i in range(N):
            for j in range(N):
                 charging_cost.append(m.getVarByName('epsilon[%s,%s,%s]' % (k,i,j) ).x)
    print(charging_cost)


    travel_cost=[]
    for k in range(K):
        for i in range(N):
            for j in range(N):
                travel_cost.append( m.getVarByName('x[%s,%s,%s]' %(k,i,j)).x  *  T[i,j] )
    print(travel_cost)


    print('MIPGAP = %s ' % m.getAttr('MIPGap'))
    obj = m.getObjective()

    sys.stdout = open('solutionTest/SL/sol_map{}_N{}_K{}_run_{}.txt'.format(map, N, K, run), 'w')
    print(varInfo)
    # m.printAttr('x')
    print('\n' + '#' * 20)
    print('Runtime = %.2f sec' % (m.Runtime))
    print('Optimal objective value = %.2f ' % (obj.getValue()))
    print('Explored nodes = %s ' % m.getAttr('NodeCount'))
    print('Total vars = %s ' % m.getAttr('NumVars'))

    #  Depend on cut added or not
    # if cut == 1:
    #     print('Added EPIs = %s ' % m._numEPI)

    print('MIPGAP = %.4f  ' % m.getAttr('MIPGap'))
    # print('Added EPIs = %s ' % m._numEPI)
    # print('Biterms = %s ' % biterm_sol)


    # "Export as csv."
    return m,charging_cost,travel_cost



#  Export csv.file (solution ) and txt.file simultaneously (solution details)

if __name__ == '__main__':
    tlimit =  3600*2
    map = 100
    vehicle = 3
    "Export into csv file "

    index = 6
    runmax = 50

    # node = 13
    # m = BVRP(map, tlimit, node, vehicle)

    for i in range(index):
        for run in range(runmax):
            run = run + 1
            node = 13 + 2 * i
            m, charging_cost, travel_cost = BVRP(map, tlimit, node, vehicle,run)
            obj = m.getObjective()
            with open(r'solutionTest/SL/Sol.csv'.format(run), mode='a') as f:
                Sol = csv.writer(f, delimiter='\t', lineterminator='\n')
                Sol.writerow(
                    ['sol_map{}_N{}_K{}_run{}'.format(map, node, vehicle,run), '%.2f' % (m.Runtime), '%.2f ' % (obj.getValue()),
                     '%i ' % m.getAttr('NodeCount'), '%.4f' % m.getAttr('MIPGap'), '%.4f' % np.sum(charging_cost),
                     '%.4f' % np.sum(travel_cost)])

