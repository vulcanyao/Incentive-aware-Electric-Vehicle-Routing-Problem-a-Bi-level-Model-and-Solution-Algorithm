# branch and cut applied in EVRP
import math
import random
from gurobipy import *
from scipy.io import loadmat
import numpy as np
import sys
from csv import writer
import time


def EPI(model, where):
    if where == GRB.Callback.MIPNODE:  #
        nodecnt = int(model.cbGet(GRB.callback.MIPNODE_NODCNT))
        nodecnt = nodecnt / model._numcheckedNode
        status = model.cbGet(GRB.Callback.MIPNODE_STATUS)  #
        # print(nodecnt)

        if status == GRB.OPTIMAL:
            if (nodecnt.is_integer()) & (model._numEPI <= model._numEPImax):  # model._numEPImax
                # if nodecnt <= 5e3:

                sol = model.cbGetNodeRel(model._vars)  # A relaxed sol. obtained at the  explored node
                SubFun = model._Gfun
                N = model._coeff[0]
                K = model._coeff[1]
                cut = np.zeros((K, N ** 2))
                Z = np.zeros(K)
                label_cut = np.zeros((K, N ** 2))
                cut_reindex = np.zeros((K, N ** 2))
                pi = np.zeros((K, N ** 2))

                for i in range(K):
                    # for j in range(N ** 2):
                    cut[i, :] = sol[i * (N ** 2): (i + 1) * (N ** 2)]

                for i in range(K):
                    Z[i] = sol[i + K * (N ** 2)]

                #  Indeices of relaxed sol in descending order (max-min).
                for i in range(K):
                    label_cut[i, :] = np.argsort(-cut[i, :])
                # label_cut.astype(int)

                for i in range(K):
                    for j in range(N ** 2):
                        pi[i, j] = SubFun[int(label_cut[i, j])]
                        cut_reindex[i, j] = cut[i, int(label_cut[i, j])]

                # EPIs construction

                # for i in range(K):
                #     if np.dot(pi[i, :], cut[i, :]) > Z[i]:
                #             model.cbCut( sum((pi[i, j] * model._vars[j + i * (N ** 2)]) for j in range(N ** 2)) <= model._vars[ i + K * (N ** 2)])
                # model._numEPI += 1

                if sum(np.dot(pi[i, :], cut_reindex[i, :]) for i in range(K)) > sum((Z[i]) for i in range(K)):
                    model.cbCut(sum(
                        sum((pi[i, j] * model._vars[int(label_cut[i, j]) + i * (N ** 2)]) for i in range(K)) for j in
                        range(N ** 2)) <= sum((model._vars[k + K * (N ** 2)]) for k in range(K)))
                    model._numEPI += 1


# print(pi)


# Initial paras


def BVRP(map, tlimit, cut, deltabar, gamma, node, vehicle, maxEPI, run):
    ###### Main function
    data = loadmat('DataMap/Vehicles_{}/RealMap_{}_test.mat'.format(run, node))

    Emax = data['Emax'][0][0]
    rmax = Emax
    # e=[]
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
    m.setParam('Heuristics', 0.00)
    # m.setParam('Threads',1)

    ini_time = time.time()

    x = m.addVars(K, N, N, vtype=GRB.BINARY, name="x")  # ****************
    Z = m.addVars(K, vtype=GRB.CONTINUOUS, name="Z")
    omega1 = m.addVars(N, vtype=GRB.BINARY, name="omega1")
    omega2 = m.addVars(N, vtype=GRB.BINARY, name="omega2")
    eta = m.addVars(K, N, N, vtype=GRB.CONTINUOUS, name="eta")  # ****************
    epsilon = m.addVars(K, N, N, vtype=GRB.CONTINUOUS, name="epsilon")  # ****************
    biterm = m.addVars(N, vtype=GRB.CONTINUOUS, name="biterm")
    r = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="r")
    E = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="E")
    E0 = m.addVars(N, K, vtype=GRB.CONTINUOUS, name="E0")
    u = m.addVars(N, vtype=GRB.CONTINUOUS, name="u")
    v = m.addVars(N, vtype=GRB.CONTINUOUS, name="v")
    q = m.addVars(N, vtype=GRB.CONTINUOUS, name="q")
    delta = m.addVars(N, vtype=GRB.CONTINUOUS, name="delta")
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
    # m.addConstrs(r[j, k] <= rmax for j in range(N) for k in range(K))
    # m.addConstrs(0 <= r[j, k] for j in range(N) for k in range(K))
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
    m.addConstrs(tau[j] - t[j]- delta[j] <= 0 for j in range(N))
    m.addConstrs(
        tau[j] >= tau[i] + T[i, j] - M * (1 - x[k, i, j]) for j in range(N) for i in range(N) for k in range(K))

    m.addConstrs(
        eta[k, i, j] >= g[i] * r[i, k] - M * (1 - x[k, i, j]) for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(
        epsilon[k, i, j] >= p[i] * r[i, k] - M * (1 - x[k, i, j]) for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(eta[k, i, j] >= 0 for i in range(N) for j in range(N) for k in range(K))
    m.addConstrs(epsilon[k, i, j] >= 0 for i in range(N) for j in range(N) for k in range(K))

    #  KKT (in)equation set
    m.addConstrs(gamma - q[j] + u[j] - v[j] == 0 for j in range(N) if 0 < j < N - 1)
    m.addConstrs(0 <= deltabar - delta[j] for j in range(N))
    m.addConstrs(deltabar - delta[j] <= M * omega1[j] for j in range(N) if 0 < j < N - 1)
    m.addConstrs(0 <= u[j] for j in range(N) if 0 < j < N - 1)
    m.addConstrs(u[j] <= M * (1 - omega1[j]) for j in range(N) if 0 < j < N - 1)
    m.addConstrs(0 <= delta[j] for j in range(N) if 0 < j < N - 1)
    m.addConstrs(delta[j] <= M * omega2[j] for j in range(N) if 0 < j < N - 1)
    m.addConstrs(0 <= v[j] for j in range(N) if 0 < j < N - 1)
    m.addConstrs(v[j] <= M * (1 - omega2[j]) for j in range(N) if 0 < j < N - 1)

    m.addConstrs(delta[j] <= deltabar for j in range(N) if 0 < j < N - 1)

    # Linearization of bi-linear
    m.addConstrs(biterm[j] >= 0 for j in range(N) if 0 < j < N - 1)
    m.addConstrs(biterm[j] >= gamma * delta[j] + deltabar * u[j] - M * (1 - x.sum('*', '*', j)) for j in range(N) if
                 0 < j < N - 1)

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
    # zz = 0
    for k in range(K):
        for i in range(N):
            for j in range(N):
                qq += eta[k, i, j] + epsilon[k, i, j]

    # m.addConstrs(    for i in range(N) for j in range(N) for k in range(K) )

    m.setObjective(qq + quicksum(biterm) + quicksum(Z))  # quicksum(Z)
    m.update()
    m._vars = m.getVars()

    # Submodular set function (quadratic form)
    GFun = []
    for i in range(N):
        for j in range(N):
            GFun.append(T[i, j] + c[i])
    m._Gfun = GFun

    m._coeff = [N, K]
    m._numEPI = 0;
    m._numEPImax = maxEPI;
    m._numcheckedNode = 1;
    m._num_inte_sol = 0;

    # m._numVehicle = K;

    solution = m._vars

    m.Params.PreCrush = 1
    A = m.getA()
    # f=m.getCoeff()
    # print("(A) Matrix coeff. %s" % A)
    # print("obj coeff.=     %s" %f)

    if cut == 0:
        print('\n' + '#' * 20)
        print('Optimize: sol_map{}_N{}_K{}_cut{}_run{}.txt'.format(map, N, K, cut, run))
        print('\n' + '#' * 20)
        m.optimize()
    else:
        print('\n' + '#' * 20)
        print('Optimize: sol_map{}_N{}_K{}_cut{}_run{}.txt'.format(map, N, K, cut, run))
        print('\n' + '#' * 20)
        m.optimize(EPI)

    " Results analysis"
    print('Added EPIs = %s ' % m._numEPI)
    # m.printAttr('x')
    # print(m.getAttr())
    # print('Runtime = %.2f sec' %(m.Runtime) )

    varInfo = [(v.varName, v.X) for v in m.getVars() if v.X > 0]
    print('MIPGAP = %s ' % m.getAttr('MIPGap'))
    obj = m.getObjective()

    print(m.getVarByName('biterm[%s]' % 4).x)

    biterm_sol = []
    for i in range(N):
        if 0 < i < N - 1:
            biterm_sol.append(m.getVarByName('biterm[%s]' % i).x)
    print(biterm_sol)

    biterm_val = 0
    if np.sum(biterm_sol) !=0:
        biterm_no = 0
        for j in range(node - 2):
            if biterm_sol[j] != 0:
                biterm_no += 1
        biterm_val = np.sum(biterm_sol) / biterm_no

    charging_cost = []
    for k in range(K):
        for i in range(N):
            for j in range(N):
                charging_cost.append(m.getVarByName('epsilon[%s,%s,%s]' % (k, i, j)).x)
    print(charging_cost)

    travel_cost = []
    for k in range(K):
        for i in range(N):
            for j in range(N):
                travel_cost.append(m.getVarByName('x[%s,%s,%s]' % (k, i, j)).x * T[i, j])
    print(travel_cost)

    if cut == 0:
        sys.stdout = open(
            'solutionTest/BVRP_gamma_{}_deltabar_{}_maxEPI_{}/sol_map{}_N{}_K{}_run{}.txt'.format(gamma, deltabar,
                                                                                                  maxEPI, map, N, K,
                                                                                                  run), 'w')

    else:
        sys.stdout = open(
            'solutionTest/BVRP_gamma_{}_deltabar_{}_maxEPI_{}/sol_EPI_map{}_N{}_K{}_run{}.txt'.format(gamma, deltabar,
                                                                                                      maxEPI, map, N, K,
                                                                                                      run), 'w')

    print(varInfo)
    # m.printAttr('x')
    print('\n' + '#' * 20)
    print('Runtime = %.2f sec' % (m.Runtime))
    print('Optimal objective value = %.2f ' % (obj.getValue()))
    print('Explored nodes = %s ' % m.getAttr('NodeCount'))
    print('Total vars = %s ' % m.getAttr('NumVars'))

    #  Depend on cut added or not
    if cut == 1:
        print('Added EPIs = %s ' % m._numEPI)

    print('MIPGAP = %.4f  ' % m.getAttr('MIPGap'))
    print('Added EPIs = %s ' % m._numEPI)
    print('Biterms = %s ' % biterm_sol)

    "Export as csv."

    return m, biterm_val, charging_cost, travel_cost


#  Export csv.file (solution ) and txt.file simultaneously (solution details)

if __name__ == '__main__':
    gamma0 = np.array([1.5, 0.01])
    deltabar0 = np.array([1.5, 1.5])
    tlimit = 3600 * 2
    map = 100
    vehicle = 3
    maxEPI0 = np.array([1, 500])
    "Export into csv file "

    cut = 0
    index = 6
    runmax = 50
    nodeI = 13

    # for k in range(run):
    for gam in range(1):
        # if gam >0 :
        gamma = gamma0[gam]
        for del1 in range(1):
            deltabar = deltabar0[del1]
            for mEPI in range(1):
                maxEPI = maxEPI0[mEPI];
                for i in range(index):
                    for run in range(runmax):
                        run = run + 1

                        # node = nodeI + 2 * i
                        # m, biterm_val, charging_cost, travel_cost = BVRP(map, tlimit, 1, deltabar, gamma, node, vehicle,
                        #                                                  maxEPI, run)
                        # obj = m.getObjective()
                        # with open(
                        #         r'solutionTest/BVRP_gamma_{}_deltabar_{}_maxEPI_{}/Sol_EPI.csv'.format((gamma),
                        #                                                                                (deltabar),
                        #                                                                                maxEPI),
                        #         mode='a') as f:
                        #     Sol = writer(f, delimiter='\t', lineterminator='\n')
                        #     Sol.writerow(
                        #         ['sol_map{}_N{}_K{}_run{}'.format(map, node, vehicle, run), '%.2f' % (m.Runtime),
                        #          '%.2f ' % (obj.getValue()),
                        #          '%i ' % m.getAttr('NodeCount'),
                        #          '%.4f' % m.getAttr('MIPGap'), '%s ' % m._numEPI, '%s ' % m._num_inte_sol,
                        #          '%s ' % int(m._numEPImax), '%s ' % int(m._numcheckedNode),
                        #          '%.4f ' % (biterm_val), '%.4f' % np.sum(charging_cost),
                        #          '%.4f' % np.sum(travel_cost)])
                        #     f.close()

                        node = nodeI + 2 * i
                        m, biterm_val, charging_cost, travel_cost = BVRP(map, tlimit, 0, deltabar, gamma, node, vehicle,
                                                                         maxEPI, run)
                        obj = m.getObjective()




                        with open(
                                r'solutionTest/BVRP_gamma_{}_deltabar_{}_maxEPI_{}/Sol.csv'.format((gamma),
                                                                                                   (deltabar),
                                                                                                   maxEPI),
                                mode='a') as f:
                            Sol = writer(f, delimiter='\t', lineterminator='\n')
                            Sol.writerow(
                                ['sol_map{}_N{}_K{}_run{}'.format(map, node, vehicle, run), '%.2f' % (m.Runtime),
                                 '%.2f ' % (obj.getValue()),
                                 '%i ' % m.getAttr('NodeCount'),
                                 '%.4f' % m.getAttr('MIPGap'), '%s ' % m._numEPI, '%.6f ' % (biterm_val),
                                 '%.4f' % np.sum(charging_cost),
                                 '%.4f' % np.sum(travel_cost)])
                            f.close()


