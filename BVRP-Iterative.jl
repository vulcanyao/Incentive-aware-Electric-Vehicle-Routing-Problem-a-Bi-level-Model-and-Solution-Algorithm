# https://atoptima.github.io/Coluna.jl/latest/start/start/

using JuMP, BlockDecomposition, Coluna, CPLEX, MAT,CSV,DataFrames,Gurobi;


# ite = 0;
# node = 29;


function BVRP_BDD(ite, node,vehicle,deltabar,gamma,wt,Index_BCP)
    coluna = optimizer_with_attributes(
    Coluna.Optimizer,
    "params" => Coluna.Params(
        solver=Coluna.Algorithm.TreeSearchAlgorithm(
            # conqueralg =  Coluna.ColCutGenConquer(colgen = Coluna.Algorithm.ColumnGeneration(opt_rtol=0.01)),
            conqueralg = Coluna.ColCutGenConquer(stages=[Coluna.Algorithm.ColumnGeneration(opt_rtol=0.1)])
        #     maxnumnodes::Int = 100000,
        #    opennodeslimit = 10
        ) # default BCP
    ),
    "default_optimizer" => CPLEX.Optimizer # CPLEX for the master & the subproblems
);
    # Input 

    # vehicle = 5
    # ite = 1;
    # node = 3;
    # wt = 1    
    # Index for whether choosing BCP or not (MIP) 
    # Index_BCP = 1
    
    
    
    # Input data preparation 
    data = matread(string("/Users/vulcanyao/OneDrive - 南方科技大学/PDF Expert/Research Yao/code_VRPPD/[J3]Bilevel EVRP with time flexibility/J3-V1/DataMap/Vehicles_",ite,"/RealMap_",node,"_test.mat"))
    # data = matread(string("/DataMap/Vehicles_",ite,"/RealMap_",node,"_test.mat"))
    
    
    T = data["T"];
    M = data["M"];
    t = data["t"];
    c = data["c"]
    Node = Integer(data["N"]);
    chi = 0.01
    
    
    
    # Define the mathematical formulation 
    
    
    # K = 1:vehicle;
    # @axis(N,1:Node);
    
    
    if Index_BCP == 1
        # BVRP = BlockModel(coluna)
        BVRP  = BlockModel(colunay,direct_model = true)s
        @axis(K,1:vehicle);
    else
        # BVRP = Model(CPLEX.Optimizer)
        BVRP = Model(Gurobi.Optimizer)
        K = 1:vehicle;
    end
    
    N = 1:Node;
    @variable(BVRP,x[k in K, i in N, j in N] ,Bin);
    @variable(BVRP,omega1[i in N] ,Bin);
    @variable(BVRP,omega2[i in N] ,Bin);
    @variable(BVRP,omega3[i in N] ,Bin);
    @variable(BVRP,omega4[i in N] ,Bin);
    @variable(BVRP,0 <= eta[   i in N]);
    @variable(BVRP, 0 <= epsilon[  i in N]);
    @variable(BVRP, 0 <= q[  i in N]);
    @variable(BVRP,0<= u[  i in N]);
    @variable(BVRP,0<= v[  i in N]);
    @variable(BVRP,0<= zeta1[  i in N]);
    @variable(BVRP,0<= zeta2[  i in N]);
    @variable(BVRP,0<= delta[  i in N]<= deltabar );
    @variable(BVRP, tau_var[  i in N]);
    
    
    # Cons: pre-processing flow
    @constraint(BVRP, a1[i=2:Node ], sum(x[ k,i,j] for j in N, k in K) <= 1 );
    @constraint(BVRP, a2[i in N, k in K], x[k,i,i] == 0 );
    @constraint(BVRP, a3[i in N,  k in K],  x[k, Node,i] == 0 );
    @constraint(BVRP, a4[i in N,  k in K],  x[ k,i,1] == 0 );
    
    # # # Cons: flow conservation
    #  Start depot and end depot constraints result in the infeasibility.
    @constraint(BVRP, [ k in K],  sum(x[k,1,i] for i in N)  - sum(x[k,i,1] for i in N) == 1  );
    @constraint(BVRP, [ k in K],  - sum(x[k,i,Node] for i in N) == -1 );
    @constraint(BVRP, [i=2:Node-1, k in K],  sum(x[k,i,j] for j in N) - sum(x[k,j,i] for j in N) == 0  );
    # time constraint 
    @constraint(BVRP,[i in N],  tau_var[i] <= t[i] + delta[i]  );
    @constraint(BVRP,[i in N],  t[i] <= tau_var[i]  );
    @constraint(BVRP, [i in N, j in N, k in K], tau_var[j] >= tau_var[i] + T[i,j] - M * (1 - x[k,i, j])   );
    # Linearization
    @constraint(BVRP,[j in N],  eta[j] >= epsilon[j] + u[j] * deltabar -zeta1[j]*chi -zeta2[j]*chi  - M * (1 - sum(x[k,i, j] for i in N for k in K))     );
    #   KKT opt. condition
    @constraint(BVRP,[i=2:length(N)-2],  -gamma * zeta1[i] +gamma * zeta2[i]  - q[i] + u[i] - v[i] == 0  );
    @constraint(BVRP,[i=2:length(N)-2],  1 - zeta1[i] - zeta2[i]   == 0  );
    
    @constraint(BVRP,[i=2:length(N)-2],  0 <= epsilon[i] + gamma * delta[i]  - chi   );
    @constraint(BVRP,[i=2:length(N)-2],  epsilon[i] + gamma * delta[i]  - chi <= M * omega1[i]  );
    @constraint(BVRP,[i=2:length(N)-2],  zeta1[i] <= M *(1-omega1[i])  )
    @constraint(BVRP,[i=2:length(N)-2],  0 <= epsilon[i] - gamma * delta[i]  + chi  );
    
    @constraint(BVRP,[i=2:length(N)-2],  zeta2[i] <= M *(1-omega2[i])  );
    @constraint(BVRP,[i=2:length(N)-2],  epsilon[i] - gamma * delta[i]  + chi <= M *(1-omega2[i])  );
    @constraint(BVRP,[i=2:length(N)-2],  u[i] <= M *(1-omega3[i])  );
    @constraint(BVRP,[i=2:length(N)-2],  deltabar -  delta[i]  <= M *(1-omega3[i])  );
    @constraint(BVRP,[i=2:length(N)-2],  v[i] <= M *(1-omega4[i])  );
    @constraint(BVRP,[i=2:length(N)-2],  delta[i]  <= M *(1-omega4[i])  );
    
    
    @objective(BVRP,Min, sum(wt * T[i, j] * x[k,i, j] + c[i] * x[k,i, j] for i in N for j in N for k in K)+ sum( eta[i] for i = 2:length(N)-1)     )
    
    
    
    
    if Index_BCP == 1
        # dantzig_wolfe_decomposition along the K axis
        @dantzig_wolfe_decomposition(BVRP, decomposition, K)
        # Define sub-problem
        master = getmaster(decomposition)
        subproblems = getsubproblems(decomposition)
        specify!.(subproblems, lower_multiplicity=1, upper_multiplicity=1)
        getsubproblems(decomposition)
    end
    optimize!(BVRP)


    # Results analysis
    solution_summary(BVRP)
    value(sum(c[i] * sum(x[k,i, j]  for j in N, k in K) for i in N))
    value(sum(T[i,j] * sum(x[k,i, j]  for k in K) for i in N, j in N))
    value(sum(eta[i]  for i in N ))
    for j in N 
        for k in K 
            for i in N
                if value(x[k,i, j]) == 1 
                    println("vehicle=",k,"-node=",i,"-ndoe=",j )
                end
                
            end    
        end
    end
    sol_time = MOI.get(BVRP, MOI.SolveTime()) 
    opt_val = objective_value(BVRP)
    
    
    
    return sol_time, opt_val
    
    
end


    
    
function main()
    index = 1 
    vehicle = 5 
    deltabar = 1.5
    gamma = 1.5
    Index_BCP = 1
    wt = .01 

    for node = 11:2:21;
        for ite = 1:1;

            time, obj_val = BVRP_BDD(ite, node,vehicle,deltabar,gamma,wt,Index_BCP)
            fields = DataFrame( run_no = ite, vehicle_no = vehicle, node_no=node, Runtime =time, obj_value=obj_val )
            fi="c:\\temp\\Boxplot_BVRP_BCP.csv"
            fi="Boxplot_BVRP_BCP.csv"
            CSV.write(fi, fields, writeheader = (index==1), append = true)  
            index = index+1
            println("The $index iteration done!!!!")
        end
    end
end



main()
    
    