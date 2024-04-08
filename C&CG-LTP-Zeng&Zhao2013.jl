using JuMP
using Gurobi
f = [400, 414, 326]
a = [18, 25, 20]
C = [[22 33 24];[33 23 30];[20 25 27]]
D = [206 + 40, 274 + 40, 220 + 40]
dl = [206, 274, 220]
du = [40, 40, 40]
k = 0

MP = Model(Gurobi.Optimizer)
set_silent(MP)
@variables(MP, begin
    x[1:3, 1:3] >= 0
    y[1:length(f)], Bin
    z[1:length(a)] >= 0
    d[1:3] >= 0
    eta >= 0
    0 <= g[1:3] <= 1
end)

@objective(MP, Min, f' * y + a' * z + eta)
@constraints(MP, begin
    [i in 1:3], z[i] <= 800 * y[i]
    sum(z) >= 772
    [i in 1:3], sum(x[i, :]) <= z[i]
    [j in 1:3], sum(x[:, j]) >= d[j]
    eta >= sum(C .* x)
    [i in 1:3], d[i] == dl[i] + du[i] * g[i]
   sum(g[i] for i in 1:3) <= 1.8
   sum(g[i] for i in 1:2) <= 1.2
end)

optimize!(MP)
LB = objective_value(MP)
println("Lower Bound: $LB")
println("d=",value.(d))
println("eta=",value.(eta))
println("z=",value.(z))

#Recourse:

z_star = value.(z)
M = 10000
SP = Model(Gurobi.Optimizer)
set_silent(SP)
@variables(SP, begin
    x_sub[1:3, 1:3] >= 0;
    d_sub[1:3] >= 0;
    0<=g_sub[1:3]<=1;
    pi[1:6] >= 0;
    v[1:6], Bin;
    w[1:3, 1:3], Bin;
end)
@objective(SP, Max, sum(C .* x_sub));
@constraints(SP, begin
    SP_Cons_1[i in 1:3], sum(x_sub[i, :]) <= z_star[i];
    [j in 1:3], sum(x_sub[:, j]) >= d_sub[j];
    [i in 1:3, j in 1:3], -pi[i] + pi[j + 3] <= C[i,j];
        
    SP_SLACK_CONS_1[i in 1:3], z_star[i] - sum(x_sub[i, :]) <= M * (1 - v[i]);
    [j in 1:3], sum(x_sub[:, j]) - d_sub[j] <= M * (1 - v[j + 3]);
    [i in 1:6], pi[i] <= M * v[i];
        
    [i in 1:3, j in 1:3], C[i,j] + pi[i] - pi[j + 3] <= M * (1 - w[i, j]);
    [i in 1:3, j in 1:3], x_sub[i, j] <= M * w[i, j];
        
    [i in 1:3], d_sub[i] == dl[i] + du[i] * g_sub[i];
   sum(g_sub[i] for i in 1:3)<=1.8;
        sum(g_sub[i] for i in 1:2)<=1.2;
end)
optimize!(SP)
println("objective:", objective_value(SP))
println("pi:", value.(pi))
println("g_sub:", value.(g_sub))


#CCG main:

UB = Inf
while abs(UB - LB) > 1e-5
    optimize!(SP)
    SP_objval = objective_value(SP)
    println("SP_objval: $SP_objval")
    UB = LB - value(eta) + SP_objval
    println("Upper Bound: $UB")
    d_sub_star = value.(d_sub)
    println("d_sub=",value.(d_sub_star))
    unregister(MP,:x_new)
    @variable(MP, x_new[1:3, 1:3] >= 0)
    k = k + 1
    @constraint(MP, [i in 1:3], sum(x_new[i, :]) <= z[i])
    @constraint(MP, [j in 1:3], sum(x_new[:, j]) >= d_sub_star[j])
    @constraint(MP, eta >= sum(C .* x_new))
    optimize!(MP)
    LB = max(LB, objective_value(MP))
    println("Lower Bound: $LB")
    z_star = value.(z)
    v_star = value.(v)
    unregister(SP,:x_new)
    delete.(SP, SP_Cons_1)
    unregister(SP,:SP_Cons_1)
    delete.(SP,SP_SLACK_CONS_1)
    unregister(SP,:SP_SLACK_CONS_1)
    @constraints(SP, begin
    SP_Cons_1[i in 1:3], sum(x_sub[i, j] for j in 1:3) <= z_star[i]
    SP_SLACK_CONS_1[i in 1:3], z_star[i] - sum(x_sub[i, j] for j in 1:3) <= M * (1 - v[i])
    end)
end
println("Iteration finished! We found the optimal solution!")
println("Final Objective:{0}",UB)
