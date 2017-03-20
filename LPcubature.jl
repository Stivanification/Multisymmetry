include("assemble_constraint_matrix.jl")

using SymPy
using JuMP
using Gurobi
using Combinatorics
using FastGaussQuadrature

function LPcube(A,moments)
  M = Model(solver=GurobiSolver())
  @variable(M, x[1:size(A,2)] >= 0)
  @constraint(M, A*x.== moments)
  @objective(M, Min, sum(0*x))
  status = JuMP.solve(M)

  weights = getvalue(x)
end

#Gudmund ist der Beste.

function return_cubature_formula(n,m,maxdegree)
  nodes,weights = try
    nodes,weights = readdlm("formulas/($n,$m,$maxdegree)-nodes.txt"),readdlm("formulas/($n,$m,$maxdegree)-weights.txt")
  catch
    println("computing formula...")
    A,moments,nodes = assemble_constraints(n,m,maxdegree)
	  tic()
    weights = LPcube(A,moments)
    t = toq()
    println("Simplex solved in $t seconds.")

    #restrict to nodes with nonzero weights
    ind = find(weights)
    weights = weights[ind]
    nodes = nodes[ind,:]

    writedlm("formulas/($n,$m,$maxdegree)-nodes.txt",nodes)
    writedlm("formulas/($n,$m,$maxdegree)-weights.txt",weights)
    nodes,weights
   end
end

## tensor product gauss formula on [0,1]^N
function gausstensor(N,maxdegree)
    degree = convert(Int64,(maxdegree+1)/2)
    nodes0,weights0 = gausslegendre(degree)

    #transform from [-1,1] to [0,1]
    nodes0 = (nodes0 + ones(length(nodes0)))/2
    weights0 = weights0/2

    weights = kron(ntuple(n->weights0,N)...)

    tmp = []
    for k in 1:N
        push!(tmp,nodes0)
    end

    tuplenodes = collect(product(tmp...))

    nodes = zeros(degree^N,N)

    for k in 1:size(nodes,1)
        for l in 1:size(nodes,2)
            nodes[k,l] = tuplenodes[k][l]
        end
    end

    nodes,weights
end
