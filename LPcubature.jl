include("symbolic.jl")
include("partition.jl")

using SymPy
using JuMP
using Clp
using PyPlot
using Combinatorics
using FastGaussQuadrature

function LPcube(A,moments)
  M = Model()
  @variable(M, x[1:size(A,2)] >= 0)
  @constraint(M, A*x.== moments)
  @objective(M, Min, sum(0*x))
  status = JuMP.solve(M)

  weights = getvalue(x)
end

# checks if a formula for the parameters (n,m,maxdegree) has already been computed.
# If not, computes and saves the formula.
function return_cubature_formula(n,m,maxdegree)
  nodes,weights = try
    nodes,weights = readdlm("formulas/($n,$m,$maxdegree)-nodes.txt"),readdlm("formulas/($n,$m,$maxdegree)-weights.txt")
  catch
    println("computing formula...")
    A,moments,nodes = assemble_constraints(n,m,maxdegree)
	println(size(A))
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

function return_cubature_formula_new(n,m,maxdegree)
  nodes,weights = try
    nodes,weights = readdlm("formulas-new/($n,$m,$maxdegree)-nodes.txt"),readdlm("formulas-new/($n,$m,$maxdegree)-weights.txt")
  catch
    println("computing formula...")
    A,moments,nodes = assemble_constraints_new(n,m,maxdegree)
	tic()
    weights = LPcube(A,moments)
    t = toq()
    println("Simplex solved in $t seconds.")

    #restrict to nodes with nonzero weights
    ind = find(weights)
    weights = weights[ind]
    nodes = nodes[ind,:]

    writedlm("formulas-new/($n,$m,$maxdegree)-nodes.txt",nodes)
    writedlm("formulas-new/($n,$m,$maxdegree)-weights.txt",weights)
    nodes,weights
   end
end



function test(n,m,maxdegree,f,mode = true)
	if mode  nodes,weights = return_cubature_formula(n,m,maxdegree)
	else nodes,weights = return_cubature_formula_new(n,m,maxdegree)
	end
  varnames = create_variables(n,m)
	boundaries =""
  lambstring =""

	for j in 1:n
      for k in 1:m
          lambstring = lambstring*varnames[k]*"$j,"
          boundaries = boundaries*"("*varnames[k]*"$j"*",0,1),"
      end
  end

  boundaries = chop(boundaries)
	lambstring = chop(lambstring)

	println("integrate($f,"*boundaries*")")
  exactf = float(eval(parse("integrate($f,"*boundaries*")")))
  ff = eval(parse("lambdify($f, ["*lambstring*"])"))
  quadf = [ff(nodes[k,:]...) for k in 1:size(nodes,1)]'weights
  error = abs(exactf-quadf)
  relerror = error/exactf
  println("exact integral: $exactf")
  println("interpolation: $quadf")
  println("error: $error")
  println("relative error: $relerror")
  println("number of nodes: $(length(weights))")
end

function testthetest()
  a1,b1,a2,b2,a3,b3,a4,b4 = Sym("a1","b1","a2","b2","a3","b3","a4","b4")
  a3,a4 = Sym("a3","a4")
  c1,d1 = Sym("c1","d1")
	f = (a1,b1,a2,b2,a3,b3,a4,b4) -> exp(b1+b2+b3+b4)
	g = (x1,x2) -> exp(x1+x2)
  h = (x1,x2,x3,x4) -> exp((x1+x2+x3+x4))
  #test(4,2,5,f(a1,b1,a2,b2,a3,b3,a4,b4))
  println("alte formeln:")
  test(4,1,7,h(a1,a2,a3,a4))
  println("neue formeln:")
  test(4,1,7,h(a1,a2,a3,a4),false)
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




