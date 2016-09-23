include("symbolic.jl")

using SymPy
using JuMP
using Clp

function LPcube(A,moments)
  M = Model()
  @variable(M, x[1:size(A,2)] >= 0)
  @constraint(M, A*x.== moments)
  @objective(M, Min, sum(0*x))
  status = JuMP.solve(M)

  weights = getvalue(x)
end

function return_cubature_formula(n,m,maxdegree)
  A,moments,nodes = assemble_constraints(n,m,maxdegree)
  tic()
  weights = LPcube(A,moments)
  t = toq()
  println("Simplex solved in $t seconds.")

  #restrict to nodes with nonzero weights
  ind = find(weights)
  weights = weights[ind]
  nodes = nodes[ind,:]

  nodes,weights
end

function test(n,m,maxdegree,f)
  ## funktioniert nur für erlesene funktionen f

  A,moments,nodes = assemble_constraints(n,m,maxdegree)
  weights = LPcube(A,moments)
  s,d = create_symmetricbasis(n,m,maxdegree)

  for i in 1:length(s)
  println("Fehler s[$i]: $([s[i](nodes[k,:]...) for k in 1:size(nodes,1)]'weights - moments[i+1])")
  end

  varnames = create_variables(n,m)
	boundaries =""

	for j in 1:n
      for k in 1:m
          boundaries = boundaries*"("*varnames[k]*"$j"*",0,1),"
      end
  end

	exactf = float(eval(parse("integrate($f,"*boundaries*")")))
  quadf = [f(nodes[k,:]...) for k in 1:size(nodes,1)]'weights
  error = abs(exactf-quadf)
  println(exactf)
  println(quadf)
  dump(error)

  weights

end
