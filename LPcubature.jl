include("symbolic.jl")

using SymPy
using JuMP
using Clp
using PyPlot

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
  A,moments,nodes = assemble_constraints(n,m,maxdegree)
  weights = LPcube(A,moments)

  ind = find(weights)
  weights = weights[ind]
  nodes = nodes[ind,:]

  s,d = create_symmetricbasis(n,m,maxdegree)
  varnames = create_variables(n,m)
	boundaries =""

	for j in 1:n
      for k in 1:m
          boundaries = boundaries*"("*varnames[k]*"$j"*",0,1),"
      end
  end

  boundaries = chop(boundaries)
  lambstring =""

	for j in 1:n
      for k in 1:m
          lambstring = lambstring*varnames[k]*"$j,"
      end
  end

	lambstring = chop(lambstring)

	fs = []
	for p in s
		push!(fs,eval(parse("lambdify($p, ["*lambstring*"])")))
	end

	println("integrate($f,"*boundaries*")")
	println(typeof(a1))
  #println(eval(parse("integrate(f"*tmp*","*boundaries*")")))
  println(boundaries)
  exactf = float(eval(parse("integrate($f,"*boundaries*")")))
  ff = eval(parse("lambdify($f, ["*lambstring*"])"))
  quadf = [ff(nodes[k,:]...) for k in 1:size(nodes,1)]'weights
  error = abs(exactf-quadf)
  relerror = error/exactf
  println("exact integral: $exactf")
  println("interpolation: $quadf")
  println("error: $error")
  println("relative error: $relerror")

  weights = weights[find(weights)]
  weights
end

function testthetest()
  a1,b1,a2,b2 = Sym("a1","b1","a2","b2")
  a3,a4 = Sym("a3","a4")
	f = (a1,b1,a2,b2) -> b1*b2*exp(a1+a2)
	g = (x1,x2) -> exp(x1+x2)
  h = (x1,x2,x3,x4) -> exp(-(x1+x2+x3+x4))
  test(2,2,7,f(a1,b1,a2,b2))
end
