using SymPy
using LegacyStrings
using Combinatorics
using Iterators
using FastGaussQuadrature

function create_symbols(n,name)
	X = name*"1"
	A = "=Sym(\""
	B = "\")"
	C = "@syms "
	for i = 2:n
		X = X*","*name*"$i"
	end
	eval(parse(C*X))
	eval(parse(X*A*X*B))
end

function check26(x)
	if maximum(x) > 25
		x[find(x.== maximum(x))] = 0
		x[find(x.== maximum(x))] += 1
		x = check26(x)
	end
	x
end

# In our case: m ... physical dimensions, ie. 2
#			   n ... amount of particles
function create_variables(n,m)
	mlog = convert(Int,floor(log(26,m))+1)
	x = zeros(Int,mlog)
	mtmp = m
	varnames = []
	name1 = ""
	name2 = ""
	while mtmp > 0
		mtmp -= 1
		x = check26(x)
		name = ""
		for i in 1:length(x)	name *= join('a' + x[i])	end
		push!(varnames, name)
		if mtmp == m-1	name1 = name end
		if mtmp == 0	name2 = name end
		create_symbols(n,name)
		x[end] += 1
	end
	println("Variablen von $(name1) bis $(name2) erzeugt, jeweils von 1 bis $n.")
	varnames
end

# Symbolic multiplication of var[i]^potenz[i]
function symbolic_prod(vars, potenz)
	if length(vars) != length(potenz) println("Nicht richtige Anzahl an Potenzen & Variablen!"); return 0 end
	myprod = Sym("prod")
	myprod = vars[1]^potenz[1]
	for i in 2:length(vars)
		myprod *= vars[i]^potenz[i]
	end
	myprod
end

# Creates all monomials in m variables of degree maxdegree
function create_monomials(m,maxdegree)
	create_symbols(m,"Y")
	Vartmp = "[Y1"
	for i in 2:m	Vartmp *= " Y$i" end
	Vartmp *= "]"
	yvars = eval(parse(Vartmp))
	tmp1 = collect(product(collect(repeated(0:maxdegree,m))...)) # erzeugt alle möglichen kombinationen von produkten
	tmp = []
	Grad = []
	for t in tmp1
		if 0 < sum(t) <= maxdegree	push!(tmp,t); push!(Grad, sum(t)) end
	end
	M = Array{Any}(length(tmp))
	for i in 1:length(tmp)
		M[i] = symbolic_prod(yvars, tmp[i])
	end
	M, Grad
end

# Creates vector of elementary polynomial applied to a vector of symbolic functions
function apply_elpol(n,m,M)
	EP1 = Array{Any}(length(M))
	varnames = create_variables(n,m)
	for i in 1:length(M)
		for j in 1:n
			ex = M[i]
			for k in 1:m
				ex = subs(ex,eval(parse("Y$k")), eval(parse(varnames[k]*"$j")))
			end
			if j == 1 EP1[i] = ex
			else EP1[i] += ex
			end
		end
	end
	EP1
end

function create_symmetricbasis(n,m, maxdegree)
	M,Grad = create_monomials(m,maxdegree)
	EP = apply_elpol(n,m,M)
	level = 0
	symbasis = EP
	tmpGrad = []
	for grad in Grad	push!(tmpGrad, grad)	end
	level = 0
	while(level < maxdegree)
		feasibleind1 = find(tmpGrad.< maxdegree)
		feasibleind2 = find(Grad.< maxdegree)
		for i in feasibleind1
			for j in feasibleind2
				if tmpGrad[i] + Grad[j] <= maxdegree
					push!(tmpGrad, tmpGrad[i] + Grad[j])
					push!(Grad, tmpGrad[i] + Grad[j])
					push!(symbasis, symbasis[i]*symbasis[j])
				end
			end
			tmpGrad[i] = maxdegree
		end
		level += 1
	end
	ls = length(symbasis)
	symbasisfinal = []
	for i in 1:ls
		if !(symbasis[i] in symbasisfinal)
			push!(symbasisfinal,symbasis[i])
		end
	end

	symbasisfinal, Grad
end

# creates vector of the 1. elementary polynomial applied to all monomials 
function create_elpolmonomials(n,m,maxdegree)
	M,Grad = create_monomials(m,maxdegree)
	EP = apply_elpol(n,m,M)
	EP, Grad
end


# from https://rosettacode.org/wiki/Combinations_with_repetitions#Julia
function combos_with_replacement(list, k)
    n = length(list)
    [[list[c[i]-i+1] for i=1:length(c)] for c in combinations(1:(n+k-1),k)]
end


function assemble_constraints(n,m,maxdegree)

	#######################################################
	#	assemble the constraint matrix
	#######################################################

	if mod(maxdegree,2) == 0
		println("Warning: maxdegree should be an odd number")
	end
	tic()
	s,d = create_symmetricbasis(n,m,maxdegree)
	t1 = toq()
	println("basis generated in $t1 seconds")
	tic()
	degree = convert(Int64,(maxdegree+1)/2)
	nodes = gausslegendre(degree)[1]

	# transform univariate gaussian quadrature rule from [-1,1] to [0,1]
	nodes = (nodes + ones(length(nodes)))/2

	# create all possible m-tuples of node values
	tmp = []
	for k in 1:m
		push!(tmp,nodes)
	end
	mnodes = collect(product(tmp...))
	# take all combinations with repetitions of the m-tuples to obtain
	# a full representative system of the nodes modulo the permutation group G
	combonodes_tmp = combos_with_replacement(mnodes,n)

	#convert the messed up type of combonodes_tmp to matrix
	combonodes = zeros(length(combonodes_tmp),n*m)
	for i in 1:length(combonodes_tmp)
		for j in 1:n
			for k in 1:m
				combonodes[i,m*(j-1)+k]=combonodes_tmp[i][j][k]
			end
		end
	end

	t2 = toq()
	println("nodes calculated in $t2 seconds")

	tic()

	# convert symbolic expressions s[i] to julia functions
	varnames = create_variables(n,m)
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

	t3 = toq()
	println("functions converted in $t3 seconds")

	tic()
	#create constraint matrix
	A = zeros(length(s)+1,length(combonodes_tmp))
	A[1,:] = ones(1,length(combonodes_tmp))
	for i in 1:length(s)
		for j in 1:length(combonodes_tmp)
			A[i+1,j] = fs[i](combonodes[j,:]...)
		end
	end

	t4 =toq()
	println("constraints assembled in $t4 seconds")
	tic()

	#########################################
	# calculate moment vector
	#########################################

	moments = zeros(length(s)+1)
	moments[1] = 1

	boundaries =""

	for j in 1:n
      for k in 1:m
          boundaries = boundaries*"("*varnames[k]*"$j"*",0,1),"
      end
  end

	boundaries = chop(boundaries)

	for i in 1:length(s)
		moments[i+1] = eval(parse("integrate($(s[i]),"*boundaries*")"))
	end
	t5 = toq()
	println("moments calculated in $t5 seconds")

	A,moments,combonodes
end
