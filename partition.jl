using Combinatorics
using Iterators
using FastGaussQuadrature
using LegacyStrings

# Creates all partitions of maximal length n of set [1,...,l]
function partitionize(n,l)
	L = [k for k in 1:l]
	B = collect(partitions(L,1))
	for i in 2:n
		append!(B,collect(partitions(L,i)))
	end
	B
end

# Computes all vectors of nonnegative integers of length(n) with sum k
function multisets(n,k)
 map(A -> [sum(A .== i) for i in 1:n],with_replacement_combinations(1:n, k))
end

# Creates vector M of m-tuples which represent the power of the monomial
function generate_monomialvec(m,maxdegree)
	Mtmp = collect(product(collect(repeated(0:maxdegree,m))...))
	M = []; D = []; Y =[]
	for j in Mtmp	if sum(j) <= maxdegree	j = collect(j); push!(M,reshape(j,m,1)); push!(D,sum(j)) end end
	M, D, Y
end

function compute_monvals(n,m,M,combonodes)
	X = zeros(length(M), size(combonodes)[1])
	for i in 1:length(M)
		tmp1 = combonodes'.^repmat(M[i],n)
		tmp2 = zeros(size(combonodes)[1],n)
		for j in 1:n	tmp2[:,j] = vec(prod(tmp1[1 + m*(j-1):j*m,:],1))	end
		X[i,:] = sum(tmp2,2)			
	end
	X
end

# Creates all possible products of monomials of maximal degree maxdegree
function generate_generatingsys(n,m,maxdegree,combonodes)
	tic()
	M, D = generate_monomialvec(m,maxdegree)
	t1 = toc()
	println("Basis creation: Step 1: $(t1)s")
	tic()
	M = M[2:end]
	D = D[2:end]
	tic()
	X = compute_monvals(n,m,M,combonodes)
	t3 = toc()
	P = [k for k in 1:length(M)]
	G = []
	GD = zeros(length(D))
	for i in 1:length(D)	GD[i] = D[i]; push!(G,M[i])	end
	count = 2
	tmp = 0
	# count ... number of products (must be smaller than maxdegree)
	while count <= maxdegree
		tmp1 = length(G)
		Ind1 = find(GD[tmp+1:end].< maxdegree)
		# Ind1 ... find products which are still feasible (degree smaller than maxdegree)
		for j in Ind1 + tmp
			Ind2 = find(D + GD[j] .<=maxdegree)
			Ind3 = find(Ind2.>=P[j])
			# Ind2 ... find monomials which are allowed to be multiplied
			# (degree smaller than maxdegree and index in M bigger than index of last factor)
			for kk in Ind3
				k = Ind2[kk]
				push!(G,[G[j] M[k]])
				push!(GD,GD[j]+D[k])
				push!(P,k)
				X = cat(1,X,(X[j,:].*X[k,:])')
			end
		end
		tmp = tmp1
		count += 1
	end
	t2 = toc()
	println("Basis creation: Step 2: $(t2)s")
	println("Basis created!")
	G,GD,X
end

# Computation of integrals (moments) according to formula (see paper)
function calculate_integral(G,n,m,maxdegree)
	B = []
	for i in 1:maxdegree	push!(B, partitionize(n,i)) end
	# Creates all partitions of the sets {1,...,maxdegree}
	I = zeros(length(G))
	for i in 1:length(I)
		for b in B[size(G[i])[2]]
			myprod = n
			for j in 1:length(b)-1	myprod *= n-j	end
			# Compute the factor n*(n-1)*...*(n-|b|+1)
			mypot = zeros(n*m)
			for l in 1:length(b) mypot[1 + (l-1)*m:l*m] = vec(sum(G[i][:,b[l]],2))	end
			mypot += 1
			# Create vector of potences
			I[i] += myprod / prod(mypot)
		end
	end
	I
end

function test(n,m,maxdegree)
	degree = convert(Int64,(maxdegree+1)/2)
	nodes = gausslegendre(degree)[1]
	tic()
	
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

	t2 = toc()
	println("nodes calculated in $t2 seconds")
	combonodes
end

function assemble_constraints_new(n,m,maxdegree)
	combonodes = test(n,m,maxdegree)
	G,D,X = generate_generatingsys(n,m,maxdegree,combonodes)
	tic()
	a,b = size(X)
	Y = ones(a+1,b)
	Y[2:end,1:end] = X
	M = calculate_integral(G,n,m,maxdegree)
	MM = ones(length(M)+1)
	MM[2:end] = M
	t1 = toc()
	println("Moments computed in $(t1) s")
	Y,MM,combonodes
end

function combos_with_replacement(list, k)
	n = length(list)
	[[list[c[i]-i+1] for i = 1:length(c)] for c in combinations(1:(n+k-1),k)]
end

function testt(n,m,maxdegree)
	a = generate_generatingsys(n,m,maxdegree,zeros(1,n*m))[1]
	Mtmp = generate_monomialvec(m,maxdegree)[1]
	vs = create_variables(n,m)
	P = Array{Any}(length(a))
	P1 = create_symmetricbasis(n,m,maxdegree)[1]
	
	create_symbols(m,"Y")
	Vartmp = "[Y1"
	for i in 2:m	Vartmp *= " Y$i" end
	Vartmp *= "]"
	yvars = eval(parse(Vartmp))
	M = Array{Any}(length(Mtmp))
	for i in 1:length(Mtmp)
		M[i] = symbolic_prod(yvars, Mtmp[i])
	end

	for l in 1:length(a)
		b = a[l]
		Ptmp2 = 1
		for i in 1:size(b)[2]
			Ptmp = 0
			count = 1
			while true
				(Mtmp[count] == reshape(b[:,i],m,1)) && break
				count += 1				
			end
			for j in 1:n
				ex = M[count]
				for k in 1:m
					ex = subs(ex,eval(parse("Y$k")), eval(parse(vs[k]*"$j")))
				end
				if j == 1 Ptmp = ex
				else Ptmp += ex
				end
			end
			Ptmp2 *= Ptmp	
		end
		P[l] = factor(Ptmp2)
	end
	temp = zeros(length(P),2)
	for i in 1:length(P)
		for j in 1:length(P)
			if P[i] == P1[j]
				temp[i,1] = i
				temp[i,2] = j
			end
			temp[i,1] != 0 && break
		end
	end
	P, P1, temp	
end















