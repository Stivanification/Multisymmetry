include("symbolic.jl")
include("LPcubature.jl")

using SymPy

function gethigh(n,m,d)

  N = n*m

  varnames = create_variables(n,m)
  nodes,weights = return_cubature_formula_new(n,m,d)

  # scale formula to [1,2]^4
  nodes = ones(size(nodes,1),size(nodes,2)) + nodes

  lambstring ="" 
  #boundaries =""

  for j in 1:n
      for k in 1:m
          lambstring = lambstring*varnames[k]*"$j,"
         # boundaries = boundaries*"("*varnames[k]*"$j"*",1,2),"
      end
  end

  lambstring = chop(lambstring); #boundaries = chop(boundaries)

  testfunctions = [x -> exp(-sum(x)), x -> N/prod(x), x -> sum(x'x/N), #x->1/(sum(x)/N) #x -> log(sum(x)/N), 
  x -> log(prod(x)), x -> cos(sum(x)), #x -> (sum(x)/N)^(0.5), #x -> (sum(x)/N)^(-3.5), 
  x -> 2^(sum(x))]

  analytic = [(exp(-1)-exp(-2))^N, N*log(2)^N, 7.0/3.0, N*(2*log(2) - 1), #
  real(0.5*((-im *(exp((2)*im)-exp(im)))^N + (im* (exp(-(2)*im) - exp(-im)))^N)),#
   (2/log(2))^N ]
  expressions = ["exp(-sum(x))", "N/prod(x)", "x'x/N", "log(prod(x))", "cos(sum(x))", "2^(sum(x))"]

  # for i in 1:length(testfunctions)
  #   push!(expressions,testfunctions[i](eval(parse("$lambstring"))))
  #   println("integrating $(expressions[end])")
  #   analytic[i]= float(eval(parse("integrate($(expressions[end]),"*boundaries*")")))
  # end

  quadf = zeros(length(testfunctions))

  quadf_tmp = [[testfunctions[i](nodes[k,:]) for k in 1:size(nodes,1)]'weights for i in 1:length(testfunctions)]

  for i in 1:length(quadf_tmp)
    quadf[i] = quadf_tmp[i][1]
  end

  println("\nExact integrals: $analytic")
  println("Degree $d Multisymmetry Rule: $quadf \n")

  symrelerror = vec([abs(analytic[i]-quadf[i])/abs(analytic[i]) for i in 1:length(testfunctions)])

  for i in 1:length(testfunctions)
    println("Function $(expressions[i])")
    println("Symmetry Rule relative error: $(symrelerror[i])\n")
  end
  println("\nnumber of nonzero nodes: $(length(find(weights)))")

end
