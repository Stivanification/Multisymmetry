include("symbolic.jl")
include("LPcubature.jl")

  using SymPy
  a1,a2,a3,a4 = Sym("a1","a2","a3","a4")
  b1,b2,b3,b4 = Sym("b1","b2","b3","b4")
  nodes5,weights5 = return_cubature_formula_new(4,2,5)
  nodes7,weights7 = return_cubature_formula_new(4,2,7)

  # scale formula to [1,2]^4
  nodes5 = ones(size(nodes5,1),size(nodes5,2)) + nodes5
  nodes7 = ones(size(nodes7,1),size(nodes7,2)) + nodes7

  testfunctions = [(a1,b1,a2,b2,a3,b3,a4,b4) -> exp(-a1-a2-a3-a4 + b1+b2+b3+b4), (a1,b1,a2,b2,a3,b3,a4,b4) -> 1+ 1/(a1+a2+a3+a4) - 1/(b1+b2+b3+b4)]
  analytics = [(1/e-1/e^2)^4*(e^2-e)^4, 1]
  quadf5 = zeros(length(testfunctions))
  quadf7 = zeros(length(testfunctions))

  quadf5_tmp = [[testfunctions[i](nodes5[k,:]...) for k in 1:size(nodes5,1)]'weights5 for i in 1:length(testfunctions)]
  quadf7_tmp = [[testfunctions[i](nodes7[k,:]...) for k in 1:size(nodes7,1)]'weights7 for i in 1:length(testfunctions)]

  for i in 1:length(quadf5_tmp)
    quadf5[i] = quadf5_tmp[i][1]
    quadf7[i] = quadf7_tmp[i][1]
  end

  println("\nExact integrals: $analytics")
  println("Degree 5 Symmetry Rule: $quadf5")
  println("Degree 7 Symmetry Rule: $quadf7")

  symrelerror5 = vec([abs(analytics[i]-quadf5[i])/abs(analytics[i]) for i in 1:length(testfunctions)])
  symrelerror7 = vec([abs(analytics[i]-quadf7[i])/abs(analytics[i]) for i in 1:length(testfunctions)])

  println("\n Degree 5 comparison \n")

  for i in 1:length(testfunctions)
    println("Function $(testfunctions[i](a1,b1,a2,b2,a3,b3,a4,b4))")
    println("Symmetry Rule relative error: $(symrelerror5[i])\n")
  end
  println("\nnumber of nonzero nodes: Symmetry Rule $(length(find(weights5)))")

  println("\n Degree 7 comparison \n")

  for i in 1:length(testfunctions)
    println("Function $(testfunctions[i](a1,b1,a2,b2,a3,b3,a4,b4))")
    println("Symmetry Rule relative error: $(symrelerror7[i])\n")
  end
  println("\nnumber of nonzero nodes: Symmetry Rule $(length(find(weights7))) \n")
