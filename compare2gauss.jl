include("symbolic.jl")
include("LPcubature.jl")

  using SymPy

  nodes5,weights5 = return_cubature_formula(4,1,5)
  nodes7,weights7 = return_cubature_formula(4,1,7)

  # scale formula to [1,2]^4
  nodes5 = ones(size(nodes5,1),size(nodes5,2)) + nodes5
  nodes7 = ones(size(nodes7,1),size(nodes7,2)) + nodes7

  testfunctions = [(a1,a2,a3,a4) -> exp(-a1-a2-a3-a4), (a1,a2,a3,a4)->1/(a1+a2+a3+a4), (a1,a2,a3,a4)->a1^2 *a2^2 * a3^2 * a4^2, (a1,a2,a3,a4) -> 1/(a1*a2*a3*a4), #
  (a1,a2,a3,a4) -> log(a1+a2+a3+a4), (a1,a2,a3,a4) -> log(a1*a2*a3*a4), (a1,a2,a3,a4) -> cos(a1+a2+a3+a4), (a1,a2,a3,a4) -> (a1+a2+a3+a4)^(0.5), #
  (a1,a2,a3,a4) -> (a1+a2+a3+a4)^(-3.5), (a1,a2,a3,a4) -> 2^(a1+a2+a3+a4)]

  analytic = [0.00292429, 0.16825, 29.64198, 0.23084, 1.78707, 1.54518, 0.81162, 2.44663, 0.00203901, 69.31355]
  gaussrelerror5 = [0.00295, 0.00058, 0.00137, 0.0125, 0.00001, 0.00157, 0.00248, 0.00008, 0.00319, 0.00020]
  gaussrelerror7 = [0.00073, 0.00001, 4.0e-06, 0.00082, 9.9e-07, 0.00004, 0.00107, 0.00009, 0.00042, 0.00003]
  quadf5 = zeros(length(testfunctions))
  quadf7 = zeros(length(testfunctions))

  quadf5_tmp = [[testfunctions[i](nodes5[k,:]...) for k in 1:size(nodes5,1)]'weights5 for i in 1:length(testfunctions)]
  quadf7_tmp = [[testfunctions[i](nodes7[k,:]...) for k in 1:size(nodes7,1)]'weights7 for i in 1:length(testfunctions)]

  for i in 1:length(quadf5_tmp)
    quadf5[i] = quadf5_tmp[i][1]
    quadf7[i] = quadf7_tmp[i][1]
  end

  println("\nExact integrals: $analytic")
  println("Degree 5 Symmetry Rule: $quadf5")
  println("Degree 7 Symmetry Rule: $quadf7")

  symrelerror5 = vec([abs(analytic[i]-quadf5[i])/abs(analytic[i]) for i in 1:length(testfunctions)])
  symrelerror7 = vec([abs(analytic[i]-quadf7[i])/abs(analytic[i]) for i in 1:length(testfunctions)])

  println("\n Degree 5 comparison \n")

  for i in 1:length(testfunctions)
    println("Function $(testfunctions[i](a1,a2,a3,a4))")
    println("Full Gauss relative error: $(gaussrelerror5[i])")
    println("Symmetry Rule relative error: $(symrelerror5[i])\n")
  end
  println("\nnumber of nodes: Full Gauss 69, Symmetry Rule $(length(weights5))")

  println("\n Degree 7 comparison \n")

  for i in 1:length(testfunctions)
    println("Function $(testfunctions[i](a1,a2,a3,a4))")
    println("Full Gauss relative error: $(gaussrelerror7[i])")
    println("Symmetry Rule relative error: $(symrelerror7[i])\n")
  end
  println("\nnumber of nodes: Full Gauss 159, Symmetry Rule $(length(weights7)) \n")
