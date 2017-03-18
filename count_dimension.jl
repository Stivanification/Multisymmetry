include("vector_partitions.jl")

using Combinatorics
using Iterators

# Computes all vectors of nonnegative integers of length(n) with sum k
function multisets(n,k)
 map(A -> [sum(A .== i) for i in 1:n],with_replacement_combinations(1:n, k))
end

function count_dimension(n,m,d)
  multidegrees = []

  for k in 1:min(n,d)
    append!(multidegrees,multisets(m,k))
  end
  count = 0

  for delta in multidegrees
    count = count + length(vector_partitions(delta))
  end

  count
end

function count_dimension2(n,m,d)
  range = vec([k for k in 1:d])

  tmp = []
  for k in 1:m
      push!(tmp,range)
  end

  multidegrees = collect(product(tmp...))
  count = 0

  for delta in multidegrees
    vp = vector_partitions(collect(delta))
    count = count + length(find(x -> length(x) <= n,vp))
  end

  count
end

function compare_dimensions(m,n,d)
  A = 0
  B = 0

  for i in 1:d
      multset = multisets(m,i)
      for mult in multset
          vp = vector_partitions(mult)
          A = A + length(find(x -> length(x) <= n,vp))
              for j in 1:length(vp)
                   tmp = find(x -> sum(x) <= n, vp[j])
                   if length(tmp) == length(vp[j])
                       B = B + 1
                   end
              end
       end
  end
  println(A)
  println(B)
end

function basis_partitions(n,m,d)
  multidegrees = []

  for k in 1:min(n,d)
    append!(multidegrees,multisets(m,k))
  end

  basis = []
  for delta in multidegrees
    append!(basis,vector_partitions(delta))
  end
  basis
end
