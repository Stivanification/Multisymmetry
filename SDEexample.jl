@everywhere using DriftDiffusionPoissonSystems
@everywhere using Sobol
@everywhere include("LPcubature.jl")

@everywhere function eval_integrand(x,y,realization)
	#(sum(realization)-length(realization)/2)*sin(x.^2).*sin(y.^2)
	hcat(mean(exp(-vec(realization)*vec(x.^2+y.^2)'),1))
end

@everywhere function integral(x,y)
	#(sum(realization)-length(realization)/2)*sin(x.^2).*sin(y.^2)
	(1-exp(-(x.^2+y.^2)))./(x.^2+y.^2)
end

@everywhere function multisym_expectation(n,d)
	m = 1
	nodes_sym, weights_sym = return_cubature_formula_new(n,m,d)
	A = [x->1;x->0;x->0;x->1]
	bddata=[1 2 3 4;'D' 'D' 'D' 'D'; (x,y)->1 (x,y)->1 (x,y)->1 (x,y)->1]
	mesh = read_mesh("mesh_s_p05.msh")
	F = [(x,y,u)->eval_integrand(x,y,nodes_sym[i,:]) for i in 1:size(nodes_sym)[1]]
	eval_tmp = pmap((j) -> solve_semlin_poisson(mesh,A,bddata,F[j],(x,y,u)->0),1:length(F))
	sum(eval_tmp.*weights_sym)
end

@everywhere function qmc_expectation(n,N)
	m = 1
	S = SobolSeq(n)
	p = hcat([next(S) for i = 1:N]...)'
	A = [x->1;x->0;x->0;x->1]
	bddata=[1 2 3 4;'D' 'D' 'D' 'D'; (x,y)->1 (x,y)->1 (x,y)->1 (x,y)->1]
	mesh = read_mesh("mesh_s_p05.msh")
	F = [(x,y,u)->eval_integrand(x,y,p[i,:]) for i in 1:N]
	eval_tmp = pmap((j) -> solve_semlin_poisson(mesh,A,bddata,F[j],(x,y,u)->0),1:length(F))
	mean(eval_tmp)
end

@everywhere function vmc_expectation(n,N)
	m = 1
	p = hcat([rand(n) for i = 1:N]...)'
	A = [x->1;x->0;x->0;x->1]
	bddata=[1 2 3 4;'D' 'D' 'D' 'D'; (x,y)->1 (x,y)->1 (x,y)->1 (x,y)->1]
	mesh = read_mesh("mesh_s_p05.msh")
	F = [(x,y,u)->eval_integrand(x,y,p[i,:]) for i in 1:N]
	eval_tmp = pmap((j) -> solve_semlin_poisson(mesh,A,bddata,F[j],(x,y,u)->0),1:length(F))
	mean(eval_tmp)
end

function exact()
	A = [x->1;x->0;x->0;x->1]
	bddata=[1 2 3 4;'D' 'D' 'D' 'D'; (x,y)->1 (x,y)->1 (x,y)->1 (x,y)->1]
	mesh = read_mesh("mesh_s_p05.msh")
	f = (x,y,u) -> (1-exp(x.^2+y.^2))./(x.^2+y.^2)
	vec(solve_semlin_poisson(mesh,A,bddata,(x,y,u)->integral(x,y),(x,y,u)->0))
end

function test(n,d,N)
	mesh = read_mesh("mesh_s_p05.msh")
	mulsym = multisym_expectation(n,d)
	qmc = qmc_expectation(n,N)
	ex = exact()
	calculate_norm(mesh,(ex-mulsym).^2), calculate_norm(mesh,(ex-qmc).^2)
end

function plotit(n,d,N)
	mesh = read_mesh("mesh_s_p05.msh")
	X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
	mulsym = multisym_expectation(n,d)
	qmc = qmc_expectation(n,N)
	ex = exact()
	clf()
	figure(1)
	surf(X,Y,abs(vec(ex)-vec(mulsym)),cmap="jet")
	figure(2)
	surf(X,Y,abs(vec(ex)-vec(qmc)),cmap="jet")
	sqrt(calculate_norm(mesh,(ex-mulsym).^2)), sqrt(calculate_norm(mesh,(ex-qmc).^2))
end

function calc_sym(n,d)
	mesh = read_mesh("mesh_s_p05.msh")
	X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
	mulsym = multisym_expectation(n,d)
	ex = exact()
	sqrt(calculate_norm(mesh,(ex-mulsym).^2)), mulsym
end

function calc_qmc(n,N)
	mesh = read_mesh("mesh_s_p05.msh")
	X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
	qmc = qmc_expectation(n,N)
	ex = exact()
	sqrt(calculate_norm(mesh,(ex-qmc).^2)), qmc
end

function calc_vmc(n,N)
	mesh = read_mesh("mesh_s_p05.msh")
	X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
	vmc = vmc_expectation(n,N)
	ex = exact()
	sqrt(calculate_norm(mesh,(ex-vmc).^2)), vmc
end

function calc_qmc_fromto(n,N1,N2,N)
	m = 1
	S = SobolSeq(n)
	p = (hcat([next(S) for i = 1:N]...)')[N1:N2]
	A = [x->1;x->0;x->0;x->1]
	bddata=[1 2 3 4;'D' 'D' 'D' 'D'; (x,y)->1 (x,y)->1 (x,y)->1 (x,y)->1]
	mesh = read_mesh("mesh_s_p05.msh")
	F = [(x,y,u)->eval_integrand(x,y,p[i,:]) for i in 1:(N2-N1)]
	eval_tmp = pmap((j) -> solve_semlin_poisson(mesh,A,bddata,F[j],(x,y,u)->0),1:length(F))
	mean(eval_tmp)
end

function calculate_norm(mesh,u)
        a1 = vec(mesh.elements[1,:])
        a2 = vec(mesh.elements[2,:])
        a3 = vec(mesh.elements[3,:])

        q1 = mesh.nodes[:,a1]
        q2 = mesh.nodes[:,a2]
        q3 = mesh.nodes[:,a3]
        uu = q2-q3
        vv = q3-q1
        ww = q1-q2
        ar = (uu[1,:].*vv[2,:]-uu[2,:].*vv[1,:])./2

        sol=((u[a1]+u[a2]+u[a3])./3)
        sol2=ones(size(sol,1),1)
        norm=sqrt(sum(abs(ar'.*sol.^2)))

        return norm
end

#=
sym_sol = []; sym_err = [];
for d in [3:2:11;]
	err, sol = calc_sym(15,d)
	push!(sym_sol,sol)
	push!(sym_err,err)
end
push!(sym_sol, vec(exact()))

qmc_err = []; qmc_sol = [];
for N in [10;100;1000;]
#	for N in [10000;]
	err, sol = calc_qmc(15,N)
	push!(qmc_sol, sol)
	push!(qmc_err, err)
end
push!(qmc_sol, vec(exact()))
=#
#=

sym_err = []; sym_sol = [];
for d in [3:2:11;]
	err, sol = calc_sym(15,d)
	push!(sym_sol,sol)
	push!(sym_err,err)
end
push!(sym_sol, vec(exact()))
mesh = read_mesh("mesh_s_p05.msh")
X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])

fig = figure()
for i in 1:2
	for j in 1:3
		ax = fig[:add_subplot](3,2,i+2*(j-1), projection = "3d")
		Z = vec(sym_sol[i+2*(j-1)])
		if i+2*(j-1) != 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(sym_sol[end]) - vec(sym_sol[i+2*(j-1)]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks

			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
	end
end

savefig("sde")
=#

function plot_sol(sol, err, ex, name)
mesh = read_mesh("mesh_s_p05.msh")
X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
fig = figure(name,figsize=(20,12),dpi=80,facecolor="w",edgecolor="k")
for i in 1:2
	for j in 1:3
		ax = fig[:add_subplot](3,2,i+2*(j-1), projection = "3d")
		Z = vec(ex)
		if i+2*(j-1) != 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(ex) - vec(sol[i+2*(j-1)]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks

			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
	end
end
savefig(name)
end

function create_plots(n)
mesh = read_mesh("mesh_s_p05.msh")
X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])

sym_err = []; sym_sol = [];
dd = [3:2:11;]

for d in dd
	tmp_err, tmp_sol = calc_sym(n,d)
	push!(sym_sol,tmp_sol)
	push!(sym_err,tmp_err)
end

ex = exact()

sol = sym_sol
fig = figure("Sym",figsize=(20,12),dpi=80,facecolor="w",edgecolor="k")
for i in 1:2
	for j in 1:3
		ax = fig[:add_subplot](3,2,i+2*(j-1), projection = "3d")
		Z = vec(ex)
		if i+2*(j-1) != 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(ex) - vec(sol[i+2*(j-1)]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks

			d = convert(Int64,dd[i+2*(j-1)])
			title("d = $d")
			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			title("Exact solution")
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
	end
end
savefig("sym.svg")

qmc_err = []; qmc_sol = [];
NN = [10;100;200;500;1000]
for N in NN
	tmp_err, tmp_sol = calc_qmc(n,N)
	push!(qmc_sol, tmp_sol)
	push!(qmc_err, tmp_err)
end

sol = qmc_sol
fig = figure("qMC",figsize=(20,12),dpi=80,facecolor="w",edgecolor="k")
for i in 1:2
	for j in 1:3
		ax = fig[:add_subplot](3,2,i+2*(j-1), projection = "3d")
		Z = vec(ex)
		if i+2*(j-1) != 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(ex) - vec(sol[i+2*(j-1)]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks
			
			N = convert(Int64,NN[i+2*(j-1)])
			title("N = $N")
			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			title("Exact solution")
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
	end
end
savefig("qmc.svg")
writedlm("qmc_sol.txt",qmc_sol)
writedlm("qmc_err.txt",qmc_err)
writedlm("sym_sol.txt",sym_sol)
writedlm("sym_err.txt",sym_err)
end


function plotdaplots(fig_x,fig_y)
mesh = read_mesh("mesh_s_p05.msh")
X = vec(mesh.nodes[1,:]); Y = vec(mesh.nodes[2,:])
ex = exact()
NN = [10;100;200;500;1000]
dd = [3:2:11;]
sol = readdlm("sym_sol.txt")
fig = figure("Sym",figsize=(fig_x,fig_y),dpi=80,facecolor="w",edgecolor="k")
for i in 1:6
		ax = fig[:add_subplot](6,1,i, projection = "3d")
		Z = vec(ex)
		if i!= 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(ex) - vec(sol[i,:]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks

			d = convert(Int64,dd[i])
			title("d = $d")
			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			title("Exact expectation")
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
end
savefig("fsc_sde_sym.pdf")

sol = readdlm("qmc_sol.txt")
fig = figure("qMC",figsize=(fig_x,fig_y),dpi=80,facecolor="w",edgecolor="k")
for i in 1:6
		ax = fig[:add_subplot](6,1,i, projection = "3d")
		Z = vec(ex)
		if i != 6 
			f = matplotlib[:ticker][:FormatStrFormatter]("%1.1f") # Define format of tick labels
			ax[:zaxis][:set_major_formatter](f) # Set format of tick labels
			tmp = abs(vec(ex) - vec(sol[i,:]))
			exponent = floor(Int64,maximum(log10(tmp)))
			Z = tmp*10^(-exponent)

			Mx = matplotlib[:ticker][:MultipleLocator](ceil(maximum(Z)/3)) # Define interval of major ticks
			ax[:zaxis][:set_major_locator](Mx) # Set interval of major ticks
			
			N = convert(Int64,NN[i])
			title("N = $N")
			zlabel("e-$(exponent)",rotation = 90)
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		else
			title("Exact expectation")
			ax[:plot_trisurf](X,Y,Z,cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
		end
end
savefig("fsc_sde_qmc.pdf")
end



























