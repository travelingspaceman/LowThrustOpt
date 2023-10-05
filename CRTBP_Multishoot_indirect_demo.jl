#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#DEMONSTRATION OF TRAJECTORY OPTIMIZATION AND NEURAL NETWORK FITTING:
#   - Initial guess generation
#   - Optimization with direct method (robust)
#   - Optimization with indirect method (sensitive but fast)
#   - Then a neural network is fit to map (δx₀, t) -> (δλ_v(t)).
#   - Neural network is tested, with results plotted at the end of the script.
#
#Note: It is normal for this to take a long time (a few minutes) to run the
#first time. This is because the "DifferentialEquations" and "ForwardDiff"
#packages are very slow to compile. Running subsequent solutions will be many
#times faster.
#
#Note: there are a handful of paths which must be set manually in this script --
#see notes in the code drawing attention to them

#On the first run, import all the packages
#Packages by Nathan, from https://github.com/travelingspaceman



#Packages from other sources:
using Plots
using DifferentialEquations
using Ipopt
using JuMP
using Interpolations
using ForwardDiff
using MATLAB
using Printf
using DelimitedFiles
using LinearAlgebra
using Test
Plots.plotlyjs()


include("GeneralCode/ode.jl")
include("src/CRTBP_prop_EP_deriv.jl")
include("src/CRTBP_stateCostate_deriv.jl")
include("src/HelperFunctions.jl")
include("src/multiShoot_CRTBP_direct.jl")
include("src/multiShoot_CRTBP_indirect.jl")
include("GeneralCode/norm_many.jl")
include("GeneralCode/sphere.jl")
println("  Done loading packages.")



#Earth-Moon system:
const global MU = 0.012150585609624037 # MU = mu_moon/(mu_moon + mu_planet)
const global DU = 384747.96285603708 #Earth-Moon distance (km) (km per distance unit)
const global TU = 375699.81732246041 #time units (seconds per time unit)
const global r_moon = 1737. #km
const global r_earth = 6378. #km
const global day = 86400. #seconds
##Calculate mu_moon based on MU for consistency with other code:
const global mu_planet = 398600.4415
const global mu_moon = (MU * mu_planet) / (1 - MU)


#Read in initial & final states:
print("Reading endpoints files...")
path_j = "C:/Users/Tyler/.julia/dev/LowThrustOpt/"
X0_states = readdlm(string(path_j, "L2_Anderson_1.txt"))
X0_times = LinRange(0, 1, size(X0_states,2))
Xf_states = readdlm(string(path_j, "L2_Anderson_2.txt"))
Xf_times = LinRange(0, 1, size(Xf_states,2))
print("  Endpoints loaded.\n")

## TRAJECTORY STACKING INITIAL GUESS ###########################################
tof1 = 10 * day / TU #TU (time in first orbit)
tof2 = 10 * day / TU #TU (time in second orbit)
tof = tof1 + tof2 #TU
n_nodes = 30
t_TU = collect(LinRange(0, tof, n_nodes))
t_TU1 = t_TU[t_TU .<  tof1]
t_TU2 = t_TU[t_TU .>= tof1]

p = 1.
rho = 1.
mass = 1e3 #kg
time_direction = -1
thrustLimit = 0.
time_direction = 1.

params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)

#Ballistically propagate first orbit:
tspan = (0., tof1)
τ1 = 0.75
state_0 = vcat(interpInitialStates(τ1, X0_times, X0_states, MU), zeros(6))
prob = ODEProblem(CRTBP_stateCostate_deriv!, state_0, tspan, params)
sol_orbit1 = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

#Find closest point on the second orbit:
τ2₀ = find_τ(Xf_times, Xf_states, sol_orbit1[1:6, end], MU)
state_f₀ = vcat(interpInitialStates(τ2₀, Xf_times, Xf_states, MU), zeros(6))

#Ballistically propagate second orbit:
tspan = (tof1, tof1 + tof2)
prob = ODEProblem(CRTBP_stateCostate_deriv!, state_f₀, tspan, params)
sol_orbit2 = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

#Fix final state:
τ2 = find_τ(Xf_times, Xf_states, sol_orbit2[1:6, end], MU)
state_f = vcat(interpInitialStates(τ2, Xf_times, Xf_states, MU), zeros(6))

#Concatenate states
XC_trajectoryStack = hcat(sol_orbit1(t_TU1)[:,:], sol_orbit2(t_TU2)[:,:] )
XC_trajectoryStack[:,1] = state_0
XC_trajectoryStack[:,end] = state_f

plot!(X0_states[1,:], X0_states[2,:], X0_states[3,:], label="Initial orbit")
plot!(Xf_states[1,:], Xf_states[2,:], Xf_states[3,:], label="Final orbit")
scatter!([state_0[1],], [state_0[2],], [state_0[3],], label="X0", color=:green)
scatter!([state_f[1],], [state_f[2],], [state_f[3],], label="Xf", color=:red)
scatter!([1-MU,], [0.,], [0.,], color=:black, label="Moon")
plot!(XC_trajectoryStack[1,:], XC_trajectoryStack[2,:], XC_trajectoryStack[3,:], color=:black, label="Guess solution")


################################################################################

## direct multiple shooting:
t_TU = collect(LinRange(0, tof, n_nodes))

#Random initial guess:
# X_all = 0.5*randn(6, n_nodes)
# X_all[1,:] += (1-MU)

#Trajectory stacking initial guess:
X_all = XC_trajectoryStack[1:6, :]

u_all = zeros(3, n_nodes)
dV1 = zeros(3)
dV2 = zeros(3)
Isp = 2000.
flagEnd = false
β = 0.
allowImpulsive = false
nsteps = 10
t0_TU = t_TU[1]
tf_TU = t_TU[end]
# thrustLimit = 10.
plot_yn = false
plotlyjs(legend=false)
maxIter = 100
@time (X_all, u_all, τ1, τ2, t_TU, dV1, dV2, defect) =
    multiShoot_CRTBP_direct(X_all, u_all, τ1, τ2, t_TU, dV1, dV2, MU, DU, TU, n_nodes,
    nsteps, mass, Isp, X0_times, X0_states, Xf_times, Xf_states,
    plot_yn, flagEnd, β, allowImpulsive, maxIter)

state_0 = interpInitialStates(τ1, X0_times, X0_states, MU)
state_f = interpInitialStates(τ2, Xf_times, Xf_states, MU)

plotTrajPlotly_direct(X_all, u_all, X0_states, Xf_states, 0.2)
gui()


################################################################################
## Indirect multiple shooting, with p = 2
#Initial guess from direct method:
XC_all = vcat(X_all, 0.1*randn(6, n_nodes))


#Fix the final state exactly:
XC_all[1:6, 1] = state_0[1:6]
XC_all[1:6, end] = state_f[1:6]

#With exact derivatives, we get issues when initializing with an exact
#numerically-integrated initial guess. Adding a tiny amount of noise makes
#it numerically stable.
XC_all[:, 2:end-1] += 1e-10 * randn(12, n_nodes-2)

p = 2.
thrustLimit = 10. #N ---> For p=2, leave thrust unconstrained. Constrain it later.

plot_yn = false
flag_adjointsOnly = true
plotlyjs(legend=false)
maxIter = 10
@time (XC_all, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
    mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)


flag_adjointsOnly = false
maxIter = 50
@time (XC_all, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
    mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)

u_all = zeros(3, n_nodes)
for ind = 1:n_nodes
    u_all[:,ind] = controlLaw_cart(XC_all[10:12,ind], thrustLimit, p, rho, mass)
end



# "densify" (add more points) for plotting
n_desired = 300
time_direction = 1
params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
(XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)
u_all = zeros(3, size(XC_dense, 2))
for ind = 1:size(XC_dense, 2)
    u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho, mass)
end


#plot control profile
plotlyjs(legend=true, aspect_ratio=:none)
plot(t_dense*TU/day, norm_many(u_all), xlabel="Time [days]", ylabel="Control [N]",
    lw=2, color=:black, label="Control Magnitude", size=(800, 400))
plot!(t_dense*TU/day, u_all[1,:], label="x", color=RGB(0,.5, .9))
plot!(t_dense*TU/day, u_all[2,:], label="y", color=RGB(.9, .4, .4))
plot!(t_dense*TU/day, u_all[3,:], label="z", color=RGB(0.4, .7, .4))


# Plot nicely
plotlyjs(legend=false)
plotTrajPlotly_indirect(XC_dense, u_all, X0_states, Xf_states, 0.2)


## Indirect multiple shooting, with p = 1
p = 1.
t_TU = collect(LinRange(0, tof, n_nodes))

#Fix the final state exactly:
XC_all[1:6, 1] = state_0[1:6]
XC_all[1:6, end] = state_f[1:6]


#With exact derivatives, we get issues when initializing with an exact
#numerically-integrated initial guess. Adding a tiny amount of noise makes
#it numerically stable.
XC_all[:, 2:end-1] += 1e-10 * randn(12, n_nodes-2)

thrustLimit = 0.05
println("thrustLimit = $thrustLimit N")
rho = 1
maxIter = 30
flag_adjointsOnly = false
plot_yn = false
@time (XC_all, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
    mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)


# "densify" (add more points) for plotting
n_desired = 300
time_direction = 1
params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
(XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)
u_all = zeros(3, size(XC_dense, 2))
for ind = 1:size(XC_dense, 2)
    u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho, mass)
end


#plot control profile
plotlyjs(legend=true, aspect_ratio=:none)
#plot(t_dense*TU/day, norm_many(u_all), xlabel="Time [days]", ylabel="Control [N]",
#   lw=2, color=:black, label="Control Magnitude", size=(800, 400))
plot!(t_dense*TU/day, u_all[1,:], label="x", color=RGB(0,.5, .9))
plot!(t_dense*TU/day, u_all[2,:], label="y", color=RGB(.9, .4, .4))
plot!(t_dense*TU/day, u_all[3,:], label="z", color=RGB(0.4, .7, .4))
plot!(title="Control profile")
gui()

# Plot nicely
plotlyjs(legend=false)
plotTrajPlotly_indirect(XC_dense, u_all, X0_states, Xf_states, 0.2)


## Reduce rho for sharper thrust on/off:
rho_current = copy(rho)
rho_target = 1e-4
plot_yn = false
@time (XC_all_new, defect, status_flag) = reduceFuel_indirect(XC_all, t_TU, MU, DU, TU, n_nodes, mass,
    thrustLimit, rho_current, rho_target)
if status_flag == 0
    println("Success!")
    XC_all = copy(XC_all_new)
    rho = copy(rho_target)
else
    warn("Failed to reduce rho... :/")
end



# "densify" (add more points) for plotting
n_desired = 300
time_direction = 1
params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
(XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)
u_all = zeros(3, size(XC_dense, 2))
for ind = 1:size(XC_dense, 2)
    u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho, mass)
end


#plot control profile
plotlyjs(legend=true, aspect_ratio=:none)
plot(t_dense*TU/day, norm_many(u_all), xlabel="Time [days]", ylabel="Control [N]",
    lw=2, color=:black, label="Control Magnitude", size=(800, 400))
plot!(t_dense*TU/day, u_all[1,:], label="x", color=RGB(0,.5, .9))
plot!(t_dense*TU/day, u_all[2,:], label="y", color=RGB(.9, .4, .4))
plot!(t_dense*TU/day, u_all[3,:], label="z", color=RGB(0.4, .7, .4))
plot!(title="Control profile")
gui()

# Plot nicely
plotlyjs(legend=false)
plotTrajPlotly_indirect(XC_dense, u_all, X0_states, Xf_states, 0.2)


println("Forcing exit.")
error("Run the next part manually, if you want to continue...")
