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
if !isdefined(:DifferentialEquations)
    println("Loading packages...")

    #Correct this path for your system:
    push!(LOAD_PATH, "/Users/nathan/Documents/GitHub")

    #Packages by Nathan, from https://github.com/travelingspaceman
    using GeneralCode
    using LowThrustOpt

    #Packages from other sources:
    using MATLAB #used for plotting
    using Plots
    using DifferentialEquations
    using Ipopt
    using JuMP
    using Interpolations
    using ForwardDiff
    plotlyjs()

    println("  Done loading packages.")
end



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
    path_j = "/Users/nathan/Dropbox/CU/Research/Julia/multipleShooting/"
    X0_states = readdlm(string(path_j, "L2_Anderson_1.txt"))
    X0_times = linspace(0, 1, size(X0_states,2))
    Xf_states = readdlm(string(path_j, "L2_Anderson_2.txt"))
    Xf_times = linspace(0, 1, size(Xf_states,2))
    print("  Endpoints loaded.\n")


########################## initial guess generation ############################
# tof = 30 * day / TU #TU
#     p = 1.
#     rho = 1.
#     mass = 1e3 #kg
#     time_direction = -1
#
#
#     #-- Propagate a single trajectory with thrust:
#     #set up the ODE problem for numerical integration
#     time_direction = 1.
#     println("----------------------------")
#     thrustLimit = 0.1
#     params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
#
#     #Run numerical integration, with optimal control policy:
#     tspan = (0., tof)
#     τ1 = 0.75
#     state_0 = vcat(interpInitialStates(τ1, X0_times, X0_states, MU), zeros(6))
#     prob = ODEProblem(CRTBP_stateCostate_deriv!, state_0, tspan, params)
#     sol_forward = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)
#
#
#     #### propagate backwards
#     time_direction = -1
#
#     τ2 = 0.5
#     state_f = interpInitialStates(τ2, Xf_times, Xf_states, MU)
#
#     xf = vcat(state_f, zeros(6))
#
#     xf[4:9] *= -1
#     params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
#     prob = ODEProblem(CRTBP_stateCostate_deriv!, xf, tspan, params)
#     sol_backward = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

################################################################################




## TRAJECTORY STACKING INITIAL GUESS ###########################################
tof1 = 10 * day / TU #TU (time in first orbit)
    tof2 = 10 * day / TU #TU (time in second orbit)
    tof = tof1 + tof2 #TU
    n_nodes = 30
    t_TU = collect(linspace(0, tof, n_nodes))
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

    plotlyjs()
    plot(X0_states[1,:], X0_states[2,:], X0_states[3,:], label="Initial orbit")
    plot!(Xf_states[1,:], Xf_states[2,:], Xf_states[3,:], label="Final orbit")
    plot!(XC_trajectoryStack[1,:], XC_trajectoryStack[2,:], XC_trajectoryStack[3,:],
        color=:black, lw=2, label="Guess solution")
    scatter!([state_0[1],], [state_0[2],], [state_0[3],], label="X0", color=:green)
    scatter!([state_f[1],], [state_f[2],], [state_f[3],], label="Xf", color=:red)
    scatter!([1-MU,], [0.,], [0.,], color=:black, label="Moon")

################################################################################

## direct multiple shooting:
t_TU = collect(linspace(0, tof, n_nodes))

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
    maxIter = 5
    @time (X_all, u_all, τ1, τ2, t_TU, dV1, dV2, defect) =
        multiShoot_CRTBP_direct(X_all, u_all, τ1, τ2, t_TU, dV1, dV2, MU, DU, TU, n_nodes,
        nsteps, mass, Isp, X0_times, X0_states, Xf_times, Xf_states,
        plot_yn, flagEnd, β, allowImpulsive, maxIter)

    state_0 = interpInitialStates(τ1, X0_times, X0_states, MU)
    state_f = interpInitialStates(τ2, Xf_times, Xf_states, MU)

    plotTrajPlotly_direct(X_all, u_all, X0_states, Xf_states, 0.2)
    gui()



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
    maxIter = 5
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
    t_TU = collect(linspace(0, tof, n_nodes))

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

## Generate Neural Net training samples ########################################

σ_r = 1e2/DU #Position 1-σ errors (DU)
σ_v = 1e-3 * TU/DU #Velocity 1-σ errors (DU/TU)

N_train = 200
    m = size(XC_all,2)

    #Normally distributed:
    δX0_all = randn(6, N_train)
    δX0_all[1:3, :] *= σ_r #Position, km
    δX0_all[4:6, :] *= σ_v #velocity, km/s

    δλ_all = zeros(6, N_train * m )

    t_all = zeros(N_train * m);

    NN_input = zeros(7, N_train * m)
    NN_output = zeros(6, N_train * m);

for ind = 1:N_train
    println("ind = $ind")

    #perturb the initial state:
    XC_all_mod = copy(XC_all)
    XC_all_mod[1:6,1] += δX0_all[:,ind]

    #add a tiny amount of noise, which can help convergence in some cases
    XC_all_mod[:, 2:end-1] += 1e-10 * randn(12, n_nodes-2)

    #Optimize transfer with perturbed initial condition:
    plot_yn = false
    (XC_all_mod, defect, status_flag) = reduceFuel_indirect(XC_all_mod, t_TU,
        MU, DU, TU, n_nodes, mass, thrustLimit, rho, rho)

    if status_flag == 0 #accept the results
        NN_input[:, (ind-1)* m + 1 : ind * m] = vcat(repmat(δX0_all[:,ind], 1, m), t_TU' )
        NN_output[:, (ind-1)* m + 1 : ind * m] = XC_all_mod[7:12, :] - XC_all[7:12, :]

    else #don't accept the results
        NN_input[:, (ind-1)* m + 1 : ind * m] = NaN * ones(7, m)
        NN_output[:, (ind-1)* m + 1 : ind * m] = NaN * ones(6, m)
    end
end



## TRAIN NEURAL NET ############################################################

#Remove any NaN's (un-converged cases)
temp = .!isnan.(NN_input[1,:])
    NN_input = NN_input[:, temp]
    NN_output = NN_output[:, temp]

    #use δλ_v only:
    NN_output = NN_output[4:6,:]

println("---------------------------------------------------------")
    println("Training NN in MATLAB...")
    warn("Make sure to correct the paths in MATLAB function and in this block of code!")

    NN_size = [10, 10, 10]

    @mput NN_size #pass into MATLAB this way to avoid weird issue with the Julia/MATLAB interface

    #Correct this path for your system:
    mat"addpath('/Users/nathan/Dropbox/CU/Research/ai/')"

    identifier = "L2L2_t0"
    fcn_name = string("NN_fcn_", identifier)

    #Call MATLAB to train network:
    scaleShift = mat"TrainNN_CreateFcn($NN_input, $NN_output, NN_size, $fcn_name)"

    #The training function creates the NN function, but we need to save the scale/shift parameters separately:
    #Correct this path for your system:
    writedlm(string("/Users/nathan/Dropbox/CU/Research/ai/NN_scaleShift_", identifier, ".txt"), scaleShift)





## CALL TRAINED NEURAL NETWORK WITH MATLAB #####################################


#first, read in the scaling parameters
println("---------------------------------------------------------")
    println("Evaluating NN on new test cases...")

    #Correct this path for your system:
    scaleShift = readdlm(string("/Users/nathan/Dropbox/CU/Research/ai/NN_scaleShift_", identifier, ".txt"))
    scale_in = scaleShift[:,1]
    shift_in = scaleShift[:,2]
    scale_out = scaleShift[:,3]
    shift_out = scaleShift[:,4]

    #remove the NaN padding:
    scale_in = scale_in[.!isnan.(scale_in)]
    shift_in = shift_in[.!isnan.(shift_in)]
    scale_out = scale_out[.!isnan.(scale_out)]
    shift_out = shift_out[.!isnan.(shift_out)]

N_val = 100
    δX0_val = randn(6, N_val)
    # δX0_val = zeros(6, N_val)
    δX0_val[1:3, :] *= σ_r #Position, DU
    δX0_val[4:6, :] *= σ_v #velocity, DU/TU

    err_xf_NoCorrect = zeros(6, N_val);
    err_xf_YesCorrect = zeros(6, N_val);

    err_r_NoCorrect = zeros(n_nodes, N_val)
    err_v_NoCorrect = zeros(n_nodes, N_val)

    err_r_YesCorrect = zeros(n_nodes, N_val)
    err_v_YesCorrect = zeros(n_nodes, N_val);



#Need to evaluate the control law at 1,000 times (nom + NN update) to build an
#interpolation object for control. So first, evaluate the nominal at 1,000 times
N_dense = 1000
    time_direction = 1
    params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
    (XC_dense_nom, t_dense) = densify(XC_all, t_TU, params, N_dense)

    nstate = 6
    λv_mag_nom = norm_many(XC_dense_nom[10:12,:]);


#Set up plot for trajectories compare:
plotlyjs(legend=false, aspect_ratio=:1)
    plotTrajPlotly_indirect(XC_dense, zeros(3, size(XC_dense,2)), X0_states, Xf_states, 0.1)
    scatter!([XC_dense[1,1],], [XC_dense[2,1],], [XC_dense[3,1],], marker=2, color=:green)
    scatter!([XC_dense[1,end],], [XC_dense[2,end],], [XC_dense[3,end],], marker=2, color=:red, camera=(45, 20))

#propagate validation samples forward in time to check error at final time.
for ind = 1:1:N_val
    #Evaluate neural network for update to λ(t)
    NN_input_val = vcat( repmat(δX0_val[:,ind], 1, N_dense), t_dense[:]') #position and velocity errors
    NN_input_val_scaled = (NN_input_val - repmat(shift_in, 1, N_dense)) ./ repmat(scale_in, 1, N_dense)

    NN_output_scaled = mat"NN_fcn_L2L2_t0($NN_input_val_scaled)"

    NN_output_val = NN_output_scaled .* repmat(scale_out, 1, N_dense) + repmat(shift_out, 1, N_dense)

    X0_val = XC_all[1:6, 1] + δX0_val[:,ind]




    ## ----- WITHOUT NN CORRECTION -----
    λv_all = XC_dense_nom[10:12,:] # --- NO UPDATE
    λv_itp = scale(interpolate(λv_all', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_dense, 1:3)

    params = (MU, DU, TU, Isp, time_direction, λv_itp, thrustLimit, p, rho)
    tspan = (t_TU[1], t_TU[end])
    prob = ODEProblem(CRTBP_prop_EP_NNControl_deriv!, X0_val, tspan, params)
    sol = DifferentialEquations.solve(prob, Tsit5(), reltol=1e-13, abstol=1e-13)
    plot!(sol[1,:], sol[2,:], sol[3,:], color=:red, label="No correction")

    #calculate final state error
    err_xf_NoCorrect[:, ind] = sol[1:6, end] - XC_all[1:6, end]

    #calculate error along path
    err_r_NoCorrect[:, ind] = norm_many(sol(t_TU)[1:3,:] - XC_all[1:3,:])
    err_v_NoCorrect[:, ind] = norm_many(sol(t_TU)[4:6,:] - XC_all[4:6,:])





    ## ----- WITH NN CORRECTION -----
    λv_all = XC_dense_nom[10:12,:] + NN_output_val # -- YES UPDATE


    ## Truth update (for debugging)
    # XC_all_mod = copy(XC_all)
    # XC_all_mod[1:6,1] += δX0_val[:,ind]
    # (XC_all_mod, defect, status_flag) = reduceFuel_indirect(XC_all_mod, t_TU, MU, DU, TU, n_nodes, mass,
    #     thrustLimit, rho, rho)
    # #densify truth:
    # params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
    # (XC_dense_truth, t_dense_truth) = densify(XC_all_mod, t_TU, params, N_dense)
    # λv_all = XC_dense_truth[10:12,:] #truth


    λv_itp = scale(interpolate(λv_all', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_dense, 1:3)

    params = (MU, DU, TU, Isp, time_direction, λv_itp, thrustLimit, p, rho)
    tspan = (t_TU[1], t_TU[end])
    prob = ODEProblem(CRTBP_prop_EP_NNControl_deriv!, X0_val, tspan, params)
    sol = DifferentialEquations.solve(prob, Tsit5(), reltol=1e-13, abstol=1e-13)
    plot!(sol[1,:], sol[2,:], sol[3,:], color=:green, label="Corrected")

    #calculate final state error
    err_xf_YesCorrect[:, ind] = sol[1:6, end] - XC_all[1:6, end]

    #calculate error along path
    err_r_YesCorrect[:, ind] = norm_many(sol(t_TU)[1:3,:] - XC_all[1:3,:])
    err_v_YesCorrect[:, ind] = norm_many(sol(t_TU)[4:6,:] - XC_all[4:6,:]);

    println("ind = $ind,  err = ", norm(err_xf_YesCorrect[1:3, ind]) * DU, " km")
end
gui()
