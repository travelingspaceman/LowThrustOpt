#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Demonstration of trajectory optimization:
#   - Initial guess generation
#   - Optimization with direct method (robust)
#   - Optimization with indirect method (sensitive but fast)

#On the first run, import all the packages
if !isdefined(:DifferentialEquations)

    using GeneralCode
    using MATLAB #used by direct MS
    using DifferentialEquations #used by indirect MS
    using Plots
    using JuMP #used by direct MS
    using Ipopt #used by direct MS
    using Interpolations
    using ForwardDiff  #used by indirect MS
    using DualNumbers #used by indirect MS

    plotlyjs()

    path_j = "/Users/nathan/Dropbox/CU/Research/Julia/multipleShooting/"

    include(string(path_j, "CRTBP_prop_EP_deriv.jl")) #used by direct MS
    include(string(path_j, "multiShoot_CRTBP_direct.jl")) #direct MS

    include("/Users/nathan/Dropbox/CU/Research/ai/CRTBP_stateCostate_deriv.jl") #used by indirect MS
    include(string(path_j, "multishoot_CRTBP_indirect.jl")) #indirect MS
end

#Simple function to calculate Jacobi constant
function jacobiConstant(state, MU, DU)
    r1 = sqrt.((state[1,:]+MU).^2 + state[2,:].^2 + state[3,:].^2)
    r2 = sqrt.((state[1,:]+MU-1).^2 + state[2,:].^2 + state[3,:].^2)
    v = norm_many(state[4:6,:])
    C = state[1,:].^2 + state[2,:].^2 + 2*(1-MU)./r1 + 2*MU./r2 - v.^2
end

#Spline interpolation for initial and final states as function of τ
function interpInitialStates(p1, X0_times, X0_states, MU)
  #check bounds on endpoint parameters (wrap to 0,1)
  while p1 > 1
    p1 -= 1
  end
  while p1 < 0
    p1 += 1
  end

  #build interpolation object:
  X0_itp = Interpolations.scale(interpolate(X0_states', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), X0_times, 1:6)

  #endpoints from fixed parameters given:
  state_0 = [X0_itp[p1, 1], X0_itp[p1, 2], X0_itp[p1, 3], X0_itp[p1, 4], X0_itp[p1, 5], X0_itp[p1, 6]]

  #outputs:
  state_0
end

#finds the value of τ that matches the X_states most closely to a new endpoint
function find_τ(X_times, X_states, state)
    τ_trial = linspace(0, 1, 1001)
    d = zeros(length(τ_trial))
    for ind = 1:length(τ_trial)
        state_trial = interpInitialStates(τ_trial[ind], X_times, X_states, MU)
        d[ind] = norm(state_trial - state)
    end

    #return close approach τ:
    τ_trial[d .== minimum(d)][1]
end

#Add more intermediate points in a solution (mainly for plotting)
function densify(XC_all, t_TU, params, n_desired)

    XC_dense = zeros(size(XC_all,1), 0)
    t_dense_db = zeros(0)
    t_dense = linspace(t_TU[1], t_TU[end], n_desired)
    for ind = 1:(size(XC_all, 2)-1)

        #propagate and save points
        tspan = (t_TU[ind], t_TU[ind+1])
        prob = ODEProblem(CRTBP_stateCostate_deriv!, XC_all[:,ind], tspan, params)
        sol_forward = DifferentialEquations.solve(prob, Vern8(),
            reltol=1e-13, abstol=1e-13)

        if any((t_dense .>= tspan[1]) .& (t_dense .< tspan[2]))
            temp = sol_forward(t_dense[(t_dense .>= tspan[1]) .& (t_dense .< tspan[2])])[:,:]
        else
            temp = zeros(size(XC_dense,1), 0)
        end

        XC_dense = hcat(XC_dense, temp)
        t_dense_db = append!(t_dense_db, t_dense[(t_dense .>= tspan[1]) .& (t_dense .< tspan[2])])

        #add on the last one:
        if ind == (size(XC_all, 2)-1)
            XC_dense = hcat(XC_dense, sol_forward[:,end])
            t_dense_db = append!(t_dense_db, t_TU[end]) #debug
        end
    end

    (XC_dense, t_dense)
end


#Reduce rho algorithmically with successive indirect multiple shooting
function reduceFuel_indirect(XC_all, t_TU, MU, DU, TU, n_nodes, mass,
    thrustLimit, rho_current, rho_target)

    #The whole point is to reduce rho down to rho_target. Here, catch mistakes:
    if rho_target > rho_current
        warn("Non-sensical rho_target. Fixing...")
        rho_target = rho_current
    end

    #initial values. p is always 1, but rho gets decreased
    p = 1.
    rho_temp = rho_current
    maxIter = 10

    #Check convergence on initial value of rho
    (XC_all_new, defect, status_flag) = multishoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
        mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho_temp)

    if (status_flag == 0) && (rho_current == rho_target)
        #We're golden!
        println("   multiShoot converged at initial rho = $rho_temp")

        return (XC_all_new, defect, status_flag)
    end

    #while not converged, increase rho
    while !(status_flag == 0) && (rho_temp < 1)
        #radius of convergence is larger for larger rho

        println("   Did not converge with rho = $rho_temp. Increasing rho...")
        rho_temp *= 5
        if rho_temp > 1.
            rho_temp = 1.
        end
        (XC_all_new, defect, status_flag) = multishoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
            mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho_temp)
    end

    if (rho_temp == 1) && !(status_flag == 0)
        warn("  Did not converge at rho = 1. Stopping...")
        return (XC_all_new, defect, status_flag)
    end
    if status_flag == 0
        println("   multiShoot converged at rho = $rho_temp")

        XC_all = copy(XC_all_new)
    end


    #Now, decrease rho again, if needed
    count = 0
    while (rho_temp > rho_target) || !(status_flag == 0)
        count += 1
        if count > 100
            status_flag = 3
            warn("  Did not converge after $count iterations!")

            # XC_all_new *= NaN
            return (XC_all_new, defect, status_flag)
        end


        #accept the new solution if it converged (need to check the defect magnitudes instead of this)
        if status_flag == 0
            XC_all = copy(XC_all_new)

            rho_temp /= 2

            if rho_temp < rho_target
                rho_temp = rho_target
            end
            println("   Reducing rho, new value:   rho = ", rho_temp)

        else
            rho_temp *= 3 * (1 + rand())
            println("   Increasing rho, new value: rho = ", rho_temp)
        end

        #Run single-shooting algorithm again:
        (XC_all_new, defect, status_flag) = multishoot_CRTBP_indirect_fixedEnds(XC_all, t_TU, MU, DU, TU, n_nodes,
            mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho_temp)
    end


    #Outputs:
    (XC_all_new, defect, status_flag)
end

#Add time to the final endpoint
function addTimeFinal(XC_all, t_TU, Δt)

    #-propagate end state forward in time.
    XC_all[7:12,end] = 0.
        tspan_edge = (t_TU[end], t_TU[end] + Δt)
        prob = ODEProblem(CRTBP_stateCostate_deriv!, XC_all[:,end], tspan, params)
        sol = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

        XC_all = hcat(XC_all, sol[:, end])
        t_TU = append!(t_TU, tspan_edge[2])

    #-densify, then interpolate on the densified points to find the new XC_all
    n_desired = 200
        (XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)

    #build interpolation object:
    XC_itp = Interpolations.scale(interpolate(XC_dense', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_dense, 1:12)

    #endpoints from fixed parameters given:
    #(This is also an opportunity to change the number / distribution of nodes.)
    # n_nodes = 10
    t_TU_new = collect(linspace(t_TU[1], t_TU[end], n_nodes))
    XC_all_new = zeros(12, n_nodes)
    for ind = 1:n_nodes
        XC_all_new[:,ind] = [XC_itp[t_TU_new[ind], 1], XC_itp[t_TU_new[ind], 2], XC_itp[t_TU_new[ind], 3],
            XC_itp[t_TU_new[ind], 4], XC_itp[t_TU_new[ind], 5], XC_itp[t_TU_new[ind], 6],
            XC_itp[t_TU_new[ind], 7], XC_itp[t_TU_new[ind], 8], XC_itp[t_TU_new[ind], 9],
            XC_itp[t_TU_new[ind], 10], XC_itp[t_TU_new[ind], 11], XC_itp[t_TU_new[ind], 12]]
    end

    #-call the function to find the value of τ that matches the Xf_states most closely to the new endpoint
    τ = find_τ(Xf_times, Xf_states, XC_all_new[1:6, end])

    #-Set XC_all[1:6, end] equal to that interpolated Xf_states value
    XC_all_new[1:6, end] = interpInitialStates(τ, Xf_times, Xf_states, MU)

    #-Re-solve indirect problem and repeat
    plot_yn = true
    (XC_all_new, defect, status_flag) = multishoot_CRTBP_indirect(XC_all_new, t_TU, MU, DU, TU, n_nodes,
        mass, thrustLimit, false, false, maxIter, p, rho)

    if status_flag == 0
        println("Successfully added time!")
        XC_all = copy(XC_all_new)
        t_TU = copy(t_TU_new)
    else
        warn("Failed! Returning original values")
        XC_all = XC_all[:, 1:n_nodes]
        t_TU = t_TU[1:n_nodes]
    end

    (XC_all, t_TU)
end


#Earth-Moon system:
MU = 0.012150585609624037 # MU = mu_moon/(mu_moon + mu_planet)
DU = 384747.96285603708 #Earth-Moon distance (km) (km per distance unit)
TU = 375699.81732246041 #time units (seconds per time unit)
r_moon = 1737. #km
r_earth = 6378. #km
day = 86400. #seconds
##Calculate mu_moon based on MU for consistency with other code:
mu_planet = 398600.4415
mu_moon = (MU * mu_planet) / (1 - MU)


#Read in initial & final states:
print("Reading endpoints files...")
    X0_states = readdlm(string(path_j, "L2_states.txt"))
    X0_times = linspace(0, 1, size(X0_states,2))
    Xf_states = readdlm(string(path_j, "NRHO.txt"))
    Xf_times = linspace(0, 1, size(Xf_states,2))
    print("  Endpoints loaded.\n")


########################## initial guess generation ############################
tof = 30 * day / TU #TU
    p = 1.
    rho = 1.
    mass = 1e3 #kg
    time_direction = -1


    #-- Propagate a single trajectory with thrust:
    #set up the ODE problem for numerical integration
    time_direction = 1.
    println("----------------------------")
    thrustLimit = 0.5
    params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)

    #Run numerical integration, with optimal control policy:
    tspan = (0., tof)
    τ1 = 0.75
    state_0 = vcat(interpInitialStates(τ1, X0_times, X0_states, MU), zeros(6))
    prob = ODEProblem(CRTBP_stateCostate_deriv!, state_0, tspan, params)
    sol_forward = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)


    #### propagate backwards
    time_direction = -1

    τ2 = 0.5
    state_f = interpInitialStates(τ2, Xf_times, Xf_states, MU)

    xf = vcat(state_f, zeros(6))

    xf[4:9] *= -1
    params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
    prob = ODEProblem(CRTBP_stateCostate_deriv!, xf, tspan, params)
    sol_backward = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

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
    τ2₀ = .75
    state_f₀ = vcat(interpInitialStates(τ2₀, Xf_times, Xf_states, MU), zeros(6))

    #Ballistically propagate second orbit:
    tspan = (tof1, tof1 + tof2)
    prob = ODEProblem(CRTBP_stateCostate_deriv!, state_f₀, tspan, params)
    sol_orbit2 = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

    #Fix final state:
    τ2 = find_τ(Xf_times, Xf_states, sol_orbit2[1:6, end])
    τ2 = 0.5
    state_f = vcat(interpInitialStates(τ2, Xf_times, Xf_states, MU), zeros(6))

    #Concatenate states
    XC_trajectoryStack = hcat(sol_orbit1(t_TU1)[:,:], sol_orbit2(t_TU2)[:,:] )
    XC_trajectoryStack[:,1] = state_0
    XC_trajectoryStack[:,end] = state_f

    plotlyjs()
    plot(X0_states[1,:], X0_states[2,:], X0_states[3,:])
    plot!(Xf_states[1,:], Xf_states[2,:], Xf_states[3,:])
    plot!(XC_trajectoryStack[1,:], XC_trajectoryStack[2,:], XC_trajectoryStack[3,:], color=:black, lw=2)
    scatter!([state_0[1],], [state_0[2],], [state_0[3],], label="X0")
    scatter!([state_f[1],], [state_f[2],], [state_f[3],], label="Xf")
    scatter!([sol_orbit1[1,end],], [sol_orbit1[2,end],], [sol_orbit1[3,end],],
        label="Final state from guess", show=true)
    scatter!([1-MU,], [0.,], [0.,], color=:black, label="Moon")

################################################################################

## direct multiple shooting:
n_nodes = 30
    t_TU = collect(linspace(0, tof, n_nodes))

    X_all = 0.5*randn(6, n_nodes)
    X_all[1,:] += (1-MU)

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
    plot_yn = true
    # gr(legend=false)
    plotlyjs(legend=false)
    maxIter = 20
    @time (X_all, u_all, τ1, τ2, t_TU, dV1, dV2, defect, Jac_full) =
        multiShoot_CRTBP_v7a(X_all, u_all, τ1, τ2, t_TU, dV1, dV2, MU, DU, TU, n_nodes,
        nsteps, mass, Isp, thrustLimit, X0_times, X0_states, Xf_times, Xf_states,
        plot_yn, flagEnd, β, allowImpulsive, maxIter)

    state_0 = interpInitialStates(τ1, X0_times, X0_states, MU)
    state_f = interpInitialStates(τ2, Xf_times, Xf_states, MU)

    plotTrajMATLAB(X_all, u_all, X0_states, Xf_states)


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
    @time (XC_all, defect, status_flag) = multishoot_CRTBP_indirect_fixedEnds(XC_all, t_TU, MU, DU, TU, n_nodes,
        mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)


    flag_adjointsOnly = false
    maxIter = 50
    @time (XC_all, defect, status_flag) = multishoot_CRTBP_indirect_fixedEnds(XC_all, t_TU, MU, DU, TU, n_nodes,
        mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)

    u_all = zeros(3, n_nodes)
    for ind = 1:n_nodes
        u_all[:,ind] = controlLaw_cart(XC_all[10:12,ind], thrustLimit, p, rho)
    end



# "densify" (add more points) for plotting
n_desired = 300
    time_direction = 1
    params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
    (XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)
    u_all = zeros(3, size(XC_dense, 2))
    for ind = 1:size(XC_dense, 2)
        u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho)
    end


#plot control profile
plotlyjs(legend=true, aspect_ratio=:none)
    plot(t_dense*TU/day, norm_many(u_all), xlabel="Time [days]", ylabel="Control [N]",
        lw=2, color=:black, label="Control Magnitude", size=(800, 400))
    plot!(t_dense*TU/day, u_all[1,:], label="x", color=RGB(0,.5, .9))
    plot!(t_dense*TU/day, u_all[2,:], label="y", color=RGB(.9, .4, .4))
    plot!(t_dense*TU/day, u_all[3,:], label="z", color=RGB(0.4, .7, .4))


# Plot nicely
mat"""
    figure

    %subplot(1,2,1)
    hold all
    plot3($X0_states(1,:), $X0_states(2,:), $X0_states(3,:), 'linewidth', 2)
    plot3($Xf_states(1,:), $Xf_states(2,:), $Xf_states(3,:), 'linewidth', 2)
    plot3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), 'k', 'linewidth', 2)
    quiver3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), $u_all(1,:), $u_all(2,:), $u_all(3,:),'r')
    axis equal
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_moon/$DU;
    x2 = earth_scale*x2 +1 - $MU; y2 = earth_scale*y2; z2 = earth_scale*z2;
    moon = surf(x2,y2,z2);
    set(moon,'FaceColor',[.0,.5,0.5],'FaceAlpha',0.7,'EdgeAlpha',0)
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_earth/$DU;
    x2 = earth_scale*x2; y2 = earth_scale*y2; z2 = earth_scale*z2;
    %earth = surf(x2,y2,z2);
    %set(earth,'FaceColor',[.2,.5,0.5],'FaceAlpha',0.7,'EdgeAlpha',0)
    rotate3d('on')
    xlabel('X [DU]')
    ylabel('Y [DU]')
    zlabel('Z [DU]')
    axis tight
    grid on
    view(3)
    title('Minimum Energy')
    """

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

    thrustLimit = 0.3
    println("thrustLimit = $thrustLimit N")
    rho = 1
    maxIter = 30
    flag_adjointsOnly = false
    plot_yn = false
    @time (XC_all, defect, status_flag) = multishoot_CRTBP_indirect_fixedEnds(XC_all, t_TU, MU, DU, TU, n_nodes,
        mass, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)


# "densify" (add more points) for plotting
n_desired = 300
    time_direction = 1
    params = (MU, DU, TU, thrustLimit, mass, time_direction, p, rho)
    (XC_dense, t_dense) = densify(XC_all, t_TU, params, n_desired)
    u_all = zeros(3, size(XC_dense, 2))
    for ind = 1:size(XC_dense, 2)
        u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho)
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
mat"""
    figure

    %subplot(1,2,1)
    hold all
    plot3($X0_states(1,:), $X0_states(2,:), $X0_states(3,:), 'linewidth', 2)
    plot3($Xf_states(1,:), $Xf_states(2,:), $Xf_states(3,:), 'linewidth', 2)
    plot3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), 'k', 'linewidth', 2)
    quiver3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), $u_all(1,:), $u_all(2,:), $u_all(3,:),'r')
    axis equal
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_moon/$DU;
    x2 = earth_scale*x2 +1 - $MU; y2 = earth_scale*y2; z2 = earth_scale*z2;
    moon = surf(x2,y2,z2);
    set(moon,'FaceColor',[.0,.5,0.5],'FaceAlpha',0.7,'EdgeAlpha',0)
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_earth/$DU;
    x2 = earth_scale*x2; y2 = earth_scale*y2; z2 = earth_scale*z2;
    %earth = surf(x2,y2,z2);
    %set(earth,'FaceColor',[.2,.5,0.5],'FaceAlpha',0.7,'EdgeAlpha',0)
    rotate3d('on')
    xlabel('X [DU]')
    ylabel('Y [DU]')
    zlabel('Z [DU]')
    axis tight
    grid on
    view(3)
    title('Smoothed Minimum Fuel')
    """

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
        u_all[:,ind] = controlLaw_cart(XC_dense[10:12,ind], thrustLimit, p, rho)
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
mat"""
    figure

    hold all
    plot3($X0_states(1,:), $X0_states(2,:), $X0_states(3,:), 'linewidth', 2)
    plot3($Xf_states(1,:), $Xf_states(2,:), $Xf_states(3,:), 'linewidth', 2)
    plot3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), 'k', 'linewidth', 2)
    quiver3($XC_dense(1,:), $XC_dense(2,:), $XC_dense(3,:), $u_all(1,:), $u_all(2,:), $u_all(3,:),'r')
    axis equal
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_moon/$DU;
    x2 = earth_scale*x2 +1 - $MU; y2 = earth_scale*y2; z2 = earth_scale*z2;
    moon = surf(x2,y2,z2);
    set(moon,'FaceColor',[.0,.5,0.5],'FaceAlpha',0.7,'EdgeAlpha',0)
    [x2,y2,z2] = sphere(20);
    earth_scale = $r_earth/$DU;
    x2 = earth_scale*x2; y2 = earth_scale*y2; z2 = earth_scale*z2;
    rotate3d('on')
    xlabel('X [DU]')
    ylabel('Y [DU]')
    zlabel('Z [DU]')
    axis tight
    grid on
    view(3)
    title('Bang-Bang Minimum Fuel')
    """
