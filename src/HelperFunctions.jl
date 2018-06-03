#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Helper functions for LowThrustOpt


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
function find_τ(X_times, X_states, state, MU)
    τ_trial = linspace(0, 1, 1001)
    d = zeros(length(τ_trial))
    for ind = 1:length(τ_trial)
        state_trial = interpInitialStates(τ_trial[ind], X_times, X_states, MU)
        d[ind] = norm(state_trial - state)
    end

    #return close approach τ:
    τ_trial[d .== minimum(d)][1]
end

#Add more intermediate points in a solution
function densify(XC_all, t_TU, params, n_desired)

    XC_dense = zeros(size(XC_all,1), 0)
    t_dense_db = zeros(0)
    t_dense = linspace(t_TU[1], t_TU[end], n_desired)
    for ind = 1:(size(XC_all, 2)-1)

        # t_save = linspace(t_TU[ind], t_TU[ind+1],
        #     Int32(round(n_desired / (size(XC_all, 2)-1))))

        #propagate and save points
        tspan = (t_TU[ind], t_TU[ind+1])
        prob = ODEProblem(CRTBP_stateCostate_deriv!, XC_all[:,ind], tspan, params)
        sol_forward = DifferentialEquations.solve(prob, Vern8(),
            reltol=1e-13, abstol=1e-13)

        # temp = t_dense[(t_dense .>= tspan[1]) .& (t_dense .< tspan[2])]
        # temp = 1
        # return temp

        if any((t_dense .>= tspan[1]) .& (t_dense .< tspan[2]))
            temp = sol_forward(t_dense[(t_dense .>= tspan[1]) .& (t_dense .< tspan[2])])[:,:]
        else
            temp = zeros(size(XC_dense,1), 0)
        end

        # if any(t_dense .== tspan[2])
        #     warn("argh")
        #     println("t overlap at $(tspan[2])")
        # end

        # display(tspan)
        # display(t_dense[(t_dense .>= tspan[1]) .& (t_dense .< tspan[2])])
        # println("---")

        # XC_dense = hcat(XC_dense, sol_forward[:,1:end-1])
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
    plot_yn = false
    flag_adjointsOnly = false

    #Check convergence on initial value of rho
    (XC_all_new, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
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
        (XC_all_new, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
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
        (XC_all_new, defect, status_flag) = multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
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
    # (XC_all_new, defect, status_flag) = reduceFuel_indirect(XC_all_new, t_TU, MU, DU, TU, n_nodes, mass,
    #     thrustLimit, rho, rho)
    (XC_all_new, defect, status_flag) = multishoot_CRTBP_indirect_fixedEnds(XC_all_new, t_TU, MU, DU, TU, n_nodes,
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
