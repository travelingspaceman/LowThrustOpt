#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Function to do direct multiple shooting in CRTBP
#
#This version:
#   - Direct method.
#   - Can optimize endpoints and time of flight.
#   - Works great for "short" transfers, but not for "long" transfers. Uses
#     linear endpoint constraints w/ quadratic endpoint penalty function
#     parameterized by β input term.
#   - Allows impulsive maneuvers at endpoints.
#   - Fixed step integrator, with mesh refinement
#
#Inputs:
#   - X_all     = [6-7 x n_nodes] Initial guess array of states. Each column is
#                 the state at time t[i]. [Nondimensional units of CRTBP]
#   - u_all     = [3 x n_nodes] Initial guess array of controls. [Newtons]
#   - τ1        = Initial guess of position on the initial orbit. In range [0-1]
#   - τ2        = Initial guess of position on the final orbit. In range [0-1]
#   - t_TU      = Guess vector of times [Nondimensional time units of CRTBP]
#   - dV1       = Initial guess of initial impulsive maneuver
#   - dV2       = Initial guess of final impulsive maneuver
#   - MU        = μ for CRTBP
#   - DU        = Distance unit for CRTBP [km]
#   - TU        = Time unit for CRTBP [sec]
#   - n_nodes   = Number of nodes
#   - nsteps    = Number of intergration steps per node
#   - mass      = Initial spacecraft mass [kg]
#   - Isp       = Specific impulse [sec]
#   - X0_times  = Vector of times on which X0_states are defined (usually
#                 linspace(0, 1, n_nodes))
#   - X0_states = [6 x n_nodes] States given to define initial orbit.
#   - Xf_times  = Vector of times on which Xf_states are defined (usually
#                 linspace(0, 1, n_nodes))
#   - Xf_states = [6 x n_nodes] States given to define final orbit.
#   - plot_yn   = Boolean, true to plot solution at every iteration.
#   - flagEnd   = Boolean, true to allow the endpoints and time of flight to
#                 vary.
#   - β         = Float to determine magnitude of quadratic endpoint cost.
#   - allowImpulsive = Boolean, true to allow impulsive maneuvers at endpoints.
#   - maxIter   = Max number of iterations.
#
#Outputs: (X_all, u_all, τ1, τ2, t_TU, dV1, dV2, defect)
#   - X_all     = [6-7 x n_nodes] Solution array of states.
#   - u_all     = [3 x n_nodes] Solution array of controls.
#   - τ1        = Solution position on the initial orbit. In range [0-1]
#   - τ2        = Solution position on the final orbit. In range [0-1]
#   - t_TU      = Solution vector of times [Nondimensional time units of CRTBP]
#   - dV1       = Initial guess of initial impulsive maneuver
#   - dV2       = Initial guess of final impulsive maneuver
#   - defect    = Differential defects at solution.



function multiShoot_CRTBP_direct(X_all, u_all, τ1, τ2, t_TU, dV1, dV2, MU, DU, TU, n_nodes,
    nsteps, mass, Isp, X0_times, X0_states, Xf_times, Xf_states, plot_yn,
    flagEnd, β, allowImpulsive, maxIter)

    global day

    ## DEFINE SUB-FUNCTIONS

    function defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
        #Calculate defect constraint violations

        #find mid-point times:
        t_mid_TU = t_TU[1:end-1] + diff(t_TU)/2

        #pre-allocate defect:
        defect1 = zeros(nstate, n_nodes-1)

        #pre-allocate integration errors:
        errors = zeros(n_nodes-1)
        for i = 1:n_nodes-1
            #propagate node i forward from time i to time between (i+1) and (i)
            #propagate node i+1 backward from time i+1 to time between (i+1) and (i)

            ########### propagate forward:
            x0 = copy(X_all[:,i])
            control = copy(u_all[:,i])
            tspan = LinRange(t_TU[i], t_mid_TU[i], nsteps)
            time_direction = 1.0 #forward
            (state_for, maxErr_for) = ode7_8(odefun, tspan, x0, MU, DU, TU, Isp, control, time_direction)

            ########### propagate backward:
            #can use the same tspan because we just need the length of time to be correct
            x0 = copy(X_all[:,i+1])
            #reverse velocity for going backward:
            x0[4:6] = -x0[4:6]
            control = copy(u_all[:,i+1]) #should this be negative? No.
            time_direction = -1.0 #backward
            (state_back, maxErr_back) = ode7_8(odefun, tspan, x0, MU, DU, TU, Isp, control, time_direction)
            stateF_back = state_back[:,end]
            #reverse velocity to be forward again:
            stateF_back[4:6] = -stateF_back[4:6]

            ########## defect between forward and backward results:
            defect1[:,i] = state_for[:,end] - stateF_back

            #save largest error:
            errors[i] = max(maxErr_for, maxErr_back)
        end

        #outputs:
        (defect1, errors)
    end

    function jacobianCalc(X_all, u_all, t_TU, defect, nstate, n_nodes, nsteps, Isp, odefun, pert)
        #Calculate Jacobian of defect constraints with respect to all states and
        #controls

        #number of variables that affect the defect between 2 nodes:
        nvar = 2*(nstate+3)

        # --- JACOBIAN (found first as a smaller, dense matrix) ---
        #The i'th defect vector is affected by X_all[:,i:i+1] and by control[:,i:i+1]

        #size of Jacobian is [# defects,    # variables]
        Jac_temp = zeros(nstate*(n_nodes-1), nvar)
        for i = 1:n_nodes-1

            XU_nom = vcat(X_all[:, i:i+1][:], u_all[:, i:i+1][:])
            for j = 1:nvar
                #create modified arrays
                XU_mod = copy(XU_nom)
                XU_mod[j] = XU_mod[j] + pert

                #reshape
                X_temp = reshape(XU_mod[1:(nstate*2)], nstate, 2)
                u_temp = reshape(XU_mod[(nstate*2)+1 : end], 3, 2)

                #compute single defect between the two affected nodes:
                (defect_temp, errors) = defectCalc(X_temp, u_temp, t_TU[i:i+1], nstate, 2, nsteps, Isp, odefun)

                #index of the rows which correspond to this defect
                idx = (i-1)*nstate+1 : i*nstate
                Jac_temp[idx,j] = (defect_temp - defect[:,i]) / pert;

            end
        end

        #re-create full Jacobian (band diagonal):
        Jac_full = zeros(nstate*(n_nodes-1), n_nodes*(nstate+3));
        for ind = 1:size(Jac_full,1) #loop through rows
            #Which node's defect are we on?
            idx = Int32(ceil(ind/(nstate)));

            #offset (in # of columns) for the state-related ones
            idx2 = (idx-1)*(nstate);

            #offset (in # of columns) for the control-related ones
            idx3 = nstate*n_nodes + (idx-1)*(3);

            #set the appropriate columns of Jac_full from Jac_temp:
            #state-related terms:
            Jac_full[ind, (idx2+1):(idx2+2*nstate)] = Jac_temp[ind,1:2*nstate];
            #control-related terms:
            Jac_full[ind, (idx3+1):(idx3+2*3)] = Jac_temp[ind,(2*nstate+1):end];
        end

        #outputs:
        Jac_full
    end

    function endpointPartials(X_inner, u_all, τ1, τ2, t_TU, dV1, dV2, Jac_XU,
        defect, X0_times, X0_states, Xf_times, Xf_states)
        #Computes partial derivatives of the defect constraints with respect to
        #the endpoint parameters (τ1 and τ2) and time of flight.

        (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states,
            Xf_times, Xf_states, MU)

        ############ partial of all defects wrt tf
        pert = 1e-5
        t0_TU = copy(t_TU[1])
        tf_TU = copy(t_TU[end])
        tf_TU_mod = t_TU[end] + pert

        tau = (t_TU - t0_TU) / (tf_TU - t0_TU) * 2 - 1 #transform t -> tau
        t_TU_mod = t0_TU + (tau+1)/2*(tf_TU_mod-t0_TU) #transform tau -> t

        (defect_mod, errors) = defectCalc(X_all, u_all, t_TU_mod, nstate,
            n_nodes, nsteps, Isp, odefun)

        ddefect_dt = (defect_mod[:] - defect[:]) / pert


        ############ partial of defect1 wrt dV1, and defectN wrt dV2
        pert = 1e-5 # (DU/TU)
        d_defect1_dV1 = zeros(nstate, 3)
        d_defectN_dV2 = zeros(nstate, 3)
        for ind = 1:3

            #apply impulsive delta-v at beginning and end
            state_0_mod = state_0
            state_0_mod[ind+3] += pert

            state_f_mod = state_f
            state_f_mod[ind+3] += pert

            ############ defect1 wrt dV1
            X_temp = hcat(state_0_mod, X_inner[:,1])
            u_temp = u_all[:, 1:2]

            (defect_temp, errors) = defectCalc(X_temp, u_temp, t_TU[1:2], nstate, 2, nsteps, Isp, odefun)

            d_defect1_dV1[:,ind] = (defect_temp - defect[:,1]) / pert


            ############ defectN wrt dV2
            X_temp = hcat(X_inner[:,end], state_f_mod)
            u_temp = u_all[:, end-1 : end]

            (defect_temp, errors) = defectCalc(X_temp, u_temp, t_TU[end-1 : end], nstate, 2, nsteps, Isp, odefun)

            d_defectN_dV2[:,ind] = (defect_temp - defect[:,end]) / pert
        end


        ############ Combine together with the rest of the Jacobian:
        n1 = (n_nodes-2)*nstate + n_nodes*3

        Jac_full = hcat(Jac_XU, zeros(nstate*(n_nodes-1), 9) )

        #partial of defect1 wrt τ1
        Jac_full[1:nstate, n1 + 1] = d_defect1_dp1

        #partial of defectN wrt τ2
        Jac_full[end-(nstate-1):end, n1 + 2] = d_defectN_dp2

        #partial of all defects wrt tf
        Jac_full[:, n1 + 3] = ddefect_dt

        #partial of defect1 wrt dV1
        Jac_full[1:nstate, (n1+4):(n1+6)] = d_defect1_dV1

        #partial of defectN wrt dV2
        Jac_full[end-(nstate-1):end, (n1+7):(n1+9)] = d_defectN_dV2


        ############ Outputs:
        Jac_full
    end

    function optimizeTraj(X_all, u_all, τ1, τ2, tf_TU, tau, dV1, dV2, defect, Jac_full,
        nstate, n_nodes, X0_times, X0_states, Xf_times, Xf_states, flagEnd, β)
        #Optimizes 1 iteration of the QP/SOCP problem.

        # global day

        #JuMP linear/quadratic constrained solver
        #
        #See link for more solvers:
        # http://www.juliaopt.org/JuMP.jl/0.18/installation.html#getting-solvers
        #Any solver capable of solving SOCP or NLP problems should work. The
        #following options have been used successfully. Gurobi is very fast, but
        #has a commercial license (free for academics). Ipopt is usable, but
        #has suffered frequently from bugs when updating Julia.
        #
        # m = Model(solver=GurobiSolver(BarConvTol=1e-10, OutputFlag=0))
        m = Model(Ipopt.Optimizer)

        #State update variables:
        @variable(m,X_jump[1:length(X_all)] )
        #constrain initial mass:
        if nstate == 7
            @constraint(m, (X_jump[7] + X_all[7, 1]) == mass)
        end

        #control update variables:
        #limits on u_jump have to be determined by u_all
        @variable(m,u_jump[1:length(u_all)])

        #endpoint parameter variables:
        if flagEnd
            d = 0.1 #unfreeze endpoints
        else
            d = 0.0 #freeze endpoints
        end
        @variable(m, (- d) <= p1_jump <= (d))
        @variable(m, (- d) <= p2_jump <= (d))

        #time of flight variable:
        if flagEnd
            d = 1 * day / TU #unfreeze endpoints
        else
            d = 0.0 #freeze endpoints
        end
        d = 0.
        @variable(m, -d <= tf_jump <= d)
        @constraint(m, (tf_TU + tf_jump) <= 40 * day/TU) #absolute max bound
        @constraint(m, (tf_TU + tf_jump) >= 0.0 * day/TU) #absolute min bound

        #impulsive maneuvers at endpoints:
        @variable(m, dV1_jump[1:3])
        @variable(m, dV2_jump[1:3])
        if !allowImpulsive
            @constraint(m, dV1_jump .== 0.) #disallow impulsive maneuvers
            @constraint(m, dV2_jump .== 0.) #
        end

        #Thrust magnitude constraints:
        #
        #Note: Thrust magnitude constraint is disabled because it was found
        #more effective to constrain thrust with the indirect method (used after
        #converging on the direct solution).
        #
        #Since u_jump is just the *update* to control, the value we want to actually
        #limit to be within constraints is (u_jump + u_all).
        #Create a vector of symbolic quadratic constraints:
        # umag = (u_jump[1:3:end] + u_all[1:3:end]).^2 +
        #     (u_jump[2:3:end] + u_all[2:3:end]).^2 +
        #     (u_jump[3:3:end] + u_all[3:3:end]).^2
        # @constraint(m, umag .<= thrustLimit^2 )


        #Need 'dt' so that unequal time steps gives the right result. If we
        #don't use 'dt' like this, then high thrust for a long segment is
        #weighted equally with high thrust for a short segment.
        t_TU_fixed = t0_TU .+ (tau.+1)./2 .*(tf_TU-t0_TU) #transform tau -> t
        dt = diff(t_TU_fixed)
        dt_temp = vcat(dt/2, dt[end]/2) + vcat(0, dt[1:end-1]/2, 0)
        dt_repeated = ones(3,1) * dt_temp' #same size as control

        #replace this by actual mass vector (repeated like dt)
        if nstate == 7
            mass_repeated = repmat(X_all[7,:], 3, 1)
        else
            mass_repeated = 1000.0 * ones(size(u_all))
        end

        #Dynamics constraints:
        #constrain bringing the defects to zero (according to linearization)
        @constraint(m, -Jac_full*vcat(X_jump, u_jump, tf_jump) .== defect[:])

        #Constraints for endpoints:

        #Finite difference derivatives of endpoints wrt parameters
        pert = 0.05 #needs to be relatively large because the endpoint orbits are only specified at 100 points (so, every 0.01 tau)
        (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)
        (state_0_mod1, state_f_mod1) = interpEndStates(τ1+pert, τ2+pert, X0_times, X0_states, Xf_times, Xf_states, MU)
        (state_0_mod2, state_f_mod2) = interpEndStates(τ1-pert, τ2-pert, X0_times, X0_states, Xf_times, Xf_states, MU)
        dstate0_dp1 = (state_0_mod1 - state_0_mod2) / (2*pert)
        dstatef_dp2 = (state_f_mod1 - state_f_mod2) / (2*pert)
        #2nd derivatives:
        ddstate0_dp1 = (state_0_mod1 - 2*state_0 + state_0_mod2) / (pert^2)
        ddstatef_dp2 = (state_f_mod1 - 2*state_f + state_f_mod2) / (pert^2)


        if flagEnd
            #linear:
            linearEnd1 = (state_0 + dstate0_dp1*(p1_jump))
            linearEnd2 = (state_f + dstatef_dp2*(p2_jump))
            @constraint(m, (X_all[1:6, 1]   + X_jump[1:6]                     + vcat(zeros(3), dV1_jump+dV1) ) - linearEnd1 .== 0 )
            @constraint(m, (X_all[1:6, end] + X_jump[(end-nstate+1):end][1:6] + vcat(zeros(3), dV2_jump+dV2) ) - linearEnd2 .== 0 )

            #quadratic:
            quadEnd1 = (state_0 + dstate0_dp1*(p1_jump) + ddstate0_dp1 / 2 * (p1_jump)^2 )
            quadEnd2 = (state_f + dstatef_dp2*(p2_jump) + ddstatef_dp2 / 2 * (p2_jump)^2 )
            c1_quad = (X_all[1:6,1]   + X_jump[1:6])                     - (state_0 + dstate0_dp1*(p1_jump) + ddstate0_dp1 / 2 * (p1_jump)^2 )
            c2_quad = (X_all[1:6,end] + X_jump[(end-nstate+1):end][1:6]) - (state_f + dstatef_dp2*(p2_jump) + ddstatef_dp2 / 2 * (p2_jump)^2 )

            #Difference from linear:
            quadCost = norm(ddstate0_dp1) * 1/2 * (p1_jump)^2 +
                       norm(ddstatef_dp2) * 1/2 * (p2_jump)^2

        else #hard fix endpoints
            quadCost = 0.
            (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)

            @constraint(m, (X_all[1:6, 1]   + X_jump[1:6]                     + vcat(zeros(3), dV1_jump+dV1) ) - state_0 .== 0 )
            @constraint(m, (X_all[1:6, end] + X_jump[(end-nstate+1):end][1:6] + vcat(zeros(3), dV2_jump+dV2) ) - state_f .== 0 )
        end

        costEnd = sum( ( (dV1 + dV1_jump) * DU/TU ).^2) + sum( ( (dV2 + dV2_jump) * DU/TU ).^2)

        #quadratic cost, with quadratic endpoint constraint violation cost added:
        cost = sum( ((u_all[:] + u_jump).^2) .* dt_repeated[:] ) + β * quadCost + costEnd

        @objective(m, Min, cost)

        #Run QP solver:
        optimize!(m)

        x_update = reshape(JuMP.value.(X_jump), nstate, n_nodes)
        u_update = reshape(JuMP.value.(u_jump), 3, n_nodes)

        p1_update = JuMP.value.(p1_jump)
        p2_update = JuMP.value.(p2_jump)
        tf_update = JuMP.value.(tf_jump)
        dV1_update = JuMP.value.(dV1_jump)
        dV2_update = JuMP.value.(dV2_jump)

        tf_TU = tf_TU + JuMP.value.(tf_jump)

        cost = JuMP.value.(cost)

        #Outputs:
        (x_update, u_update, p1_update, p2_update, tf_update, dV1_update, dV2_update, cost, dstate0_dp1, dstatef_dp2, ddstate0_dp1, ddstatef_dp2)
    end

    function lineSearch(X_all, x_update, u_all, u_update, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
        #Very simple line search function. Seeks to minimize the sum of squared
        #defect constraint violations, while forcing some step to be taken.

        alpha = 1.0

        alpha_all = LinRange(0.1, 1, 10)

        er = zeros(size(alpha_all))

        for ind = 1:length(alpha_all)
            alpha = alpha_all[ind]

            X_all_trial = X_all + x_update*alpha
            u_all_trial = u_all + u_update *alpha

            #check defects
            (defect, errors) = defectCalc(X_all_trial, u_all_trial, t_TU, nstate, n_nodes, nsteps, Isp, odefun)

            er[ind] = sum(defect[:].^2)
        end

        #Output:
        alpha = alpha_all[ er .== minimum(er) ]
        alpha = alpha[1]
    end



    function interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)
        #Interpolate the initial and final states from the given set of states

        #check bounds on endpoint parameters (wrap to 0,1)
        while τ1 > 1
            τ1 -= 1
        end
        while τ1 < 0
            τ1 += 1
        end
        while τ2 > 1
            τ2 -= 1
        end
        while τ2 < 0
            τ2 += 1
        end

        #build interpolation object:
        X0_itp = Interpolations.scale(interpolate(X0_states', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), X0_times, 1:6)
        Xf_itp = Interpolations.scale(interpolate(Xf_states', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), Xf_times, 1:6)

        #endpoints from fixed parameters given:
        state_0 = [X0_itp[τ1, 1], X0_itp[τ1, 2], X0_itp[τ1, 3], X0_itp[τ1, 4], X0_itp[τ1, 5], X0_itp[τ1, 6]]
        state_f = [Xf_itp[τ2, 1], Xf_itp[τ2, 2], Xf_itp[τ2, 3], Xf_itp[τ2, 4], Xf_itp[τ2, 5], Xf_itp[τ2, 6]]

        #outputs:
        (state_0, state_f)
    end



################################################################################
    ## Done defining sub-functions. Now, run the algorithm.

    #number of state elements (should be 6 or 7)
    nstate = size(X_all, 1)

    #derivatives function to use
    odefun = CRTBP_prop_EP_deriv

    pert = 1e-8;

    #find mid-point times:
    t_mid_TU = t_TU[1:end-1] + diff(t_TU)/2

    t0_TU = copy(t_TU[1])
    tf_TU = copy(t_TU[end])
    tau = (t_TU .- t0_TU) ./ (tf_TU .- t0_TU) .* 2 .- 1; #transform t -> tau

    (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)

    ############ First nominal run:
    (defect, errors) = defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)

    #iterate!
    iterCount = 0;
    er = 1.0 #initialize error
    while er > 1e-6
        iterCount += 1
        if iterCount > maxIter
            print("Reached max iteration count at ", iterCount," iterations")
            break
        end


        ############ Compute Jacobian for dynamics constraints
        Jac_full = jacobianCalc(X_all, u_all, t_TU, defect, nstate, n_nodes, nsteps, Isp, odefun, pert)


        #Partial wrt tf
        pert_tf = 1e-3
        t0_TU = copy(t_TU[1])
        tf_TU = copy(t_TU[end])
        tf_TU_mod1 = tf_TU + pert_tf
        tf_TU_mod2 = tf_TU - pert_tf
        t_TU_mod1 = t0_TU .+ (tau.+1)./2 .*(tf_TU_mod1-t0_TU) #transform tau -> t
        t_TU_mod2 = t0_TU .+ (tau.+1)./2 .*(tf_TU_mod2-t0_TU) #transform tau -> t
        (defect_mod1, errors) = defectCalc(X_all, u_all, t_TU_mod1, nstate, n_nodes, nsteps, Isp, odefun)
        (defect_mod2, errors) = defectCalc(X_all, u_all, t_TU_mod2, nstate, n_nodes, nsteps, Isp, odefun)

        ddefect_dt = (defect_mod1 - defect_mod2) / (2*pert_tf)

        Jac_full = hcat(Jac_full, ddefect_dt[:])



        ############ OPTIMIZE QP SUBPROBLEM

        #every other iteration, allow the endpoints to vary
        flagEnd_temp = copy(flagEnd)
        if mod(iterCount,2) == 0
            flagEnd_temp = false * flagEnd
        end

        (x_update, u_update, p1_update, p2_update, tf_update, dV1_update, dV2_update,
             cost, dstate0_dp1, dstatef_dp2, ddstate0_dp1, ddstatef_dp2) =
            optimizeTraj(X_all, u_all, τ1, τ2, tf_TU, tau, dV1, dV2, defect, Jac_full,
            nstate, n_nodes, X0_times, X0_states, Xf_times, Xf_states, flagEnd_temp, β)


        #Plot endpoint constraint lines:
        #check the linear/quadratic approximations of the initial and final orbits:
        (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)

        num_test = 30
        x1_1 = zeros(6, num_test)
        x2_1 = zeros(6, num_test)
        x1_2 = zeros(6, num_test)
        x2_2 = zeros(6, num_test)
        temp = LinRange(-1, 1, num_test)
        for ind = 1:num_test
            #linear:
            x1_1[:,ind] = (state_0 + dstate0_dp1*(p1_update)*temp[ind])
            x2_1[:,ind] = (state_f + dstatef_dp2*(p2_update)*temp[ind])

            #quadratic:
            x1_2[:,ind] = (state_0 + dstate0_dp1*(p1_update)*temp[ind] + ddstate0_dp1/2*( (p1_update)*temp[ind] )^2)
            x2_2[:,ind] = (state_f + dstatef_dp2*(p2_update)*temp[ind] + ddstatef_dp2/2*( (p2_update)*temp[ind] )^2)
        end



        ############ LINE SEARCH
        alpha = 1.0
        #Call line search:
        if iterCount > 10
            alpha = lineSearch(X_all, x_update, u_all, u_update, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
        end

        X_all = X_all + x_update * alpha
        u_all = u_all + u_update * alpha
        τ1 = τ1 + p1_update * alpha
        τ2 = τ2 + p2_update * alpha
        tf_TU = tf_TU + tf_update * alpha
        dV1 = dV1 + dV1_update * alpha
        dV2 = dV2 + dV2_update * alpha


        ############ PLOT UPDATE
        if plot_yn
            #Can plot either with Julia (slow) or with MATLAB (less slow)

            #Plot with Julia:
            plotTrajPlotly_direct(X_all, zeros(size(u_all)), X0_states, Xf_states, 0.2)
            gui()
        end


        t_TU = t0_TU .+ (tau.+1)./2 .*(tf_TU-t0_TU) #transform tau -> t

        ############ CHECK UPDATE
        (defect, errors) = defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)

        #print progress report
        er = norm(defect[:], Inf)
        @printf "Iter %d. Max defect = %.2e. Cost = %.5f. tf = %.2f days. α = %.3f.\n" iterCount er cost tf_TU*TU/day alpha
    end

    #outputs:
    (X_all, u_all, τ1, τ2, t_TU, dV1, dV2, defect)
end


function meshRefine_direct(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
    ############ MESH REFINEMENT

    (defect, errors) = defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)

    #These can be adjusted up or down depending on the problem.
    tol_min = 1e-20
    tol_max = 1e-18

    println("Starting with ",n_nodes," nodes.")
    n_nodes_old = n_nodes

    # First, remove extra points
    while minimum(errors) < tol_min
        #find the index of the minimum & maximum errors
        idx_min = find(errors .== minimum(errors))[1]

        #ensure that we dont remove the endpoints with this.
        if idx_min == 1
            idx_min = 2
        end

        #Remove the node because it's just making things messy:
        X_all = hcat(X_all[:, 1:(idx_min-1)], X_all[:, (idx_min+1):end])
        u_all = hcat(u_all[:, 1:(idx_min-1)], u_all[:, (idx_min+1):end])
        t_TU = vcat(t_TU[1:(idx_min-1)], t_TU[(idx_min+1):end])
        n_nodes = n_nodes - 1;

        #Update the defects & error estimates:
        (defect, errors) = defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
    end

    if ~(n_nodes == n_nodes_old)
        println("Removed nodes, now n_nodes = ", n_nodes)
    end

    n_nodes_old2 = n_nodes

    #Add more points where needed:
    while maximum(errors) > tol_max

        #find the index of the minimum & maximum errors
        idx_max = find(errors .== maximum(errors))[1]


        #find mid-point times:
        t_mid_TU = t_TU[1:end-1] + diff(t_TU)/2

        #new time to add:
        t_new = t_TU[idx_max] + diff( t_TU[idx_max:idx_max+1] )[1]/2

        #Propagate to the midpoint:
        x0 = X_all[:,idx_max]
        control = u_all[:,idx_max]
        tspan = [t_TU[idx_max],    t_new ]
        time_direction = 1.0 #forward
        state_new = ode7(odefun, tspan, x0, MU, DU, TU, Isp, control, time_direction)

        #new control is the average of the control on either side:
        u_new = (u_all[:, idx_max] + u_all[:, idx_max+1]) / 2

        #add a node because it's not good enough
        X_all = hcat(X_all[:, 1:idx_max], state_new[:,2], X_all[:,(idx_max+1:end)])
        u_all = hcat(u_all[:, 1:idx_max], u_new, u_all[:,(idx_max+1:end)])
        t_TU    = vcat(t_TU[1:idx_max], t_new, t_TU[idx_max+1:end])
        n_nodes = n_nodes + 1;

        #Update the defects & error estimates:
        (defect, errors) = defectCalc(X_all, u_all, t_TU, nstate, n_nodes, nsteps, Isp, odefun)
    end

    if !(n_nodes == n_nodes_old2)
        println("Added nodes, now n_nodes = ", n_nodes)
    end
    if n_nodes == n_nodes_old
        println("Did not need to add or remove nodes.")
    else
        println("Refined the mesh. Now have ",n_nodes," nodes.")
    end


    #outputs:
    (X_all, u_all, t_TU, n_nodes)
end

function plotTrajPlotly_direct(X_all, u_all, X0_states, Xf_states, u_scale)
    #Plot trajectory with PlotlyJS package

    Plots.plot(X0_states[1,:], X0_states[2,:], X0_states[3,:], w=3,
        xlabel="X (DU)", ylabel="Y (DU)", zlabel="Z (DU)",
        label="X₀")
    Plots.plot!(Xf_states[1,:], Xf_states[2,:], Xf_states[3,:], w=2, label="Xf")
    Plots.plot!([1-MU], [0.], [0.], m=(0,2,:black) , label="Moon")
    Plots.plot!(X_all[1,:],X_all[2,:],X_all[3,:], w=3, color=:black, label="Transfer trajectory")

    arrows_start = X_all[1:3,:]
    arrows_end = X_all[1:3,:] + u_all*u_scale

    #remove zero thrust vectors:
    arrows_start = arrows_start[:, norm_many(u_all) .!== 0.]
    arrows_end   = arrows_end[:, norm_many(u_all) .!== 0.]

    arrows_x = vcat(arrows_start[1,:]', arrows_end[1,:]')
    arrows_y = vcat(arrows_start[2,:]', arrows_end[2,:]')
    arrows_z = vcat(arrows_start[3,:]', arrows_end[3,:]')


    Plots.plot!(arrows_x, arrows_y, arrows_z, color=:red, show=true)



    #create data points for a sphere surface:
    (x,y,z) = sphere(32)

    #scale and move the sphere:
    r_moon = 1737. #km
    x = (x .* r_moon./DU) .+ (1-MU)
    y = y .* r_moon./DU
    z = z .* r_moon./DU

    #plot the sphere as a surface:
    Plots.surface!(x,y,z, c=:blues)
end
