#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Function to do indirect multiple shooting in CRTBP
#
#This version:
#   - Fixed endpoints.
#   - Impulsive maneuvers not allowed.
#   - Uses sigmoid control smoothing
#   - Variable step integration.
#   - Uses dual numbers for automatic differentiation of the Jacobian (only the
#     non-zero elements are computed). ~Same speed as finite differences,
#     but accuracy is perfect.
#   - Exploits sparsity of Jacobian for least squares solution. (speeds up
#     matrix inversion by ~10-100x).
#   - Mesh refinement not implemented yet.
#   - When close to converging, uses second-order-correction step to get a
#     "free" iteration
#
#Inputs:
#   - XC_all    = [12 x n_nodes] Initial guess array of states and costates.
#                 Each column is the state and costate at time t[i].
#                 [Nondimensional units of CRTBP]
#   - t_TU      = Guess vector of times [Nondimensional time units of CRTBP]
#   - MU        = μ for CRTBP
#   - DU        = Distance unit for CRTBP [km]
#   - TU        = Time unit for CRTBP [sec]
#   - n_nodes   = Number of nodes
#   - mass      = Initial spacecraft mass [kg]
#   - thrustLimit = Thrust limit [Newtons]
#   - plot_yn   = Boolean, true to plot solution at every iteration.
#   - flag_adjointsOnly = Boolean, true to only change the adjoints (not the
#                 states). This is helpful to run for a few iterations when we
#                 have a good guess of the states, but a poor guess of the
#                 adjoints.
#   - maxIter   = Max number of iterations.
#   - p         = Chooses the control law.
#                 With p=2, minimum energy objective.
#                 With p=1, minimum fuel objective.
#   - rho       = Smoothing term on sigmoid-smoothed control. Usually start with
#                 rho=1 or rho=0.5. Then use continuation to get down to
#                 rho=1e-3 or rho=1e-4.
#
#Outputs: (XC_all, defect, status_flag)
#   - XC_all    = [12 x n_nodes] Solution array of states and costates.
#                 Each column is the state and costate at time t[i].
#                 [Nondimensional units of CRTBP]
#   - defect    = Array of defect constraint violations.
#   - status_flag =
#                 0 for successful convergence.
#                 1 for failure.

using Printf


function multiShoot_CRTBP_indirect(XC_all, t_TU, MU, DU, TU, n_nodes,
    mass0, thrustLimit, plot_yn, flag_adjointsOnly, maxIter, p, rho)

    ## DEFINE SUB-FUNCTIONS

    function defectCalc(XC_all, t_TU, nstate, n_nodes, odefun, params)
        #Calculate defect constraint violations

        #pre-allocate defect:
        defect1 = zeros(nstate*2, n_nodes-1)

        #pre-allocate integration errors:
        errors = zeros(n_nodes-1)
        for i = 1:n_nodes-1
            #propagate node i forward from time (i) to time (i+1)

            ########### propagate forward:
            x0 = copy(XC_all[:,i])
            tspan = (t_TU[i], t_TU[i+1])

            prob = ODEProblem(odefun, x0, tspan, params)
            sol_for = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

            ########## defect calc
            defect1[:,i] = sol_for[:,end] - XC_all[:,i+1]

            #save largest error:
            errors[i] = 0.
        end

        #outputs:
        (defect1, errors)
    end


    function jacobianCalc(XC_all, t_TU, nstate, n_nodes, odefun, params)
        #Compute the Jacobian sparsely, with automatic differentation

        #Notes:
        #- The i'th defect vector is affected by XC_all[:,i]
        #- Size/shape of Jacobian is [# defects,  # variables]
        #- uses Automatic Differentation (takes a while to compile, but runtime
        #  performance is good)

        #simple function used by automatic differentiation:
        function f(XC0)
            #need to use 'eltype' because 'ForwardDiff' will call this function with
            #dual numbers and with regular numbers, and it needs to work with either
            #one.
            prob = ODEProblem(odefun, eltype(XC0).(XC0), eltype(XC0).(tspan), eltype(XC0).(params))

            #propagate forward:
            sol_for = DifferentialEquations.solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)[:,end]
        end

        #Build up non-zero elements of sparse Jacobian:
        Jac_temp = zeros(nstate*2*(n_nodes-1), 2*nstate*2)
        tspan = (t_TU[1], t_TU[2]) #elevate scope
        for i = 1:n_nodes-1
            x0 = copy(XC_all[:,i])
            tspan = (t_TU[i], t_TU[i+1])

            #automatic differentiation with dual numbers:
            Jac_AD_ForwardDiff_prop = ForwardDiff.jacobian(f,x0)

            Jac_temp[(i-1)*nstate*2 + 1 : i*nstate*2, :] = hcat(Jac_AD_ForwardDiff_prop, -1 .* I(nstate*2))
        end


        #re-create full Jacobian (band diagonal):
        Jac_full = zeros( nstate*2*(n_nodes-1), n_nodes*nstate*2 )
        for ind = 1:size(Jac_full,1) #loop through rows
            #Which node's defect are we on?
            idx = Int32(ceil(ind/(nstate*2)));

            #offset (in # of columns) for the state-related ones
            idx2 = (idx-1)*(2*nstate);

            #set the appropriate columns of Jac_full from Jac_temp:
            Jac_full[ind, (idx2+1):(idx2+2*2*nstate)] = Jac_temp[ind,:]
        end

        #Remove the first and last states:
        Jac_full[:, 1:nstate] .= 0
        Jac_full[:, (end-2*nstate+1):(end-nstate)] .= 0

        #outputs:
        Jac_full
    end


    function optimizeTraj_OLS(XC_all, t_TU, defect, Jac_full, nstate, n_nodes, odefun, params)
        #Solve the linearized problem with Ordinary Least Squares (makes use of
        #sparse Jacobian)

        #use Ordinary Least Squares to find a solution that satisfies the constraints
        #Remove the endpoints from the problem. Also remove the initial mass from
        #optimization (just leave it fixed at mass0).

        #concatenate all the defect constraints into one vector:
        if nstate == 6
            defect_vec = defect[:]
        elseif nstate == 7
            defect_vec = defect[:][1:end-1] #leave off final mass defect
        end

        #If true, we only adjust the adjoints. This is helpful in the case where we
        #do not have a good guess for the adjoints, but we do have a good guess for
        #the states. If false, adjust the states & adjoints.
        (r,c) = size(Jac_full)
        temp = trues(c)
        if flag_adjointsOnly
            #Remove the columns of Jac_full corresponding to the all states
            idx = 1:nstate
            for ind = 1:(n_nodes-1)
                temp[idx] .= false
                idx = idx .+ (2*nstate)
            end

            Jac_full = Jac_full[:, temp]
        end

        #output:
        Jac_sparse = MATLAB.sparse(Jac_full) #create sparse Jacobian (~2-3 ms)
        xc_update = -Jac_sparse \ defect_vec #sparse matrix pseudo-inverse
        xc_update2 = zeros(size(temp))
        xc_update2[temp] = xc_update
        xc_update = reshape(xc_update2, 2*nstate, n_nodes)

        #Second Order Correction step:
        #We only do the SOC when the update is small (necessary for assumptions
        #to be valid)
        if norm(xc_update, Inf) < 1e-1
            #this trick approximates using the second derivatives, but in reality
            #just uses the same first derivatives again. In some cases, it
            #essentially gives us another iteration for ~zero computational cost.
            XC_all_soc = XC_all + xc_update

            #calculate the defects at the first order step:
            (defect, errors) = defectCalc(XC_all_soc, t_TU, nstate, n_nodes, odefun, params)

            #concatenate all the defect constraints into one vector:
            if nstate == 6
                defect_vec = defect[:]
            elseif nstate == 7
                defect_vec = defect[:][1:end-1] #leave off final mass defect
            end

            #solve for SOC step:
            xc_update_soc = -Jac_sparse \ defect_vec #sparse matrix pseudo-inverse

            #apply SOC step:
            xc_update_soc2 = zeros(size(temp))
            xc_update_soc2[temp] = xc_update_soc
            xc_update_soc = reshape(xc_update_soc2, 2*nstate, n_nodes)
            xc_update += xc_update_soc
        end

        #Output:
        xc_update
    end


    function lineSearch(XC_all, xc_update, t_TU, nstate, n_nodes, odefun, params)
        #Very simple line search function. Seeks to minimize the sum of squared
        #defect constraint violations, while forcing some step to be taken.

        alpha = 1.0

        alpha_all = LinRange(0.1, 1, 20)

        er = zeros(size(alpha_all))
        cost = zeros(size(alpha_all))

        for ind = 1:length(alpha_all)
            alpha = alpha_all[ind]

            XC_all_trial = XC_all + xc_update*alpha

            #check defects
            (defect, errors) = defectCalc(XC_all_trial, t_TU, nstate, n_nodes, odefun, params)

            er[ind] = sum(defect[:].^2)
        end

        #Output:
        alpha = alpha_all[ er .== minimum(er) ]
        alpha = alpha[1]
    end



################################################################################



    #number of state elements (should be 6 or 7)
    nstate = Int32(size(XC_all, 1)/2)

    #derivatives function to use
    odefun = CRTBP_stateCostate_deriv!
    time_direction = 1.0 #forward
    params = (MU, DU, TU, thrustLimit, mass0, time_direction, p, rho)


    #default to success:
    status_flag = 0

    t0_TU = copy(t_TU[1])
    tf_TU = copy(t_TU[end])
    tau = (t_TU .- t0_TU) ./ (tf_TU - t0_TU) .* 2 .- 1 #transform t -> tau

    state_0 = XC_all[1:nstate, 1]
    state_f = XC_all[1:nstate, end]

    ############ First nominal run:
    (defect, errors) = defectCalc(XC_all, t_TU, nstate, n_nodes, odefun, params)


    #Iterate until converged
    iterCount = 0;
    er = 1.0 #initialize error
    while er > 1e-10
        iterCount += 1
        if iterCount > maxIter
            print("Reached max iteration count at ", iterCount," iterations")
            status_flag = 1
            break
        end


        ############ Compute Jacobian for dynamics constraints
        Jac_full = jacobianCalc(XC_all, t_TU, nstate, n_nodes, odefun, params)


        ############ OPTIMIZE LINEAR SUBPROBLEM
        xc_update = optimizeTraj_OLS(XC_all, t_TU, defect, Jac_full, nstate, n_nodes, odefun, params)


        ############ LINE SEARCH
        alpha = 1.0
        #Call line search:
        if iterCount > 3
            alpha = lineSearch(XC_all, xc_update, t_TU, nstate, n_nodes, odefun, params)
        end

        XC_all = XC_all + xc_update * alpha


        ############ PLOT UPDATE
        if plot_yn
            (XC_dense, t_dense) = densify(XC_all, t_TU, params, 100)

            #compute control, or just use zeros (plotting is faster when not
            #drawing control arrows)
            u_all = zeros(3, size(XC_dense,2))
            # for ind = 1:size(XC_dense,2)
            #     u_all[:,ind] = controlLaw_cart(XC_dense[nstate+4:nstate+6,ind], thrustLimit, p, rho, mass0)
            # end

            plotlyjs(legend=false)
            plotTrajPlotly_indirect(XC_dense, u_all, X0_states, Xf_states, 0.2)
            gui()
        end

        #Make sure the end states didn't change:
        XC_all[1:6, 1] = state_0
        XC_all[1:6, end] = state_f

        ############ CHECK UPDATE
        (defect, errors) = defectCalc(XC_all, t_TU, nstate, n_nodes, odefun, params)

        #print progress report
        er = norm(defect[:], Inf)
        @printf "Iter %d. Max defect = %.2e. tf = %.2f days. α = %.3f.\n" iterCount er tf_TU*TU/day alpha
        if er > 1e3
            warn("Not likely to converge. Aborting.")
            iterCount += 100
        end
    end

    if isnan(XC_all[1])
        status_flag = 2
    end

    #outputs:
    (XC_all, defect, status_flag)
end


function plotTrajPlotly_indirect(XC_all, u_all, X0_states, Xf_states, u_scale)
    #Plot in Julia with PlotlyJS package

    Plots.plot(X0_states[1,:], X0_states[2,:], X0_states[3,:], w=3,
        xlabel="X (DU)", ylabel="Y (DU)", zlabel="Z (DU)")
    Plots.plot!(Xf_states[1,:], Xf_states[2,:], Xf_states[3,:], w=3)

    Plots.plot!(XC_all[1,:],XC_all[2,:],XC_all[3,:], w=3, color=:black)

    arrows_start = XC_all[1:3,:]
    arrows_end = XC_all[1:3,:] + u_all*u_scale

    #remove zero thrust vectors:
    arrows_start = arrows_start[:, norm_many(u_all) .!== 0.]
    arrows_end   = arrows_end[:, norm_many(u_all) .!== 0.]

    arrows_x = vcat(arrows_start[1,:]', arrows_end[1,:]')
    arrows_y = vcat(arrows_start[2,:]', arrows_end[2,:]')
    arrows_z = vcat(arrows_start[3,:]', arrows_end[3,:]')

    Plots.plot!(arrows_x, arrows_y, arrows_z, color=:red)

    #invisible points to get axis scaling right
    plot!([1-MU,], [0,], [.1,], marker=1, markeralpha=0.)
    plot!([1-MU,], [0,], [-.1,], marker=1, markeralpha=0.)

    #create data points for a sphere surface:
    (x,y,z) = sphere(32)

    #scale and move the sphere:
    r_moon = 1737. #km
    x = (x .* r_moon ./ DU) .+ (1-MU)
    y = y .* r_moon./DU
    z = z .* r_moon./DU

    #plot the sphere as a surface:
    #clibrary(:cmocean)
    Plots.surface!(x,y,z, c=:gray)
end


function controlLaw_cart(lambda_v, thrustLimit, p, rho, mass)
    #Calculate control (only used for plotting)

    #Control law for indirect multiple shooting (compute the control as a function
    #of the costates).
    #
    #Valid for any Cartesian representation (inertial or rotating frame).
    #
    #Inputs:
    #     lambda_v = [3 x 1]
    #     thrustLimit
    #     p = exponent to determine which control law to use
    #             (p=2 for min energy, p=1 for min fuel)
    #Outputs:
    #     controls = [3 x 1]

    #magnitude of primer vector:
    lambda_v_mag = norm(lambda_v);

    #direction of thrust, based on primer vector:
    uhat = - lambda_v ./ lambda_v_mag

    accelLimit = thrustLimit / mass / 1e3 * TU^2 / DU #(kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration

    #Thrust magnitude from control law:
    if p == 0
      #Force thrust to always be on, in the optimal direction. This limits how big
      #of a search we need to make for the neural net stuff.
      umag = accelLimit #DU/TU^2

    elseif p == 1 #minimum mass problem
        umag = 1/2 * (1 + tanh( (lambda_v_mag - 1) / (2*rho) )) * accelLimit #DU/TU^2

    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(lambda_v))^(1 / (p-1)) #DU/TU^2

        if umag > accelLimit
          umag = accelLimit #DU/TU^2
        end
    else
      error("Invalid value of p!")
    end

    if isnan(umag)
        umag = 0.
    end

    #units of umag are acceleration (DU/TU^2)

    #Control vector, scaled properly:
    control = -umag * lambda_v ./ norm(lambda_v) * mass * DU * 1e3 / TU^2 #(DU/TU^2) -> (km/s^2) -> (kg*m/s^2) force
end
