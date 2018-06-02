#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#A set of Ordinary Differential Equations solvers of varying order. These are
#mostly made obsolete by 'DifferentialEquations.jl', but they may be useful in
#the future.
#
#ode7_8 is still used in my multiple shooting implementation.
#
#Included:
#   ode4            - fixed step, 4th order
#   ode5            - fixed step, 5th order
#   ode7            - fixed step, 7th order
#   ode78           - variable step, 7th order with 8th order checking
#   ode78_events    - variable step, 7th order with 8th order checking, and
#                     events functions allowed
#   ode7_8          - fixed step, 7th order with 8th order checking

function ode4(odefun,tspan,y0, varargin...)
  #ODE4  Solve differential equations with a non-adaptive method of order 4.
  #   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates
  #   the system of differential equations y' = f(t,y) by stepping from T0 to
  #   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
  #   The vector Y0 is the initial conditions at T0. Each row in the solution
  #   array Y corresponds to a time specified in TSPAN.
  #
  #   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters
  #   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...).
  #
  #   This is a non-adaptive solver. The step sequence is determined by TSPAN
  #   but the derivative function ODEFUN is evaluated multiple times per step.
  #   The solver implements the classical Runge-Kutta method of order 4.
  #
  #   Example
  #         tspan = 0:0.1:20;
  #         y = ode4(@vdp1,tspan,[2 0]);
  #         plot(tspan,y(:,1));
  #     solves the system y' = vdp1(t,y) with a constant step size of 0.1,
  #     and plots the first component of the solution.
  #


  h = diff(tspan);

  f0 = odefun(tspan[1], y0, varargin...)

  # y0 = y0(:);   # Make a column vector.
  # if ~isequal(size(y0),size(f0))
  #   error('Inconsistent sizes of Y0 and f(t0,y0).');
  # end

  neq = length(y0);
  N = length(tspan);
  Y = zeros(neq,N);
  F = zeros(neq,4);

  Y[:,1] = y0;
  for i = 2:N
    ti = tspan[i-1];
    hi = h[i-1];
    yi = Y[:,i-1];
    F[:,1] = odefun(ti,yi,varargin...);
    F[:,2] = odefun(ti+0.5*hi,yi+0.5*hi*F[:,1],varargin...);
    F[:,3] = odefun(ti+0.5*hi,yi+0.5*hi*F[:,2],varargin...);
    F[:,4] = odefun(tspan[i],yi+hi*F[:,3],varargin...);
    Y[:,i] = yi + (hi/6)*(F[:,1] + 2*F[:,2] + 2*F[:,3] + F[:,4]);
  end
  # Y = Y.';

  Y
end


function ode5(odefun,tspan,y0, varargin...)
  #ODE5  Solve differential equations with a non-adaptive method of order 5.
  #   Y = ODE5(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates
  #   the system of differential equations y' = f(t,y) by stepping from T0 to
  #   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
  #   The vector Y0 is the initial conditions at T0. Each row in the solution
  #   array Y corresponds to a time specified in TSPAN.
  #
  #   Y = ODE5(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters
  #   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...).
  #
  #   This is a non-adaptive solver. The step sequence is determined by TSPAN
  #   but the derivative function ODEFUN is evaluated multiple times per step.
  #   The solver implements the Dormand-Prince method of order 5 in a general
  #   framework of explicit Runge-Kutta methods.
  #
  #   Example
  #         tspan = 0:0.1:20;
  #         y = ode5(@vdp1,tspan,[2 0]);
  #         plot(tspan,y(:,1));
  #     solves the system y' = vdp1(t,y) with a constant step size of 0.1,
  #     and plots the first component of the solution.

  #There are no checks on the inputs! So be careful.

  h = diff(tspan);

  #Evaluate at initial time
  f0 = odefun(tspan[1], y0, varargin...)
  # f0 = feval(odefun,tspan(1),y0,varargin{:});

  neq = length(y0);
  N = length(tspan);
  Y = zeros(neq,N);

  # Method coefficients -- Butcher's tableau
  #
  #   C | A
  #   --+---
  #     | B

  C = [1/5; 3/10; 4/5; 8/9; 1.];

  A = zeros(5,5);
  A[:,1] = [1/5,          0,           0,            0,         0]
  A[:,2] = [3/40,         9/40,        0,            0,         0]
  A[:,3] = [44/45,        -56/15,       32/9,         0,         0]
  A[:,4] = [19372/6561,  -25360/2187,  64448/6561,  -212/729,   0]
  A[:,5] = [9017/3168,   -355/33,      46732/5247,   49/176,   -5103/18656];


  B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];


  nstages = length(B);
  F = zeros(neq,nstages);

  Y[:,1] = y0;
  for i = 2:N
    ti = tspan[i-1];
    hi = h[i-1];
    yi = Y[:,i-1];

    # General explicit Runge-Kutta framework
    F[:,1] = odefun(ti,yi,varargin...);
    for stage = 2:nstages
      tstage = ti + C[(stage-1)]*hi;
      ystage = yi + F[:,1:stage-1]*(hi*A[1:stage-1,stage-1]);
      F[:,stage] = odefun(tstage,ystage,varargin...);
    end
    Y[:,i] = yi + F*(hi*B);

  end

  Y
end


function ode7(odefun, tspan, x0, varargin...)
  # Copyright (C) 2001, 2000 Marc Compere
  # This file is intended for use with Octave.
  # ode78.m is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2, or (at your option)
  # any later version.
  #
  # ode78.m is distributed in the hope that it will be useful, but
  # WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  # General Public License for more details at www.gnu.org/copyleft/gpl.html.
  #
  # --------------------------------------------------------------------
  #
  # ode78 (v1.14) Integrates a system of ordinary differential equations using
  # 7th order formulas.
  #
  # This is a 7th-order accurate integrator therefore the local error normally
  # expected is O(h^8).  However, because this particular implementation
  # uses the 8th-order estimate for xout (i.e. local extrapolation) moving
  # forward with the 8th-order estimate will yield errors on the order of O(h^9).
  #
  # The order of the RK method is the order of the local *truncation* error, d,
  # which is the principle error term in the portion of the Taylor series
  # expansion that gets dropped, or intentionally truncated.  This is different
  # from the local error which is the difference between the estimated solution
  # and the actual, or true solution.  The local error is used in stepsize
  # selection and may be approximated by the difference between two estimates of
  # different order, l(h) = x_(O(h+1)) - x_(O(h)).  With this definition, the
  # local error will be as large as the error in the lower order method.
  # The local truncation error is within the group of terms that gets multipled
  # by h when solving for a solution from the general RK method.  Therefore, the
  # order-p solution created by the RK method will be roughly accurate to O(h^(p+1))
  # since the local truncation error shows up in the solution as h*d, which is
  # h times an O(h^(p)) term, or rather O(h^(p+1)).
  # Summary:   For an order-p accurate RK method,
  #            - the local truncation error is O(h^p)
  #            - the local error used for stepsize adjustment and that
  #              is actually realized in a solution is O(h^(p+1))
  #
  # This requires 13 function evaluations per integration step.
  #
  # Relevant discussion on step size choice can be found on pp.90,91 in
  # U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
  # and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
  # (SIAM), Philadelphia, 1998
  #
  # More may be found in the original author's text containing numerous
  # applications on ordinary and partial differential equations using Matlab:
  #
  #     Howard Wilson and Louis Turcotte, 'Advanced Mathematics and
  #     Mechanics Applications Using MATLAB', 2nd Ed, CRC Press, 1997
  #
  #
  # [tout, xout] = ode78(FUN,tspan,x0,ode_fcn_format,tol,trace,count,hmax)
  #
  # INPUT:
  # FUN   - String containing name of user-supplied problem description.
  #         Call: xprime = fun(t,x) where FUN = 'fun'.
  #         t      - Time (scalar).
  #         x      - Solution column-vector.
  #         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
  # tspan - [ tstart, tfinal ]
  # x0    - Initial value COLUMN-vector.
  # ode_fcn_format - this specifies if the user-defined ode function is in
  #         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
  #         or:           xprime = fun(x,t)   (ode_fcn_format=1)
  #         Matlab's solvers comply with ode_fcn_format=0 while
  #         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
  # tol   - The desired accuracy. (optional, default: tol = 1.e-6).
  # trace - If nonzero, each step is printed. (optional, default: trace = 0).
  # count - if nonzero, variable 'rhs_counter' is initalized, made global
  #         and counts the number of state-dot function evaluations
  #         'rhs_counter' is incremented in here, not in the state-dot file
  #         simply make 'rhs_counter' global in the file that calls ode78
  # hmax  - limit the maximum stepsize to be less than or equal to hmax
  #
  # OUTPUT:
  # tout  - Returned integration time points (row-vector).
  # xout  - Returned solution, one solution column-vector per tout-value.
  #
  # The result can be displayed by: plot(tout, xout).

  #   Daljeet Singh & Howard Wilson
  #   Dept. Of Electrical Engg., The University of Alabama.
  #   11-24-1988.
  #
  # modified by:
  # Marc Compere
  # CompereM@asme.org
  # created : 06 October 1999
  # modified: 19 May 2001
  #
  #Julia version by Nathan Parrish, 2017

  # The Fehlberg coefficients:
  alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];

  beta_ = zeros(13,12); #pre-allocate (MUCH faster this way)
  beta_[:,1] = [2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,2] = [1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,3] = [1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,4] = [5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,5] = [0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,6] = [-25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,7] = [31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0]
  beta_[:,8] = [2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0]
  beta_[:,9] = [-91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0 ]
  beta_[:,10] = [2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0]
  beta_[:,11] = [3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0 ]
  beta_[:,12] = [-1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0 ]

  chi_  = [ 0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840];
  psi_  = [ 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ];

  pow = 1/8; # see p.91 in the Ascher & Petzold reference for more infomation.

  # hmax = (tspan[2] - tspan[1])/2.5;

  # Initialization
  # t = tspan[1];
  # tfinal = tspan[2];
  # hmin = (tfinal - t)/1e7; #could also make this an input...
  # h = (tfinal - t)/50;

  h = diff(tspan)

  neq = length(x0);
  N = length(tspan)
  Xout = zeros(neq, N)

  nstages = length(psi_)

  # xout = copy(x0);
  # x = copy(x0); #copy() so we don't modify x0 accidentally
  f = zeros(neq,nstages);  # f needs to be an neq x 13 matrix where neq=number of rows in x
  # tout = [t]; #initialize as an array so it can grow later
  # tau = tol * max(norm(x, Inf), 1);

  Xout[:,1] = x0

  # The main loop
  for ind = 2:N
    ti = tspan[ind-1]
    hi = h[ind-1];
    xi = Xout[:,ind-1]

    # if t + h > tfinal
    #   h = tfinal - t;
    # end


    # Compute the slopes
    f[:,1] = odefun(ti,xi,varargin...)
    for j = 1:12
      # println("\n")
      # println(size(xi))
      # println(size(hi))
      # println(size(f))
      # println(size(beta_))
      # println(size(alpha_))
      # println(size(ti))

      # temp1 = ti+alpha_[j]*hi;
      # temp2 = xi+hi*f*beta_[:,j];

      f[:,j+1] = odefun(ti+alpha_[j]*hi, xi+hi*f*beta_[:,j], varargin...)
    end

    Xout[:,ind] = xi + hi*f*chi_

    # Truncation error term
    # gamma1 = h*41/840*f*psi_;
    # Estimate the error and the acceptable error
    # delta = norm(gamma1, Inf);
    # tau = tol*max(norm(x, Inf),1);
    # Update the solution only if the error is acceptable
    # if delta <= tau
    #   # t += h;
    #   # x += h*f*chi_;  # this integrator uses local extrapolation
    #
    #   #Efficiently grow the size of arrays. The 'push!' operator requires xout
    #   #to be a vector (not array). It will be reshaped later.
    #   push!(tout, t)
    #   for i = 1:neq
    #     push!(xout, x[i])
    #   end
    # end

    # Update the step size
    # if delta == 0.0
    #   delta = 1e-16;
    # end
    # h = min(hmax, 0.8*h*(tau/delta)^pow);

  end;

  # if (t < tfinal)
  #   warn("Encountered a singularity at t = ", t)
  # end

  # Xout = reshape(xout, neq, length(tout))

  #outputs:
  # (tout, Xout)
  Xout
end


function ode78(odefun, tspan, x0, tol, varargin...)
  # Copyright (C) 2001, 2000 Marc Compere
  # This file is intended for use with Octave.
  # ode78.m is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2, or (at your option)
  # any later version.
  #
  # ode78.m is distributed in the hope that it will be useful, but
  # WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  # General Public License for more details at www.gnu.org/copyleft/gpl.html.
  #
  # --------------------------------------------------------------------
  #
  # ode78 (v1.14) Integrates a system of ordinary differential equations using
  # 7th order formulas.
  #
  # This is a 7th-order accurate integrator therefore the local error normally
  # expected is O(h^8).  However, because this particular implementation
  # uses the 8th-order estimate for xout (i.e. local extrapolation) moving
  # forward with the 8th-order estimate will yield errors on the order of O(h^9).
  #
  # The order of the RK method is the order of the local *truncation* error, d,
  # which is the principle error term in the portion of the Taylor series
  # expansion that gets dropped, or intentionally truncated.  This is different
  # from the local error which is the difference between the estimated solution
  # and the actual, or true solution.  The local error is used in stepsize
  # selection and may be approximated by the difference between two estimates of
  # different order, l(h) = x_(O(h+1)) - x_(O(h)).  With this definition, the
  # local error will be as large as the error in the lower order method.
  # The local truncation error is within the group of terms that gets multipled
  # by h when solving for a solution from the general RK method.  Therefore, the
  # order-p solution created by the RK method will be roughly accurate to O(h^(p+1))
  # since the local truncation error shows up in the solution as h*d, which is
  # h times an O(h^(p)) term, or rather O(h^(p+1)).
  # Summary:   For an order-p accurate RK method,
  #            - the local truncation error is O(h^p)
  #            - the local error used for stepsize adjustment and that
  #              is actually realized in a solution is O(h^(p+1))
  #
  # This requires 13 function evaluations per integration step.
  #
  # Relevant discussion on step size choice can be found on pp.90,91 in
  # U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
  # and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
  # (SIAM), Philadelphia, 1998
  #
  # More may be found in the original author's text containing numerous
  # applications on ordinary and partial differential equations using Matlab:
  #
  #     Howard Wilson and Louis Turcotte, 'Advanced Mathematics and
  #     Mechanics Applications Using MATLAB', 2nd Ed, CRC Press, 1997
  #
  #
  # [tout, xout] = ode78(FUN,tspan,x0,ode_fcn_format,tol,trace,count,hmax)
  #
  # INPUT:
  # FUN   - String containing name of user-supplied problem description.
  #         Call: xprime = fun(t,x) where FUN = 'fun'.
  #         t      - Time (scalar).
  #         x      - Solution column-vector.
  #         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
  # tspan - [ tstart, tfinal ]
  # x0    - Initial value COLUMN-vector.
  # ode_fcn_format - this specifies if the user-defined ode function is in
  #         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
  #         or:           xprime = fun(x,t)   (ode_fcn_format=1)
  #         Matlab's solvers comply with ode_fcn_format=0 while
  #         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
  # tol   - The desired accuracy. (optional, default: tol = 1.e-6).
  # trace - If nonzero, each step is printed. (optional, default: trace = 0).
  # count - if nonzero, variable 'rhs_counter' is initalized, made global
  #         and counts the number of state-dot function evaluations
  #         'rhs_counter' is incremented in here, not in the state-dot file
  #         simply make 'rhs_counter' global in the file that calls ode78
  # hmax  - limit the maximum stepsize to be less than or equal to hmax
  #
  # OUTPUT:
  # tout  - Returned integration time points (row-vector).
  # xout  - Returned solution, one solution column-vector per tout-value.
  #
  # The result can be displayed by: plot(tout, xout).

  #   Daljeet Singh & Howard Wilson
  #   Dept. Of Electrical Engg., The University of Alabama.
  #   11-24-1988.
  #
  # modified by:
  # Marc Compere
  # CompereM@asme.org
  # created : 06 October 1999
  # modified: 19 May 2001
  #
  #Julia version by Nathan Parrish, 2017

  # The Fehlberg coefficients:
  alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];
  beta_  = [ 2/27  0  0  0  0  0  0  0  0  0  0  0  0 ;
            1/36  1/12  0  0  0  0  0  0  0  0  0  0  0 ;
            1/24  0  1/8  0  0  0  0  0  0  0  0  0  0 ;
            5/12  0  -25/16  25/16  0  0  0  0  0  0  0  0  0 ;
            0.05  0  0  0.25  0.2  0  0  0  0  0  0  0  0 ;
            -25/108  0  0  125/108  -65/27  125/54  0  0  0  0  0  0  0 ;
            31/300  0  0  0  61/225  -2/9  13/900  0  0  0  0  0  0 ;
            2  0  0  -53/6  704/45  -107/9  67/90  3  0  0  0  0  0 ;
            -91/108  0  0  23/108  -976/135  311/54  -19/60  17/6  -1/12  0  0  0  0 ;
            2383/4100  0  0  -341/164  4496/1025  -301/82  2133/4100  45/82  45/164  18/41  0  0  0 ;
            3/205  0  0  0  0  -6/41  -3/205  -3/41  3/41  6/41  0  0  0 ;
            -1777/4100  0  0  -341/164  4496/1025  -289/82  2193/4100  51/82  33/164  12/41  0  1  0 ]';
  chi_  = [ 0  0  0  0  0  34/105  9/35  9/35  9/280  9/280  0  41/840  41/840]';
  psi_  = [ 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ];

  pow = 1/8; # see p.91 in the Ascher & Petzold reference for more infomation.

  hmax = (tspan[2] - tspan[1])/2.5;

  # Initialization
  t = tspan[1];
  tfinal = tspan[2];
  hmin = (tfinal - t)/1e7; #could also make this an input...
  h = (tfinal - t)/50;

  N = length(x0);

  xout = copy(x0);
  x = copy(x0); #copy() so we don't modify x0 accidentally
  f = x*zeros(1,13);  # f needs to be an Nx13 matrix where N=number of rows in x
  tout = [t]; #initialize as an array so it can grow later
  tau = tol * max(norm(x, Inf), 1);


  # The main loop
  while (t < tfinal) & (h >= hmin)
    if t + h > tfinal
      h = tfinal - t;
    end

    # Compute the slopes
    f[:,1] = odefun(t,x,varargin...)
    for j = 1:12
      f[:,j+1] = odefun(t+alpha_[j]*h, x+h*f*beta_[:,j], varargin...)
    end

    # Truncation error term
    gamma1 = h*41/840*f*psi_;

    # Estimate the error and the acceptable error
    delta = norm(gamma1, Inf);
    tau = tol*max(norm(x, Inf),1);

    # Update the solution only if the error is acceptable
    if delta <= tau
      t += h;
      x += h*f*chi_;  # this integrator uses local extrapolation

      #Efficiently grow the size of arrays. The 'push!' operator requires xout
      #to be a vector (not array). It will be reshaped later.
      push!(tout, t)
      for i = 1:N
        push!(xout, x[i])
      end
    end

    # Update the step size
    if delta == 0.0
      delta = 1e-16;
    end
    h = min(hmax, 0.8*h*(tau/delta)^pow);

  end;

  if (t < tfinal)
    warn("Encountered a singularity at t = ", t)
  end

  Xout = reshape(xout, N, length(tout))

  #outputs:
  (tout, Xout)
end


function ode78_events(odefun, evfun, tspan, x0, tol, etol, varargin...)
  # Copyright (C) 2001, 2000 Marc Compere
  # This file is intended for use with Octave.
  # ode78.m is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2, or (at your option)
  # any later version.
  #
  # ode78.m is distributed in the hope that it will be useful, but
  # WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  # General Public License for more details at www.gnu.org/copyleft/gpl.html.
  #
  # --------------------------------------------------------------------
  #
  # ode78 (v1.14) Integrates a system of ordinary differential equations using
  # 7th order formulas.
  #
  # This is a 7th-order accurate integrator therefore the local error normally
  # expected is O(h^8).  However, because this particular implementation
  # uses the 8th-order estimate for xout (i.e. local extrapolation) moving
  # forward with the 8th-order estimate will yield errors on the order of O(h^9).
  #
  # The order of the RK method is the order of the local *truncation* error, d,
  # which is the principle error term in the portion of the Taylor series
  # expansion that gets dropped, or intentionally truncated.  This is different
  # from the local error which is the difference between the estimated solution
  # and the actual, or true solution.  The local error is used in stepsize
  # selection and may be approximated by the difference between two estimates of
  # different order, l(h) = x_(O(h+1)) - x_(O(h)).  With this definition, the
  # local error will be as large as the error in the lower order method.
  # The local truncation error is within the group of terms that gets multipled
  # by h when solving for a solution from the general RK method.  Therefore, the
  # order-p solution created by the RK method will be roughly accurate to O(h^(p+1))
  # since the local truncation error shows up in the solution as h*d, which is
  # h times an O(h^(p)) term, or rather O(h^(p+1)).
  # Summary:   For an order-p accurate RK method,
  #            - the local truncation error is O(h^p)
  #            - the local error used for stepsize adjustment and that
  #              is actually realized in a solution is O(h^(p+1))
  #
  # This requires 13 function evaluations per integration step.
  #
  # Relevant discussion on step size choice can be found on pp.90,91 in
  # U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
  # and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
  # (SIAM), Philadelphia, 1998
  #
  # More may be found in the original author's text containing numerous
  # applications on ordinary and partial differential equations using Matlab:
  #
  #     Howard Wilson and Louis Turcotte, 'Advanced Mathematics and
  #     Mechanics Applications Using MATLAB', 2nd Ed, CRC Press, 1997
  #
  #
  # [tout, xout] = ode78(FUN,tspan,x0,ode_fcn_format,tol,trace,count,hmax)
  #
  # INPUT:
  # FUN   - String containing name of user-supplied problem description.
  #         Call: xprime = fun(t,x) where FUN = 'fun'.
  #         t      - Time (scalar).
  #         x      - Solution column-vector.
  #         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
  # tspan - [ tstart, tfinal ]
  # x0    - Initial value COLUMN-vector.
  # ode_fcn_format - this specifies if the user-defined ode function is in
  #         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
  #         or:           xprime = fun(x,t)   (ode_fcn_format=1)
  #         Matlab's solvers comply with ode_fcn_format=0 while
  #         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
  # tol   - The desired accuracy. (optional, default: tol = 1.e-6).
  # trace - If nonzero, each step is printed. (optional, default: trace = 0).
  # count - if nonzero, variable 'rhs_counter' is initalized, made global
  #         and counts the number of state-dot function evaluations
  #         'rhs_counter' is incremented in here, not in the state-dot file
  #         simply make 'rhs_counter' global in the file that calls ode78
  # hmax  - limit the maximum stepsize to be less than or equal to hmax
  #
  # OUTPUT:
  # tout  - Returned integration time points (row-vector).
  # xout  - Returned solution, one solution column-vector per tout-value.
  #
  # The result can be displayed by: plot(tout, xout).

  #   Daljeet Singh & Howard Wilson
  #   Dept. Of Electrical Engg., The University of Alabama.
  #   11-24-1988.
  #
  # modified by:
  # Marc Compere
  # CompereM@asme.org
  # created : 06 October 1999
  # modified: 19 May 2001
  #
  #Julia version by Nathan Parrish, 2017
  #Events detection added based on code by Travis Swenson and Brian Anderson.

  # The Fehlberg coefficients:
  alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];
  beta_  = [ 2/27  0  0  0  0  0  0  0  0  0  0  0  0 ;
            1/36  1/12  0  0  0  0  0  0  0  0  0  0  0 ;
            1/24  0  1/8  0  0  0  0  0  0  0  0  0  0 ;
            5/12  0  -25/16  25/16  0  0  0  0  0  0  0  0  0 ;
            0.05  0  0  0.25  0.2  0  0  0  0  0  0  0  0 ;
            -25/108  0  0  125/108  -65/27  125/54  0  0  0  0  0  0  0 ;
            31/300  0  0  0  61/225  -2/9  13/900  0  0  0  0  0  0 ;
            2  0  0  -53/6  704/45  -107/9  67/90  3  0  0  0  0  0 ;
            -91/108  0  0  23/108  -976/135  311/54  -19/60  17/6  -1/12  0  0  0  0 ;
            2383/4100  0  0  -341/164  4496/1025  -301/82  2133/4100  45/82  45/164  18/41  0  0  0 ;
            3/205  0  0  0  0  -6/41  -3/205  -3/41  3/41  6/41  0  0  0 ;
            -1777/4100  0  0  -341/164  4496/1025  -289/82  2193/4100  51/82  33/164  12/41  0  1  0 ]';
  chi_  = [ 0  0  0  0  0  34/105  9/35  9/35  9/280  9/280  0  41/840  41/840]';
  psi_  = [ 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ];

  pow = 1/8; # see p.91 in the Ascher & Petzold reference for more infomation.

  hmax = (tspan[2] - tspan[1])/2.5;

  # Initialization
  t = tspan[1];
  tfinal = tspan[2];
  hmin = (tfinal - t)/1e10; #could also make this an input...
  h = (tfinal - t)/50;

  N = length(x0);

  xout = copy(x0);
  x = copy(x0); #copy() so we don't modify x0 accidentally
  f = x*zeros(1,13);  # f needs to be an Nx13 matrix where N=number of rows in x
  tout = [t]; #initialize as an array so it can grow later
  tau = tol * max(norm(x, Inf), 1);

  #initialize event checking
  integrating = true #event has not happened
  (evf, dirn) = evfun(t,x,varargin...)

  # The main loop
  while (t < tfinal) & (h >= hmin) & integrating
    if t + h > tfinal
      h = tfinal - t;
    end

    # Compute the slopes
    f[:,1] = odefun(t,x,varargin...)
    for j = 1:12
      f[:,j+1] = odefun(t+alpha_[j]*h, x+h*f*beta_[:,j], varargin...)
    end

    # Truncation error term
    gamma1 = h*41/840*f*psi_;

    # Estimate the error and the acceptable error
    delta = norm(gamma1, Inf);
    tau = tol*max(norm(x, Inf),1);

    # Update the solution only if the error is acceptable
    if delta <= tau
      t += h;
      x += h*f*chi_;  # this integrator uses local extrapolation

      #Efficiently grow the size of arrays. The 'push!' operator requires xout
      #to be a vector (not array). It will be reshaped later.
      push!(tout, t)
      for i = 1:N
        push!(xout, x[i])
      end
    end





    # determine event trigger
    restep = false
    # compute event function, save previous value
    evf_prev    = evf;
    (evf, dirn) = evfun(t,x,varargin...)

    # determine sign change in event function
    cond1   = ~(sign(evf) == sign(evf_prev))    #sign switch
    cond2   = (sign(evf) == dirn) || (dirn == 0) #switch direction
    cond3   = ~(sign(evf_prev) == 0.0)            #avoid trigger at prev. exact event
    if cond1 && cond2 && cond3
      # If we need to update algorithm
      if abs(evf) > etol
        # Delete last step to zoom in and refine
        tout = tout[1:end-1]
        xout = xout[1:end-N]
        t           = t - h
        x           = x - h * f * chi_
        evf         = evf_prev
        restep      = true
      else
        integrating = false
      end
    end


    # Update the step size
    if delta == 0.0
      delta = 1e-16;
    end

    if ~(delta == 0.0) && ~restep #avoid adaptive step size when refining event
      h = min(hmax, 0.8*h*(tau/delta)^pow);
    end

    # reduce step size when backtracking to find exact event trigger
    if restep
      hmax    = h / 5
      h       = hmax
    end

  end

  if (t < tfinal) & (integrating)
    warn("Encountered a singularity at t = ", t)
  end

  Xout = reshape(xout, N, length(tout))

  #outputs:
  (tout, Xout)
end


function ode7_8(odefun, tspan, x0, varargin...)
  #7th order fixed step numerical integration, that also outputs the 8th order
  #error estimate for use outside of this function. maxErr is the estimated
  #maximum error at any timestep


  # Copyright (C) 2001, 2000 Marc Compere
  # This file is intended for use with Octave.
  # ode78.m is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2, or (at your option)
  # any later version.
  #
  # ode78.m is distributed in the hope that it will be useful, but
  # WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  # General Public License for more details at www.gnu.org/copyleft/gpl.html.
  #
  # --------------------------------------------------------------------
  #
  # ode78 (v1.14) Integrates a system of ordinary differential equations using
  # 7th order formulas.
  #
  # This is a 7th-order accurate integrator therefore the local error normally
  # expected is O(h^8).  However, because this particular implementation
  # uses the 8th-order estimate for xout (i.e. local extrapolation) moving
  # forward with the 8th-order estimate will yield errors on the order of O(h^9).
  #
  # The order of the RK method is the order of the local *truncation* error, d,
  # which is the principle error term in the portion of the Taylor series
  # expansion that gets dropped, or intentionally truncated.  This is different
  # from the local error which is the difference between the estimated solution
  # and the actual, or true solution.  The local error is used in stepsize
  # selection and may be approximated by the difference between two estimates of
  # different order, l(h) = x_(O(h+1)) - x_(O(h)).  With this definition, the
  # local error will be as large as the error in the lower order method.
  # The local truncation error is within the group of terms that gets multipled
  # by h when solving for a solution from the general RK method.  Therefore, the
  # order-p solution created by the RK method will be roughly accurate to O(h^(p+1))
  # since the local truncation error shows up in the solution as h*d, which is
  # h times an O(h^(p)) term, or rather O(h^(p+1)).
  # Summary:   For an order-p accurate RK method,
  #            - the local truncation error is O(h^p)
  #            - the local error used for stepsize adjustment and that
  #              is actually realized in a solution is O(h^(p+1))
  #
  # This requires 13 function evaluations per integration step.
  #
  # Relevant discussion on step size choice can be found on pp.90,91 in
  # U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
  # and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
  # (SIAM), Philadelphia, 1998
  #
  # More may be found in the original author's text containing numerous
  # applications on ordinary and partial differential equations using Matlab:
  #
  #     Howard Wilson and Louis Turcotte, 'Advanced Mathematics and
  #     Mechanics Applications Using MATLAB', 2nd Ed, CRC Press, 1997
  #
  #
  # [tout, xout] = ode78(FUN,tspan,x0,ode_fcn_format,tol,trace,count,hmax)
  #
  # INPUT:
  # FUN   - String containing name of user-supplied problem description.
  #         Call: xprime = fun(t,x) where FUN = 'fun'.
  #         t      - Time (scalar).
  #         x      - Solution column-vector.
  #         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
  # tspan - [ tstart, tfinal ]
  # x0    - Initial value COLUMN-vector.
  # ode_fcn_format - this specifies if the user-defined ode function is in
  #         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
  #         or:           xprime = fun(x,t)   (ode_fcn_format=1)
  #         Matlab's solvers comply with ode_fcn_format=0 while
  #         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
  # tol   - The desired accuracy. (optional, default: tol = 1.e-6).
  # trace - If nonzero, each step is printed. (optional, default: trace = 0).
  # count - if nonzero, variable 'rhs_counter' is initalized, made global
  #         and counts the number of state-dot function evaluations
  #         'rhs_counter' is incremented in here, not in the state-dot file
  #         simply make 'rhs_counter' global in the file that calls ode78
  # hmax  - limit the maximum stepsize to be less than or equal to hmax
  #
  # OUTPUT:
  # tout  - Returned integration time points (row-vector).
  # xout  - Returned solution, one solution column-vector per tout-value.
  #
  # The result can be displayed by: plot(tout, xout).

  #   Daljeet Singh & Howard Wilson
  #   Dept. Of Electrical Engg., The University of Alabama.
  #   11-24-1988.
  #
  # modified by:
  # Marc Compere
  # CompereM@asme.org
  # created : 06 October 1999
  # modified: 19 May 2001
  #
  #Julia version by Nathan Parrish, 2017

  # The Fehlberg coefficients:
  alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];

  beta_ = zeros(13,12); #pre-allocate (MUCH faster this way)
  beta_[:,1] = [2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,2] = [1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,3] = [1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,4] = [5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,5] = [0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,6] = [-25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,7] = [31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0]
  beta_[:,8] = [2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0]
  beta_[:,9] = [-91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0 ]
  beta_[:,10] = [2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0]
  beta_[:,11] = [3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0 ]
  beta_[:,12] = [-1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0 ]

  chi_  = [ 0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840];
  psi_  = [ 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ];

  pow = 1/8; # see p.91 in the Ascher & Petzold reference for more infomation.

  # hmax = (tspan[2] - tspan[1])/2.5;

  # Initialization
  # t = tspan[1];
  # tfinal = tspan[2];
  # hmin = (tfinal - t)/1e7; #could also make this an input...
  # h = (tfinal - t)/50;

  h = diff(tspan)

  neq = length(x0);
  N = length(tspan)
  Xout = zeros(neq, N)

  nstages = length(psi_)

  # xout = copy(x0);
  # x = copy(x0); #copy() so we don't modify x0 accidentally
  f = zeros(neq,nstages);  # f needs to be an neq x 13 matrix where neq=number of rows in x
  # tout = [t]; #initialize as an array so it can grow later
  # tau = tol * max(norm(x, Inf), 1);

  Xout[:,1] = x0

  #initialize error estimate:
  maxErr = 0.0

  # The main loop
  for ind = 2:N
    ti = tspan[ind-1]
    hi = h[ind-1];
    xi = Xout[:,ind-1]


    # Compute the slopes
    f[:,1] = odefun(ti,xi,varargin...)
    for j = 1:12

      f[:,j+1] = odefun(ti+alpha_[j]*hi, xi+hi*f*beta_[:,j], varargin...)
    end

    Xout[:,ind] = xi + hi*f*chi_

    # Truncation error term
    gamma1 = hi*41/840*f*psi_;

    # Estimate the error and the acceptable error
    delta = norm(gamma1, Inf);

    #save worst case error
    if delta > maxErr
      maxErr = copy(delta)
    end
  end

  #outputs:
  (Xout, maxErr)
end


function ode7_8!(odefun, tspan, x0, params)
  #7th order fixed step numerical integration, that also outputs the 8th order
  #error estimate for use outside of this function. maxErr is the estimated
  #maximum error at any timestep


  # Copyright (C) 2001, 2000 Marc Compere
  # This file is intended for use with Octave.
  # ode78.m is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2, or (at your option)
  # any later version.
  #
  # ode78.m is distributed in the hope that it will be useful, but
  # WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  # General Public License for more details at www.gnu.org/copyleft/gpl.html.
  #
  # --------------------------------------------------------------------
  #
  # ode78 (v1.14) Integrates a system of ordinary differential equations using
  # 7th order formulas.
  #
  # This is a 7th-order accurate integrator therefore the local error normally
  # expected is O(h^8).  However, because this particular implementation
  # uses the 8th-order estimate for xout (i.e. local extrapolation) moving
  # forward with the 8th-order estimate will yield errors on the order of O(h^9).
  #
  # The order of the RK method is the order of the local *truncation* error, d,
  # which is the principle error term in the portion of the Taylor series
  # expansion that gets dropped, or intentionally truncated.  This is different
  # from the local error which is the difference between the estimated solution
  # and the actual, or true solution.  The local error is used in stepsize
  # selection and may be approximated by the difference between two estimates of
  # different order, l(h) = x_(O(h+1)) - x_(O(h)).  With this definition, the
  # local error will be as large as the error in the lower order method.
  # The local truncation error is within the group of terms that gets multipled
  # by h when solving for a solution from the general RK method.  Therefore, the
  # order-p solution created by the RK method will be roughly accurate to O(h^(p+1))
  # since the local truncation error shows up in the solution as h*d, which is
  # h times an O(h^(p)) term, or rather O(h^(p+1)).
  # Summary:   For an order-p accurate RK method,
  #            - the local truncation error is O(h^p)
  #            - the local error used for stepsize adjustment and that
  #              is actually realized in a solution is O(h^(p+1))
  #
  # This requires 13 function evaluations per integration step.
  #
  # Relevant discussion on step size choice can be found on pp.90,91 in
  # U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
  # and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
  # (SIAM), Philadelphia, 1998
  #
  # More may be found in the original author's text containing numerous
  # applications on ordinary and partial differential equations using Matlab:
  #
  #     Howard Wilson and Louis Turcotte, 'Advanced Mathematics and
  #     Mechanics Applications Using MATLAB', 2nd Ed, CRC Press, 1997
  #
  #
  # [tout, xout] = ode78(FUN,tspan,x0,ode_fcn_format,tol,trace,count,hmax)
  #
  # INPUT:
  # FUN   - String containing name of user-supplied problem description.
  #         Call: xprime = fun(t,x) where FUN = 'fun'.
  #         t      - Time (scalar).
  #         x      - Solution column-vector.
  #         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
  # tspan - [ tstart, tfinal ]
  # x0    - Initial value COLUMN-vector.
  # ode_fcn_format - this specifies if the user-defined ode function is in
  #         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
  #         or:           xprime = fun(x,t)   (ode_fcn_format=1)
  #         Matlab's solvers comply with ode_fcn_format=0 while
  #         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
  # tol   - The desired accuracy. (optional, default: tol = 1.e-6).
  # trace - If nonzero, each step is printed. (optional, default: trace = 0).
  # count - if nonzero, variable 'rhs_counter' is initalized, made global
  #         and counts the number of state-dot function evaluations
  #         'rhs_counter' is incremented in here, not in the state-dot file
  #         simply make 'rhs_counter' global in the file that calls ode78
  # hmax  - limit the maximum stepsize to be less than or equal to hmax
  #
  # OUTPUT:
  # tout  - Returned integration time points (row-vector).
  # xout  - Returned solution, one solution column-vector per tout-value.
  #
  # The result can be displayed by: plot(tout, xout).

  #   Daljeet Singh & Howard Wilson
  #   Dept. Of Electrical Engg., The University of Alabama.
  #   11-24-1988.
  #
  # modified by:
  # Marc Compere
  # CompereM@asme.org
  # created : 06 October 1999
  # modified: 19 May 2001
  #
  #Julia version by Nathan Parrish, 2017

  # The Fehlberg coefficients:
  alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];

  beta_ = zeros(13,12); #pre-allocate (MUCH faster this way)
  beta_[:,1] = [2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,2] = [1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,3] = [1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  beta_[:,4] = [5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,5] = [0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,6] = [-25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0]
  beta_[:,7] = [31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0]
  beta_[:,8] = [2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0]
  beta_[:,9] = [-91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0 ]
  beta_[:,10] = [2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0]
  beta_[:,11] = [3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0 ]
  beta_[:,12] = [-1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0 ]

  chi_  = [ 0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840];
  psi_  = [ 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ];

  pow = 1/8; # see p.91 in the Ascher & Petzold reference for more infomation.

  h = diff(tspan)

  neq = length(x0);
  N = length(tspan)
  Xout = zeros(neq, N)

  nstages = length(psi_)

  f = zeros(neq,nstages);  # f needs to be an neq x 13 matrix where neq=number of rows in x

  Xout[:,1] = x0

  #initialize error estimate:
  maxErr = 0.0

  # The main loop
  for ind = 2:N
    ti = tspan[ind-1]
    hi = h[ind-1];
    xi = Xout[:,ind-1]


    # Compute the slopes
    dstate = zeros(size(x0))
    odefun(ti, xi, params, dstate)
    f[:,1] = dstate
    for j = 1:12

        # f[:,j+1] = odefun(ti+alpha_[j]*hi, xi+hi*f*beta_[:,j], varargin...)
        odefun(ti+alpha_[j]*hi, xi+hi*f*beta_[:,j], params, dstate)
        f[:,j+1] = dstate
    end

    Xout[:,ind] = xi + hi*f*chi_

    # Truncation error term
    gamma1 = hi*41/840*f*psi_;

    # Estimate the error and the acceptable error
    delta = norm(gamma1 ./  Xout[:,ind], Inf);

    #save worst case error
    if delta > maxErr
      maxErr = copy(delta)
    end
  end

  #outputs:
  (Xout, maxErr)
end
