#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Derivatives for states & costates -- for CRTBP

#For use with NEW version of DifferentialEquations.jl
function CRTBP_stateCostate_deriv!(dstate_costate, statesCostates, params, t)
  #Constant mass.

  #pull out parameters
  (MU, DU, TU, thrustLimit, mass, time_direction, p, rho) = params

  λv = statesCostates[10:12]

  # States
  X1 = statesCostates[1]
  X2 = statesCostates[2]
  X3 = statesCostates[3]
  X4 = statesCostates[4]
  X5 = statesCostates[5]
  X6 = statesCostates[6]

  #Costates:
  L1 = statesCostates[7]
  L2 = statesCostates[8]
  L3 = statesCostates[9]
  L4 = statesCostates[10]
  L5 = statesCostates[11]
  L6 = statesCostates[12]

  accelLimit = thrustLimit / mass / 1e3 * TU^2 / DU #(kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration

  #Thrust magnitude from control law:
  if p == 0
    #Force thrust to always be on, in the optimal direction. This limits how big
    #of a search we need to make for the neural net stuff.
    umag = accelLimit #N

  elseif p == 1 #minimum mass problem
      g = norm(λv) - 1
      umag = 1/2 * (1 + tanh( g / (2*rho))) * accelLimit

  elseif p > 1 #p is in the range (1, 2]
      umag = (1/p * norm(λv))^(1 / (p-1));

      if umag > accelLimit
        umag = accelLimit
      end
  else
    error("Invalid value of p!")
  end


  #Control vector, scaled properly:
  control_accel = -umag * λv ./ norm(λv) #DU/TU^2

  if isnan(control_accel[1])
    for i = 1:length(control_accel)
      control_accel[i] = 0.
    end 
    
  end


  #calculate some stuff that shows up a lot:
  #calculate r1^3 and r2^3 once each
  r1_3 = ((X1+MU)^2 + X2^2 + X3^2)^(3/2);
  r2_3 = ((X1+MU-1)^2 + X2^2 + X3^2)^(3/2);

  temp1 = ((MU + X1 - 1)^2 + X2^2 + X3^2);
  temp2 = ((MU + X1)^2 + X2^2 + X3^2);
  temp3 = (2*MU + 2*X1 - 2);

  #################################################### Derivatives:
  #Output in place:
  dstate_costate[1:3] = statesCostates[4:6]
  dstate_costate[4] = -(1-MU)*(X1+MU)/r1_3 - MU*(X1-1+MU)/r2_3 + 2*time_direction*X5 + X1 + control_accel[1];
  dstate_costate[5] = -(1-MU)*X2/r1_3 - MU*X2/r2_3 - 2*time_direction*X4 + X2 + control_accel[2];
  dstate_costate[6] = -(1-MU)*X3/r1_3 - MU*X3/r2_3 + control_accel[3];

  dstate_costate[7] = - L5*((3*MU*X2*temp3)/(2*temp1^(5/2)) - (3*X2*(MU - 1)*(2*MU + 2*X1))/(2*temp2^(5/2))) - L6*((3*MU*X3*temp3)/(2*temp1^(5/2)) - (3*X3*(MU - 1)*(2*MU + 2*X1))/(2*temp2^(5/2))) - L4*((MU - 1)/temp2^(3/2) - MU/temp1^(3/2) + (3*MU*(MU + X1 - 1)*temp3)/(2*temp1^(5/2)) - (3*(MU + X1)*(MU - 1)*(2*MU + 2*X1))/(2*temp2^(5/2)) + 1);
  dstate_costate[8] = L6*((3*X2*X3*(MU - 1))/temp2^(5/2) - (3*MU*X2*X3)/temp1^(5/2)) - L5*((MU - 1)/temp2^(3/2) - MU/temp1^(3/2) - (3*X2^2*(MU - 1))/temp2^(5/2) + (3*MU*X2^2)/temp1^(5/2) + 1) - L4*((3*MU*X2*(MU + X1 - 1))/temp1^(5/2) - (3*X2*(MU + X1)*(MU - 1))/temp2^(5/2));
  dstate_costate[9] = L6*(MU/temp1^(3/2) - (MU - 1)/temp2^(3/2) + (3*X3^2*(MU - 1))/temp2^(5/2) - (3*MU*X3^2)/temp1^(5/2)) + L5*((3*X2*X3*(MU - 1))/temp2^(5/2) - (3*MU*X2*X3)/temp1^(5/2)) - L4*((3*MU*X3*(MU + X1 - 1))/temp1^(5/2) - (3*X3*(MU + X1)*(MU - 1))/temp2^(5/2));
  dstate_costate[10] = 2*L5*time_direction - L1;
  dstate_costate[11] = - L2 - 2*L4*time_direction;
  dstate_costate[12] = -L3

end
