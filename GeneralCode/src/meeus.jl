#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Approximate planetary ephemerides, from:
#J. Meeus, Astronomical Algorithms. Richmond, VA: William-Bell, Inc., First
#English ed., 1991

function meeus_ephem(T_JD,planet_num)

  #
  #Inputs:
  #   T_JD = julian date
  #       should be a vector or a scalar
  #   planet_num = Number of planet (i.e. Earth=3, Mars=4)
  #       should be a scalar (so all the calls are for the same planet)
  #Outputs:
  #   R       = km
  #   V       = km/s
  #   sma     = km
  #   ecc     = --
  #   inc     = rad
  #   RAAN    = rad
  #   AOP     = rad
  #   M       = rad
  #
  #All vectorized operations (Jan 30, 2017)

  T = (T_JD - 2451545)./36525; #JD -> Julian centuries since J2000.0

  mu_sun = 1.32712440018e11; #km3./s2
  AU = 149599650; #1 AU in km

  ## find values of constants
  if planet_num == 1
      ## Mercury
      a0_L = 252.250906;
      a1_L = 149472.6746358;
      a2_L = -0.00000535;
      a3_L = 0.000000002;

      a0_a = 0.387098310;
      a1_a = 0;
      a2_a = 0;

      a0_e = 0.20563175;
      a1_e = 0.000020406;
      a2_e = -0.0000000284;
      a3_e = -0.00000000017;

      a0_i = 7.004986;
      a1_i = -0.0059516;
      a2_i = 0.00000081;
      a3_i = 0.000000041;

      a0_RAAN = 48.330893;
      a1_RAAN = -0.1254229;
      a2_RAAN = -0.00008833;
      a3_RAAN = -0.000000196;

      a0_Pi = 77.456119;
      a1_Pi = 0.1588643;
      a2_Pi = -0.00001343;
      a3_Pi = 0.000000039;

  elseif planet_num == 2
      ## Venus
      a0_L = 181.979801;
      a1_L = 58517.8156760;
      a2_L = 0.00000165;
      a3_L = -0.000000002;

      a0_a = 0.72332982;
      a1_a = 0;
      a2_a = 0;

      a0_e = 0.00677188;
      a1_e = -0.000047766;
      a2_e = 0.0000000975;
      a3_e = 0.00000000044;

      a0_i = 3.394662;
      a1_i = -0.0008568;
      a2_i = -0.00003244;
      a3_i = 0.000000010;

      a0_RAAN = 76.679920;
      a1_RAAN = -0.2780080;
      a2_RAAN = -0.00014256;
      a3_RAAN = -0.000000198;

      a0_Pi = 131.563707;
      a1_Pi = 0.0048646;
      a2_Pi = -0.00138232;
      a3_Pi = -0.000005332;

  elseif planet_num == 3
      ## Earth
      a0_L = 100.466449;
      a1_L = 35999.3728519;
      a2_L = - 0.00000568;
      a3_L = 0;

      a0_a = 1.000001018;
      a1_a = 0;
      a2_a = 0;

      a0_e = 0.01670862;
      a1_e = - 0.000042037;
      a2_e = - 0.0000001236;
      a3_e = 0.00000000004;

      a0_i = 0.0;
      a1_i = 0.0130546;
      a2_i = - 0.00000931;
      a3_i = - 0.000000034;

      a0_RAAN = 174.873174;
      a1_RAAN = -0.2410908;
      a2_RAAN = 0.00004067;
      a3_RAAN = -0.000001327;

      a0_Pi = 102.937348;
      a1_Pi = 0.3225557;
      a2_Pi = 0.00015026;
      a3_Pi = 0.000000478;

  elseif planet_num == 4
      ## Mars
      a0_L =355.433275;
      a1_L =19140.2993313;
      a2_L =0.00000261;
      a3_L =- 0.000000003;

      a0_a =1.523679342;
      a1_a = 0;
      a2_a = 0;

      a0_e =0.09340062;
      a1_e =0.000090483;
      a2_e =- 0.0000000806;
      a3_e =- 0.00000000035;

      a0_i =1.849726;
      a1_i =- 0.0081479;
      a2_i =- 0.00002255;
      a3_i =- 0.000000027;

      a0_RAAN =49.558093;
      a1_RAAN =- 0.2949846;
      a2_RAAN =- 0.00063993;
      a3_RAAN =- 0.000002143;

      a0_Pi =336.060234;
      a1_Pi =0.4438898;
      a2_Pi =- 0.00017321;
      a3_Pi =0.000000300;

  elseif planet_num == 5
      ## Jupiter
      a0_L =34.351484;
      a1_L =3034.9056746;
      a2_L =-0.00008501;
      a3_L =0.000000004;

      a0_a =5.202603191;
      a1_a = 0.0000001913;
      a2_a = 0;

      a0_e =0.04849485;
      a1_e =0.000163244;
      a2_e =-0.0000004719;
      a3_e =-0.00000000197;

      a0_i =1.303270;
      a1_i =-0.0019872;
      a2_i =0.00003318;
      a3_i =0.000000092;

      a0_RAAN =100.464441;
      a1_RAAN =0.1766828;
      a2_RAAN =0.00090387;
      a3_RAAN =-0.000007032;

      a0_Pi =14.331309;
      a1_Pi =0.2155525;
      a2_Pi =0.00072252;
      a3_Pi =-0.000004590;

  elseif planet_num == 6
      ## Saturn
      a0_L = 50.077471;
      a1_L = 1222.1137943;
      a2_L = 0.00021004;
      a3_L = -0.000000019;

      a0_a = 9.554909596;
      a1_a = -0.0000021389;
      a2_a = 0;

      a0_e = 0.05550862;
      a1_e = -0.000346818;
      a2_e = -0.0000006456;
      a3_e = 0.00000000338;

      a0_i = 2.488878;
      a1_i = 0.0025515;
      a2_i = -0.00004903;
      a3_i = 0.000000018;

      a0_RAAN = 113.665524;
      a1_RAAN = -0.2566649;
      a2_RAAN = -0.00018345;
      a3_RAAN = 0.000000357;

      a0_Pi = 93.056787;
      a1_Pi = 0.5665496;
      a2_Pi = 0.00052809;
      a3_Pi = 0.000004882;

  elseif planet_num == 7
      ## Uranus
      a0_L = 314.055005;
      a1_L = 428.4669983;
      a2_L = -0.00000486;
      a3_L = 0.000000006;

      a0_a = 19.218446062;
      a1_a = -0.0000000372;
      a2_a = 0.00000000098;

      a0_e = 0.04629590;
      a1_e = -0.000027337;
      a2_e = 0.0000000790;
      a3_e = 0.00000000025;

      a0_i = 0.773196;
      a1_i = -0.0016869;
      a2_i = 0.00000349;
      a3_i = 0.000000016;

      a0_RAAN = 74.005947;
      a1_RAAN = 0.0741461;
      a2_RAAN = 0.00040540;
      a3_RAAN = 0.000000104;

      a0_Pi = 173.005159;
      a1_Pi = 0.0893206;
      a2_Pi = -0.00009470;
      a3_Pi = 0.000000413;

  elseif planet_num == 8
      ## Neptune
      a0_L = 304.348665;
      a1_L = 218.4862002;
      a2_L = 0.00000059;
      a3_L = -0.000000002;

      a0_a = 30.110386869;
      a1_a = -0.0000001663;
      a2_a = 0.00000000069;

      a0_e = 0.00898809;
      a1_e = 0.000006408;
      a2_e = -0.0000000008;
      a3_e = -0.00000000005;

      a0_i = 1.769952;
      a1_i = 0.0002257;
      a2_i = 0.00000023;
      a3_i = 0.0;

      a0_RAAN = 131.784057;
      a1_RAAN = -0.0061651;
      a2_RAAN = -0.00000219;
      a3_RAAN = -0.000000078;

      a0_Pi = 48.123691;
      a1_Pi = 0.0291587;
      a2_Pi = 0.00007051;
      a3_Pi = -0.000000023;

  elseif planet_num == 9
      ## Pluto
      a0_L = 238.92903833;
      a1_L = 145.20780515;
      a2_L = 0;
      a3_L = 0;

      a0_a = 39.48211675;
      a1_a = -0.00031596;
      a2_a = 0;

      a0_e = 0.24882730;
      a1_e = 0.00005170;
      a2_e = 0;
      a3_e = 0;

      a0_i = 17.14001206;
      a1_i = 0.00004818;
      a2_i = 0;
      a3_i = 0;

      a0_RAAN = 110.30393684;
      a1_RAAN = -0.01183482;
      a2_RAAN = 0;
      a3_RAAN = 0;

      a0_Pi = 224.06891629;
      a1_Pi = -0.04062942;
      a2_Pi = 0;
      a3_Pi = 0;

  end

  ## find orbital elements
  deg2rad = pi/180;

  L = a0_L + a1_L.*T + a2_L.*T.*T + a3_L.*T.*T.*T; # Mean longitude (deg)
  sma = a0_a + a1_a.*T + a2_a.*T.*T; # Semimajor axis (AU)
  sma = sma.*AU; # convert AU -> km
  ecc = a0_e + a1_e.*T + a2_e.*T.*T + a3_e.*T.*T.*T; # eccentricity
  inc = a0_i + a1_i.*T + a2_i.*T.*T + a3_i.*T.*T.*T; # inclination (deg)
  inc = inc.*deg2rad;
  RAAN = a0_RAAN + a1_RAAN.*T + a2_RAAN.*T.*T + a3_RAAN.*T.*T.*T; # RAAN (deg)
  longitude_peri = a0_Pi + a1_Pi.*T + a2_Pi.*T.*T + a3_Pi.*T.*T.*T; # Longitude of the Perihelion (deg)

  #convert to radians:
  L = deg2rad.*L;
  longitude_peri = deg2rad.*longitude_peri;
  RAAN = deg2rad.*RAAN;

  AOP = longitude_peri- RAAN; # argument of perihelion (rad)
  M = L -longitude_peri; # mean anomaly (rad)

  ecc2 = ecc.*ecc;
  ecc3 = ecc2.*ecc;
  ecc4 = ecc3.*ecc;
  ecc5 = ecc4.*ecc;

  C_cen = (2.*ecc - ecc3./4 + 5./96.*ecc5).*sin(M) + (5./4.*ecc2 - 11./24.*ecc4).*sin(2.*M) + (13./12.*ecc3 - 43./64.*ecc5).*sin(3.*M) + 103./96.*ecc4.*sin(4.*M) + 1097./960.*ecc5.*sin(5.*M);
  nu = M + C_cen; # true anomaly (rad)



  ## ---------- convert COEs to r, v (Code taken from Jeff Parker) ----------
  p = sma.*(1 - ecc.*ecc);                    # semiparameter (km)

  ##  Position Coordinates in Perifocal Coordinate System
  x  = (p.*cos(nu)) ./ (1.0 + ecc.*cos(nu)); # x-coordinate (km)
  y  = (p.*sin(nu)) ./ (1.0 + ecc.*cos(nu)); # y-coordinate (km)
  vx = -sqrt(mu_sun./p) .* sin(nu);      # velocity in x (km./s)
  vy = sqrt(mu_sun./p) .* (ecc + cos(nu));  # velocity in y (km./s)

  r_1 = x.*(cos(RAAN).*cos(AOP) - sin(RAAN).*cos(inc).*sin(AOP)) - y.*(cos(RAAN).*sin(AOP) + sin(RAAN).*cos(inc).*cos(AOP));
  r_2 = x.*(sin(RAAN).*cos(AOP) + cos(RAAN).*cos(inc).*sin(AOP)) - y.*(sin(RAAN).*sin(AOP) - cos(RAAN).*cos(inc).*cos(AOP));
  r_3 = y.*cos(AOP).*sin(inc) + x.*sin(inc).*sin(AOP);

  v_1 = vx.*(cos(RAAN).*cos(AOP) - sin(RAAN).*cos(inc).*sin(AOP)) - vy.*(cos(RAAN).*sin(AOP) + sin(RAAN).*cos(inc).*cos(AOP));
  v_2 = vx.*(sin(RAAN).*cos(AOP) + cos(RAAN).*cos(inc).*sin(AOP)) - vy.*(sin(RAAN).*sin(AOP) - cos(RAAN).*cos(inc).*cos(AOP));
  v_3 = vy.*cos(AOP).*sin(inc) + vx.*sin(inc).*sin(AOP);

  R = vcat(r_1', r_2', r_3')
  V = vcat(v_1', v_2', v_3')

  [R,V, sma,ecc,inc,RAAN,AOP,M]
end
