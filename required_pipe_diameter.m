function safety_diam = required_pipe_diameter(max_flow)
% Calculate the required diameter of a pipe carrying liquid CO2, with inlet
% flow rate given by max_flow (in tCO2/h). Assumes a pressure of 12 MPa,
% temperature of 15 C, and velocity of 2 m/s.

% Set parameters.
mass_rate_kg = max_flow*1000/3600; % kgCO2/s
velocity_in = 2; % m/s

density_init = 907.8021; % kg/m3

% Calculate pipe diameter.
pipe_diam = sqrt((4*mass_rate_kg)/(pi*density_init*velocity_in)); % m

% Increase by 3% to take account of pressure drop (found to be
% approximately correct in a numerical study).
safety_diam = 1.03*pipe_diam;
