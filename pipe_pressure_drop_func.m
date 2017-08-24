% MECH5825 Pipeline repressurisation stations function
% Created 15 June 2016 by Reyes Gonzalez Ferreras
% Density and viscosity obtained from the CoolProp Database and code
% Function calculates the outlet conditions of CO2 in a pipe, the
% number conditions and location of booster station
% pipe length must be entered in m and max_flow in tCO2/h

function [outlet_cond,booster_station] = ...
          pipe_pressure_drop_func(pipe_length,max_flow)

% Add Coolprop database
 addpath(genpath('C:\Users\preapi\Dropbox\GATEWAY\Hierarchical Clustering Method 6\CO2datamain'));

% Set Parameters
pipe_length_m = pipe_length*1000; % km to m
mass_rate_kg = max_flow*1000/3600; %kgCO2/s
recommended_press = 12000000; %Pa
pressure_min = 8500000; %Pa
pipe_segment_length = 1000; %m (1km)
pipe_segments = round(pipe_length_m/pipe_segment_length,0);
temp_cel = 15; % celcius, max temp of UK ground below 1m
temp_kel = 273.15 + temp_cel; %K
velocity_in = 2; %m/s 
pipe_roughness = 0.045*10^-3; %m for a steel pipe (pipelines.com)

% Calculate density and viscosity using the CoolProp database
density_init = CoolProp.PropsSI('D', 'T', temp_kel, 'P',...
                               recommended_press,'CarbonDioxide');  %kg/m3
dynamic_visc_init = CoolProp.PropsSI('V', 'T', temp_kel, 'P',...
                               recommended_press,'CarbonDioxide');  %Pa.s
kinematic_visc_init = dynamic_visc_init/density_init;   %m2/s

% Calculate pipe diameter, Reynolds number and Darcy
% friction factor. These are assumed constant for the entire pipe.
pipe_diam = sqrt((4*mass_rate_kg)/(pi*density_init*velocity_in));   %m
safety_diam = 1.03*pipe_diam;
Reynolds_num = velocity_in*safety_diam/kinematic_visc_init;
% Swamee-Jain equation for a full-flowing circular pipe
f_D = 1.325/((log((pipe_roughness/(3.7*safety_diam))+...
                   (5.74/(Reynolds_num.^(0.9)))))^2);

% Setting up values for loops
kilometer_count = 0;
repressurisation_count = 0;
pressure_in = zeros(1,pipe_segments);
density_in = zeros(1,pipe_segments);
pressure_change = zeros(1,pipe_segments);
pressure_out = zeros(1,pipe_segments);
velocity_out = zeros(1,pipe_segments);

% Initial pressure and density values 
pressure_in(1) = recommended_press;
density_in(1) = density_init;
pressure_change(1) = pipe_segment_length*f_D*(density_in(1)/2)*...
                    ((velocity_in^2)/safety_diam);
pressure_out(1) = recommended_press - pressure_change(1);

% Calculating the number of repressurisation stations needed, the
% velocity, pressure drop, pressure in and pressure out per km
    for i = 1:pipe_segments
        if pressure_out(i) > pressure_min
            kilometer_count = kilometer_count+1;
            pressure_in(i+1) = pressure_out(i);
            density_in(i+1) = CoolProp.PropsSI('D', 'T', temp_kel,...
                             'P',pressure_out(i),'CarbonDioxide');   %kg/m3
            velocity_out(i) = (4*mass_rate_kg)/(density_in(i+1)*pi*...
                                                (safety_diam^2));             
            pressure_change(i+1) = pipe_segment_length*f_D*...
                     (density_in(i+1)/2)*((velocity_out(i)^2)/safety_diam);
            pressure_out(i+1) = pressure_in(i+1) - pressure_change(i+1);
            
        else
            kilometer_count = kilometer_count+1;
            repressurisation_count = repressurisation_count+1;
            if repressurisation_count == 1
            repress_location = kilometer_count; % km
            end
            pressure_in(i+1) = recommended_press;
            density_in(i+1) = CoolProp.PropsSI('D', 'T', temp_kel,...
                           'P',pressure_in(i+1),'CarbonDioxide');   %kg/m3
            velocity_out(i) = (4*mass_rate_kg)/(density_in(i+1)*pi*...
                                                (safety_diam^2));           
            pressure_change(i+1) = pipe_segment_length*f_D*...
                     (density_in(i+1)/2)*((velocity_out(i)^2)/safety_diam);
            pressure_out(i+1) = pressure_in(i+1) - pressure_change(i+1);
                       
        end
    end

    if isempty(velocity_out)
        velocity_out = velocity_in;
    end
 
% Outputs
    if (repressurisation_count > 1) || (repressurisation_count == 1)
    booster_station = [repressurisation_count,repress_location,...
                      min(pressure_out),max(velocity_out),min(density_in)];
    else
    booster_station = [0,0,...
                      min(pressure_out),max(velocity_out),min(density_in)];
    end
    
outlet_cond = [pressure_out(end),velocity_out(end),density_in(end)];    
  
end