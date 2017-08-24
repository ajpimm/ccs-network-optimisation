function pipe_cost = calc_pipe_costs(pipe_diameter,pipe_length,scaling)
% Calculates the capital costs of a set of pipes given their diameters and
% lengths. Cost equations are taken from the IEA Pipeline Infrastructure
% 2013 report, based on the National Energy Technology Laboratory 2013. All
% costs are given in USD and based on 2011 costings.

% Convert to inches and miles.
pipe_diameter_inch = 39.37*pipe_diameter;
pipe_length_miles = 0.6214*pipe_length;

% Value by which to scale the cost overheads for each pipe. Half-sigmoid
% function (tanh) used.
switch scaling
    case 'half-sigmoid'
        overhead_scaler = tanh(pipe_length_miles/0.25);
    case '1'
        overhead_scaler = 1;
    otherwise
        error('Overhead scaling factor not understood')
end

materials_cost = 70350*overhead_scaler + 2.01*pipe_length_miles...
    .*(330.5*(pipe_diameter_inch.^2)+686.7*pipe_diameter_inch+26960);
labour_cost = 371850*overhead_scaler + 2.01*pipe_length_miles...
    .*(343.2*(pipe_diameter_inch.^2)+2074*pipe_diameter_inch+170013);
miscellaneous_costs = 147250*overhead_scaler...
    + 1.55*pipe_length_miles.*(8471*pipe_diameter_inch+7234);
right_of_way_cost = 51200*overhead_scaler...
    + 1.28*pipe_length_miles.*(577*pipe_diameter_inch+29788);
                
pipe_cost = materials_cost + labour_cost + miscellaneous_costs + right_of_way_cost;
