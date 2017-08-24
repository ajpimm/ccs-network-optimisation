function [total_cost,pipe_cost,pipe_length,pipe_diameter] = cost_fun_w_mass_flows(junction_coords,edges,source_coords,reservoir_coords,pipe_flows,scaling)
% Calculate the cost of the topology defined within the edges matrix.

% Form a vector of node positions, then reshape it into an array.
all_coords = [source_coords;junction_coords;reservoir_coords];
all_coords_mat = reshape(all_coords,2,[])';

% Count the number of edges and initialise results vectors.
num_edges = size(edges,1);
pipe_length = zeros(size(edges,1),1);
pipe_diameter = zeros(size(edges,1),1);
pipe_cost = zeros(size(edges,1),1);

% Step through each edge and add the cost of that edge to the total cost.
total_cost = 0;
for edge_num = 1:num_edges
    
    % Extract the latitude and longitude of the current pipe.
    n1 = edges(edge_num,1);
    n2 = edges(edge_num,2);
    lat1 = all_coords_mat(n1,2);
    lon1 = all_coords_mat(n1,1);
    lat2 = all_coords_mat(n2,2);
    lon2 = all_coords_mat(n2,1);
    
    % Calculate the length and required diameter of the current pipe.
    pipe_length(edge_num,1) = pos2dist(lat1,lon1,lat2,lon2,2); % km
    pipe_diameter(edge_num,1) = required_pipe_diameter(pipe_flows(edge_num)); % m
    
    % Add the cost of the current pipe to the total cost.
    pipe_cost(edge_num,1) = calc_pipe_costs(pipe_diameter(edge_num,1),pipe_length(edge_num,1),scaling);
    total_cost = total_cost + pipe_cost(edge_num,1);
end
