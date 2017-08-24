function dist = dist_fun(XI,XJ,reservoir_coords)
% Calculate the 'length reduction potential' between pairs of sources,
% which is defined as the difference between the total length of the
% optimal route between the two sources and the reservoir and the total
% length of the direct route between the two sources and the reservoir.
% 
% XI contains the coordinates of a single source node.
% XJ contains the coordinates of any number of nodes. The length reduction
% potential is found for each node in XJ when paired up with the node in
% XI, and output as dist.
% reservoir_coords contains the coordinates of a single reservoir node.

% Transpose the reservoir coordinates.
reservoir_coords = reservoir_coords';

% Set node 1 to be the node in XI.
n1 = XI';

% Create an edge list, connecting the two sources (1 & 2) and reservoir (4)
% to a junction node (3).
edges = [1,3;2,3;3,4];

% Initialise the distance vector.
dist = zeros(size(XJ,1),1);

% Step through each node (row) in XJ and find the total length of the
% network that optimally connects that node and node 1 to the reservoir
% with a single junction, as well as the total length when directly
% connecting those two nodes to the reservoir.
for row_num = 1:size(XJ,1)
    % Set node 2 to be the current row of XJ.
    n2 = XJ(row_num,:)';
    
    % Create the source coordinates vector.
    source_coords = [n1;n2];
    
    % Create the fixed coordinates vector.
    fixed_coords = [source_coords;reservoir_coords];
    fixed_coords_x = fixed_coords(1:2:end);
    fixed_coords_y = fixed_coords(2:2:end);
    
    % Create an initial guess for the optimal location of the junction.
    j0 = (n1+n2)/2;
    
    % Set the lower and upper bounds for the position of the junction as
    % the extents of a box bounding the source nodes and the reservoir.
    lb = [min(fixed_coords_x); min(fixed_coords_y)];
    ub = [max(fixed_coords_x); max(fixed_coords_y)];
    
    % Optimise the position of the junction and find the network length.
    options = optimset('Display', 'off');
    [~,total_length_opt] = fmincon(@(x)total_network_length(x,edges,source_coords,...
        reservoir_coords),j0,[],[],[],[],lb,ub,[],options);
    
    %% Calculate the total direct length from the two nodes to the reservoir.
    total_length_direct = total_network_length([],[1,2],source_coords(1:2,:),reservoir_coords)...
        + total_network_length([],[1,2],source_coords(3:4),reservoir_coords);
    
    % Calculate the 'distance' between node 1 and node 2 as the length
    % reduction potential.
    dist(row_num,1) = total_length_opt - total_length_direct;
end
