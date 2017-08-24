function total_length = total_network_length(junction_coords,edges,source_coords,reservoir_coords)
% Calculate the total length of the topology defined within the edges matrix.

% Reshape the node position vector into an array.
all_coords = [source_coords;junction_coords;reservoir_coords];
all_coords_mat = reshape(all_coords,2,[])';

% Count the number of edges.
num_edges = size(edges,1);

% Step through each edge and add the pipe length.
total_length = 0;
for edge_num = 1:num_edges
    n1 = edges(edge_num,1);
    n2 = edges(edge_num,2);
    total_length = total_length + link_length(all_coords_mat(n1,:),all_coords_mat(n2,:));
end
