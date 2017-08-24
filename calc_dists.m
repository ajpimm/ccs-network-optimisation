function Y = calc_dists(source_coords,reservoir_coords)
% Calculate the distance vector to be used by the linkage function that
% performs the hierarchical clustering. Here, 'distance' between a pair of
% sources is defined not as Euclidean distance, but as the 'length
% reduction potential' for the pair of sources.

% Find the number of sources.
num_sources = size(source_coords,1);

% Evaluate the length reduction potential for each pair of nodes.
k = 1;
for i = 1:num_sources-1
    Y(k:(k+num_sources-i-1)) = dist_fun(source_coords(i,:),source_coords((i+1):num_sources,:),...
        reservoir_coords);
    k = k + (num_sources-i);
end
