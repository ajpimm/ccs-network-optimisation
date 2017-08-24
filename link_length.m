function link_length = link_length(pos1,pos2)
% Calculate the straight-line distance between position pos1 and position
% pos2. pos1 and pos2 can be multiple rows long to deal with multiple pairs
% of nodes.

% Initialise the vector of link lengths.
link_length = zeros(size(pos1,1),1);

% Find the straight-line distance between each pair of nodes.
for i1 = 1:size(pos1,1)
    
    % Extract the latitudes and longitudes of the two nodes.
    lat1 = pos1(i1,2);
    lon1 = pos1(i1,1);
    lat2 = pos2(i1,2);
    lon2 = pos2(i1,1);

    % Calculate the distance between the current pair of nodes.
    link_length(i1) = pos2dist(lat1,lon1,lat2,lon2,2); % km
end