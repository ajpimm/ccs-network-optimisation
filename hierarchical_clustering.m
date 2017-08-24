function [edges,all_coords,cost] = hierarchical_clustering(source_coords,reservoir_coords,source_flows,method)
% Carries out hierarchical clustering of the sources and finds the optimal
% junction positions, ultimately giving a set of edges and node
% coordinates, as well as the network cost.
% 
% The positions of the sources are defined within source_coords and the CO2
% emissions from each source are defined within source_flows. The
% hierarchical clustering takes into account the position of the reservoir
% as defined within reservoir_coords. method is a string which specifies
% the method of evaluating the inter-cluster distances.

% Determine the number of sources and unpack the coordinates matrix into
% vectors containing the x (longitude) and y (latitude) components.
num_sources = size(source_coords,1);
source_x = source_coords(:,1);
source_y = source_coords(:,2);

% Calculate the 'distances' between all pairs of sources, where 'distance'
% is defined as the 'length reduction potential' of the pair of sources.
% The 'length reduction potential' between a pair of sources is the total
% length of the optimal route from the two sources to the reservoir
% (allowing a single junction) minus the total length of the direct routes
% from the two sources to the reservoir.
Y = calc_dists(source_coords,reservoir_coords);

% Cluster the sources through agglomerative hierarchical clustering, using
% 'length reduction potential' as the distance metric. This gives the
% topology to be used.
Z = linkage(Y,method);

% Find how many junction nodes are in the topology, and randomly position
% the junctions within the extents of the sources.
num_junctions = size(Z,1);
xmin = min(source_x); xmax = max(source_x); xrng = xmax - xmin;
ymin = min(source_y); ymax = max(source_y); yrng = ymax - ymin;
junction_x = xmin + xrng*rand(num_junctions,1);
junction_y = ymin + yrng*rand(num_junctions,1);
junction_coords = [junction_x, junction_y];

% Create an array containing the coordinates of the sources, randomly
% positioned junctions, and reservoir, and create separate variables
% containing the x and y components.
all_coords = [source_coords; junction_coords; reservoir_coords];
all_coords_x = all_coords(:,1); all_coords_y = all_coords(:,2);

% Create the edges array for the topology defined in Z. Each row of edges
% corresponds to a different edge. In each row, column 1 defines the child
% node in the tree and column 2 defines the parent node.
child = [Z(:,1); Z(:,2); num_sources + num_junctions];
parent = [(num_sources+1:num_sources+num_junctions)';...
    (num_sources+1:num_sources+num_junctions)';...
    1 + num_sources + num_junctions];
edges = [child,parent];

% Sort the edges so that the child nodes are in ascending order.
[~,idx] = sort(edges(:,1),'ascend');
edges = edges(idx,:);

% Calculate the flow through each node specified in all_coords.
node_flows = 0*all_coords_x;
node_flows(1:num_sources) = source_flows;
for junction_num = 1:num_junctions+1
    edges_idx = edges(:,2) == junction_num + num_sources;
    child_idx = edges(edges_idx,1);
    node_flows(junction_num + num_sources) = sum(node_flows(child_idx));
end

% Find the pipe flows.
pipe_flows = get_pipe_flows(edges,source_flows);

% Create lower and upper bounds for the junction position optimisation, as
% the extents of the sources and reservoir.
lb = [0*junction_x' + min(all_coords_x);...
    0*junction_y' + min(all_coords_y)];
lb = lb(:);
ub = [0*junction_x' + max(all_coords_x);...
    0*junction_y' + max(all_coords_y)];
ub = ub(:);

% Set the initial guess at the junction positions (using the randomly
% chosen positions).
x0 = junction_coords'; x0 = x0(:);

% Run the optimisation routine.
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
[junction_coords,cost] = fmincon(@(x)cost_fun_w_mass_flows(x,edges,...
    reshape(source_coords',[],1),reshape(reservoir_coords',[],1),...
    pipe_flows,'half-sigmoid'),x0,[],[],[],[],lb,ub,[],options);

% Add the fixed network costs.
cost = cost + 1244724 + 111907;

% Prepare arrays storing the coordinates of all the nodes.
junction_x = junction_coords(1:2:end);
junction_y = junction_coords(2:2:end);
all_coords = [source_coords; [junction_x,junction_y]; reservoir_coords];

%% Prepare for trying to remove particularly short pipes from the network.
% Find the shortest straight-line distance between all pairs of fixed
% nodes (used to check if an edge is particularly short).
fixed_coords = [all_coords(1:num_sources,:);all_coords(end,:)];
all_dists = zeros(((length(fixed_coords)-1)^2 + (length(fixed_coords)-1))/2,1);
all_dists_nodes = zeros(((length(fixed_coords)-1)^2 + (length(fixed_coords)-1))/2,2);
k = 1;
for i1 = 1:length(fixed_coords)
    for i2 = i1+1:length(fixed_coords)
        all_dists(k,1) = pos2dist(fixed_coords(i1,2),fixed_coords(i1,1),fixed_coords(i2,2),fixed_coords(i2,1),2);
        all_dists_nodes(k,:) = [i1,i2];
        k = k+1;
    end
end

% Calculate the length of each pipe.
junction_coords = all_coords(num_sources+1:end-1,:);
[old_cost,~,pipe_lengths,~]...
    = cost_fun_w_mass_flows(reshape(junction_coords',[],1),edges,reshape(source_coords',[],1),...
    reshape(reservoir_coords',[],1),pipe_flows,'1');

% Add the fixed network costs.
old_cost = old_cost + 1244724 + 111907;

%% Remove particularly short edges arising from the optimisation.
% Find the pipe flows.
pipe_flows = get_pipe_flows(edges,source_flows);

% Find the length of each pipe.
junction_coords = all_coords(num_sources+1:end-1,:);
[~,~,pipe_lengths,~] = cost_fun_w_mass_flows(...
    reshape(junction_coords',[],1),edges,reshape(...
    source_coords',[],1),reshape(reservoir_coords',[],1),pipe_flows,'half-sigmoid');

% Sort the edges, pipe lengths and flows to put the shortest pipes first.
[~,idx] = sort(pipe_lengths,'ascend');
edges = edges(idx,:);
pipe_lengths = pipe_lengths(idx,:);
pipe_flows = pipe_flows(idx,:);

pipe_num = 1;
while pipe_num <= size(edges,1)
    % Find out which nodes can be moved (i.e. which are junctions).
    can_be_moved = edges > num_sources & edges < size(all_coords,1);

    % Check if the current pipe can be removed, and try removing it if so.
    if sum(can_be_moved(pipe_num,:)) ~= 0
        % Both nodes can be moved.
        if sum(can_be_moved(pipe_num,:)) == 2
            % Arbitrarily pick one of the nodes to remove.
            remaining_node = edges(pipe_num,1);
            node_to_remove = edges(pipe_num,2);

        % Only one node (node 1) can be moved.
        elseif can_be_moved(pipe_num,1)
            remaining_node = edges(pipe_num,2);
            node_to_remove = edges(pipe_num,1);

        % Only one node (node 2) can be moved.
        elseif can_be_moved(pipe_num,2)
            remaining_node = edges(pipe_num,1);
            node_to_remove = edges(pipe_num,2);
        end

        % Remove one of the nodes.
        all_coords_new = all_coords;
        % If both nodes can be moved, place the single remaining node at
        % the midpoint of the pipe that is being removed.
        if sum(can_be_moved(pipe_num,:)) == 2
            all_coords_new(remaining_node,:) = 0.5*(all_coords_new(node_to_remove,:) + all_coords_new(remaining_node,:));
        end
        all_coords_new(node_to_remove,:) = [];
        
        % Connect all broken paths to the remaining node, and modify the
        % node numbering in the edges matrix accordingly.
        edges_new = edges;
        edges_new(edges_new == node_to_remove) = remaining_node;
        edges_new(edges_new > node_to_remove) = edges_new(edges_new > node_to_remove) - 1;

        % Remove the current edge.
        edges_new(pipe_num,:) = [];

        %% Re-run the optimisation on the new topology.
        all_coords_x = all_coords_new(:,1);
        all_coords_y = all_coords_new(:,2);

        % Sort the edges so that the child nodes are in ascending
        % order.
        [~,idx] = sort(edges_new(:,1),'ascend');
        edges_new = edges_new(idx,:);

        % Find the pipe flows.
        pipe_flows = get_pipe_flows(edges_new,source_flows);

        % Only try to optimise junction positions if there are junctions.
        num_junctions = size(all_coords_new,1) - num_sources - 1;
        if num_junctions ~= 0
            % Create lower and upper bounds for the junction position
            % optimisation, as the extents of the sources and reservoir.
            lb = [zeros(1,num_junctions) + min(all_coords_x);...
                zeros(1,num_junctions) + min(all_coords_y)];
            lb = lb(:);
            ub = [zeros(1,num_junctions) + max(all_coords_x);...
                zeros(1,num_junctions) + max(all_coords_y)];
            ub = ub(:);

            % Run the optimisation routine.
            x0 = all_coords_new(num_sources+1:end-1,:)'; x0 = x0(:);
            options = optimoptions('fmincon', 'Display', 'off','Algorithm', 'sqp');
            [junction_coords,~] = fmincon(@(x)cost_fun_w_mass_flows(x,edges_new,...
                reshape(source_coords',[],1),reshape(reservoir_coords',[],1),...
                pipe_flows,'half-sigmoid'),x0,[],[],[],[],lb,ub,[],options);
        else
            % Otherwise, set junction_coords to an empty array.
            junction_coords = [];
        end
            
        % Calculate the actual costs (i.e. with no scaling of fixed
        % pipe costs).
        [cost,~,pipe_lengths,~] = cost_fun_w_mass_flows(reshape(junction_coords',[],1),...
            edges_new,reshape(source_coords',[],1),reshape(reservoir_coords',[],1),pipe_flows,'1');

        % Add the fixed network costs.
        cost = cost + 1244724 + 111907;
        
        % Check if there has been an improvement in cost. If so, go back to
        % the first (i.e. shortest) pipe again, store the coordinates, set
        % old cost to be the best cost so far, and sort the pipes by length
        % again, otherwise move onto the next pipe.
        improvement = cost < old_cost;
        if improvement
            pipe_num = 1;
            edges = edges_new;
            
            % Prepare arrays storing coordinates.
            junction_x = junction_coords(1:2:end);
            junction_y = junction_coords(2:2:end);
            all_coords = [source_coords; [junction_x,junction_y];...
                reservoir_coords];
            
            % Set the old cost to be the current cost.
            old_cost = cost;
            
            % Sort the edges, pipe lengths and flows to put the shortest
            % pipes first.
            [~,idx] = sort(pipe_lengths,'ascend');
            edges = edges(idx,:);
            pipe_lengths = pipe_lengths(idx,:);
            pipe_flows = get_pipe_flows(edges,source_flows);
        else
            pipe_num = pipe_num + 1;
        end
    
    % Move onto the next pipe if the current pipe couldn't be removed (i.e.
    % because it was between two fixed nodes).
    else
        pipe_num = pipe_num + 1;
    end
end
    
% Find the pipe flows.
pipe_flows = get_pipe_flows(edges,source_flows);

% Run the optimisation routine once more on the final topology, with no
% overhead scaling.
num_junctions = size(all_coords,1) - num_sources - 1;
if num_junctions ~= 0
    % Create lower and upper bounds for the junction position optimisation,
    % as the extents of the sources and reservoir.
    all_coords_x = all_coords(:,1);
    all_coords_y = all_coords(:,2);
    lb = [zeros(1,num_junctions) + min(all_coords_x);...
        zeros(1,num_junctions) + min(all_coords_y)];
    lb = lb(:);
    ub = [zeros(1,num_junctions) + max(all_coords_x);...
        zeros(1,num_junctions) + max(all_coords_y)];
    ub = ub(:);
    x0 = all_coords(num_sources+1:end-1,:)'; x0 = x0(:);
    options = optimoptions('fmincon', 'Display', 'off','Algorithm', 'sqp');
    [junction_coords,cost] = fmincon(@(x)cost_fun_w_mass_flows(x,edges,...
        reshape(source_coords',[],1),reshape(reservoir_coords',[],1),...
        pipe_flows,'1'),x0,[],[],[],[],lb,ub,[],options);

    cost = cost + 1244724 + 111907;

    junction_x = junction_coords(1:2:end);
    junction_y = junction_coords(2:2:end);
    all_coords = [source_coords; [junction_x,junction_y];...
        reservoir_coords];
end

%% Fixing step for sub-optimal local minima for pairs of sources.
% Either both sources are connected together with a pipe, in which case one
% (and only one) of the sources will have other pipes attached to it, or
% the sources aren't connected together, in which case they will both be
% connected to the same node.

% If they're connected together, we could try moving all the pipes attached
% to one of the sources onto the other source. However, for now we'll do
% nothing.

% If they're not connected together, try connecting them together in the
% two ways possible.
connected_pairs = Z(:,1) <= num_sources & Z(:,2) <= num_sources;
num_connected_pairs = sum(connected_pairs);
for pair_num = 1:size(Z,1)
    if connected_pairs(pair_num)
        n1 = Z(pair_num,1);
        n2 = Z(pair_num,2);
        
        % Set connected_already back to 0.
        connected_already = 0;
        
        % Find the edges with n1 on, and all connecting nodes.
        n1_edges = find(edges(:,1) == n1);
        conn_nodes = edges(n1_edges,2);
        if isempty(n1_edges)
            n1_edges = find(edges(:,2) == n1);
            conn_nodes = edges(n1_edges,1);
        end
        % Check that n1 and n2 aren't connected together.
        if max(conn_nodes == n2)
            connected_already = 1;
        end
        
        n2_edges = find(edges(:,1) == n2);
        conn_nodes = edges(n2_edges,2);
        if isempty(n2_edges)
            n2_edges = find(edges(:,2) == n2);
            conn_nodes = edges(n2_edges,1);
        end
        % Check that n1 and n2 aren't connected together.
        if max(conn_nodes == n1)
            connected_already = 1;
        end
        
        % Only continue with this pair of nodes if they're not already
        % connected together.
        if ~connected_already
            % Connect n1 to n2, remove n1 to conn_node and calculate cost.
            edges_new1 = edges;
            edges_new1 = [edges_new1;n1,n2];
            edges_new1(n1_edges,:) = [];
            pipe_flows1 = get_pipe_flows(edges_new1,source_flows);
            [cost_new1,~,~,~] = cost_fun_w_mass_flows(reshape(junction_coords',[],1),...
                edges_new1,reshape(source_coords',[],1),reshape(reservoir_coords',[],1),pipe_flows1,'1');
            cost_new1 = cost_new1 + 1244724 + 111907;

            % Connect n1 to n2, remove n2 to conn_node and calculate cost.
            edges_new2 = edges;
            edges_new2 = [edges_new2;n1,n2];
            edges_new2(n2_edges,:) = [];
            try
                pipe_flows2 = get_pipe_flows(edges_new2,source_flows);
            catch
                true;
            end
            [cost_new2,~,~,~] = cost_fun_w_mass_flows(reshape(junction_coords',[],1),...
                edges_new2,reshape(source_coords',[],1),reshape(reservoir_coords',[],1),pipe_flows2,'1');
            cost_new2 = cost_new2 + 1244724 + 111907;

            % Find the lowest cost option.
            if cost_new1 < cost_new2
                if cost_new1 < cost
                    edges = edges_new1;
                    cost = cost_new1;
                end
            else
                if cost_new2 < cost
                    edges = edges_new2;
                    cost = cost_new2;
                end
            end
        end
    end
end
