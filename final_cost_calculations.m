
edge_coords = zeros(size(edges,1),4);

for i1 = 1:size(edges,1)
    edge_coords(i1,:) = [all_coords(edges(i1,1),2),all_coords(edges(i1,1),1),all_coords(edges(i1,2),2),all_coords(edges(i1,2),1)];
end

num_sources = size(source_coords,1);
num_junctions = size(all_coords,1) - num_sources - 1;

% Find the pipe flows.
pipe_flows = get_pipe_flows(edges,source_flows);

% Calculate the cost of each pipe and total cost.
junction_coords = all_coords(num_sources+1:end-1,:);
[total_cost,pipe_costs,pipe_lengths,pipe_diameters]...
    = cost_fun_w_mass_flows(reshape(junction_coords',[],1),edges,reshape(source_coords',[],1),...
    reshape(reservoir_coords',[],1),pipe_flows,'1');

total_cost = total_cost + 1244724 + 111907;

% Calculate how many booster stations are needed along each pipe.
for i1 = 1:length(pipe_lengths)
    [outlet_cond(i1,:),booster_station(i1,:)] = pipe_pressure_drop_func(pipe_lengths(i1),pipe_flows(i1));
end

table_for_spreadsheet = [edge_coords,pipe_lengths,pipe_flows,pipe_diameters,pipe_costs];

type_markers = [ones(num_sources,1);2*ones(size(all_coords,1)-num_sources-1,1);3];
table_for_GIS = [(1:size(all_coords,1))',all_coords(:,2),all_coords(:,1),type_markers];
