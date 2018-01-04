
clear
close all
clc

tic

% Create a random set of source positions and set the reservoir position.
% num_sources = 15;
% source_coords = 10*rand(num_sources,2);
% source_flows = 10*rand(num_sources,1); clear num_sources
% reservoir_coords = [5,5];

% Load a saved set of source and reservoir coordinates.
cd('Data')
loaded_data = importdata('GATEWAY_Case_A.txt','\t');
cd('..')
source_coords = loaded_data.data(1:end-1,1:2); % Long-lat
source_flows = loaded_data.data(1:end-1,3)*0.9/(365*24); % tCO2/hr
reservoir_coords = loaded_data.data(end,1:2); % Long-lat
clear loaded_data

% Carry out optimisation using the hierarchical clustering method, trying
% out multiple methods of measuring distance between clusters.
method = 'single'; % Single linkage clustering.
[edges_single,all_coords_single,cost_single] = ...
    hierarchical_clustering(source_coords,reservoir_coords,source_flows,method);

method = 'complete'; % Complete linkage clustering.
[edges_complete,all_coords_complete,cost_complete] = ...
    hierarchical_clustering(source_coords,reservoir_coords,source_flows,method);

method = 'weighted'; % Weighted average clustering (WPGMA).
[edges_weighted,all_coords_weighted,cost_weighted] = ...
    hierarchical_clustering(source_coords,reservoir_coords,source_flows,method);

method = 'average'; % Unweighted average clustering (UPGMA).
[edges_average,all_coords_average,cost_average] = ...
    hierarchical_clustering(source_coords,reservoir_coords,source_flows,method);

% Find which inter-cluster distance measure gives the lowest cost network.
[cost_min,min_idx] = min([cost_single,cost_complete,cost_weighted,...
    cost_average]);
switch min_idx
    case 1
        edges = edges_single;
        all_coords = all_coords_single;
        method_used = 'single';
    case 2
        edges = edges_complete;
        all_coords = all_coords_complete;
        method_used = 'complete';
    case 3
        edges = edges_weighted;
        all_coords = all_coords_weighted;
        method_used = 'weighted';
    case 4
        edges = edges_average;
        all_coords = all_coords_average;
        method_used = 'average';
end

% Create a graph object defining the lowest-cost network.
all_coords_x = all_coords(:,1); all_coords_y = all_coords(:,2);
EdgeTable = table(edges,'VariableNames',{'EndNodes'});
NodeTable = table(all_coords,'VariableNames',{'Coords'});
G = graph(EdgeTable,NodeTable);

% Plot the lowest-cost network.
figure
hold on
plot(G,'XData',all_coords_x,'YData',all_coords_y,'NodeLabel',[])
scatter(reservoir_coords(1), reservoir_coords(2), 'r+', 'LineWidth', 1.5)
hold off
xlabel('Longitude')
ylabel('Latitude')
grid on
axis equal
source_x = source_coords(:,1); source_y = source_coords(:,2);
text(source_x, source_y, cellstr(num2str((1:length(source_x))')))
title(sprintf('Cost = $%.2fm. Method used = %s',cost_min/1e6,method_used))

% Plot a dendrogram of the lowest-cost network.
% figure
% Y = calc_dists(source_coords,reservoir_coords);
% Z = linkage(Y,method);
% dendrogram(Z)
% grid on
% xlabel('Source node #')
% ylabel('Length reduction potential')

% Calculate the total cost of the network and cost breakdown.
final_cost_calculations

% Show the user the total runtime.
fprintf('Total runtime = %.1f s\n',toc)
