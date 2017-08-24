function pipe_flows = get_pipe_flows(edges,source_flows)
% Runs a tree search algorithm to calculate the flows in each pipe of the
% tree network.

[pipe_flows,checked_nodes,downstream_edges] = ...
    treesearch(max(max(edges)),edges,[],[],0*edges(:,1),source_flows);