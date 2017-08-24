function [pipe_flows,checked_nodes,downstream_edges] = treesearch(curr_node,edges,downstream_edges,checked_nodes,pipe_flows,source_flows)
% Process the current node, i.e. check if it is a source node and, if so,
% add its flow rate to all the downstream connected edges.
num_sources = length(source_flows);
if curr_node <= num_sources
    pipe_flows(downstream_edges) = pipe_flows(downstream_edges) + source_flows(curr_node);
end

% Add the current node to the list of checked nodes.
checked_nodes = [checked_nodes;curr_node];

% Find all connected upstream nodes.
% First find all connected nodes.
connected_nodes = edges(edges(:,1)==curr_node,2);
connected_nodes = [connected_nodes;edges(edges(:,2)==curr_node,1)];
% Then remove any nodes that have already been checked (i.e. aren't
% upstream).
upstream_nodes = connected_nodes;
for i1 = length(upstream_nodes):-1:1
    if max(checked_nodes==upstream_nodes(i1)) ~= 0
        upstream_nodes(i1) = [];
    end
end

% If there are any connected upstream nodes, carry out the tree search on
% the connected upstream nodes.
if ~isempty(upstream_nodes)
    downstream_edges = [downstream_edges;0];
    for i1 = 1:length(upstream_nodes)
        curr_edge = find((edges(:,1)==curr_node & edges(:,2)==upstream_nodes(i1))...
            | (edges(:,2)==curr_node & edges(:,1)==upstream_nodes(i1)));
        downstream_edges(end) = curr_edge;
        [pipe_flows,checked_nodes,downstream_edges] = ...
            treesearch(upstream_nodes(i1),edges,downstream_edges,...
            checked_nodes,pipe_flows,source_flows);
    end
    % Remove from the list of downstream edges.
    downstream_edges(end) = [];
end
