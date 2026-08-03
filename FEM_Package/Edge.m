classdef Edge
    % Edge: edge.
    properties
        node (:, :); % Indices of vertices of edge (by column).
        % For boundary edge, nodes are counter-clockwise in boundary element.
        % For interior edge, nodes are in arbitrary order.
        elem (:, :); % Indices of connected elements of edge (by column).
        % elem > 0: edge is counter-clockwise in element.
        % elem < 0: edge is clockwise in element.
        % elem = 0: no connected element.
        type (1, :); % Type of edge.
        % For example: for 2D rectangular domain,
        % type = 0: interior edge.
        % type = 1/2/3/4: lower/right/upper/left boundary edge.
    end
    properties (Dependent)
        nEdge; % Number of edges.
        nNode; % Number of vertices of edge.
        nElem; % Number of connected elements of edge.
    end
    methods
        % Constructor.
        function edge = Edge(node, elem, type)
            arguments
                node (:, :);
                elem (:, :) = [];
                type (1, :) = [];
            end
            if ~isempty(elem)
                assert(size(node, 2) == size(elem, 2));
            end
            if ~isempty(type)
                assert(size(node, 2) == size(type, 2));
            end
            edge.node = node;
            edge.elem = elem;
            edge.type = type;
        end
        % Get functions.
        function nEdge = get.nEdge(edge)
            nEdge = size(edge.node, 2);
        end
        function nNode = get.nNode(edge)
            nNode = size(edge.node, 1);
        end
        function nElem = get.nElem(edge)
            nElem = size(edge.elem, 1);
        end
    end
end
