classdef Elem
    % Elem: element.
    properties
        node (:, :); % Indices of vertices of element (by column).
        % Nodes are counter-clockwise in element.
        edge (:, :); % Indices of edges of element (by column).
        % edge > 0: edge is counter-clockwise in element.
        % edge < 0: edge is clockwise in element.
        % Edge order matches node order: [node1, edge1, node2, edge2, ...].
        type (1, :); % Type of element.
        % For example: for 2D rectangular domain,
        % type = 0: interior element.
        % type = + 1/2/3/4: lower/right/upper/left boundary element.
        % type = - 1/2/3/4: lower-left/lower-right/upper-right/upper-left corner element.
    end
    properties (Dependent)
        nElem; % Number of elements.
        nNode; % Number of vertices of element.
        nEdge; % Number of edges of element.
    end
    methods
        % Constructor.
        function elem = Elem(node, edge, type)
            arguments
                node (:, :);
                edge (:, :) = [];
                type (1, :) = [];
            end
            if ~isempty(edge)
                assert(size(node, 2) == size(edge, 2));
            end
            if ~isempty(type)
                assert(size(node, 2) == size(type, 2));
            end
            elem.node = node;
            elem.edge = edge;
            elem.type = type;
        end
        % Get functions.
        function nElem = get.nElem(elem)
            nElem = size(elem.node, 2);
        end
        function nNode = get.nNode(elem)
            nNode = size(elem.node, 1);
        end
        function nEdge = get.nEdge(elem)
            nEdge = size(elem.edge, 1);
        end
    end
end
