classdef Node
    % Node: node.
    properties
        coord (:, :); % Coordinate of node (by column).
        type (1, :); % Type of node.
        % For example: for 2D rectangular domain,
        % type = 0: interior node.
        % type = +1/2/3/4: lower/right/upper/left boundary node.
        % type = -1/2/3/4: lower-left/lower-right/upper-right/upper-left corner node.
    end
    properties (Dependent)
        dim; % Dimension.
        nNode; % Number of nodes.
    end
    methods
        % Constructor.
        function node = Node(coord, type)
            arguments
                coord (:, :);
                type (1, :) = [];
            end
            if ~isempty(type)
                assert(size(coord, 2) == size(type, 2));
            end
            node.coord = coord;
            node.type = type;
        end
        % Get functions.
        function dim = get.dim(node)
            dim = size(node.coord, 1);
        end
        function nNode = get.nNode(node)
            nNode = size(node.coord, 2);
        end
    end
end
