classdef Node
    % Node: node.
    properties
        coord (:, :); % Coordinate of node (by column).
        type (1, :); % Type of node.
        % type = 0: interior node (inside domain).
        % type = +1/2/3/4: boundary node (on specific boundary part).
        % type = -1/2/3/4: intersection node (at intersection of boundary parts).
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
