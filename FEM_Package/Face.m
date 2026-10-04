classdef Face
    % Face: face.
    properties
        node (:, :); % Indices of vertices of face (by column).
        % Nodes are counter-clockwise viewed from outside of the first connected element (normal is outward of it).
        % For boundary face, normal is outward of the domain.
        elem (:, :); % Indices of connected elements of face (by column).
        % elem > 0: normal of face is outward of element.
        % elem < 0: normal of face is inward of element.
        % elem = 0: no connected element.
        edge (:, :); % Indices of edges of face (by column).
        % Edges of face: [node1, node2], [node2, node3], [node3, node1].
        % edge > 0: edge has the same orientation as in face.
        % edge < 0: edge has the opposite orientation.
        type (1, :); % Type of face.
        % For example: for 3D cuboid domain,
        % type = 0: interior face.
        % type = 1/2/3/4/5/6: x-/x+/y-/y+/z-/z+ boundary face.
    end
    properties (Dependent)
        nFace; % Number of faces.
        nNode; % Number of vertices of face.
        nElem; % Number of connected elements of face.
        nEdge; % Number of edges of face.
    end
    methods
        % Constructor.
        function face = Face(node, elem, edge, type)
            arguments
                node (:, :);
                elem (:, :) = [];
                edge (:, :) = [];
                type (1, :) = [];
            end
            if ~isempty(elem)
                assert(size(node, 2) == size(elem, 2));
            end
            if ~isempty(edge)
                assert(size(node, 2) == size(edge, 2));
            end
            if ~isempty(type)
                assert(size(node, 2) == size(type, 2));
            end
            face.node = node;
            face.elem = elem;
            face.edge = edge;
            face.type = type;
        end
        % Get functions.
        function nFace = get.nFace(face)
            nFace = size(face.node, 2);
        end
        function nNode = get.nNode(face)
            nNode = size(face.node, 1);
        end
        function nElem = get.nElem(face)
            nElem = size(face.elem, 1);
        end
        function nEdge = get.nEdge(face)
            nEdge = size(face.edge, 1);
        end
    end
end
