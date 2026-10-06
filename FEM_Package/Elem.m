classdef Elem
    % Elem: element.
    properties
        node (:, :); % Indices of vertices of element (by column).
        % 2D: nodes are counter-clockwise in element.
        % 3D: element is positively oriented, i.e. det([v2 - v1, v3 - v1, v4 - v1]) > 0.
        edge (:, :); % Indices of edges of element (by column).
        % 2D: edge > 0 (< 0): edge is counter-clockwise (clockwise) in element.
        % 2D: edge order matches node order: [node1, edge1, node2, edge2, ...].
        % 3D: edge > 0 (< 0): local edge [a, b] (see MshEnt) has the same (opposite) orientation as global edge.
        type (1, :); % Type of element.
        % For example: for 2D rectangular domain,
        % type = 0: interior element.
        % type = + 1/2/3/4: lower/right/upper/left boundary element.
        % type = - 1/2/3/4: lower-left/lower-right/upper-right/upper-left corner element.
        face (:, :); % Indices of faces of element (by column, 3D only).
        % face > 0 (< 0): normal of global face is outward (inward) of element.
        facePerm (:, :); % Permutation code of faces of element (by column, 3D only).
        % Local face nodes L and global face nodes G satisfy L = G(P(c, :)), where c is the code and
        % P = [1, 2, 3; 2, 3, 1; 3, 1, 2; 1, 3, 2; 3, 2, 1; 2, 1, 3] (c <= 3: even, face > 0; c >= 4: odd, face < 0).
    end
    properties (Dependent)
        nElem; % Number of elements.
        nNode; % Number of vertices of element.
        nEdge; % Number of edges of element.
        nFace; % Number of faces of element.
    end
    methods
        % Constructor.
        function elem = Elem(node, edge, type, face, facePerm)
            arguments
                node (:, :);
                edge (:, :) = [];
                type (1, :) = [];
                face (:, :) = [];
                facePerm (:, :) = [];
            end
            if ~isempty(edge)
                assert(size(node, 2) == size(edge, 2));
            end
            if ~isempty(type)
                assert(size(node, 2) == size(type, 2));
            end
            if ~isempty(face)
                assert(size(node, 2) == size(face, 2));
            end
            if ~isempty(facePerm)
                assert(isequal(size(facePerm), size(face)));
            end
            elem.node = node;
            elem.edge = edge;
            elem.type = type;
            elem.face = face;
            elem.facePerm = facePerm;
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
        function nFace = get.nFace(elem)
            nFace = size(elem.face, 1);
        end
    end
end
