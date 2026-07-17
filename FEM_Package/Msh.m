classdef Msh
    % Msh: mesh.

    % Dim |  Point  |  Line  |  Face  | Volume
    %-----+---------+--------+--------+-------
    %  1  |  Node   |  Elem  |        |
    %  2  |  Node   |  Edge  |  Elem  |
    %  3  |  Node   |  Edge  |  Face  | Elem

    properties
        type {mustBeMember(type, ["VOID", "D1", "D2T"])} = "VOID"; % Type of mesh.
        % D: dimension.
        % T: triangle.
        node Node; % Node.
        elem Elem; % Element.
        edge Edge; % Edge.
    end
    properties (Dependent)
        dim; % Dimension.
        nNode; % Number of nodes.
        nElem; % Number of elements.
        nEdge; % Number of edges.
    end
    methods
        %% Constructor.
        function msh = Msh(type, node, elem, edge)
            arguments
                type {mustBeMember(type, ["D1", "D2T"])};
                node Node;
                elem Elem;
                edge Edge = Edge.empty;
            end
            switch type
                case "D1"
                    assert(node.dim == 1);
                    assert(elem.nNode == 2);
                    msh.type = type;
                    msh.node = node;
                    msh.elem = elem;
                case "D2T"
                    assert(node.dim == 2);
                    assert(elem.nNode == 3);
                    assert(isempty(elem.edge) || elem.nEdge == 3);
                    if ~isempty(edge)
                        assert(edge.nNode == 2);
                        assert(isempty(edge.elem) || edge.nElem == 2);
                    end
                    msh.type = type;
                    msh.node = node;
                    msh.elem = elem;
                    msh.edge = edge;
            end
        end
        % Get functions.
        function dim = get.dim(msh)
            dim = msh.node.dim;
        end
        function nNode = get.nNode(msh)
            nNode = msh.node.nNode;
        end
        function nElem = get.nElem(msh)
            nElem = msh.elem.nElem;
        end
        function nEdge = get.nEdge(msh)
            nEdge = msh.edge.nEdge;
        end
        %% Public functions.
        function nEnt = nEnt(msh, EntDim)
            % Msh.nEnt: get number of mesh entities of specified dimension.
            % EntDim = 0: point.
            % EntDim = 1: line.
            % EntDim = 2: face.
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2])};
            end
            switch msh.dim
                case 1
                    switch EntDim
                        case 0
                            nEnt = msh.node.nNode;
                        case 1
                            nEnt = msh.elem.nElem;
                        otherwise
                            error('Invalid EntDim.');
                    end
                case 2
                    switch EntDim
                        case 0
                            nEnt = msh.node.nNode;
                        case 1
                            nEnt = msh.edge.nEdge;
                        case 2
                            nEnt = msh.elem.nElem;
                        otherwise
                            error('Invalid EntDim.');
                    end
            end
        end
        function figure(msh)
            % Msh.figure: plot mesh.
            switch msh.type
                case "D2T"
                    figure;
                    hold on;
                    for iElem = 1:msh.nElem
                        vert_x = msh.node.coord(1, msh.elem.node(:, iElem));
                        vert_y = msh.node.coord(2, msh.elem.node(:, iElem));
                        patch('XData', vert_x, 'YData', vert_y, 'FaceColor', 'none', 'LineWidth', 0.5, 'EdgeAlpha', 0.25);
                        text(mean(vert_x), mean(vert_y), num2str(iElem), 'Color', "#0072BD", 'FontSize', 8, 'HorizontalAlignment', 'center');
                    end
                    for iEdge = 1:msh.nEdge
                        vert_x = msh.node.coord(1, msh.edge.node(:, iEdge));
                        vert_y = msh.node.coord(2, msh.edge.node(:, iEdge));
                        arrow_dx = (vert_x(2) - vert_x(1)) / 5;
                        arrow_dy = (vert_y(2) - vert_y(1)) / 5;
                        quiver(mean(vert_x), mean(vert_y), arrow_dx, arrow_dy, 'Color', "#D95319", 'MaxHeadSize', 1, 'LineWidth', 0.5);
                        text(mean(vert_x), mean(vert_y), num2str(iEdge), 'Color', "#D95319", 'FontSize', 8, 'HorizontalAlignment', 'center');
                    end
                    for iNode = 1:msh.nNode
                        text(msh.node.coord(1, iNode), msh.node.coord(2, iNode), num2str(iNode), 'Color', "#77AC30", 'FontSize', 8, 'HorizontalAlignment', 'center');
                    end
                    axis equal;
                    hold off;
            end
        end
    end
    %% Static functions.
    methods (Static)
        function msh = auto(type, node, ElNode)
            % Msh.auto: auto generate mesh information.
            arguments (Input)
                type {mustBeMember(type, "D2T")};
                node Node;
                ElNode (:, :); % Indices of vertices of element (by column).
                % Nodes are counter-clockwise in element.
            end
            arguments (Output)
                msh Msh;
            end
            switch type
                case "D2T"
                    % Nodes are counter-clockwise in element.
                    for iElem = 1:size(ElNode, 2)
                        ElVert = node.coord(:, ElNode(:, iElem));
                        area = det(ElVert(:, 2:3) - ElVert(:, 1));
                        if area < 0
                            ElNode([2,3], iElem) = ElNode([3,2], iElem);
                        end
                    end
                    % Edge.node.
                    ElEgNode = [ElNode([1, 2], :), ElNode([2, 3], :), ElNode([3, 1], :)];
                    ElEgSgn = ones(1, size(ElEgNode, 2));
                    [ElEgNode, sortIdx] = sort(ElEgNode, 1);
                    ElEgSgn(sortIdx(1, :) > sortIdx(2, :)) = -1;
                    [EgNode, ~, ElEg2Eg] = unique(ElEgNode', 'rows');
                    EgNode = EgNode'; nEdge = size(EgNode, 2);
                    % Edge.type.
                    EgType = zeros(1, nEdge);
                    for iEdge = 1:nEdge
                        EgNdType = node.type(EgNode(:, iEdge));
                        if all(EgNdType ~= 0)
                            if all(EgNdType > 0)
                                if EgNdType(1) == EgNdType(2)
                                    EgType(iEdge) = EgNdType(1);
                                end
                            elseif any(EgNdType > 0)
                                EgType(iEdge) = EgNdType(EgNdType > 0);
                            else
                                error('Invalid edge type.');
                            end
                        end
                    end
                    % Edge is counter-clockwise in boundary element.
                    for iEdge = 1:nEdge
                        if EgType(iEdge) ~= 0
                            iElEg = find(ElEg2Eg == iEdge);
                            if ElEgSgn(iElEg) == -1
                                ElEgNode([1, 2], iElEg) = ElEgNode([2, 1], iElEg);
                                ElEgSgn(iElEg) = 1;
                                EgNode([1, 2], iEdge) = EgNode([2, 1], iEdge);
                            end
                        end
                    end
                    % Elem.edge.
                    nElem = size(ElNode, 2);
                    ElEdge = reshape(ElEg2Eg' .* ElEgSgn, nElem, 3)';
                    % Edge.elem.
                    EgElem = zeros(2, nEdge);
                    for iElem = 1:nElem
                        for iElEg = 1:3
                            iEdge = abs(ElEdge(iElEg, iElem));
                            if EgElem(1, iEdge) == 0
                                EgElem(1, iEdge) = iElem * sign(ElEdge(iElEg, iElem));
                            else
                                EgElem(2, iEdge) = iElem * sign(ElEdge(iElEg, iElem));
                            end
                        end
                    end
                    % Elem.type.
                    ElType = zeros(1, nElem);
                    nBdType = length(unique(node.type(node.type > 0)));
                    for iElem = 1:nElem
                        ElNdType = node.type(ElNode(:, iElem));
                        if any(ElNdType < 0)
                            assert(sum(ElNdType < 0) == 1);
                            ElType(iElem) = ElNdType(ElNdType < 0);
                        elseif any(ElNdType > 0)
                            ElNdType = unique(ElNdType(ElNdType > 0));
                            if isscalar(ElNdType)
                                ElType(iElem) = ElNdType(1);
                            elseif length(ElNdType) == 2
                                if ElNdType(2) - ElNdType(1) == 1
                                    ElType(iElem) =- ElNdType(2);
                                elseif ElNdType(2) - ElNdType(1) == nBdType - 1
                                    ElType(iElem) =- ElNdType(1);
                                else
                                    error('Invalid element type.');
                                end
                            else
                                error('Invalid element type.');
                            end
                        end
                    end
                    elem = Elem(ElNode, ElEdge, ElType);
                    edge = Edge(EgNode, EgElem, EgType);
                    msh = Msh(type, node, elem, edge);
            end
        end
    end
end
