classdef Msh
    % Msh: mesh.

    % Dim |  Point  |  Line  |  Face  | Volume
    %-----+---------+--------+--------+-------
    %  1  |  Node   |  Elem  |        |
    %  2  |  Node   |  Edge  |  Elem  |
    %  3  |  Node   |  Edge  |  Face  | Elem

    properties
        type {mustBeMember(type, ["VOID", "D1", "D2T", "D3T"])} = "VOID"; % Type of mesh.
        % D: dimension.
        % T: triangle (2D), tetrahedron (3D).
        node Node; % Node.
        elem Elem; % Element.
        edge Edge; % Edge.
        face Face; % Face.
    end
    properties (Dependent)
        dim; % Dimension.
        nNode; % Number of nodes.
        nElem; % Number of elements.
        nEdge; % Number of edges.
        nFace; % Number of faces.
        domn; % Domain of function on mesh.
        ElDomn; % Domain of element.
        FtDomn; % Domain of facet.
        FtRefDomn; % Reference domain of facet.
        facet; % Facet: mesh entity of codimension 1 (edge in 2D, face in 3D).
    end
    methods
        %% Constructor.
        function msh = Msh(type, node, elem, edge, face)
            arguments
                type {mustBeMember(type, ["D1", "D2T", "D3T"])};
                node Node;
                elem Elem;
                edge Edge = Edge.empty;
                face Face = Face.empty;
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
                case "D3T"
                    assert(node.dim == 3);
                    assert(elem.nNode == 4);
                    assert(isempty(elem.edge) || elem.nEdge == 6);
                    assert(isempty(elem.face) || elem.nFace == 4);
                    if ~isempty(edge)
                        assert(edge.nNode == 2);
                    end
                    if ~isempty(face)
                        assert(face.nNode == 3);
                        assert(isempty(face.elem) || face.nElem == 2);
                        assert(isempty(face.edge) || face.nEdge == 3);
                    end
                    msh.type = type;
                    msh.node = node;
                    msh.elem = elem;
                    msh.edge = edge;
                    msh.face = face;
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
        function nFace = get.nFace(msh)
            nFace = msh.face.nFace;
        end
        function domn = get.domn(msh)
            switch msh.type
                case "D2T"
                    domn = "D2";
                case "D3T"
                    domn = "D3";
                otherwise
                    domn = "VOID";
            end
        end
        function ElDomn = get.ElDomn(msh)
            switch msh.type
                case "D2T"
                    ElDomn = "D2T";
                case "D3T"
                    ElDomn = "D3T";
                otherwise
                    ElDomn = "VOID";
            end
        end
        function FtDomn = get.FtDomn(msh)
            switch msh.type
                case "D2T"
                    FtDomn = "D2L";
                case "D3T"
                    FtDomn = "D3F";
                otherwise
                    FtDomn = "VOID";
            end
        end
        function FtRefDomn = get.FtRefDomn(msh)
            FtRefDomn = MshEnt(msh.FtDomn).dual.type;
        end
        function facet = get.facet(msh)
            switch msh.type
                case "D2T"
                    facet = msh.edge;
                case "D3T"
                    facet = msh.face;
                otherwise
                    facet = [];
            end
        end
        %% Public functions.
        function nEnt = nEnt(msh, EntDim)
            % Msh.nEnt: get number of mesh entities in specified dimension.
            % EntDim = 0: point.
            % EntDim = 1: line.
            % EntDim = 2: face.
            % EntDim = 3: volume.
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2, 3])};
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
                case 3
                    switch EntDim
                        case 0
                            nEnt = msh.node.nNode;
                        case 1
                            nEnt = msh.edge.nEdge;
                        case 2
                            nEnt = msh.face.nFace;
                        case 3
                            nEnt = msh.elem.nElem;
                    end
            end
        end
        function ent = ent(msh, EntDim)
            % Msh.ent: get mesh entities in specified dimension.
            % EntDim = 0: point.
            % EntDim = 1: line.
            % EntDim = 2: face.
            % EntDim = 3: volume.
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2, 3])};
            end
            switch msh.dim
                case 1
                    switch EntDim
                        case 0
                            ent = msh.node;
                        case 1
                            ent = msh.elem;
                        otherwise
                            error('Invalid EntDim.');
                    end
                case 2
                    switch EntDim
                        case 0
                            ent = msh.node;
                        case 1
                            ent = msh.edge;
                        case 2
                            ent = msh.elem;
                        otherwise
                            error('Invalid EntDim.');
                    end
                case 3
                    switch EntDim
                        case 0
                            ent = msh.node;
                        case 1
                            ent = msh.edge;
                        case 2
                            ent = msh.face;
                        case 3
                            ent = msh.elem;
                    end
            end
        end
        function EntIdx = bdEnt(msh, EntDim, type)
            % Msh.bdEnt: get indices of mesh entities on facets with specified types.
            % Types are labeled on facets (edges in 2D, faces in 3D): entities of lower dimension are those on
            % these facets, and elements are those connected to these facets.
            % EntDim = 0: point.
            % EntDim = 1: line.
            % EntDim = 2: face.
            % EntDim = 3: volume.
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2, 3])};
                type (1, :); % Types of facets.
            end
            assert(ismember(msh.type, ["D2T", "D3T"]));
            facet = msh.facet;
            FtIdx = find(ismember(facet.type, type));
            switch EntDim
                case msh.dim - 1
                    EntIdx = FtIdx;
                case 0
                    EntIdx = unique(facet.node(:, FtIdx))';
                case msh.dim
                    CnElem = facet.elem(:, FtIdx);
                    EntIdx = unique(abs(CnElem(CnElem ~= 0)))';
                case 1
                    % Edges of faces (3D).
                    EntIdx = unique(abs(facet.edge(:, FtIdx)))';
                otherwise
                    error('Invalid EntDim.');
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
                case "D3T"
                    % Boundary faces colored by face type.
                    figure;
                    isBd = sum(msh.face.elem ~= 0, 1) == 1;
                    patch('Faces', msh.face.node(:, isBd)', 'Vertices', msh.node.coord', 'FaceVertexCData', msh.face.type(isBd)', ...
                        'FaceColor', 'flat', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 0.5);
                    axis equal;
                    view(3);
            end
        end
    end
    %% Static functions.
    methods (Static)
        function msh = auto(type, node, ElNode, options)
            % Msh.auto: auto generate mesh information.
            arguments (Input)
                type {mustBeMember(type, ["D2T", "D3T"])};
                node Node;
                ElNode (:, :); % Indices of vertices of element (by column).
                % Orientation of element is corrected automatically.
                options.FcType = []; % Type of boundary face (3D): function handle @(cen, nor) of centroids and
                % outward unit normals of boundary faces (by column), returning types (by row). Default type is 1.
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
                                    ElType(iElem) = -ElNdType(2);
                                elseif ElNdType(2) - ElNdType(1) == nBdType - 1
                                    ElType(iElem) = -ElNdType(1);
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
                case "D3T"
                    nElem = size(ElNode, 2);
                    X = node.coord;
                    % Element is positively oriented.
                    v1 = X(:, ElNode(1, :));
                    vol = dot(X(:, ElNode(2, :)) - v1, cross(X(:, ElNode(3, :)) - v1, X(:, ElNode(4, :)) - v1));
                    isNeg = vol < 0;
                    ElNode([3, 4], isNeg) = ElNode([4, 3], isNeg);
                    % Edge: global edge is oriented from smaller to larger node index.
                    LcEgNode = MshEnt("D3T").edge.node;
                    ElEgNode = reshape(ElNode(LcEgNode(:), :), 2, 6 * nElem);
                    [ElEgNode, sortIdx] = sort(ElEgNode, 1);
                    ElEgSgn = ones(1, 6 * nElem);
                    ElEgSgn(sortIdx(1, :) > sortIdx(2, :)) = -1;
                    [EgNode, ~, ElEg2Eg] = unique(ElEgNode', 'rows');
                    EgNode = EgNode'; nEdge = size(EgNode, 2);
                    ElEdge = reshape(ElEg2Eg' .* ElEgSgn, 6, nElem);
                    % Face: global face nodes follow local face nodes in the first connected element,
                    % i.e. normal of global face is outward of the first connected element.
                    LcFcNode = MshEnt("D3T").face.node;
                    ElFcNode = reshape(ElNode(LcFcNode(:), :), 3, 4 * nElem);
                    [~, iFirst, ElFc2Fc] = unique(sort(ElFcNode, 1)', 'rows');
                    FcNode = ElFcNode(:, iFirst); nFace = size(FcNode, 2);
                    assert(all(accumarray(ElFc2Fc, 1) <= 2), 'Face is shared by more than two elements.');
                    % Elem.face and Elem.facePerm: local face nodes L = G(perm(c, :)) for global face nodes G.
                    perm = [1, 2, 3; 2, 3, 1; 3, 1, 2; 1, 3, 2; 3, 2, 1; 2, 1, 3];
                    GlFcNode = FcNode(:, ElFc2Fc);
                    ElFcPerm = zeros(1, 4 * nElem);
                    for iPerm = 1:6
                        ElFcPerm(all(ElFcNode == GlFcNode(perm(iPerm, :), :), 1)) = iPerm;
                    end
                    assert(all(ElFcPerm > 0));
                    ElFcSgn = 1 - 2 * (ElFcPerm > 3);
                    ElFace = reshape(ElFc2Fc' .* ElFcSgn, 4, nElem);
                    ElFcPerm = reshape(ElFcPerm, 4, nElem);
                    % Face.elem.
                    ElIdx = repelem(1:nElem, 4);
                    isFirst = false(1, 4 * nElem);
                    isFirst(iFirst) = true;
                    FcElem = zeros(2, nFace);
                    FcElem(1, ElFc2Fc(isFirst)) = ElIdx(isFirst) .* ElFcSgn(isFirst);
                    FcElem(2, ElFc2Fc(~isFirst)) = ElIdx(~isFirst) .* ElFcSgn(~isFirst);
                    % Face.edge: edges [node1, node2], [node2, node3], [node3, node1] of face.
                    FcEgNode = reshape(FcNode([1, 2, 2, 3, 3, 1], :), 2, 3 * nFace);
                    [~, FcEg2Eg] = ismember(sort(FcEgNode, 1)', EgNode', 'rows');
                    FcEgSgn = 1 - 2 * (FcEgNode(1, :) > FcEgNode(2, :));
                    FcEdge = reshape(FcEg2Eg' .* FcEgSgn, 3, nFace);
                    % Face.type: interior face is 0, boundary face is given by `FcType`.
                    FcType = zeros(1, nFace);
                    isBd = FcElem(2, :) == 0;
                    if isempty(options.FcType)
                        FcType(isBd) = 1;
                    else
                        BdFcNd1 = X(:, FcNode(1, isBd));
                        BdFcNd2 = X(:, FcNode(2, isBd));
                        BdFcNd3 = X(:, FcNode(3, isBd));
                        nor = cross(BdFcNd2 - BdFcNd1, BdFcNd3 - BdFcNd1);
                        nor = nor ./ vecnorm(nor);
                        FcType(isBd) = options.FcType((BdFcNd1 + BdFcNd2 + BdFcNd3) / 3, nor);
                    end
                    % Edge.type is not used in 3D (types are labeled on faces).
                    elem = Elem(ElNode, ElEdge, [], ElFace, ElFcPerm);
                    edge = Edge(EgNode, [], zeros(1, nEdge));
                    face = Face(FcNode, FcElem, FcEdge, FcType);
                    msh = Msh(type, node, elem, edge, face);
            end
        end
    end
end
