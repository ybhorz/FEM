function msh = mshSplit(msh)
    % mshSplit: refine mesh using Alfeld splitting, i.e. each element is split into sub-elements by its barycenter.
    % Newly generated nodes, edges and faces are labeled as "1i" type; facets of the original mesh keep their types.

    % Mesh | Sub-element of element K (local numbering)                     | Facet of K (primal facet)
    % ---- | -------------------------------------------------------------- | --------------------------
    % D2T  | [v_i, v_i+1, c], i = 1, 2, 3                                   | edge 1 of sub-element
    % D3T  | [face nodes of i-th face of K (reversed), c], i = 1, 2, 3, 4   | face 4 of sub-element
    % c: barycenter of K. Sub-elements of K are numbered consecutively; other facets of sub-elements contain c (dual facets).

    arguments (Input)
        msh Msh;
    end
    arguments (Output)
        msh Msh;
    end
    switch msh.type
        case "D2T"
            msh = split2D(msh);
        case "D3T"
            msh = split3D(msh);
        otherwise
            error("Unsupported mesh type.");
    end
end
%% Local functions.
function msh = split2D(msh)
    % split2D: Alfeld splitting of 2D triangular mesh.

    % Node.
    CenCrd = zeros(2, msh.nElem);
    for iElem = 1:msh.nElem
        CenCrd(:, iElem) = mean(msh.node.coord(:, msh.elem.node(:, iElem)), 2);
    end
    NdCrd = [msh.node.coord, CenCrd];
    NdType = [msh.node.type, repmat(1i, [1, msh.nElem])];

    % Element.
    nSplit = msh.elem.nNode;
    nElem = msh.nElem * nSplit;
    ElNode = zeros(3, nElem);
    ElType = zeros(1, nElem);
    for iElem = 1:msh.nElem
        cumIdx = (iElem - 1) * nSplit;
        ElNode(1, cumIdx + (1:nSplit)) = msh.elem.node(:, iElem)';
        ElNode(2, cumIdx + (1:nSplit)) = circshift(msh.elem.node(:, iElem)', -1);
        ElNode(3, cumIdx + (1:nSplit)) = msh.nNode + iElem;
        ElType(cumIdx + (1:nSplit)) = msh.elem.type(iElem);
    end

    % Edge.
    nEdge = msh.nEdge + msh.nElem * nSplit;
    EgNode = zeros(2, nEdge);
    EgType = zeros(1, nEdge);
    EgNode(:, 1:msh.nEdge) = msh.edge.node;
    EgType(1:msh.nEdge) = msh.edge.type;
    EgType(msh.nEdge + 1:end) = 1i;
    for iElem = 1:msh.nElem
        cumIdx = msh.nEdge + (iElem - 1) * nSplit;
        EgNode(1, cumIdx + (1:nSplit)) = msh.nNode + iElem;
        EgNode(2, cumIdx + (1:nSplit)) = msh.elem.node(:, iElem)';
    end

    % Element - edge.
    ElEdge = zeros(3, nElem);
    for iElem = 1:msh.nElem
        cumElIdx = (iElem - 1) * nSplit;
        cumEgIdx = msh.nEdge + (iElem - 1) * nSplit;
        ElEdge(1, cumElIdx + (1:nSplit)) = msh.elem.edge(:, iElem)';
        ElEdge(2, cumElIdx + (1:nSplit)) = - (cumEgIdx + circshift(1:nSplit, -1));
        ElEdge(3, cumElIdx + (1:nSplit)) = cumEgIdx + (1:nSplit);
    end

    % Edge - element.
    EgElem = zeros(2, nEdge);
    for iElem = 1:msh.nElem
        cumElIdx = (iElem - 1) * nSplit;
        cumEgIdx = msh.nEdge + (iElem - 1) * nSplit;
        EgElem(1, cumEgIdx + (1:nSplit)) = cumElIdx + (1:nSplit);
        EgElem(2, cumEgIdx + (1:nSplit)) = - (cumElIdx + circshift(1:nSplit, 1));
        for iElEg = 1:msh.elem.nEdge
            EgIdx = abs(msh.elem.edge(iElEg, iElem));
            sgn = sign(msh.elem.edge(iElEg, iElem));
            EgElem(msh.edge.elem(:, EgIdx) == iElem * sgn, EgIdx) = (cumElIdx + iElEg) * sgn;
        end
    end

    node = Node(NdCrd, NdType);
    elem = Elem(ElNode, ElEdge, ElType);
    edge = Edge(EgNode, EgElem, EgType);
    msh = Msh("D2T", node, elem, edge);
end
function msh = split3D(msh)
    % split3D: Alfeld splitting of 3D tetrahedral mesh.
    nNode = msh.nNode;
    nElem = msh.nElem;
    % Node: barycenters of elements are appended.
    CenCrd = reshape(mean(reshape(msh.node.coord(:, msh.elem.node), 3, 4, nElem), 2), 3, nElem);
    NdCrd = [msh.node.coord, CenCrd];
    NdType = msh.node.type;
    if isempty(NdType)
        NdType = zeros(1, nNode);
    end
    NdType = [NdType, repmat(1i, [1, nElem])];
    % Element: i-th sub-element of K consists of the nodes of i-th face of K in reversed order (normal points to c)
    % and c, hence is positively oriented and its 4th face is the i-th face of K.
    LcFcNode = MshEnt("D3T").face.node;
    ElNode = zeros(4, 4, nElem);
    for iFace = 1:4
        ElNode(1:3, iFace, :) = reshape(msh.elem.node(LcFcNode([1, 3, 2], iFace), :), 3, 1, nElem);
        ElNode(4, iFace, :) = reshape(nNode + (1:nElem), 1, 1, nElem);
    end
    ElNode = reshape(ElNode, 4, 4 * nElem);
    sub = Msh.auto("D3T", Node(NdCrd, NdType), ElNode);
    % Orientation is not corrected by `Msh.auto`, i.e. barycenter remains the 4th node.
    assert(isequal(sub.elem.node, ElNode));
    % Face type: dual face (containing barycenter) is "1i", primal face keeps type of face of original mesh.
    isDual = any(sub.face.node > nNode, 1);
    [isPrimal, FcIdx] = ismember(sort(sub.face.node, 1)', sort(msh.face.node, 1)', "rows");
    assert(isequal(isPrimal', ~isDual));
    FcType = zeros(1, sub.nFace);
    FcType(isDual) = 1i;
    FcType(~isDual) = msh.face.type(FcIdx(~isDual));
    % Edge type: dual edge (containing barycenter) is "1i".
    EgType = zeros(1, sub.nEdge);
    EgType(any(sub.edge.node > nNode, 1)) = 1i;
    % Element type is inherited from original element.
    ElType = msh.elem.type;
    if ~isempty(ElType)
        ElType = repelem(ElType, 4);
    end
    elem = Elem(sub.elem.node, sub.elem.edge, ElType, sub.elem.face, sub.elem.facePerm);
    edge = Edge(sub.edge.node, sub.edge.elem, EgType);
    face = Face(sub.face.node, sub.face.elem, sub.face.edge, FcType);
    msh = Msh("D3T", sub.node, elem, edge, face);
end
