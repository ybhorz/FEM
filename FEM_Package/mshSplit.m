function msh = mshSplit(msh)
    % mshSplit: refine 2D mesh using Alfeld splitting.
    % Newly generated nodes and edges are labeled as "1i" type.
    arguments (Input)
        msh Msh;
    end
    arguments (Output)
        msh Msh;
    end
    assert(msh.type == "D2T");

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
