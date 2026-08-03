function msh = mshD2TS(domn, nSub)
    % mshD2TS: generate a structured 2D triangular mesh over a rectangular domain.

    % Mesh entity type:
    % Type | Location
    % ---- | -------------------
    %   0  | interior
    %  +1  | lower boundary
    %  +2  | right boundary
    %  +3  | upper boundary
    %  +4  | left boundary
    %  -1  | lower left corner
    %  -2  | lower right corner
    %  -3  | upper right corner
    %  -4  | upper left corner

    arguments
        domn (1, 4); % Rectangular domain: [xmin, xmax, ymin, ymax].
        nSub (1, :); % Number of subdivisions in each direction.
        % If nSub is scalar, nx = ny = nSub, else nx = nSub(1), ny = nSub(2).
    end
    a1 = domn(1); a2 = domn(2); b1 = domn(3); b2 = domn(4);
    if isscalar(nSub)
        nx = nSub;
        ny = nSub;
    else
        nx = nSub(1);
        ny = nSub(2);
    end
    
    % Node.
    nNode = (nx + 1) * (ny + 1);
    xI = linspace(a1, a2, nx + 1);
    yJ = linspace(b1, b2, ny + 1);
    [XI, YJ] = meshgrid(xI, yJ);
    XI = XI'; YJ = YJ';
    NdIdx = reshape(1:nNode, [nx + 1, ny + 1]);

    NdCrd = [XI(:)'; YJ(:)'];
    NdType = zeros(1, nNode);
    NdType(NdIdx(2:nx, 1)) = 1;
    NdType(NdIdx(nx + 1, 2:ny)) = 2;
    NdType(NdIdx(2:nx, ny + 1)) = 3;
    NdType(NdIdx(1, 2:ny)) = 4;
    NdType(NdIdx(1, 1)) = -1;
    NdType(NdIdx(nx + 1, 1)) = -2;
    NdType(NdIdx(nx + 1, ny + 1)) = -3;
    NdType(NdIdx(1, ny + 1)) = -4;

    % Element.
    nElem = 2 * nx * ny;
    El1Idx = reshape(1:2:nElem, [nx, ny]);
    El2Idx = reshape(2:2:nElem, [nx, ny]);

    ElNode = zeros(3, nElem);
    ElNode(:, 1:2:nElem) = [reshape(NdIdx(1:nx, 1:ny), [1, nx * ny]); ...
                            reshape(NdIdx(2:nx + 1, 1:ny), [1, nx * ny]); ...
                            reshape(NdIdx(1:nx, 2:ny + 1), [1, nx * ny])];
    ElNode(:, 2:2:nElem) = [reshape(NdIdx(2:nx + 1, 1:ny), [1, nx * ny]); ...
                            reshape(NdIdx(2:nx + 1, 2:ny + 1), [1, nx * ny]); ...
                            reshape(NdIdx(1:nx, 2:ny + 1), [1, nx * ny])];
    ElType = zeros(1, nElem);
    ElType(El1Idx(2:nx - 1, 1)) = 1;
    ElType(El2Idx(2:nx - 1, 1)) = 1;
    ElType(El1Idx(nx, 2:ny - 1)) = 2;
    ElType(El2Idx(nx, 2:ny - 1)) = 2;
    ElType(El1Idx(2:nx - 1, ny)) = 3;
    ElType(El2Idx(2:nx - 1, ny)) = 3;
    ElType(El1Idx(1, 2:ny - 1)) = 4;
    ElType(El2Idx(1, 2:ny - 1)) = 4;
    ElType(El1Idx(1, 1)) = -1;
    ElType(El2Idx(1, 1)) = -1;
    ElType(El1Idx(nx, 1)) = -2;
    ElType(El2Idx(nx, 1)) = -2;
    ElType(El1Idx(nx, ny)) = -3;
    ElType(El2Idx(nx, ny)) = -3;
    ElType(El1Idx(1, ny)) = -4;
    ElType(El2Idx(1, ny)) = -4;

    % Edge.
    nEdge1 = nx * (ny + 1); nEdge2 = nx * ny; nEdge3 = (nx + 1) * ny;
    Eg1Idx = reshape(1:nEdge1, [nx, ny + 1]);
    Eg2Idx = nEdge1 + reshape(1:nEdge2, [nx, ny]);
    Eg3Idx = nEdge1 + nEdge2 + reshape(1:nEdge3, [nx + 1, ny]);
    nEdge = nEdge1 + nEdge2 + nEdge3;

    EgNode = zeros(2, nEdge);
    EgNode(:, 1:nEdge1) = [reshape(NdIdx(1:nx, 1:ny + 1), [1, nx * (ny + 1)]); ...
                            reshape(NdIdx(2:nx + 1, 1:ny + 1), [1, nx * (ny + 1)])];
    EgNode(:, nEdge1 + (1:nEdge2)) = [reshape(NdIdx(2:nx + 1, 1:ny), [1, nx * ny]); ...
                                        reshape(NdIdx(1:nx, 2:ny + 1), [1, nx * ny])];
    EgNode(:, nEdge1 + nEdge2 + (1:nEdge3)) = [reshape(NdIdx(1:nx + 1, 2:ny + 1), [1, (nx + 1) * ny]); ...
                                                reshape(NdIdx(1:nx + 1, 1:ny), [1, (nx + 1) * ny])];
    EgType = zeros(1, nEdge);
    EgType(Eg1Idx(1:nx, 1)) = 1;
    EgType(Eg3Idx(nx + 1, 1:ny)) = 2;
    EgType(Eg1Idx(1:nx, ny + 1)) = 3;
    EgType(Eg3Idx(1, 1:ny)) = 4;

    % Element - Edge.
    ElEdge = zeros(3, nElem);
    ElEdge(:, 1:2:nElem) = [reshape(Eg1Idx(1:nx, 1:ny), [1, nx * ny]); ...
                            reshape(Eg2Idx(1:nx, 1:ny), [1, nx * ny]); ...
                            reshape(Eg3Idx(1:nx, 1:ny), [1, nx * ny])];
    ElEdge(:, 2:2:nElem) = -[reshape(Eg3Idx(2:nx + 1, 1:ny), [1, nx * ny]); ...
                            reshape(Eg1Idx(1:nx, 2:ny + 1), [1, nx * ny]); ...
                            reshape(Eg2Idx(1:nx, 1:ny), [1, nx * ny])];

    % Edge - Element.
    EgElem = zeros(2, nEdge);
    EgElem(1, 1:nx * ny) = El1Idx(:)';
    EgElem(1, nEdge1 + (1:nEdge2)) = El1Idx(:)';
    EgElem(1, reshape(Eg3Idx(1:nx, 1:ny), [1, nx * ny])) = El1Idx(:)';
    EgElem(2, reshape(Eg3Idx(2:nx + 1, 1:ny), [1, nx * ny])) = -El2Idx(:)';
    EgElem(2, nx + 1:nEdge1) = -El2Idx(:)';
    EgElem(2, nEdge1 + (1:nEdge2)) = -El2Idx(:)';

    % Edge is counter-clockwise in boundary element.
    EgNeg = [Eg1Idx(1:nx, ny + 1)', Eg3Idx(nx + 1, 1:ny)];
    EgNode([1, 2], EgNeg) = EgNode([2, 1], EgNeg);
    EgElem(:, EgNeg) = -EgElem(:, EgNeg);

    ElNeg1 = El2Idx(1:nx, ny);
    ElEdge(2, ElNeg1) = -ElEdge(2, ElNeg1);
    ElNeg2 = El2Idx(nx, 1:ny);
    ElEdge(1, ElNeg2) = -ElEdge(1, ElNeg2);

    node = Node(NdCrd, NdType);
    elem = Elem(ElNode, ElEdge, ElType);
    edge = Edge(EgNode, EgElem, EgType);
    msh = Msh("D2T", node, elem, edge);
end
