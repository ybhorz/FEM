function msh = mshD3TS(domn, nSub)
    % mshD3TS: generate a structured 3D tetrahedral mesh over a cuboid domain.
    % Each cube is split into 6 tetrahedra sharing its main diagonal (Kuhn subdivision), which is conforming.

    % Face type:
    % Type | Location
    % ---- | -------------
    %   0  | interior
    %  +1  | x = xmin
    %  +2  | x = xmax
    %  +3  | y = ymin
    %  +4  | y = ymax
    %  +5  | z = zmin
    %  +6  | z = zmax

    arguments
        domn (1, 6); % Cuboid domain: [xmin, xmax, ymin, ymax, zmin, zmax].
        nSub (1, :); % Number of subdivisions in each direction.
        % If nSub is scalar, nx = ny = nz = nSub, else nx = nSub(1), ny = nSub(2), nz = nSub(3).
    end
    if isscalar(nSub)
        nx = nSub; ny = nSub; nz = nSub;
    else
        nx = nSub(1); ny = nSub(2); nz = nSub(3);
    end

    % Node.
    sNd = [nx + 1, ny + 1, nz + 1];
    [XI, YJ, ZK] = ndgrid(linspace(domn(1), domn(2), nx + 1), linspace(domn(3), domn(4), ny + 1), linspace(domn(5), domn(6), nz + 1));
    NdCrd = [XI(:)'; YJ(:)'; ZK(:)'];

    % Element: vertices v0, v0 + e_p1, v0 + e_p1 + e_p2, v0 + e1 + e2 + e3 of cube with lower corner v0,
    % for each permutation p of [1, 2, 3].
    [I, J, K] = ndgrid(1:nx, 1:ny, 1:nz);
    nCube = nx * ny * nz;
    prm = perms(1:3);
    ElNode = zeros(4, 6 * nCube);
    for iPrm = 1:6
        IJK = [I(:)'; J(:)'; K(:)'];
        TetNode = zeros(4, nCube);
        TetNode(1, :) = sub2ind(sNd, IJK(1, :), IJK(2, :), IJK(3, :));
        for iStep = 1:3
            IJK(prm(iPrm, iStep), :) = IJK(prm(iPrm, iStep), :) + 1;
            TetNode(iStep + 1, :) = sub2ind(sNd, IJK(1, :), IJK(2, :), IJK(3, :));
        end
        ElNode(:, (iPrm - 1) * nCube + (1:nCube)) = TetNode;
    end

    msh = Msh.auto("D3T", Node(NdCrd), ElNode, "FcType", @(cen, nor) fcType(nor));
end
%% Local functions.
function type = fcType(nor)
    % fcType: type of boundary face from its outward unit normal.
    type = zeros(1, size(nor, 2));
    for iDim = 1:3
        type(nor(iDim, :) < -0.5) = 2 * iDim - 1;
        type(nor(iDim, :) > 0.5) = 2 * iDim;
    end
end
