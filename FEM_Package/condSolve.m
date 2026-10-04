function sol = condSolve(Stiff, Load, FESs, iCond)
    % condSolve: solve linear system by static condensation of element-local DoFs.
    % DoFs of FE spaces `FESs(iCond)` (e.g. discontinuous spaces of HDG) belong to single elements, so that the
    % stiffness matrix restricted to them is block diagonal (one block per element). They are eliminated element by
    % element, the Schur complement system for the other DoFs (e.g. traces) is solved, and the eliminated DoFs are
    % recovered. The result equals `Stiff \ Load` up to rounding.
    arguments (Input)
        Stiff (:, :); % Stiffness matrix (boundary conditions imposed).
        Load (:, 1); % Load vector.
        FESs (1, :) FES; % FE spaces of the system.
        iCond (1, :); % Indices of FE spaces whose DoFs are eliminated.
    end
    arguments (Output)
        sol (:, 1); % Solution.
    end
    msh = FESs.getMsh;
    nDoF = size(Stiff, 1);
    assert(nDoF == FESs.cumDoF && length(Load) == nDoF);
    % Eliminated DoFs grouped by element: ElDoF(:, K) are the DoFs of element K.
    ElDoF = zeros(0, msh.nElem);
    for iFES = iCond
        fES = FESs(iFES);
        assert(ismember(fES.elem, msh.ElDomn));
        assert(isempty(fES.BC));
        L2G = abs(fES.Lc2Gl);
        % Element-local: each DoF belongs to exactly one element.
        assert(all(accumarray(L2G(:), 1, [fES.nGlDoF, 1]) == 1));
        ElDoF = [ElDoF; FESs.cumDoF(iFES - 1) + L2G];
    end
    [nBlk, nElem] = size(ElDoF);
    I = ElDoF(:);
    B = setdiff((1:nDoF)', I);
    % Inverse of block-diagonal matrix A_II (consecutive blocks of size nBlk).
    [r, c, v] = find(Stiff(I, I));
    iBlk = ceil(r / nBlk);
    assert(all(iBlk == ceil(c / nBlk)));
    Blk = zeros(nBlk, nBlk, nElem);
    Blk(sub2ind(size(Blk), r - (iBlk - 1) * nBlk, c - (iBlk - 1) * nBlk, iBlk)) = v;
    BlkInv = pageinv(Blk);
    [rLc, cLc, K] = ndgrid(1:nBlk, 1:nBlk, 1:nElem);
    AIIInv = sparse(rLc(:) + (K(:) - 1) * nBlk, cLc(:) + (K(:) - 1) * nBlk, BlkInv(:), length(I), length(I));
    % Schur complement system for the other DoFs.
    ABI = Stiff(B, I);
    AIB = Stiff(I, B);
    S = Stiff(B, B) - ABI * (AIIInv * AIB);
    xB = S \ (Load(B) - ABI * (AIIInv * Load(I)));
    sol = zeros(nDoF, 1);
    sol(B) = xB;
    sol(I) = AIIInv * (Load(I) - AIB * xB);
end
