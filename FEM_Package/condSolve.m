function sol = condSolve(Stiff, Load, FESs, iCond)
    % condSolve: solve linear system by static condensation of local DoFs.
    % DoFs of FE spaces `FESs(iCond)` are local: elements connected by shared DoFs of these spaces form small groups
    % (one element for discontinuous spaces of HDG; one macro element of Alfeld splitting for velocity of SDG, whose
    % DoFs on dual facets are shared only inside the macro element), and the stiffness matrix restricted to these DoFs
    % is block diagonal (one block per group). They are eliminated group by group, the Schur complement system for the
    % other DoFs is solved, and the eliminated DoFs are recovered. The result equals `Stiff \ Load` up to rounding.
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
    % Incidence of eliminated DoFs and elements.
    DoFIdx = zeros(0, 1);
    ElIdx = zeros(0, 1);
    for iFES = iCond
        fES = FESs(iFES);
        assert(ismember(fES.elem, msh.ElDomn));
        assert(isempty(fES.BC));
        L2G = abs(fES.Lc2Gl);
        DoFIdx = [DoFIdx; FESs.cumDoF(iFES - 1) + L2G(:)];
        ElIdx = [ElIdx; reshape(repmat(1:msh.nElem, size(L2G, 1), 1), [], 1)];
    end
    % Groups: connected components of elements sharing eliminated DoFs.
    Inc = sparse(DoFIdx, ElIdx, 1, nDoF, msh.nElem);
    ElGrp = conncomp(graph(double(Inc' * Inc ~= 0)));
    nGrp = max(ElGrp);
    % Eliminated DoFs must split into several groups (shared DoFs of a continuous space form one global group).
    assert(nGrp > 1);
    DoFGrp = zeros(nDoF, 1);
    DoFGrp(DoFIdx) = ElGrp(ElIdx);
    % Eliminated DoFs sorted by group; blocks of A_II.
    [grp, I] = sort(DoFGrp(DoFGrp > 0));
    I0 = find(DoFGrp > 0);
    I = I0(I);
    nI = length(I);
    AII = Stiff(I, I);
    [r, c] = find(AII);
    assert(all(grp(r) == grp(c)));
    % Inverse of A_II: blocks of the same size are inverted together.
    sGrp = accumarray(grp, 1, [nGrp, 1]);
    GrpStart = cumsum([0; sGrp(1:end - 1)]);
    rInv = zeros(0, 1); cInv = zeros(0, 1); vInv = zeros(0, 1);
    for sBlk = unique(sGrp)'
        iGrp = find(sGrp == sBlk)';
        nBlk = length(iGrp);
        % Pos(:, k): positions (in I) of DoFs of k-th group of this size.
        Pos = GrpStart(iGrp)' + (1:sBlk)';
        rPos = reshape(repmat(reshape(Pos, sBlk, 1, nBlk), 1, sBlk, 1), [], 1);
        cPos = reshape(repmat(reshape(Pos, 1, sBlk, nBlk), sBlk, 1, 1), [], 1);
        Blk = reshape(full(AII(sub2ind([nI, nI], rPos, cPos))), sBlk, sBlk, nBlk);
        BlkInv = pageinv(Blk);
        rInv = [rInv; rPos];
        cInv = [cInv; cPos];
        vInv = [vInv; BlkInv(:)];
    end
    AIIInv = sparse(rInv, cInv, vInv, nI, nI);
    % Schur complement system for the other DoFs.
    B = setdiff((1:nDoF)', I);
    ABI = Stiff(B, I);
    AIB = Stiff(I, B);
    S = Stiff(B, B) - ABI * (AIIInv * AIB);
    xB = S \ (Load(B) - ABI * (AIIInv * Load(I)));
    sol = zeros(nDoF, 1);
    sol(B) = xB;
    sol(I) = AIIInv * (Load(I) - AIB * xB);
end
