function sol = condSolve(Stiff, Load, FESs, iCond)
    % condSolve: solve linear system by static condensation of local DoFs.
    % DoFs of FE spaces `FESs(iCond)` are local: elements connected by shared DoFs of these spaces form small groups
    % (one element for discontinuous spaces of HDG; one macro element of Alfeld splitting for velocity of SDG, whose
    % DoFs on dual facets are shared only inside the macro element), and the stiffness matrix restricted to these DoFs
    % is block diagonal (one block per group). They are eliminated group by group, the Schur complement system for the
    % other DoFs is solved, and the eliminated DoFs are recovered. The result equals `Stiff \ Load` up to rounding.
    % Elimination in several levels: `iCond` is a cell array, e.g. {[1], [2, 3]} eliminates spaces 1 by their groups,
    % then spaces 2 and 3 by their groups from the Schur complement system.
    % DoFs of a trace space (facets) belong to the elements connected to their facets. Boundary DoFs (`FES.BC`) are not
    % eliminated (their rows are identity rows).
    arguments (Input)
        Stiff (:, :); % Stiffness matrix (boundary conditions imposed).
        Load (:, 1); % Load vector.
        FESs (1, :) FES; % FE spaces of the system.
        iCond; % Indices of FE spaces whose DoFs are eliminated (numeric, or cell array of levels).
    end
    arguments (Output)
        sol (:, 1); % Solution.
    end
    msh = FESs.getMsh;
    nDoF = size(Stiff, 1);
    assert(nDoF == FESs.cumDoF && length(Load) == nDoF);
    if ~iscell(iCond)
        iCond = {iCond};
    end
    % Group of each DoF (0 if not eliminated) for each level.
    DoFGrps = cell(1, length(iCond));
    for iLev = 1:length(iCond)
        DoFGrps{iLev} = grpDoF(msh, FESs, iCond{iLev}, nDoF);
    end
    sol = condLevel(Stiff, Load, DoFGrps);
end
%% Local functions.
function DoFGrp = grpDoF(msh, FESs, iCond, nDoF)
    % grpDoF: groups of eliminated DoFs: connected components of elements sharing eliminated DoFs.
    DoFIdx = zeros(0, 1);
    ElIdx = zeros(0, 1);
    for iFES = iCond
        fES = FESs(iFES);
        L2G = abs(fES.Lc2Gl);
        isBd = false(fES.nGlDoF, 1);
        if ~isempty(fES.BC)
            isBd(fES.BC.DoFIdx) = true;
        end
        switch fES.elem
            case msh.ElDomn
                El = repmat(1:msh.nElem, size(L2G, 1), 1);
                Dof = L2G;
            case msh.FtRefDomn
                % Trace space: DoFs of a facet belong to its connected elements.
                CnElem = abs(msh.facet.elem);
                El = [repmat(CnElem(1, :), size(L2G, 1), 1), repmat(CnElem(2, :), size(L2G, 1), 1)];
                Dof = [L2G, L2G];
            otherwise
                error("Unsupported FE space.");
        end
        isVal = El > 0 & ~isBd(Dof);
        DoFIdx = [DoFIdx; FESs.cumDoF(iFES - 1) + Dof(isVal)];
        ElIdx = [ElIdx; El(isVal)];
    end
    Inc = sparse(DoFIdx, ElIdx, 1, nDoF, msh.nElem);
    ElGrp = conncomp(graph(double(Inc' * Inc ~= 0)));
    % Eliminated DoFs must split into several groups (shared DoFs of a continuous space form one global group).
    assert(numel(unique(ElGrp(ElIdx))) > 1);
    DoFGrp = zeros(nDoF, 1);
    DoFGrp(DoFIdx) = ElGrp(ElIdx);
end
function sol = condLevel(Stiff, Load, DoFGrps)
    % condLevel: eliminate DoFs grouped by DoFGrps{1}, then solve the Schur complement system (recursively for further
    % levels) and recover the eliminated DoFs.
    nDoF = size(Stiff, 1);
    DoFGrp = DoFGrps{1};
    [grp, I] = sort(DoFGrp(DoFGrp > 0));
    I0 = find(DoFGrp > 0);
    I = I0(I);
    % Renumber groups consecutively.
    [~, ~, grp] = unique(grp);
    nGrp = max(grp);
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
    rhs = Load(B) - ABI * (AIIInv * Load(I));
    if length(DoFGrps) > 1
        xB = condLevel(S, rhs, cellfun(@(g) g(B), DoFGrps(2:end), "UniformOutput", false));
    else
        xB = S \ rhs;
    end
    sol = zeros(nDoF, 1);
    sol(B) = xB;
    sol(I) = AIIInv * (Load(I) - AIB * xB);
end
