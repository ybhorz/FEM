function sys = brinkmanSDG2D(nSub, useEntIdx)
    % brinkmanSDG2D: system of the hybridized 2D Brinkman SDG scheme (same as FEM_Example/Brinkman_SDG, nu = 1,
    % alpha = 10) on mshSplit(mshD2TS(unit square, nSub)). Unknowns (FE spaces `trls`): velocity gradient s
    % (discontinuous), velocity u, tangential multiplier gamma on dual edges, multiplier lambda on interior primal edges,
    % pressure p.
    % useEntIdx = true: multipliers have DoFs only on their edges (option "EntIdx" of FES);
    % useEntIdx = false: multipliers on all edges, set to zero on the other edges by boundary conditions.
    % sys: msh, trls, Stiff, Load, edge sets DlEdge, PrOEdge.
    v = 1; a = 10;
    u = Fcn("D2", "[sin(pi*x)^2*sin(2*pi*y); - sin(2*pi*x)*sin(pi*y)^2]");
    p = Fcn("D2", "10*(x-1/2)*(y-1/2)");
    d0_u = [0, 0; 0, 0]'; div_u = [1, 0; 0, 1]; grad_u = cat(3, [1, 0; 1, 0]', [0, 1; 0, 1]');
    s = dif(u, grad_u);
    d0_s = [0, 0; 0, 0; 0, 0; 0, 0]'; div_s = [1, 0; 1, 0; 0, 1; 0 ,1]';
    d0_p = [0, 0]'; grad_p = cat(3, [1, 0]', [0, 1]');
    f = - sum(dif(s, div_s), 2) * v + u * a + dif(p, grad_p);
    d0_g = 0; d0_l = [0, 0];
    UNV = MshEnt("D2L").UNV; UTV = MshEnt("D2L").UTV;
    msh = mshSplit(mshD2TS([0, 1, 0, 1], nSub));
    PrEdge = find(ismember(msh.edge.type, [0, 1, 2, 3, 4]));
    PrOEdge = find(ismember(msh.edge.type, 0));
    DlEdge = find(ismember(msh.edge.type, 1i));
    BdEdge = find(ismember(msh.edge.type, [1, 2, 3, 4]));
    D2TElem = MshEnt("D2T").msh;
    SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", D2TElem, 1, [0, 1], d0_p, "EntIdx", 1), ...
        NdDoF("D2", D2TElem, 0, [], d0_p, "EntIdx", 3, "share", false)]);
    SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
        [NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true), ...
        NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", 1, "orien", true, "share", false)]);
    DG_1M_DoFs = NdDoF.empty;
    for iComp = 1:4
        d0_si = nan(2, 4); d0_si(:, iComp) = 0;
        DG_1M_DoFs(end + 1) = NdDoF("D2", D2TElem, 0, [], d0_si, "share", false);
    end
    DG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), DG_1M_DoFs, "map", "affine");
    TP_1S_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_g, "share", false));
    TP_1V_FE = FE("D2LR", FE.repFS("[1,s]", [2, 1]), ...
        [NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_l, "coef", Fcn.cst([1; 0]), "share", false), ...
        NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_l, "coef", Fcn.cst([0; 1]), "share", false)]);
    if useEntIdx
        Gh = FES(msh, TP_1S_FE, BC.empty, "EntIdx", DlEdge);
        Lh = FES(msh, TP_1V_FE, BC.empty, "EntIdx", PrOEdge);
    else
        Gh = FES(msh, TP_1S_FE, BC(Fcn("D2", "0"), "edge", PrEdge));
        Lh = FES(msh, TP_1V_FE, BC(Fcn("D2", "0"), "edge", [DlEdge, BdEdge]));
    end
    trls = [FES(msh, DG_1M_FE), FES(msh, SDG_1V_FE), Gh, Lh, FES(msh, SDG_1S_FE, BC(p, "DoF", 1))];
    is = 1; iu = 2; ig = 3; il = 4; ip = 5; it = 1; iv = 2; ih = 3; im = 4; iq = 5;
    d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; d0_h = d0_g; d0_m = d0_l; d0_q = d0_p; grad_q = grad_p;
    Auv = DLF(msh, 2, Fcn.cst(a), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", 2));
    Ast = DLF(msh, 2, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D2T", 2));
    But = [DLF(msh, 2, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, [-UNV, UNV, UNV], d0_u, d0_t, "EntIdx", DlEdge, "iTrl", iu, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl, coef(1)) .* dot(tst * coef(2), coef(3)), "GInt", GInt("D2L", 2))];
    Bgt = DLF.interface(msh, 1, [-UNV, UTV], d0_g, d0_t, "EntIdx", DlEdge, "iTrl", ig, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) trl .* dot(tst * coef(1), coef(2)), "GInt", GInt("D2L", 2));
    Blt = DLF.interface(msh, 1, -UNV, d0_l, d0_t, "EntIdx", PrOEdge, "iTrl", il, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D2L", 2));
    Bsv = [DLF(msh, 2, Fcn.cst(-v), div_s, d0_v, "iTrl", is, "iTst", iv, "form", @(coef, trl, tst) coef .* dot(sum(trl, 2), tst), "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, [UNV * v, UNV, UNV], d0_s, d0_v, "EntIdx", DlEdge, "iTrl", is, "iTst", iv, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(tst, coef(1)) .* dot(trl * coef(2), coef(3)), "GInt", GInt("D2L", 2))];
    Bsh = DLF.interface(msh, 1, [UNV * v, UTV], d0_s, d0_h, "EntIdx", DlEdge, "iTrl", is, "iTst", ih, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl * coef(1), coef(2)) .* tst, "GInt", GInt("D2L", 2));
    Bsm = DLF.interface(msh, 1, UNV * v, d0_s, d0_m, "EntIdx", PrOEdge, "iTrl", is, "iTst", im, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D2L", 2));
    Bpv = [DLF(msh, 2, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, UNV, d0_p, d0_v, "EntIdx", PrEdge, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L", 2))];
    Buq = [DLF(msh, 2, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, UNV, d0_u, d0_q, "EntIdx", DlEdge, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", 2))];
    Fv = SLF(msh, 2, f, d0_v, "iTst", iv, "GInt", GInt("D2T", 3));
    [Stiff, Load] = assemble(msh, trls, trls, [Auv, Ast, But, Bgt, Blt, Bsv, Bsh, Bsm, Bpv, Buq], Fv);
    sys = struct("msh", msh, "trls", trls, "Stiff", Stiff, "Load", Load, "DlEdge", DlEdge, "PrOEdge", PrOEdge);
end
