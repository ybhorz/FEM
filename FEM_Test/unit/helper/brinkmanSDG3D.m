function sys = brinkmanSDG3D(nSub, v)
    % brinkmanSDG3D: system of the hybridized 3D Brinkman SDG scheme (same as FEM_Example/Brinkman_SDG_3D, alpha = 10)
    % on mshSplit(mshD3TS(unit cube, nSub)) with viscosity v. Unknowns (FE spaces `trls`): velocity gradient s
    % (discontinuous), velocity u, tangential multiplier on dual faces, multiplier on interior primal faces, pressure p.
    % sys: msh, trls, Stiff, Load, and exact solutions s, u, p.
    a = 10;
    d0_p = [0; 0; 0]; grad_p = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
    phi = Fcn("D3", "sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2").dif(grad_p).fun;
    u = Fcn("D3", [phi(2) - phi(3); phi(3) - phi(1); phi(1) - phi(2)]);
    p = Fcn("D3", "10*(x-1/2)*(y-1/2)*(z-1/2)");
    d0_u = zeros(3); div_u = eye(3);
    grad_u = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
    s = dif(u, grad_u);
    d0_s = zeros(3, 9); div_s = repelem(eye(3), 1, 3);
    f = - sum(dif(s, div_s), 2) * v + u * a + dif(p, grad_p);
    d0_g = zeros(2, 2); d0_l = zeros(2, 3);
    UNV = MshEnt("D3F").UNV;
    FcParm = MshEnt("D3F").parm;
    Tmat = Fcn("D3F", [FcParm(:, 2) - FcParm(:, 1), FcParm(:, 3) - FcParm(:, 1)]);
    msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], nSub));
    PrFace = find(ismember(msh.face.type, 0:6));
    PrOFace = find(ismember(msh.face.type, 0));
    DlFace = find(ismember(msh.face.type, 1i));
    BdFace = msh.bdEnt(2, 1:6);
    trls = [FES(msh, stdFE("DG1M")), FES(msh, stdFE("SDG1V")), ...
        FES(msh, stdFE("TP1V2"), BC(Fcn("D3", "0"), "face", PrFace)), ...
        FES(msh, stdFE("TP1V3"), BC(Fcn("D3", "0"), "face", [DlFace, BdFace])), ...
        FES(msh, stdFE("SDG1S"), BC(p, "DoF", 1))];
    is = 1; iu = 2; ig = 3; il = 4; ip = 5; it = 1; iv = 2; ih = 3; im = 4; iq = 5;
    d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; d0_h = d0_g; d0_m = d0_l; d0_q = d0_p; grad_q = grad_p;
    Auv = DLF(msh, 3, Fcn.cst(a), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D3T", 2));
    Ast = DLF(msh, 3, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D3T", 2));
    But = [DLF(msh, 3, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D3T", 1)), ...
        DLF.interface(msh, 2, [-UNV, UNV, UNV], d0_u, d0_t, "EntIdx", DlFace, "iTrl", iu, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl, coef(1)) .* dot(tst * coef(2), coef(3)), "GInt", GInt("D3F", 2))];
    Bgt = DLF.interface(msh, 2, [-UNV, Tmat], d0_g, d0_t, "EntIdx", DlFace, "iTrl", ig, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(coef(2) * trl, tst * coef(1)), "GInt", GInt("D3F", 2));
    Blt = DLF.interface(msh, 2, -UNV, d0_l, d0_t, "EntIdx", PrOFace, "iTrl", il, "iTst", it, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D3F", 2));
    Bsv = [DLF(msh, 3, Fcn.cst(-v), div_s, d0_v, "iTrl", is, "iTst", iv, "form", @(coef, trl, tst) coef .* dot(sum(trl, 2), tst), "GInt", GInt("D3T", 1)), ...
        DLF.interface(msh, 2, [UNV * v, UNV, UNV], d0_s, d0_v, "EntIdx", DlFace, "iTrl", is, "iTst", iv, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(tst, coef(1)) .* dot(trl * coef(2), coef(3)), "GInt", GInt("D3F", 2))];
    Bsh = DLF.interface(msh, 2, [UNV * v, Tmat], d0_s, d0_h, "EntIdx", DlFace, "iTrl", is, "iTst", ih, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl * coef(1), coef(2) * tst), "GInt", GInt("D3F", 2));
    Bsm = DLF.interface(msh, 2, UNV * v, d0_s, d0_m, "EntIdx", PrOFace, "iTrl", is, "iTst", im, "trlOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D3F", 2));
    Bpv = [DLF(msh, 3, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", 1)), ...
        DLF.interface(msh, 2, UNV, d0_p, d0_v, "EntIdx", PrFace, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D3F", 2))];
    Buq = [DLF(msh, 3, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", 1)), ...
        DLF.interface(msh, 2, UNV, d0_u, d0_q, "EntIdx", DlFace, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D3F", 2))];
    Fv = SLF(msh, 3, f, d0_v, "iTst", iv, "GInt", GInt("D3T", 3));
    [Stiff, Load] = assemble(msh, trls, trls, [Auv, Ast, But, Bgt, Blt, Bsv, Bsh, Bsm, Bpv, Buq], Fv);
    sys = struct("msh", msh, "trls", trls, "Stiff", Stiff, "Load", Load, "s", s, "u", u, "p", p);
end
