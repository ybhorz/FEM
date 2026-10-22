function fE = stdFE(name)
    % stdFE: standard 3D finite elements used in tests, cached within one run of `runTests` (cleared there).
    % Construction of 3D elements (e.g. P2, P3) takes seconds; the same elements are used by several test classes.

    % Name   | Element
    % ------ | --------------------------------------------------------------
    % P0     | constant, DoF at barycenter of element
    % P1     | linear, DoFs at vertices
    % P2     | quadratic, DoFs at vertices and midpoints of edges
    % P3     | cubic, DoFs at vertices, two points on each edge, barycenters of faces
    % DG1    | linear, DoFs at vertices, not shared
    % CR     | linear, DoFs at barycenters of faces
    % RT0    | lowest order Raviart-Thomas, DoFs are normal fluxes through faces
    % BDM1   | first order Brezzi-Douglas-Marini, DoFs are moments of normal component against barycentric coordinates of faces
    % NED1   | lowest order Nedelec (first kind), DoFs are tangential moments on edges
    % TP1    | linear trace on faces (reference face "D3FR"), DoFs at face vertices, not shared
    % SDG0S  | scalar SDG_0 on Alfeld split (face 4 primal): constant, DoF at barycenter of face 4
    % SDG0V  | vector SDG_0: constant, DoFs are normal fluxes through faces 1, 2, 3 (dual faces)
    % SDG1S  | scalar SDG_1: linear, DoFs at vertices of face 4 (shared) and at vertex 4 (not shared)
    % SDG1V  | vector SDG_1: linear, BDM1 moments on faces 1, 2, 3 (shared) and on face 4 (not shared)
    % StSDG1M| matrix SDG_1 of Stokes: linear, moments of each component of sigma * n against barycentric coordinates
    %        | on faces 1, 2, 3 (shared) and on face 4 (not shared), rows mapped by contravariant Piola transformation
    % StSDG1V| vector SDG_1 of Stokes: linear, each component at vertices of face 4 (shared) and at vertex 4 (not shared)
    % StSDG1S| scalar SDG_1 of Stokes: linear, value at vertex 4 and at midpoints of edges [1, 4], [2, 4], [3, 4]
    % DG1M   | matrix linear, each component at vertices, not shared (velocity gradient of hybridized Brinkman SDG)
    % TP1V2  | linear trace with 2 components on faces (reference face "D3FR"), each component at face vertices
    % TP1V3  | linear trace with 3 components on faces, each component at face vertices

    arguments (Input)
        name (1, 1) string {mustBeMember(name, ["P0", "P1", "P2", "P3", "DG1", "CR", "RT0", "BDM1", "NED1", "TP1", ...
            "SDG0S", "SDG0V", "SDG1S", "SDG1V", "StSDG1M", "StSDG1V", "StSDG1S", "DG1M", "TP1V2", "TP1V3"])};
    end
    arguments (Output)
        fE FE;
    end
    persistent cache
    if isempty(cache)
        cache = containers.Map();
    end
    if isKey(cache, char(name))
        fE = cache(char(name));
        return;
    end
    D3TElem = MshEnt("D3T").msh;
    d0 = [0; 0; 0];
    P2FS = "[1,x,y,z,x^2,y^2,z^2,x*y,y*z,z*x]";
    P3FS = "[1,x,y,z,x^2,y^2,z^2,x*y,y*z,z*x,x^3,y^3,z^3,x^2*y,x^2*z,y^2*x,y^2*z,z^2*x,z^2*y,x*y*z]";
    UNV = MshEnt("D3F").UNV;
    FcBar = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")];
    switch name
        case "P0"
            fE = FE("D3T", "1", NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], d0), "map", "affine");
        case "P1"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 0, [], d0), "map", "affine");
        case "P2"
            fE = FE("D3T", P2FS, [NdDoF("D3", D3TElem, 0, [], d0), NdDoF("D3", D3TElem, 1, 1/2, d0)], "map", "affine");
        case "P3"
            fE = FE("D3T", P3FS, [NdDoF("D3", D3TElem, 0, [], d0), NdDoF("D3", D3TElem, 1, [1/3, 2/3], d0), ...
                NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0)], "map", "affine");
        case "DG1"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 0, [], d0, "share", false), "map", "affine");
        case "CR"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0), "map", "affine");
        case "RT0"
            fE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", MoDoF("D3", D3TElem, 2, Fcn.cst(1), zeros(3), ...
                "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 1)), "map", "piolaDiv");
        case "BDM1"
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), MoDoF("D3", D3TElem, 2, ...
                [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")], zeros(3), ...
                "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 2)), "map", "piolaDiv");
        case "NED1"
            fE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; 0,-z,y; z,0,-x; -y,x,0].'", MoDoF("D3", D3TElem, 1, Fcn.cst(1), zeros(3), ...
                "coef", MshEnt("D3L").UTV, "orien", true, "GInt", GInt("D3L", 1)), "map", "piolaCurl");
        case "SDG0S"
            fE = FE("D3T", "1", NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0, "EntIdx", 4), "map", "affine");
        case "SDG0V"
            fE = FE("D3T", FE.repFS("1", [3, 1]), MoDoF("D3", D3TElem, 2, Fcn.cst(1), zeros(3), "coef", UNV, ...
                "EntIdx", [1, 2, 3], "orien", true, "GInt", GInt("D3F", 1)), "map", "piolaDiv");
        case "SDG1S"
            fE = FE("D3T", "[1,x,y,z]", [NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], d0, "EntIdx", 4), ...
                NdDoF("D3", D3TElem, 0, [], d0, "EntIdx", 4, "share", false)], "map", "affine");
        case "SDG1V"
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
                [MoDoF("D3", D3TElem, 2, FcBar, zeros(3), "coef", UNV, "EntIdx", [1, 2, 3], "orien", true, "GInt", GInt("D3F", 2)), ...
                MoDoF("D3", D3TElem, 2, FcBar, zeros(3), "coef", UNV, "EntIdx", 4, "orien", true, "share", false, ...
                "GInt", GInt("D3F", 2))], "map", "piolaDiv");
        case "StSDG1M"
            % Row i of sigma: components i, i + 3, i + 6 (column-major); other components are disabled (NaN order).
            DoFs = MoDoF.empty;
            for EntIdx = {[1, 2, 3], 4}
                for iRow = 1:3
                    ord = nan(3, 9);
                    ord(:, iRow:3:9) = 0;
                    DoFs(end + 1) = MoDoF("D3", D3TElem, 2, FcBar, ord, "coef", UNV, "EntIdx", EntIdx{1}, "orien", true, ...
                        "share", isequal(EntIdx{1}, [1, 2, 3]), "form", @(coef, fcn, tst) sum(fcn * coef(1)) .* tst, ...
                        "GInt", GInt("D3F", 2));
                end
            end
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 3]), DoFs, "map", "piolaDiv");
        case "StSDG1V"
            DoFs = NdDoF.empty;
            for iComp = 1:3
                ord = nan(3, 3);
                ord(:, iComp) = 0;
                DoFs(end + 1) = NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], ord, "EntIdx", 4);
            end
            for iComp = 1:3
                ord = nan(3, 3);
                ord(:, iComp) = 0;
                DoFs(end + 1) = NdDoF("D3", D3TElem, 0, [], ord, "EntIdx", 4, "share", false);
            end
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), DoFs, "map", "affine");
        case "StSDG1S"
            fE = FE("D3T", "[1,x,y,z]", [NdDoF("D3", D3TElem, 0, [], d0, "EntIdx", 4), ...
                NdDoF("D3", D3TElem, 1, 1/2, d0, "EntIdx", [3, 5, 6])], "map", "affine");
        case "DG1M"
            DoFs = NdDoF.empty;
            for iComp = 1:9
                ord = nan(3, 9);
                ord(:, iComp) = 0;
                DoFs(end + 1) = NdDoF("D3", D3TElem, 0, [], ord, "share", false);
            end
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 3]), DoFs, "map", "affine");
        case {"TP1V2", "TP1V3"}
            % Component selected by constant coefficient (derivative order of trace DoF must be zero).
            nComp = str2double(extractAfter(name, "TP1V"));
            DoFs = NdDoF.empty;
            for iComp = 1:nComp
                DoFs(end + 1) = NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], zeros(2, nComp), ...
                    "coef", Fcn.cst(double(1:nComp == iComp)'), "share", false);
            end
            fE = FE("D3FR", FE.repFS("[1,s,t]", [nComp, 1]), DoFs);
        case "TP1"
            fE = FE("D3FR", "[1,s,t]", NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], [0; 0], "share", false));
    end
    cache(char(name)) = fE;
end
