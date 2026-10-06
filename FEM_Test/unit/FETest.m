classdef FETest < matlab.unittest.TestCase
    % FETest: unit tests of FE.
    % Base functions are checked to be dual to degrees of freedom: DoF_i(base_j) = delta_ij.
    properties
        P2T = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1]; % Vertices of a counter-clockwise triangle (by column).
        P3T = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3]; % Vertices of a positively oriented tetrahedron.
    end
    methods (Test, TestTags = {'Fast'})
        function P0(tc)
            fE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0]));
            tc.verifyEqual(fE.nDoF, 1);
            tc.verifyEqual(fE.base.domn, "D2T");
            tc.verifyTrue(symEq(fE.base.fun, sym(1)));
        end
        function P1(tc)
            fE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
            tc.verifyEqual(fE.nDoF, 3);
            pnts = ["[x1; y1]", "[x2; y2]", "[x3; y3]"];
            tc.verifyTrue(symEq(evalPnt(fE.base, pnts), sym(eye(3))));
            tc.verifyTrue(symEq(sum([fE.base.fun]), sym(1)));
        end
        function P2(tc)
            fE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
                NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]);
            tc.verifyEqual(fE.nDoF, 6);
            pnts = ["[x1; y1]", "[x2; y2]", "[x3; y3]", ...
                "[(x1 + x2)/2; (y1 + y2)/2]", "[(x2 + x3)/2; (y2 + y3)/2]", "[(x3 + x1)/2; (y3 + y1)/2]"];
            tc.verifyTrue(symEq(evalPnt(fE.base, pnts), sym(eye(6))));
        end
        function CR(tc)
            fE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0]));
            tc.verifyEqual(fE.nDoF, 3);
            pnts = ["[(x1 + x2)/2; (y1 + y2)/2]", "[(x2 + x3)/2; (y2 + y3)/2]", "[(x3 + x1)/2; (y3 + y1)/2]"];
            tc.verifyTrue(symEq(evalPnt(fE.base, pnts), sym(eye(3))));
        end
        function VP1(tc)
            fE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
                NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]);
            tc.verifyEqual(fE.nDoF, 6);
            pnts = ["[x1; y1]", "[x2; y2]", "[x3; y3]"];
            val = sym(zeros(6));
            for iBase = 1:6
                for iPnt = 1:3
                    vec = fE.base(iBase).eval(pnts(iPnt));
                    val(iPnt, iBase) = vec(1);
                    val(3 + iPnt, iBase) = vec(2);
                end
            end
            tc.verifyTrue(symEq(val, sym(eye(6))));
        end
        function RT0(tc)
            UNV = MshEnt("D2T").UNV;
            fE = FE("D2T", "[1,0; 0,1; x,y].'", ...
                MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4)));
            tc.verifyEqual(fE.nDoF, 3);
            val = sym(zeros(3));
            for iBase = 1:3
                for iEdge = 1:3
                    val(iEdge, iBase) = int(dot(fE.base(iBase), UNV(iEdge)), "D2L", iEdge).fun;
                end
            end
            tc.verifyEqual(numSym(val, "D2T", tc.P2T), eye(3), "AbsTol", 1e-10);
        end
        function BDM1(tc)
            UNV = MshEnt("D2T").UNV;
            fE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4)), ...
                MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4))]);
            tc.verifyEqual(fE.nDoF, 6);
            % Parameter s of i-th edge (from i-th vertex to next vertex) minus 1/2.
            tsts = [Fcn("D2T", "((x - x1)*(x2 - x1) + (y - y1)*(y2 - y1))/((x2 - x1)^2 + (y2 - y1)^2) - 1/2"), ...
                Fcn("D2T", "((x - x2)*(x3 - x2) + (y - y2)*(y3 - y2))/((x3 - x2)^2 + (y3 - y2)^2) - 1/2"), ...
                Fcn("D2T", "((x - x3)*(x1 - x3) + (y - y3)*(y1 - y3))/((x1 - x3)^2 + (y1 - y3)^2) - 1/2")];
            val = sym(zeros(6));
            for iBase = 1:6
                for iEdge = 1:3
                    val(iEdge, iBase) = int(dot(fE.base(iBase), UNV(iEdge)), "D2L", iEdge).fun;
                    val(3 + iEdge, iBase) = int(dot(fE.base(iBase), UNV(iEdge)) .* tsts(iEdge), "D2L", iEdge).fun;
                end
            end
            tc.verifyEqual(numSym(val, "D2T", tc.P2T), eye(6), "AbsTol", 1e-10);
        end
        function TP1(tc)
            fE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
            tc.verifyEqual(fE.nDoF, 2);
            tc.verifyEqual(fE.base(1).domn, "D2LR");
            val = sym(zeros(2));
            for iBase = 1:2
                val(:, iBase) = [fE.base(iBase).eval(0); fE.base(iBase).eval(1)];
            end
            tc.verifyTrue(symEq(val, sym(eye(2))));
        end
        function P0_3D(tc)
            fE = FE("D3T", "1", NdDoF("D3", MshEnt("D3T").msh, 3, [1/4; 1/4; 1/4], [0; 0; 0]));
            tc.verifyEqual(fE.nDoF, 1);
            tc.verifyEqual(fE.base.domn, "D3T");
            tc.verifyTrue(symEq(fE.base.fun, sym(1)));
            % Moment DoF: base function is 1 / volume.
            fE = FE("D3T", "1", MoDoF("D3", MshEnt("D3T").msh, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)));
            P = tc.P3T;
            tc.verifyEqual(numSym(fE.base.fun, "D3T", P) * det(P(:, 2:4) - P(:, 1)) / 6, 1, "AbsTol", 1e-12);
        end
        function P1_3D(tc)
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 0, [], [0; 0; 0]));
            tc.verifyEqual(fE.nDoF, 4);
            tc.verifyEqual(fE.base(1).domn, "D3T");
            pnts = ["[x1; y1; z1]", "[x2; y2; z2]", "[x3; y3; z3]", "[x4; y4; z4]"];
            tc.verifyTrue(symEq(evalPnt(fE.base, pnts), sym(eye(4))));
            tc.verifyTrue(symEq(sum([fE.base.fun]), sym(1)));
            % Linear reproduction: sum_i x_i * base_i = x.
            tc.verifyTrue(symEq([fE.base.fun] * str2sym("[x1; x2; x3; x4]"), str2sym("x")));
        end
        function CR_3D(tc)
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 2, [1/3; 1/3], [0; 0; 0]));
            tc.verifyEqual(fE.nDoF, 4);
            % Barycenter of face i (opposite to vertex i).
            pnts = ["[x2 + x3 + x4; y2 + y3 + y4; z2 + z3 + z4]/3", "[x1 + x3 + x4; y1 + y3 + y4; z1 + z3 + z4]/3", ...
                "[x1 + x2 + x4; y1 + y2 + y4; z1 + z2 + z4]/3", "[x1 + x2 + x3; y1 + y2 + y3; z1 + z2 + z3]/3"];
            tc.verifyTrue(symEq(evalPnt(fE.base, pnts), sym(eye(4))));
            tc.verifyTrue(symEq(sum([fE.base.fun]), sym(1)));
        end
        function VP1_3D(tc)
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
                [NdDoF("D3", MshEnt("D3T").msh, 0, [], [0, nan, nan; 0, nan, nan; 0, nan, nan]), ...
                NdDoF("D3", MshEnt("D3T").msh, 0, [], [nan, 0, nan; nan, 0, nan; nan, 0, nan]), ...
                NdDoF("D3", MshEnt("D3T").msh, 0, [], [nan, nan, 0; nan, nan, 0; nan, nan, 0])]);
            tc.verifyEqual(fE.nDoF, 12);
            pnts = ["[x1; y1; z1]", "[x2; y2; z2]", "[x3; y3; z3]", "[x4; y4; z4]"];
            val = sym(zeros(12));
            for iBase = 1:12
                for iPnt = 1:4
                    vec = fE.base(iBase).eval(pnts(iPnt));
                    val(iPnt + (0:4:8), iBase) = vec;
                end
            end
            tc.verifyTrue(symEq(val, sym(eye(12))));
        end
        function P2_3D(tc)
            % With map = "none", symbolic inversion takes about 20 minutes (R2025a).
            fE = stdFE("P2");
            tc.verifyEqual(fE.nDoF, 10);
            tc.verifyEqual(fE.base(1).domn, "D3T");
            % Vertices and edge midpoints of a fixed tetrahedron.
            P = tc.P3T;
            EgNode = MshEnt("D3T").edge.node;
            tc.verifyEqual(evalNum(fE.base, [P, (P(:, EgNode(1, :)) + P(:, EgNode(2, :))) / 2], P), eye(10), "AbsTol", 1e-10);
        end
        function P3_3D(tc)
            fE = stdFE("P3");
            tc.verifyEqual(fE.nDoF, 20);
            % DoF order: vertices, 1st points on edges, 2nd points on edges, face barycenters.
            P = tc.P3T;
            EgNode = MshEnt("D3T").edge.node;
            FcNode = MshEnt("D3T").face.node;
            pnts = [P, (2 * P(:, EgNode(1, :)) + P(:, EgNode(2, :))) / 3, (P(:, EgNode(1, :)) + 2 * P(:, EgNode(2, :))) / 3, ...
                (P(:, FcNode(1, :)) + P(:, FcNode(2, :)) + P(:, FcNode(3, :))) / 3];
            tc.verifyEqual(evalNum(fE.base, pnts, P), eye(20), "AbsTol", 1e-10);
        end
        function affineMap(tc)
            % Base functions with map = "affine" equal those with map = "none".
            D2TElem = MshEnt("D2T").msh;
            D3TElem = MshEnt("D3T").msh;
            args = {{"D2T", "[1,x,y]", NdDoF("D2", D2TElem, 0, [], [0; 0])}, ...
                {"D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", D2TElem, 0, [], [0; 0]), NdDoF("D2", D2TElem, 1, 1/2, [0; 0])]}, ...
                {"D2T", "[1,x,y]", NdDoF("D2", D2TElem, 1, 1/2, [0; 0])}, ...
                {"D3T", "1", NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], [0; 0; 0])}, ...
                {"D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 0, [], [0; 0; 0])}, ...
                {"D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 2, [1/3; 1/3], [0; 0; 0])}};
            Ps = {tc.P2T, tc.P2T, tc.P2T, tc.P3T, tc.P3T, tc.P3T};
            rng(1);
            for iArg = 1:length(args)
                fE1 = FE(args{iArg}{:});
                fE2 = FE(args{iArg}{:}, "map", "affine");
                tc.verifyEqual(fE2.map, "affine");
                tc.verifyEqual([fE2.base.domn], [fE1.base.domn]);
                P = Ps{iArg};
                pnts = P * rand(size(P, 2), 5);
                tc.verifyEqual(evalNum(fE2.base, pnts, P), evalNum(fE1.base, pnts, P), "AbsTol", 1e-12);
            end
            % Vector-valued base functions.
            fE1 = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", D2TElem, 0, [], [0, nan; 0, nan]), NdDoF("D2", D2TElem, 0, [], [nan, 0; nan, 0])]);
            fE2 = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", D2TElem, 0, [], [0, nan; 0, nan]), NdDoF("D2", D2TElem, 0, [], [nan, 0; nan, 0])], "map", "affine");
            pnt = tc.P2T * [0.2; 0.3; 0.5];
            for iBase = 1:6
                fun1 = fE1.base(iBase).getFun;
                fun2 = fE2.base(iBase).getFun;
                tc.verifyEqual(fun2(pnt, tc.P2T), fun1(pnt, tc.P2T), "AbsTol", 1e-12);
            end
            % Invalid: moment DoF, derivative DoF, coefficient, trace element.
            tc.verifyError(@() FE("D3T", "1", MoDoF("D3", D3TElem, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)), ...
                "map", "affine"), "MATLAB:assertion:failed");
            tc.verifyError(@() FE("D2T", "[1,x,y]", [NdDoF("D2", D2TElem, 0, [], [0; 0], "EntIdx", 1), ...
                NdDoF("D2", D2TElem, 0, [], [1; 0], "EntIdx", 1), NdDoF("D2", D2TElem, 0, [], [0; 1], "EntIdx", 1)], ...
                "map", "affine"), "MATLAB:assertion:failed");
            tc.verifyError(@() FE("D2T", "[1,x,y]", NdDoF("D2", D2TElem, 1, 1/2, [0; 0], "coef", MshEnt("D2L").len), ...
                "map", "affine"), "MATLAB:assertion:failed");
            tc.verifyError(@() FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0), "map", "affine"), ...
                "MATLAB:assertion:failed");
        end
        function piolaDiv(tc)
            % 2D: base functions with map = "piolaDiv" equal those with map = "none".
            D2TElem = MshEnt("D2T").msh;
            UNV2 = MshEnt("D2L").UNV;
            dofs = {MoDoF("D2", D2TElem, 1, Fcn.cst(1), zeros(2), "coef", UNV2, "orien", true, "GInt", GInt("D2L", 4)), ...
                [MoDoF("D2", D2TElem, 1, Fcn.cst(1), zeros(2), "coef", UNV2, "orien", true, "GInt", GInt("D2L", 4)), ...
                MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s-1/2"), zeros(2), "coef", UNV2, "orien", true, "GInt", GInt("D2L", 4))]};
            FSs = {"[1,0; 0,1; x,y].'", FE.repFS("[1,x,y]", [2, 1])};
            pnts = tc.P2T * [0.2, 0.5, 0.1; 0.3, 0.1, 0.6; 0.5, 0.4, 0.3];
            for iFE = 1:2
                fE1 = FE("D2T", FSs{iFE}, dofs{iFE});
                fE2 = FE("D2T", FSs{iFE}, dofs{iFE}, "map", "piolaDiv");
                tc.verifyEqual(fE2.map, "piolaDiv");
                for iBase = 1:fE1.nDoF
                    fun1 = fE1.base(iBase).getFun;
                    fun2 = fE2.base(iBase).getFun;
                    for iPnt = 1:size(pnts, 2)
                        tc.verifyEqual(fun2(pnts(:, iPnt), tc.P2T), fun1(pnts(:, iPnt), tc.P2T), "AbsTol", 1e-12);
                    end
                end
            end
            % 3D RT0: normal flux of j-th base function through i-th face is delta_ij.
            D3TElem = MshEnt("D3T").msh;
            fE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", ...
                MoDoF("D3", D3TElem, 2, Fcn.cst(1), zeros(3), "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 1)), ...
                "map", "piolaDiv");
            tc.verifyEqual(fE.nDoF, 4);
            P = tc.P3T;
            flux = zeros(4);
            for iBase = 1:4
                fun = fE.base(iBase).getFun;
                for iFace = 1:4
                    FcNd = P(:, D3TElem.face.node(:, iFace));
                    % Normal component of RT0 function is constant on face: flux = v(barycenter) . (area * unit normal).
                    nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1)) / 2;
                    flux(iFace, iBase) = dot(fun(mean(FcNd, 2), P), nor);
                end
            end
            tc.verifyEqual(flux, eye(4), "AbsTol", 1e-12);
            % 3D BDM1: moments of normal component against barycentric coordinates (1 - s - t, s, t) of each face.
            % DoF index: iFace + 4 * (k - 1) for k-th test function on iFace-th face.
            fE = stdFE("BDM1");
            tc.verifyEqual(fE.nDoF, 12);
            gInt = GInt("D3F", 2);
            mom = zeros(12);
            for iBase = 1:12
                fun = fE.base(iBase).getFun;
                for iFace = 1:4
                    FcNd = P(:, D3TElem.face.node(:, iFace));
                    nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1)); % Unit normal times 2 * area.
                    for iPnt = 1:size(gInt.pnt, 2)
                        st = gInt.pnt(:, iPnt);
                        vn = dot(fun(FcNd * [1 - sum(st); st], P), nor);
                        mom(iFace + [0, 4, 8], iBase) = mom(iFace + [0, 4, 8], iBase) + gInt.wgt(iPnt) * vn * [1 - sum(st); st];
                    end
                end
            end
            tc.verifyEqual(mom, eye(12), "AbsTol", 1e-12);
            % 3D matrix-valued (rows mapped, Stokes SDG_1 stress): moments of (sigma n)_i against barycentric coordinates
            % of faces, DoF groups by face set (faces 1, 2, 3; face 4) and row i.
            % Values of base functions are evaluated numerically (`mapVal`, checked against symbolic base functions in
            % `mapValues`).
            fE = stdFE("StSDG1M");
            tc.verifyEqual(fE.nDoF, 36);
            mom = zeros(36);
            iDoF0 = 0;
            for iGrp = 1:length(fE.DoFs)
                dof = fE.DoFs(iGrp);
                iRow = find(~isnan(dof.ord(1, :)), 1);
                for iEnt = 1:dof.nEnt
                    FcNd = P(:, D3TElem.face.node(:, dof.EntIdx(iEnt)));
                    nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
                    lam = [1 - sum(gInt.pnt, 1); gInt.pnt];
                    % S(k, i, p): k-th component (column-major 3 x 3) of i-th base function at p-th Gauss point.
                    S = mapVal(fE.RefBase, fE.map, zeros(3, 9), FcNd * lam, P, fE.RefKey);
                    Sn = reshape(sum(reshape(S(iRow:3:9, :, :), 3, 36, []) .* nor, 1), 36, []);
                    idx = iDoF0 + iEnt + dof.nEnt * (0:2);
                    mom(idx, :) = mom(idx, :) + (lam .* gInt.wgt) * Sn.';
                end
                iDoF0 = iDoF0 + dof.nDoF;
            end
            tc.verifyEqual(mom, eye(36), "AbsTol", 1e-10);
            % Invalid: nodal DoF, coefficient not on facet, scalar function space.
            tc.verifyError(@() FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", ...
                NdDoF("D3", D3TElem, 2, [1/3; 1/3], zeros(3), "coef", MshEnt("D3F").UNV, "orien", true), "map", "piolaDiv"), ...
                "MATLAB:assertion:failed");
            tc.verifyError(@() FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", ...
                MoDoF("D3", D3TElem, 2, Fcn.cst(1), [0, nan, nan; 0, nan, nan; 0, nan, nan], "GInt", GInt("D3F", 1)), ...
                "map", "piolaDiv"), "MATLAB:assertion:failed");
            tc.verifyError(@() FE("D3T", "[1,x,y,z]", MoDoF("D3", D3TElem, 2, Fcn.cst(1), [0; 0; 0], "GInt", GInt("D3F", 1)), ...
                "map", "piolaDiv"), "MATLAB:assertion:failed");
        end
        function mapValues(tc)
            % Numerical evaluation of mapped base functions (mapVal) equals derivatives of symbolic base functions.
            grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
            gradV = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
            UNV2 = MshEnt("D2L").UNV;
            RT0_2D = FE("D2T", "[1,0; 0,1; x,y].'", MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), zeros(2), "coef", UNV2, ...
                "orien", true, "GInt", GInt("D2L", 4)), "map", "piolaDiv");
            cases = {stdFE("P2"), [0; 0; 0]; stdFE("P2"), grad; stdFE("P2"), [2; 0; 0]; stdFE("P2"), [0; 1; 1]; ...
                stdFE("BDM1"), zeros(3); stdFE("BDM1"), eye(3); stdFE("BDM1"), gradV; ...
                stdFE("NED1"), zeros(3); stdFE("NED1"), gradV; RT0_2D, zeros(2); RT0_2D, eye(2); ...
                stdFE("StSDG1M"), zeros(3, 9); stdFE("StSDG1M"), repelem(eye(3), 1, 3)};
            for iCase = 1:size(cases, 1)
                fE = cases{iCase, 1};
                ord = cases{iCase, 2};
                if isequal(fE.elem, "D3T")
                    P = tc.P3T;
                    bary = [0.1, 0.4, 0.25; 0.2, 0.1, 0.25; 0.3, 0.2, 0.25; 0.4, 0.3, 0.25];
                else
                    P = tc.P2T;
                    bary = [0.2, 0.5, 0.1; 0.3, 0.1, 0.6; 0.5, 0.4, 0.3];
                end
                X = P * bary;
                val = mapVal(fE.RefBase, fE.map, ord, X, P);
                % Symbolic derivatives of mapped base functions are expensive: check at most 6 base functions.
                iBases = unique(round(linspace(1, fE.nDoF, min(fE.nDoF, 6))));
                % Cached function handles (by key) give the same values.
                tc.verifyNotEqual(fE.RefKey, "");
                for iCall = 1:2
                    tc.verifyEqual(mapVal(fE.RefBase, fE.map, ord, X, P, fE.RefKey), val, "AbsTol", 1e-14);
                end
                for iBase = iBases
                    fun = fE.base(iBase).dif(ord).getFun;
                    for iPnt = 1:size(X, 2)
                        ref = fun(X(:, iPnt), P);
                        tc.verifyEqual(val(:, iBase, iPnt), ref(:), "AbsTol", 1e-10, sprintf("Case %d, base %d.", iCase, iBase));
                    end
                end
                % Batch of two elements (second one reflected and translated) equals elements one by one.
                P2 = P(:, [1, 3, 2, 4:end]) + 0.3;
                X2 = P2 * bary;
                val2 = mapVal(fE.RefBase, fE.map, ord, cat(3, X, X2), cat(3, P, P2));
                tc.verifyEqual(val2(:, :, :, 1), val, "AbsTol", 1e-12);
                tc.verifyEqual(val2(:, :, :, 2), mapVal(fE.RefBase, fE.map, ord, X2, P2), "AbsTol", 1e-12);
            end
            % Base functions without map have no reference base functions.
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 0, [], [0; 0; 0]));
            tc.verifyEmpty(fE.RefBase);
            tc.verifyEqual(fE.RefKey, "");
        end
        function piolaCurl(tc)
            % 2D lowest order Nedelec: base functions with map = "piolaCurl" equal those with map = "none".
            D2TElem = MshEnt("D2T").msh;
            dof = MoDoF("D2", D2TElem, 1, Fcn.cst(1), zeros(2), "coef", MshEnt("D2L").UTV, "orien", true, "GInt", GInt("D2L", 4));
            fE1 = FE("D2T", "[1,0; 0,1; -y,x].'", dof);
            fE2 = FE("D2T", "[1,0; 0,1; -y,x].'", dof, "map", "piolaCurl");
            pnts = tc.P2T * [0.2, 0.5, 0.1; 0.3, 0.1, 0.6; 0.5, 0.4, 0.3];
            for iBase = 1:3
                fun1 = fE1.base(iBase).getFun;
                fun2 = fE2.base(iBase).getFun;
                for iPnt = 1:3
                    tc.verifyEqual(fun2(pnts(:, iPnt), tc.P2T), fun1(pnts(:, iPnt), tc.P2T), "AbsTol", 1e-12);
                end
            end
            % 3D lowest order Nedelec: tangential moment of j-th base function on i-th edge is delta_ij.
            D3TElem = MshEnt("D3T").msh;
            fE = stdFE("NED1");
            tc.verifyEqual(fE.nDoF, 6);
            P = tc.P3T;
            mom = zeros(6);
            for iBase = 1:6
                fun = fE.base(iBase).getFun;
                for iEdge = 1:6
                    EgNd = P(:, D3TElem.edge.node(:, iEdge));
                    tan = EgNd(:, 2) - EgNd(:, 1); % Unit tangent times length.
                    % Tangential component is linear on edge: 2-point Gauss rule is exact.
                    for s = (1 + [-1, 1] / sqrt(3)) / 2
                        mom(iEdge, iBase) = mom(iEdge, iBase) + dot(fun(EgNd * [1 - s; s], P), tan) / 2;
                    end
                end
            end
            tc.verifyEqual(mom, eye(6), "AbsTol", 1e-12);
            % Invalid: moments on faces, normal coefficient.
            tc.verifyError(@() FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", MoDoF("D3", D3TElem, 2, Fcn.cst(1), zeros(3), ...
                "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 1)), "map", "piolaCurl"), "MATLAB:assertion:failed");
        end
        function TP1_3D(tc)
            fE = stdFE("TP1");
            tc.verifyEqual(fE.nDoF, 3);
            tc.verifyEqual(fE.base(1).domn, "D3FR");
            val = sym(zeros(3));
            for iBase = 1:3
                val(:, iBase) = [fE.base(iBase).eval([0; 0]); fE.base(iBase).eval([1; 0]); fE.base(iBase).eval([0; 1])];
            end
            tc.verifyTrue(symEq(val, sym(eye(3))));
        end
        function repFS(tc)
            tc.verifyTrue(symEq(FE.repFS("[1,x,y]", [2, 1]), str2sym("[1, x, y, 0, 0, 0; 0, 0, 0, 1, x, y]")));
            FS = FE.repFS("[1,x,y]", [2, 2]);
            tc.verifyEqual(size(FS), [2, 2, 12]);
            tc.verifyTrue(symEq(FS(:, :, 1), sym([1, 0; 0, 0])));
            tc.verifyTrue(symEq(FS(:, :, 5), str2sym("[0, 0; x, 0]")));
            tc.verifyTrue(symEq(FS(:, :, 9), str2sym("[0, y; 0, 0]")));
            tc.verifyTrue(symEq(FS(:, :, 12), str2sym("[0, 0; 0, y]")));
        end
    end
end
%% Local functions.
function val = evalPnt(bases, pnts)
    % evalPnt: evaluate scalar base functions at points, val(i, j) = bases(j)(pnts(i)).
    val = sym(zeros(length(pnts), length(bases)));
    for iBase = 1:length(bases)
        for iPnt = 1:length(pnts)
            val(iPnt, iBase) = bases(iBase).eval(pnts(iPnt));
        end
    end
end
function val = evalNum(bases, pnts, parm)
    % evalNum: evaluate scalar base functions at points numerically, val(i, j) = bases(j)(pnts(:, i)).
    val = zeros(size(pnts, 2), length(bases));
    for iBase = 1:length(bases)
        fun = bases(iBase).getFun;
        val(:, iBase) = fun(pnts, parm);
    end
end
