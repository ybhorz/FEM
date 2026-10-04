classdef FESTest < matlab.unittest.TestCase
    % FESTest: unit tests of FES and FEF.
    % Mesh: mshD2TS([0, 1, 0, 1], 2), 9 nodes, 16 edges, 8 elements.
    % 3D mesh: mshD3TS([0, 1, 0, 1, 0, 1], 2), 27 nodes, 98 edges, 120 faces, 48 elements.
    % Local-to-global maps are checked against mesh data, and projection must reproduce functions in the space.
    properties
        msh; % Mesh.
        msh3; % 3D mesh.
    end
    methods (TestClassSetup)
        function setMsh(tc)
            tc.msh = mshD2TS([0, 1, 0, 1], 2);
            tc.msh3 = mshD3TS([0, 1, 0, 1, 0, 1], 2);
        end
    end
    methods (Test, TestTags = {'Fast'})
        function P0(tc)
            fES = FES(tc.msh, FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0])));
            tc.verifyEqual(fES.nGlDoF, 8);
            tc.verifyEqual(fES.Lc2Gl, 1:8);
            verifyProj(tc, fES, Fcn("D2", "3"));
        end
        function P1(tc)
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])));
            tc.verifyEqual(fES.nGlDoF, 9);
            tc.verifyEqual(fES.Lc2Gl, tc.msh.elem.node);
            tc.verifyEqual(fES.ElParm, reshape(tc.msh.node.coord(:, tc.msh.elem.node), 2, 3, 8));
            verifyProj(tc, fES, Fcn("D2", "x + 2*y - 1"));
        end
        function P2(tc)
            fES = FES(tc.msh, FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
                NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]));
            tc.verifyEqual(fES.nGlDoF, 25);
            tc.verifyEqual(fES.Lc2Gl, [tc.msh.elem.node; 9 + abs(tc.msh.elem.edge)]);
            verifyProj(tc, fES, Fcn("D2", "x^2 + x*y - 2*y^2 + x"));
        end
        function DG1(tc)
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false)));
            tc.verifyEqual(fES.nGlDoF, 24);
            tc.verifyEqual(fES.Lc2Gl, reshape(1:24, 3, 8));
            verifyProj(tc, fES, Fcn("D2", "x + 2*y - 1"));
        end
        function CR(tc)
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])));
            tc.verifyEqual(fES.nGlDoF, 16);
            tc.verifyEqual(fES.Lc2Gl, abs(tc.msh.elem.edge));
            verifyProj(tc, fES, Fcn("D2", "x + 2*y - 1"));
        end
        function VP1(tc)
            fES = FES(tc.msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
                NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]));
            tc.verifyEqual(fES.nGlDoF, 18);
            tc.verifyEqual(fES.Lc2Gl, [tc.msh.elem.node; 9 + tc.msh.elem.node]);
            verifyProj(tc, fES, Fcn("D2", "[x + y; x - 2*y]"));
        end
        function RT0(tc)
            fES = FES(tc.msh, FE("D2T", "[1,0; 0,1; x,y].'", ...
                MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4))));
            tc.verifyEqual(fES.nGlDoF, 16);
            % Sign of base function follows orientation of edge in element.
            tc.verifyEqual(fES.Lc2Gl, tc.msh.elem.edge);
            verifyProj(tc, fES, Fcn("D2", "[1 + 2*x; -1 + 2*y]"));
        end
        function BDM1(tc)
            fES = FES(tc.msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4)), ...
                MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4))]));
            tc.verifyEqual(fES.nGlDoF, 32);
            % Second moment: sign of normal and sign of odd test function cancel.
            tc.verifyEqual(fES.Lc2Gl, [tc.msh.elem.edge; 16 + abs(tc.msh.elem.edge)]);
            verifyProj(tc, fES, Fcn("D2", "[x + 2*y; 3*x - y]"));
        end
        function TP1(tc)
            fES = FES(tc.msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0)));
            tc.verifyEqual(fES.nGlDoF, 32);
            tc.verifyEqual(fES.Lc2Gl, [1:16; 17:32]);
            tc.verifyEqual(fES.ElParm, reshape(tc.msh.node.coord(:, tc.msh.edge.node), 2, 2, 16));
            % Projection reproduces linear function on each edge.
            fcn = Fcn("D2", "x + 2*y - 1");
            fEF = fES.proj(fcn);
            fun = fEF.getFun; exFun = fcn.getFun;
            for iEdge = 1:tc.msh.nEdge
                EgParm = fEF.ElParm(:, :, iEdge);
                for s = [0.25, 0.6]
                    tc.verifyEqual(fun(s, EgParm, fEF.ElCoef(:, iEdge)), exFun((1 - s) * EgParm(:, 1) + s * EgParm(:, 2)), ...
                        "AbsTol", 1e-12);
                end
            end
        end
        function P0_3D(tc)
            fE = stdFE("P0");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 48);
            tc.verifyEqual(fES.Lc2Gl, 1:48);
            verifyProj(tc, fES, Fcn("D3", "3"));
        end
        function P1_3D(tc)
            fE = stdFE("P1");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 27);
            tc.verifyEqual(fES.Lc2Gl, tc.msh3.elem.node);
            tc.verifyEqual(fES.ElParm, reshape(tc.msh3.node.coord(:, tc.msh3.elem.node), 3, 4, 48));
            verifyShare(tc, fES, fE);
            verifyProj(tc, fES, Fcn("D3", "x + 2*y - 3*z + 1"));
        end
        function P2_3D(tc)
            fE = stdFE("P2");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 125);
            tc.verifyEqual(fES.Lc2Gl, [tc.msh3.elem.node; 27 + abs(tc.msh3.elem.edge)]);
            verifyShare(tc, fES, fE);
            verifyProj(tc, fES, Fcn("D3", "x^2 - y*z + 2*z^2 + x*y - y + 1"));
        end
        function P3_3D(tc)
            % Two samples on each edge: reversed if edge orientation in element is opposite to global orientation.
            fE = stdFE("P3");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 27 + 2 * 98 + 120);
            tc.verifyTrue(all(fES.Lc2Gl > 0, "all"));
            verifyShare(tc, fES, fE);
            verifyProj(tc, fES, Fcn("D3", "x^3 - x*y*z + y^2*z - z^3 + x^2 - y"));
        end
        function DG1_3D(tc)
            fE = stdFE("DG1");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 192);
            tc.verifyEqual(fES.Lc2Gl, reshape(1:192, 4, 48));
            verifyShare(tc, fES, fE);
            verifyProj(tc, fES, Fcn("D3", "x + 2*y - 3*z + 1"));
        end
        function CR_3D(tc)
            fE = stdFE("CR");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.nGlDoF, 120);
            tc.verifyEqual(fES.Lc2Gl, abs(tc.msh3.elem.face));
            verifyShare(tc, fES, fE);
            verifyProj(tc, fES, Fcn("D3", "x + 2*y - 3*z + 1"));
        end
        function RT0_3D(tc)
            msh = tc.msh3;
            fE = stdFE("RT0");
            fES = FES(msh, fE);
            tc.verifyEqual(fES.nGlDoF, 120);
            % Sign of base function follows orientation of face in element.
            tc.verifyEqual(fES.Lc2Gl, msh.elem.face);
            verifyProj(tc, fES, Fcn("D3", "[1 + 2*x; -1 + 2*y; 3 + 2*z]"));
            % Normal continuity: for random DoF values, v.n from both sides of each interior face agree.
            rng(4);
            fEF = FEF(fES, rand(fES.nGlDoF, 1));
            fun = fEF.getFun;
            IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, 1:6));
            for iFace = IntFace
                FcNd = msh.node.coord(:, msh.face.node(:, iFace));
                nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
                pnt = FcNd * [0.2; 0.3; 0.5];
                K = abs(msh.face.elem(:, iFace));
                v1 = fun(pnt, fEF.ElParm(:, :, K(1)), fEF.ElCoef(:, K(1)));
                v2 = fun(pnt, fEF.ElParm(:, :, K(2)), fEF.ElCoef(:, K(2)));
                tc.verifyEqual(dot(v1, nor), dot(v2, nor), "AbsTol", 1e-12, sprintf("Face %d.", iFace));
            end
        end
        function facePermTab(tc)
            % Permutation of samples on face: local sample i of a face with code c is global sample permTab(i, c).
            P = [1, 2, 3; 2, 3, 1; 3, 1, 2; 1, 3, 2; 3, 2, 1; 2, 1, 3];
            D3TElem = MshEnt("D3T").msh;
            G = [0.1, 1.2, 0.3; -0.2, 0.4, 1.1; 0.3, 0.2, 0.9]; % Global face nodes (by column).
            % Barycenter, orbit of 3 points, orbit of 6 points.
            crds = {[1/3; 1/3], [1/4, 1/2, 1/4; 1/4, 1/4, 1/2], ...
                [0.3, 0.6, 0.1, 0.6, 0.1, 0.3; 0.6, 0.3, 0.6, 0.1, 0.3, 0.1]};
            for iCrd = 1:length(crds)
                crd = crds{iCrd};
                permTab = FES.facePerm(NdDoF("D3", D3TElem, 2, crd, [0; 0; 0]));
                for c = 1:6
                    L = G(:, P(c, :));
                    for iSamp = 1:size(crd, 2)
                        jSamp = permTab(iSamp, c);
                        tc.verifyEqual(L * [1 - sum(crd(:, iSamp)); crd(:, iSamp)], G * [1 - sum(crd(:, jSamp)); crd(:, jSamp)], ...
                            "AbsTol", 1e-14, sprintf("Samples %d, code %d.", iCrd, c));
                    end
                end
            end
            % Not invariant under permutation.
            tc.verifyError(@() FES.facePerm(NdDoF("D3", D3TElem, 2, [1/6; 2/3], [0; 0; 0])), "MATLAB:assertion:failed");
            % Moments against barycentric coordinates of face: lambda_L(i) = lambda_G(P(c, i)).
            tst = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")];
            [permTab, sgnTab] = FES.facePerm(MoDoF("D3", D3TElem, 2, tst, [0; 0; 0], "GInt", GInt("D3F", 2)));
            tc.verifyEqual(permTab, P');
            tc.verifyEqual(sgnTab, ones(3, 6));
            % Odd test functions change sign: s - t under codes with swapped 2nd and 3rd nodes.
            tst = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s - t")];
            tc.verifyError(@() FES.facePerm(MoDoF("D3", D3TElem, 2, tst, [0; 0; 0], "GInt", GInt("D3F", 2))), ...
                "MATLAB:assertion:failed");
        end
        function BDM1_3D(tc)
            msh = tc.msh3;
            fES = FES(msh, stdFE("BDM1"));
            tc.verifyEqual(fES.nGlDoF, 360);
            verifyProj(tc, fES, Fcn("D3", "[1 + 2*x - y; z - 3*x + 1; x + y + 2*z]"));
            % Normal continuity: for random DoF values, v.n from both sides of each interior face agree.
            rng(5);
            fEF = FEF(fES, rand(fES.nGlDoF, 1));
            fun = fEF.getFun;
            IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, 1:6));
            for iFace = IntFace
                FcNd = msh.node.coord(:, msh.face.node(:, iFace));
                nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
                K = abs(msh.face.elem(:, iFace));
                for bary = [0.2, 0.7, 0.1; 0.3, 0.1, 0.6; 0.5, 0.2, 0.3]
                    pnt = FcNd * bary;
                    v1 = fun(pnt, fEF.ElParm(:, :, K(1)), fEF.ElCoef(:, K(1)));
                    v2 = fun(pnt, fEF.ElParm(:, :, K(2)), fEF.ElCoef(:, K(2)));
                    tc.verifyEqual(dot(v1, nor), dot(v2, nor), "AbsTol", 1e-12, sprintf("Face %d.", iFace));
                end
            end
        end
        function NED1_3D(tc)
            msh = tc.msh3;
            fES = FES(msh, stdFE("NED1"));
            tc.verifyEqual(fES.nGlDoF, 98);
            % Sign of base function follows orientation of edge in element.
            tc.verifyEqual(fES.Lc2Gl, msh.elem.edge);
            % a + b x x with a = [1; 2; 3], b = [1; -1; 2].
            verifyProj(tc, fES, Fcn("D3", "[1 - z - 2*y; 2 + 2*x - z; 3 + y + x]"));
            % Tangential continuity: for random DoF values, v.t agrees in all elements sharing an edge.
            rng(6);
            fEF = FEF(fES, rand(fES.nGlDoF, 1));
            fun = fEF.getFun;
            for iEdge = 1:msh.nEdge
                EgNd = msh.node.coord(:, msh.edge.node(:, iEdge));
                tan = EgNd(:, 2) - EgNd(:, 1);
                [~, K] = find(abs(msh.elem.edge) == iEdge);
                for s = [0.3, 0.8]
                    pnt = EgNd * [1 - s; s];
                    vt = zeros(1, length(K));
                    for iK = 1:length(K)
                        vt(iK) = dot(fun(pnt, fEF.ElParm(:, :, K(iK)), fEF.ElCoef(:, K(iK))), tan);
                    end
                    tc.verifyEqual(vt, repmat(vt(1), 1, length(K)), "AbsTol", 1e-12, sprintf("Edge %d.", iEdge));
                end
            end
        end
        function TP1_3D(tc)
            msh = tc.msh3;
            fES = FES(msh, stdFE("TP1"));
            tc.verifyEqual(fES.nGlDoF, 360);
            tc.verifyEqual(fES.Lc2Gl, [1:120; 121:240; 241:360]);
            tc.verifyEqual(fES.ElParm, reshape(msh.node.coord(:, msh.face.node), 3, 3, 120));
            % Projection reproduces linear function on each face.
            fcn = Fcn("D3", "x + 2*y - z + 1");
            fEF = fES.proj(fcn);
            fun = fEF.getFun; exFun = fcn.getFun;
            for iFace = 1:msh.nFace
                FcParm = fEF.ElParm(:, :, iFace);
                for st = [0.2, 0.5; 0.3, 0.1]
                    tc.verifyEqual(fun(st, FcParm, fEF.ElCoef(:, iFace)), exFun(FcParm * [1 - sum(st); st]), "AbsTol", 1e-12);
                end
            end
            % Boundary condition on trace space: DoFs on given faces.
            BdFace = msh.bdEnt(2, 1:6);
            fES = FES(msh, stdFE("TP1"), BC(fcn, "face", BdFace));
            tc.verifyEqual(fES.BC.DoFIdx, [BdFace, 120 + BdFace, 240 + BdFace]);
        end
        function SDG_3D(tc)
            % SDG elements on Alfeld split mesh: scalar continuous on interior primal faces, normal component of vector
            % continuous on dual faces.
            msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], 1));
            PrOFace = find(msh.face.type == 0);
            DlFace = find(msh.face.type == 1i);
            nPr = nnz(msh.face.type ~= 1i);
            nDl = length(DlFace);
            % SDG_0.
            fES = FES(msh, stdFE("SDG0S"));
            tc.verifyEqual(fES.nGlDoF, nPr);
            verifyShare(tc, fES, stdFE("SDG0S"));
            verifyProj(tc, fES, Fcn("D3", "2"));
            fES = FES(msh, stdFE("SDG0V"));
            tc.verifyEqual(fES.nGlDoF, nDl);
            verifyProj(tc, fES, Fcn("D3", "[1; -2; 3]"));
            verifyFaceCont(tc, fES, DlFace, "normal");
            % SDG_1.
            fES = FES(msh, stdFE("SDG1S"));
            tc.verifyEqual(fES.nGlDoF, 3 * nPr + msh.nElem);
            verifyShare(tc, fES, stdFE("SDG1S"));
            verifyProj(tc, fES, Fcn("D3", "1 + 2*x - y + 3*z"));
            verifyFaceCont(tc, fES, PrOFace, "value");
            fES = FES(msh, stdFE("SDG1V"));
            tc.verifyEqual(fES.nGlDoF, 3 * nDl + 3 * msh.nElem);
            verifyProj(tc, fES, Fcn("D3", "[1 + 2*x - y; z - 3*x + 1; x + y + 2*z]"));
            verifyFaceCont(tc, fES, DlFace, "normal");
            % Not continuous on the other faces: scalar jumps on some dual face, normal component on some primal face.
            rng(7);
            fEF = FEF(FES(msh, stdFE("SDG1S")), rand(3 * nPr + msh.nElem, 1));
            tc.verifyGreaterThan(maxJump(fEF, DlFace, "value"), 1e-3);
            fEF = FEF(FES(msh, stdFE("SDG1V")), rand(3 * nDl + 3 * msh.nElem, 1));
            tc.verifyGreaterThan(maxJump(fEF, PrOFace, "normal"), 1e-3);
        end
        function StokesSDG_3D(tc)
            % SDG elements of Stokes on Alfeld split mesh: traction of matrix continuous on dual faces, vector continuous
            % on interior primal faces, scalar continuous on dual faces.
            msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], 1));
            PrOFace = find(msh.face.type == 0);
            DlFace = find(msh.face.type == 1i);
            nPr = nnz(msh.face.type ~= 1i);
            nDl = length(DlFace);
            fES = FES(msh, stdFE("StSDG1M"));
            tc.verifyEqual(fES.nGlDoF, 9 * nDl + 9 * msh.nElem);
            verifyProj(tc, fES, Fcn("D3", "[1 + x, y - z, 2*z; x - y, 3, z + 2*x; y, x + y + z, 1 - x]"));
            verifyFaceCont(tc, fES, DlFace, "traction");
            rng(7);
            tc.verifyGreaterThan(maxJump(FEF(fES, rand(fES.nGlDoF, 1)), PrOFace, "traction"), 1e-3);
            fES = FES(msh, stdFE("StSDG1V"));
            tc.verifyEqual(fES.nGlDoF, 9 * nPr + 3 * msh.nElem);
            verifyShare(tc, fES, stdFE("StSDG1V"));
            verifyProj(tc, fES, Fcn("D3", "[1 + 2*x - y; z - 3*x + 1; x + y + 2*z]"));
            verifyFaceCont(tc, fES, PrOFace, "value");
            tc.verifyGreaterThan(maxJump(FEF(fES, rand(fES.nGlDoF, 1)), DlFace, "value"), 1e-3);
            fES = FES(msh, stdFE("StSDG1S"));
            tc.verifyEqual(fES.nGlDoF, msh.nElem / 4 + nnz(msh.edge.type == 1i));
            verifyShare(tc, fES, stdFE("StSDG1S"));
            verifyProj(tc, fES, Fcn("D3", "1 + 2*x - y + 3*z"));
            verifyFaceCont(tc, fES, DlFace, "value");
            tc.verifyGreaterThan(maxJump(FEF(fES, rand(fES.nGlDoF, 1)), PrOFace, "value"), 1e-3);
        end
        function faceSign3D(tc)
            % Oriented DoF on face: sign of base function follows orientation of face in element.
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 2, [1/3; 1/3], [0; 0; 0], "orien", true), ...
                "map", "affine");
            fES = FES(tc.msh3, fE);
            tc.verifyEqual(fES.Lc2Gl, tc.msh3.elem.face);
            % Each interior face appears once with each sign, each boundary face once with positive sign.
            cnt = accumarray(abs(fES.Lc2Gl(:)), sign(fES.Lc2Gl(:)));
            BdFace = tc.msh3.bdEnt(2, 1:6);
            tc.verifyEqual(cnt(BdFace), ones(length(BdFace), 1));
            tc.verifyEqual(nnz(cnt), length(BdFace));
            % Samples on face must be invariant under permutation of face vertices.
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 2, [1/6; 2/3], [0; 0; 0]), "map", "affine");
            tc.verifyError(@() FES(tc.msh3, fE), "MATLAB:assertion:failed");
            % 2D element on 3D mesh.
            tc.verifyError(@() FES(tc.msh3, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]))), ...
                "MATLAB:assertion:failed");
        end
        function spaceMethod(tc)
            P1_FES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])));
            P0_FES = FES(tc.msh, FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0])));
            FESs = [P1_FES, P0_FES];
            tc.verifyEqual(FESs.cumDoF, 17);
            tc.verifyEqual(FESs.cumDoF(1), 9);
            tc.verifyEqual(FESs.getMsh, tc.msh);
            % FEF.multi splits DoF values.
            [fEF1, fEF2] = FEF.multi(FESs, (1:17)');
            tc.verifyEqual(fEF1.ElCoef, tc.msh.elem.node);
            tc.verifyEqual(fEF2.ElCoef, 10:17);
            % Transformation of FEF gives plain Fcn.
            fcn = tfm(fEF1, "D2TR");
            tc.verifyClass(fcn, "Fcn");
            tc.verifyEqual(fcn.domn, "D2TR");
            tc.verifyTrue(symEq(fcn.coef, fEF1.coef));
            % FEF.subt.
            fEF = subt(FEF(P1_FES, 3 * ones(9, 1)), FEF(P1_FES, ones(9, 1)));
            tc.verifyEqual(fEF.ElCoef, 2 * ones(3, 8));
        end
    end
end
%% Local functions.
function verifyProj(tc, fES, fcn)
    % verifyProj: verify that projection of function in FE space reproduces the function.
    fEF = fES.proj(fcn);
    fun = fEF.getFun; exFun = fcn.getFun;
    msh = fES.msh;
    % Barycentric coordinates of sample points (by column).
    switch msh.dim
        case 2
            bary = [1/3, 0.6, 0.2; 1/3, 0.3, 0.1; 1/3, 0.1, 0.7];
        case 3
            bary = [1/4, 0.1, 0.4; 1/4, 0.2, 0.1; 1/4, 0.3, 0.2; 1/4, 0.4, 0.3];
    end
    for iElem = 1:msh.nElem
        pnt = msh.node.coord(:, msh.elem.node(:, iElem)) * bary;
        for iPnt = 1:size(pnt, 2)
            val = fun(pnt(:, iPnt), fEF.ElParm(:, :, iElem), fEF.ElCoef(:, iElem));
            exVal = exFun(pnt(:, iPnt));
            tc.verifyEqual(val(:), exVal(:), "AbsTol", 1e-10, sprintf("Element %d, point %d.", iElem, iPnt));
        end
    end
end
function verifyFaceCont(tc, fES, FcIdx, mode)
    % verifyFaceCont: verify that for random DoF values, function (see `maxJump` for `mode`) from both sides of given
    % interior faces agrees.
    rng(5);
    fEF = FEF(fES, rand(fES.nGlDoF, 1));
    tc.verifyLessThan(maxJump(fEF, FcIdx, mode), 1e-12);
end
function jump = maxJump(fEF, FcIdx, mode)
    % maxJump: maximal jump across given interior faces at three points per face of
    % "value": function; "normal": normal component of vector; "traction": matrix times normal.
    msh = fEF.msh;
    fun = fEF.getFun;
    jump = 0;
    for iFace = FcIdx
        FcNd = msh.node.coord(:, msh.face.node(:, iFace));
        nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
        nor = nor / norm(nor);
        K = abs(msh.face.elem(:, iFace));
        for bary = [0.2, 0.7, 0.1; 0.3, 0.1, 0.6; 0.5, 0.2, 0.3]
            pnt = FcNd * bary;
            v1 = fun(pnt, fEF.ElParm(:, :, K(1)), fEF.ElCoef(:, K(1)));
            v2 = fun(pnt, fEF.ElParm(:, :, K(2)), fEF.ElCoef(:, K(2)));
            switch mode
                case "value"
                    jump = max(jump, max(abs(v1 - v2), [], "all"));
                case "normal"
                    jump = max(jump, abs(dot(v1 - v2, nor)));
                case "traction"
                    jump = max(jump, max(abs((v1 - v2) * nor)));
            end
        end
    end
end
function verifyShare(tc, fES, fE)
    % verifyShare: verify that each local nodal DoF and its global DoF sample the same physical point.
    msh = fES.msh;
    LcMsh = MshEnt(fE.elem).msh;
    for iDoF = 1:length(fE.DoFs)
        LcDoF = fE.DoFs(iDoF);
        GlDoF = fES.GlDoFs(iDoF);
        if LcDoF.EntDim == 0
            LcEntNode = 1:LcMsh.nNode;
        else
            LcEntNode = LcMsh.ent(LcDoF.EntDim).node;
            GlEntNode = msh.ent(LcDoF.EntDim).node;
        end
        LcCrd = double(LcDoF.coord);
        GlCrd = double(GlDoF.coord);
        for iElem = 1:msh.nElem
            ElNode = msh.elem.node(:, iElem);
            for iEnt = 1:LcDoF.nEnt
                for iSamp = 1:LcDoF.nSamp
                    GlIdx = abs(fES.Lc2Gl(fE.DoFs.sub2ind(iDoF, iEnt, iSamp), iElem));
                    [jDoF, jEnt, jSamp] = ind2sub(fES.GlDoFs, GlIdx);
                    tc.assertEqual(jDoF, iDoF);
                    LcVert = msh.node.coord(:, ElNode(LcEntNode(:, LcDoF.EntIdx(iEnt))));
                    if LcDoF.EntDim == 0
                        GlVert = msh.node.coord(:, GlDoF.EntIdx(jEnt));
                        LcPnt = LcVert;
                        GlPnt = GlVert;
                    else
                        GlVert = msh.node.coord(:, GlEntNode(:, GlDoF.EntIdx(jEnt)));
                        LcPnt = LcVert * [1 - sum(LcCrd(:, iSamp)); LcCrd(:, iSamp)];
                        GlPnt = GlVert * [1 - sum(GlCrd(:, jSamp)); GlCrd(:, jSamp)];
                    end
                    tc.verifyEqual(LcPnt, GlPnt, "AbsTol", 1e-14, ...
                        sprintf("DoF group %d, element %d, entity %d, sample %d.", iDoF, iElem, iEnt, iSamp));
                end
            end
        end
    end
end
