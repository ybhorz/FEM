classdef NormTest < matlab.unittest.TestCase
    % NormTest: unit tests of Norm and eNorm.
    % Mesh: mshD2TS([0, 1, 0, 1], 4), unit square.
    % 3D mesh: mshD3TS([0, 1, 0, 1, 0, 1], 2), unit cube.
    properties
        msh; % Mesh.
        IntEdge; % Indices of interior edges.
        P1_FES; % Continuous P1 space.
    end
    methods (TestClassSetup)
        function setMsh(tc)
            tc.msh = mshD2TS([0, 1, 0, 1], 4);
            tc.IntEdge = find(tc.msh.edge.type == 0);
            tc.P1_FES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])));
        end
    end
    methods (Test, TestTags = {'Fast'})
        function elemNorm(tc)
            Uh = tc.P1_FES;
            d0 = [0; 0]; grad = cat(3, [1; 0], [0; 1]);
            L2Norm = Norm(tc.msh, 2, d0, "GInt", GInt("D2T", 2));
            H1Norm = Norm(tc.msh, 2, grad, "GInt", GInt("D2T", 0));
            % Function in FE space has zero error.
            u = Fcn("D2", "x + 2*y - 1");
            tc.verifyEqual(eNorm(tc.msh, u, Uh.proj(u), L2Norm), 0, "AbsTol", 1e-12);
            tc.verifyEqual(eNorm(tc.msh, u, Uh.proj(u), H1Norm), 0, "AbsTol", 1e-12);
            % |1|_L2 = 1 and |x|_H1 = 1 on unit square.
            one = FEF(Uh, ones(Uh.nGlDoF, 1));
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), one, L2Norm), 1, "AbsTol", 1e-12);
            xh = FEF(Uh, tc.msh.node.coord(1, :)');
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), xh, H1Norm), 1, "AbsTol", 1e-12);
            % |x|_L2 = sqrt(1/3).
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), xh, L2Norm), sqrt(1/3), "AbsTol", 1e-12);
            % Combined norm: (|x|_L2^2 + |grad x|_L2^2)^(1/2).
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), xh, [L2Norm, H1Norm]), sqrt(1/3 + 1), "AbsTol", 1e-12);
            % Norms with different powers cannot be combined.
            L1Norm = Norm(tc.msh, 2, 0, "pow", 1, "GInt", GInt("D2T", 2));
            tc.verifyError(@() eNorm(tc.msh, Fcn.cst(0), xh, [L2Norm, L1Norm]), "MATLAB:assertion:failed");
        end
        function jumpNorm(tc)
            DG1_FES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false)));
            JNorm = Norm(tc.msh, 1, [0; 0], "EntIdx", tc.IntEdge, "fcnOpr", "jump", "GInt", GInt("D2L", 2));
            % Jump of continuous function vanishes.
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), DG1_FES.proj(Fcn("D2", "x + 2*y")), JNorm), 0, "AbsTol", 1e-12);
            % Jump of indicator of element 1: sqrt of total length of its interior edges.
            c = zeros(DG1_FES.nGlDoF, 1);
            c(abs(DG1_FES.Lc2Gl(:, 1))) = 1;
            EgLen = vecnorm(tc.msh.node.coord(:, tc.msh.edge.node(2, :)) - tc.msh.node.coord(:, tc.msh.edge.node(1, :)));
            EgIdx = tc.IntEdge(any(abs(tc.msh.edge.elem(:, tc.IntEdge)) == 1, 1));
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(0), FEF(DG1_FES, c), JNorm), sqrt(sum(EgLen(EgIdx))), "AbsTol", 1e-12);
        end
        function norm3D(tc)
            msh = mshD3TS([0, 1, 0, 1, 0, 1], 2);
            Uh = FES(msh, stdFE("P1"));
            d0 = [0; 0; 0]; grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
            L2Norm = Norm(msh, 3, d0, "GInt", GInt("D3T", 2));
            H1Norm = Norm(msh, 3, grad, "GInt", GInt("D3T", 0));
            u = Fcn("D3", "x + 2*y - 3*z + 1");
            tc.verifyEqual(eNorm(msh, u, Uh.proj(u), [L2Norm, H1Norm]), 0, "AbsTol", 1e-12);
            % |1|_L2 = 1, |x|_L2 = sqrt(1/3), |x + y + z|_H1 = sqrt(3) on unit cube.
            tc.verifyEqual(eNorm(msh, Fcn.cst(0), FEF(Uh, ones(27, 1)), L2Norm), 1, "AbsTol", 1e-12);
            tc.verifyEqual(eNorm(msh, Fcn.cst(0), FEF(Uh, msh.node.coord(1, :)'), L2Norm), sqrt(1/3), "AbsTol", 1e-12);
            tc.verifyEqual(eNorm(msh, Fcn.cst(0), FEF(Uh, sum(msh.node.coord, 1)'), H1Norm), sqrt(3), "AbsTol", 1e-12);
            % Error of P1 interpolant of x^2 in L2 norm is positive.
            tc.verifyGreaterThan(eNorm(msh, Fcn("D3", "x^2"), Uh.proj(Fcn("D3", "x^2")), L2Norm), 1e-3);
            % Jump norm on interior faces.
            DG1_FES = FES(msh, stdFE("DG1"));
            IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, 1:6));
            JNorm = Norm(msh, 2, d0, "EntIdx", IntFace, "fcnOpr", "jump", "GInt", GInt("D3F", 2));
            tc.verifyEqual(eNorm(msh, Fcn.cst(0), DG1_FES.proj(u), JNorm), 0, "AbsTol", 1e-12);
            c = zeros(DG1_FES.nGlDoF, 1);
            c(abs(DG1_FES.Lc2Gl(:, 1))) = 1;
            FcNd = @(i) msh.node.coord(:, msh.face.node(i, :));
            FcArea = vecnorm(cross(FcNd(2) - FcNd(1), FcNd(3) - FcNd(1))) / 2;
            FcIdx = IntFace(any(abs(msh.face.elem(:, IntFace)) == 1, 1));
            tc.verifyEqual(eNorm(msh, Fcn.cst(0), FEF(DG1_FES, c), JNorm), sqrt(sum(FcArea(FcIdx))), "AbsTol", 1e-12);
        end
        function traceNorm(tc)
            TP1_FES = FES(tc.msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0)));
            L2Norm = Norm(tc.msh, 1, 0, "GInt", GInt("D2L", 2));
            u = Fcn("D2", "x + 2*y");
            tc.verifyEqual(eNorm(tc.msh, u, TP1_FES.proj(u), L2Norm), 0, "AbsTol", 1e-12);
            % |1|_L2 on all edges = sqrt of total length of edges.
            EgLen = vecnorm(tc.msh.node.coord(:, tc.msh.edge.node(2, :)) - tc.msh.node.coord(:, tc.msh.edge.node(1, :)));
            tc.verifyEqual(eNorm(tc.msh, Fcn.cst(1), FEF(TP1_FES, zeros(TP1_FES.nGlDoF, 1)), L2Norm), sqrt(sum(EgLen)), "AbsTol", 1e-12);
        end
    end
end
