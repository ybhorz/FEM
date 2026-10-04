classdef AssembleTest < matlab.unittest.TestCase
    % AssembleTest: unit tests of assemble (DLF, SLF, LDLF, LSLF).
    % Mesh: mshD2TS([0, 1, 0, 1], 4), unit square.
    % 3D mesh: mshD3TS([0, 1, 0, 1, 0, 1], 2), unit cube, 27 nodes, 98 edges, 120 faces (48 on boundary), 48 elements.
    properties
        msh; % Mesh.
        BdNode; % Indices of boundary nodes.
        BdEdge; % Indices of boundary edges.
        IntEdge; % Indices of interior edges.
        P1_FES; % Continuous P1 space.
        DG1_FES; % Discontinuous P1 space.
        msh3; % 3D mesh.
        BdFace3; % Indices of boundary faces of 3D mesh.
        IntFace3; % Indices of interior faces of 3D mesh.
        P1_FES3; % Continuous P1 space on 3D mesh.
        DG1_FES3; % Discontinuous P1 space on 3D mesh.
    end
    methods (TestClassSetup)
        function setMsh(tc)
            tc.msh = mshD2TS([0, 1, 0, 1], 4);
            tc.BdNode = find(tc.msh.node.type ~= 0);
            tc.BdEdge = find(tc.msh.edge.type ~= 0);
            tc.IntEdge = find(tc.msh.edge.type == 0);
            tc.P1_FES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])));
            tc.DG1_FES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false)));
            tc.msh3 = mshD3TS([0, 1, 0, 1, 0, 1], 2);
            tc.BdFace3 = tc.msh3.bdEnt(2, 1:6);
            tc.IntFace3 = setdiff(1:tc.msh3.nFace, tc.BdFace3);
            tc.P1_FES3 = FES(tc.msh3, stdFE("P1"));
            tc.DG1_FES3 = FES(tc.msh3, stdFE("DG1"));
        end
    end
    methods (Test, TestTags = {'Fast'})
        function formProp(tc)
            % Default Gauss integration of forms follows mesh.
            d0 = [0; 0];
            warnState = warning("off");
            cleanup = onCleanup(@() warning(warnState));
            tc.verifyEqual(DLF(tc.msh, 2, Fcn.cst(1), d0, d0).GInt.domn, "D2T");
            tc.verifyEqual(DLF(tc.msh, 1, Fcn.cst(1), d0, d0).GInt.domn, "D2L");
            tc.verifyEqual(SLF(tc.msh, 2, Fcn.cst(1), d0).GInt.domn, "D2T");
            tc.verifyEqual(SLF(tc.msh, 1, Fcn.cst(1), d0).GInt.domn, "D2L");
            tc.verifyEqual(Norm(tc.msh, 2, d0).GInt.domn, "D2T");
            tc.verifyEqual(Norm(tc.msh, 1, d0).GInt.domn, "D2L");
            clear cleanup;
            % Invalid combinations of entity dimension, domain and Gauss integration.
            tc.verifyError(@() DLF(tc.msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2L", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() DLF(tc.msh, 2, Fcn("D2L", "x"), d0, d0, "GInt", GInt("D2T", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() SLF(tc.msh, 1, Fcn("D2T", "x"), d0, "GInt", GInt("D2L", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() Norm(tc.msh, 2, d0, "fcnOpr", "jump", "GInt", GInt("D2T", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() DLF.interface(tc.msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2T", 1)), "MATLAB:assertion:failed");
        end
        function mass(tc)
            Uh = tc.P1_FES;
            d0 = [0; 0];
            M = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2T", 2)), SLF.empty);
            tc.verifyEqual(size(M), [Uh.nGlDoF, Uh.nGlDoF]);
            tc.verifyTrue(issparse(M));
            tc.verifyEqual(full(sum(M, "all")), 1, "AbsTol", 1e-12);
            tc.verifyEqual(full(max(abs(M - M'), [], "all")), 0, "AbsTol", 1e-14);
            % Full matrix type gives the same result.
            MFull = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2T", 2)), SLF.empty, "matType", "full");
            tc.verifyEqual(MFull, full(M), "AbsTol", 1e-14);
        end
        function stiff(tc)
            Uh = tc.P1_FES;
            grad = cat(3, [1; 0], [0; 1]);
            K = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 0)), SLF.empty);
            tc.verifyEqual(full(max(abs(K - K'), [], "all")), 0, "AbsTol", 1e-13);
            % Constants are in the kernel.
            tc.verifyEqual(full(K * ones(Uh.nGlDoF, 1)), zeros(Uh.nGlDoF, 1), "AbsTol", 1e-12);
            % Energy of u = x is |Omega| = 1.
            u = tc.msh.node.coord(1, :)';
            tc.verifyEqual(u' * K * u, 1, "AbsTol", 1e-12);
        end
        function loadVec(tc)
            Uh = tc.P1_FES;
            d0 = [0; 0];
            [~, F] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF(tc.msh, 2, Fcn.cst(1), d0, "GInt", GInt("D2T", 1)));
            tc.verifyEqual(sum(F), 1, "AbsTol", 1e-12);
            [~, F] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF(tc.msh, 2, Fcn("D2", "x"), d0, "GInt", GInt("D2T", 2)));
            tc.verifyEqual(sum(F), 1/2, "AbsTol", 1e-12);
            % Boundary load: perimeter of unit square is 4.
            [~, F] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF(tc.msh, 1, Fcn.cst(1), d0, "EntIdx", tc.BdEdge, "GInt", GInt("D2L", 1)));
            tc.verifyEqual(sum(F), 4, "AbsTol", 1e-12);
            % Flux of [x; 0] through boundary is 1.
            [~, F] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF(tc.msh, 1, dot(Fcn("D2", "[x; 0]"), MshEnt("D2L").UNV), d0, ...
                "EntIdx", tc.BdEdge, "GInt", GInt("D2L", 2)));
            tc.verifyEqual(sum(F), 1, "AbsTol", 1e-12);
        end
        function interfaceForm(tc)
            Uh = tc.DG1_FES;
            d0 = [0; 0];
            J = assemble(tc.msh, Uh, Uh, DLF.interface(tc.msh, 1, Fcn.cst(1), d0, d0, "EntIdx", tc.IntEdge, ...
                "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D2L", 2)), SLF.empty);
            tc.verifyEqual(full(max(abs(J - J'), [], "all")), 0, "AbsTol", 1e-14);
            % Jump of continuous function vanishes.
            u = evalDoF(Uh, Fcn("D2", "x + 2*y - 1"));
            tc.verifyEqual(full(J * u), zeros(Uh.nGlDoF, 1), "AbsTol", 1e-12);
            % Jump penalty is positive semi-definite.
            tc.verifyGreaterThan(min(eig(full(J))), -1e-12);
            % Jump of elementwise constant function: sum over interior edges of |e| * jump^2.
            c = zeros(Uh.nGlDoF, 1);
            c(abs(Uh.Lc2Gl(:, 1))) = 1;
            EgLen = vecnorm(tc.msh.node.coord(:, tc.msh.edge.node(2, :)) - tc.msh.node.coord(:, tc.msh.edge.node(1, :)));
            EgIdx = tc.IntEdge(any(abs(tc.msh.edge.elem(:, tc.IntEdge)) == 1, 1));
            tc.verifyEqual(c' * J * c, sum(EgLen(EgIdx)), "AbsTol", 1e-12);
            % Consistency term of average and jump: -int {du/dn} [v] on continuous u = x.
            grad = cat(3, [1; 0], [0; 1]);
            A = assemble(tc.msh, Uh, Uh, DLF.interface(tc.msh, 1, MshEnt("D2L").UNV, grad, d0, "EntIdx", tc.IntEdge, ...
                "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D2L", 1)), SLF.empty);
            % For any test function continuous across edges, the jump vanishes.
            tc.verifyEqual(u' * A * evalDoF(Uh, Fcn("D2", "x")), 0, "AbsTol", 1e-12);
        end
        function traceSpace(tc)
            % Coupling between element space and trace space on edges.
            TP1_FES = FES(tc.msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0)));
            trls = [tc.DG1_FES, TP1_FES];
            d0 = [0; 0];
            % Mass matrix on edges of trace space: sum is total length of edges.
            M = assemble(tc.msh, trls, trls, DLF(tc.msh, 1, Fcn.cst(1), 0, 0, "iTrl", 2, "iTst", 2, "GInt", GInt("D2L", 2)), SLF.empty);
            EgLen = vecnorm(tc.msh.node.coord(:, tc.msh.edge.node(2, :)) - tc.msh.node.coord(:, tc.msh.edge.node(1, :)));
            tc.verifyEqual(full(sum(M, "all")), sum(EgLen), "AbsTol", 1e-12);
            % int_e u * mu over all edges (from positive side) with u = 1 and mu = 1.
            B = assemble(tc.msh, trls, trls, DLF(tc.msh, 1, Fcn.cst(1), d0, 0, "iTrl", 1, "iTst", 2, "GInt", GInt("D2L", 2)), SLF.empty);
            nDoF1 = tc.DG1_FES.nGlDoF;
            tc.verifyEqual(full(sum(B(nDoF1 + 1:end, 1:nDoF1), "all")), sum(EgLen), "AbsTol", 1e-12);
            % Load on trace space.
            [~, F] = assemble(tc.msh, trls, trls, DLF.empty, SLF(tc.msh, 1, Fcn("D2", "x"), 0, "iTst", 2, "GInt", GInt("D2L", 2)));
            EgMid = (tc.msh.node.coord(:, tc.msh.edge.node(1, :)) + tc.msh.node.coord(:, tc.msh.edge.node(2, :))) / 2;
            tc.verifyEqual(sum(F(nDoF1 + 1:end)), sum(EgLen .* EgMid(1, :)), "AbsTol", 1e-12);
        end
        function patchP1(tc)
            % Harmonic linear solution is reproduced exactly.
            u = Fcn("D2", "x + 2*y - 1");
            grad = cat(3, [1; 0], [0; 1]);
            Uh = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])), BC(u, "node", tc.BdNode));
            [K, F] = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 0)), SLF.empty);
            tc.verifyEqual(K \ F, evalDoF(Uh, u), "AbsTol", 1e-12);
        end
        function patchP2(tc)
            % Harmonic quadratic solution is reproduced exactly.
            u = Fcn("D2", "x^2 - y^2 + x*y");
            grad = cat(3, [1; 0], [0; 1]);
            Uh = FES(tc.msh, FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
                NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]), BC(u, "node", tc.BdNode, "edge", tc.BdEdge));
            [K, F] = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 2)), SLF.empty);
            tc.verifyEqual(K \ F, evalDoF(Uh, u), "AbsTol", 1e-12);
        end
        function formProp3D(tc)
            d0 = [0; 0; 0];
            warnState = warning("off");
            cleanup = onCleanup(@() warning(warnState));
            tc.verifyEqual(DLF(tc.msh3, 3, Fcn.cst(1), d0, d0).GInt.domn, "D3T");
            tc.verifyEqual(DLF(tc.msh3, 2, Fcn.cst(1), d0, d0).GInt.domn, "D3F");
            tc.verifyEqual(SLF(tc.msh3, 3, Fcn.cst(1), d0).GInt.domn, "D3T");
            tc.verifyEqual(SLF(tc.msh3, 2, Fcn.cst(1), d0).GInt.domn, "D3F");
            tc.verifyEqual(Norm(tc.msh3, 3, d0).GInt.domn, "D3T");
            tc.verifyEqual(Norm(tc.msh3, 2, d0).GInt.domn, "D3F");
            clear cleanup;
            % Forms on edges of 3D mesh or on volume of 2D mesh are invalid.
            tc.verifyError(@() DLF(tc.msh3, 1, Fcn.cst(1), d0, d0, "GInt", GInt("D3L", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() SLF(tc.msh3, 1, Fcn.cst(1), d0, "GInt", GInt("D3L", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() Norm(tc.msh3, 1, d0, "GInt", GInt("D3L", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() DLF(tc.msh, 3, Fcn.cst(1), [0; 0], [0; 0], "GInt", GInt("D3T", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() DLF(tc.msh3, 3, Fcn("D3F", "x"), d0, d0, "GInt", GInt("D3T", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() SLF(tc.msh3, 2, Fcn.cst(1), d0, "GInt", GInt("D3T", 1)), "MATLAB:assertion:failed");
            tc.verifyError(@() DLF.interface(tc.msh3, 3, Fcn.cst(1), d0, d0, "GInt", GInt("D3T", 1)), "MATLAB:assertion:failed");
        end
        function massStiff3D(tc)
            Uh = tc.P1_FES3;
            d0 = [0; 0; 0];
            grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
            M = assemble(tc.msh3, Uh, Uh, DLF(tc.msh3, 3, Fcn.cst(1), d0, d0, "GInt", GInt("D3T", 2)), SLF.empty);
            tc.verifyEqual(size(M), [27, 27]);
            tc.verifyEqual(full(sum(M, "all")), 1, "AbsTol", 1e-12);
            tc.verifyEqual(full(max(abs(M - M'), [], "all")), 0, "AbsTol", 1e-14);
            % Element mass matrix of P1: |K| (1 + delta_ij) / 20.
            M1 = assemble(tc.msh3, Uh, Uh, DLF(tc.msh3, 3, Fcn.cst(1), d0, d0, "EntIdx", 1, "GInt", GInt("D3T", 2)), SLF.empty);
            ElNode = tc.msh3.elem.node(:, 1);
            tc.verifyEqual(full(M1(ElNode, ElNode)), (ones(4) + eye(4)) / 20 / 48, "AbsTol", 1e-14);
            K = assemble(tc.msh3, Uh, Uh, DLF(tc.msh3, 3, Fcn.cst(1), grad, grad, "GInt", GInt("D3T", 0)), SLF.empty);
            tc.verifyEqual(full(max(abs(K - K'), [], "all")), 0, "AbsTol", 1e-13);
            tc.verifyEqual(full(K * ones(27, 1)), zeros(27, 1), "AbsTol", 1e-12);
            % Energy of u = x is 1; energy of u = x + y + z is 3.
            u = tc.msh3.node.coord(1, :)';
            tc.verifyEqual(u' * K * u, 1, "AbsTol", 1e-12);
            u = sum(tc.msh3.node.coord, 1)';
            tc.verifyEqual(u' * K * u, 3, "AbsTol", 1e-12);
        end
        function loadVec3D(tc)
            Uh = tc.P1_FES3;
            d0 = [0; 0; 0];
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 3, Fcn.cst(1), d0, "GInt", GInt("D3T", 1)));
            tc.verifyEqual(sum(F), 1, "AbsTol", 1e-12);
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 3, Fcn("D3", "x*y"), d0, "GInt", GInt("D3T", 3)));
            tc.verifyEqual(sum(F), 1/4, "AbsTol", 1e-12);
            % Boundary load: surface area of unit cube is 6.
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 2, Fcn.cst(1), d0, "EntIdx", tc.BdFace3, "GInt", GInt("D3F", 1)));
            tc.verifyEqual(sum(F), 6, "AbsTol", 1e-12);
            % Flux through boundary (normal of boundary face is outward): div [x; 0; 0] = 1, div [x; y; z] = 3.
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 2, dot(Fcn("D3", "[x; 0; 0]"), MshEnt("D3F").UNV), d0, ...
                "EntIdx", tc.BdFace3, "GInt", GInt("D3F", 1)));
            tc.verifyEqual(sum(F), 1, "AbsTol", 1e-12);
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 2, dot(Fcn("D3", "[x; y; z]"), MshEnt("D3F").UNV), d0, ...
                "EntIdx", tc.BdFace3, "GInt", GInt("D3F", 1)));
            tc.verifyEqual(sum(F), 3, "AbsTol", 1e-12);
            % Flux of [x^2; 0; 0] weighted by P1 test functions: int_{x=1} v = load on face x = 1 only.
            [~, F] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 2, dot(Fcn("D3", "[x^2; 0; 0]"), MshEnt("D3F").UNV), d0, ...
                "EntIdx", tc.BdFace3, "GInt", GInt("D3F", 2)));
            [~, FRef] = assemble(tc.msh3, Uh, Uh, DLF.empty, SLF(tc.msh3, 2, Fcn.cst(1), d0, "EntIdx", tc.msh3.bdEnt(2, 2), ...
                "GInt", GInt("D3F", 1)));
            tc.verifyEqual(F, FRef, "AbsTol", 1e-12);
        end
        function interfaceForm3D(tc)
            msh = tc.msh3;
            Uh = tc.DG1_FES3;
            d0 = [0; 0; 0];
            J = assemble(msh, Uh, Uh, DLF.interface(msh, 2, Fcn.cst(1), d0, d0, "EntIdx", tc.IntFace3, ...
                "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D3F", 2)), SLF.empty);
            tc.verifyEqual(full(max(abs(J - J'), [], "all")), 0, "AbsTol", 1e-14);
            % Jump of continuous function vanishes.
            u = evalDoF(Uh, Fcn("D3", "x + 2*y - 3*z + 1"));
            tc.verifyEqual(full(J * u), zeros(Uh.nGlDoF, 1), "AbsTol", 1e-12);
            tc.verifyGreaterThan(min(eig(full(J))), -1e-12);
            % Jump of indicator of element 1: total area of its interior faces.
            c = zeros(Uh.nGlDoF, 1);
            c(abs(Uh.Lc2Gl(:, 1))) = 1;
            FcNd = @(i) msh.node.coord(:, msh.face.node(i, :));
            FcArea = vecnorm(cross(FcNd(2) - FcNd(1), FcNd(3) - FcNd(1))) / 2;
            FcIdx = tc.IntFace3(any(abs(msh.face.elem(:, tc.IntFace3)) == 1, 1));
            tc.verifyEqual(c' * J * c, sum(FcArea(FcIdx)), "AbsTol", 1e-12);
            % Consistency term int {du/dn} [v] (trial u, test v): vanishes for test function continuous across faces.
            grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
            A = assemble(msh, Uh, Uh, DLF.interface(msh, 2, MshEnt("D3F").UNV, grad, d0, "EntIdx", tc.IntFace3, ...
                "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D3F", 1)), SLF.empty);
            tc.verifyEqual(evalDoF(Uh, Fcn("D3", "x*y + z"))' * A * u, 0, "AbsTol", 1e-12);
            % Test function c (indicator of element K with a boundary face): sum of int_F {du/dn} [c] over interior faces of K
            % is the outward flux of grad u through them. Total outward flux of constant field grad u = [1; 2; -3] through the
            % surface of K is 0, so it equals minus the outward flux through boundary faces of K.
            K = abs(msh.face.elem(1, tc.BdFace3(1)));
            c = zeros(Uh.nGlDoF, 1);
            c(abs(Uh.Lc2Gl(:, K))) = 1;
            BdIdx = tc.BdFace3(any(abs(msh.face.elem(:, tc.BdFace3)) == K, 1));
            [~, FBd] = assemble(msh, Uh, Uh, DLF.empty, SLF(msh, 2, dot(Fcn("D3", "[1; 2; -3]"), MshEnt("D3F").UNV), d0, ...
                "EntIdx", BdIdx, "GInt", GInt("D3F", 1)));
            tc.verifyEqual(c' * A * u, -c' * FBd, "AbsTol", 1e-12);
            tc.verifyGreaterThan(abs(c' * FBd), 1e-3);
        end
        function traceSpace3D(tc)
            % Coupling between element space and trace space on faces.
            msh = tc.msh3;
            TP1_FES = FES(msh, stdFE("TP1"));
            trls = [tc.DG1_FES3, TP1_FES];
            d0 = [0; 0; 0];
            FcNd = @(i) msh.node.coord(:, msh.face.node(i, :));
            FcArea = vecnorm(cross(FcNd(2) - FcNd(1), FcNd(3) - FcNd(1))) / 2;
            % Mass matrix of trace space: sum is total area of faces.
            M = assemble(msh, trls, trls, DLF(msh, 2, Fcn.cst(1), [0; 0], [0; 0], "iTrl", 2, "iTst", 2, "GInt", GInt("D3F", 2)), SLF.empty);
            tc.verifyEqual(full(sum(M, "all")), sum(FcArea), "AbsTol", 1e-12);
            % int_F u * mu over all faces (from positive side) with u = 1 and mu = 1.
            B = assemble(msh, trls, trls, DLF(msh, 2, Fcn.cst(1), d0, [0; 0], "iTrl", 1, "iTst", 2, "GInt", GInt("D3F", 2)), SLF.empty);
            nDoF1 = tc.DG1_FES3.nGlDoF;
            tc.verifyEqual(full(sum(B(nDoF1 + 1:end, 1:nDoF1), "all")), sum(FcArea), "AbsTol", 1e-12);
            % Load on trace space.
            [~, F] = assemble(msh, trls, trls, DLF.empty, SLF(msh, 2, Fcn("D3", "x"), [0; 0], "iTst", 2, "GInt", GInt("D3F", 2)));
            FcCen = (FcNd(1) + FcNd(2) + FcNd(3)) / 3;
            tc.verifyEqual(sum(F(nDoF1 + 1:end)), sum(FcArea .* FcCen(1, :)), "AbsTol", 1e-12);
            % Error norm on trace space.
            u = Fcn("D3", "x + 2*y");
            tc.verifyEqual(eNorm(msh, u, TP1_FES.proj(u), Norm(msh, 2, [0; 0], "GInt", GInt("D3F", 2))), 0, "AbsTol", 1e-12);
        end
        function condSolveHDG(tc)
            % Static condensation of discontinuous spaces gives the same solution as the full system (2D HDG).
            msh = tc.msh;
            d0 = [0; 0]; d0V = [0, 0; 0, 0]; divV = [1, 0; 0, 1]; grad = cat(3, [1; 0], [0; 1]);
            UNV = MshEnt("D2L").UNV;
            p = Fcn("D2", "sin(pi*x)*y + x");
            VP1_FES = FES(msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan], "share", false), ...
                NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0], "share", false)]));
            TP1_FES = FES(msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0, "share", false)), ...
                BC(p, "edge", tc.BdEdge));
            trls = [VP1_FES, tc.DG1_FES, TP1_FES];
            Auv = [DLF(msh, 2, Fcn.cst(1), d0V, d0V, "iTrl", 1, "iTst", 1, "GInt", GInt("D2T", 2)), ...
                DLF(msh, 2, Fcn.cst(1), d0, divV, "iTrl", 2, "iTst", 1, "GInt", GInt("D2T", 1)), ...
                DLF.interface(msh, 1, -UNV, 0, d0V, "iTrl", 3, "iTst", 1, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
                DLF(msh, 2, Fcn.cst(1), d0V, grad, "iTrl", 1, "iTst", 2, "GInt", GInt("D2T", 1)), ...
                DLF(msh, 1, -UNV, d0V, d0, "iTrl", 1, "iTst", 2, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, UNV, d0V, d0, "iTrl", -1, "iTst", -2, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(1), d0, d0, "iTrl", 2, "iTst", 2, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(1), d0, d0, "iTrl", -2, "iTst", -2, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(-1), 0, d0, "iTrl", 3, "iTst", 2, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(-1), 0, d0, "iTrl", 3, "iTst", -2, "GInt", GInt("D2L", 2)), ...
                DLF.interface(msh, 1, UNV, d0V, 0, "iTrl", 1, "iTst", 3, "trlOpr", "jump", "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(-1), d0, 0, "iTrl", 2, "iTst", 3, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(-1), d0, 0, "iTrl", -2, "iTst", 3, "GInt", GInt("D2L", 2)), ...
                DLF(msh, 1, Fcn.cst(2), 0, 0, "iTrl", 3, "iTst", 3, "EntIdx", tc.IntEdge, "GInt", GInt("D2L", 2))];
            f = -dif(p, [2; 0]) - dif(p, [0; 2]);
            [Stiff, Load] = assemble(msh, trls, trls, Auv, SLF(msh, 2, f, d0, "iTst", 2, "GInt", GInt("D2T", 2)));
            sol = condSolve(Stiff, Load, trls, [1, 2]);
            tc.verifyEqual(sol, Stiff \ Load, "AbsTol", 1e-10);
            % Shared (continuous) DoFs cannot be eliminated element by element.
            tc.verifyError(@() condSolve(Stiff, Load, [tc.P1_FES, trls(2:3)], 1), "MATLAB:assertion:failed");
        end
        function patch3D(tc)
            msh = tc.msh3;
            grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
            BdNode = msh.bdEnt(0, 1:6);
            % Harmonic linear solution is reproduced exactly by P1.
            u = Fcn("D3", "x + 2*y - 3*z + 1");
            Uh = FES(msh, stdFE("P1"), BC(u, "node", BdNode));
            [K, F] = assemble(msh, Uh, Uh, DLF(msh, 3, Fcn.cst(1), grad, grad, "GInt", GInt("D3T", 0)), SLF.empty);
            tc.verifyEqual(K \ F, evalDoF(Uh, u), "AbsTol", 1e-12);
            % Harmonic quadratic solution is reproduced exactly by P2.
            u = Fcn("D3", "x^2 - y^2 + x*z + y*z");
            Uh = FES(msh, stdFE("P2"), BC(u, "node", BdNode, "edge", msh.bdEnt(1, 1:6)));
            [K, F] = assemble(msh, Uh, Uh, DLF(msh, 3, Fcn.cst(1), grad, grad, "GInt", GInt("D3T", 2)), SLF.empty);
            tc.verifyEqual(K \ F, evalDoF(Uh, u), "AbsTol", 1e-11);
        end
        function linearized(tc)
            % With previous solution w = 1, linearized forms reduce to ordinary forms.
            Uh = tc.P1_FES;
            d0 = [0; 0];
            w = FEF(Uh, ones(Uh.nGlDoF, 1));
            [M, F] = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2T", 2)), ...
                SLF(tc.msh, 2, Fcn("D2", "x"), d0, "GInt", GInt("D2T", 2)));
            [MLin, FLin] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF.empty, "preSol", w, ...
                "Awuv", LDLF(tc.msh, 2, Fcn.cst(1), d0, d0, d0, "GInt", GInt("D2T", 3)), ...
                "Fwv", LSLF(tc.msh, 2, Fcn("D2", "x"), d0, d0, "GInt", GInt("D2T", 3)));
            tc.verifyEqual(full(MLin), full(M), "AbsTol", 1e-13);
            tc.verifyEqual(FLin, F, "AbsTol", 1e-13);
        end
        function linearizedFacet(tc)
            % On facets, linearized forms with previous solution w = x reduce to ordinary forms with coefficient x.
            % Previous solution is taken from positive side (iPre = 1), negative side (iPre = -1) or trace space (iPre = 2).
            Uh = tc.DG1_FES;
            TP1_FES = FES(tc.msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0)));
            d0 = [0; 0];
            x = Fcn("D2", "x");
            ws = [tc.P1_FES.proj(x), TP1_FES.proj(x)];
            M = assemble(tc.msh, Uh, Uh, DLF(tc.msh, 1, x, d0, d0, "EntIdx", tc.IntEdge, "iTrl", 1, "iTst", -1, ...
                "GInt", GInt("D2L", 3)), SLF.empty);
            [~, F] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF(tc.msh, 1, x, d0, "EntIdx", tc.BdEdge, "GInt", GInt("D2L", 2)));
            iPres = [1, -1, 2];
            preOrds = {d0, d0, 0};
            for i = 1:length(iPres)
                MLin = assemble(tc.msh, Uh, Uh, DLF.empty, SLF.empty, "preSol", ws, "Awuv", LDLF(tc.msh, 1, Fcn.cst(1), ...
                    preOrds{i}, d0, d0, "EntIdx", tc.IntEdge, "iPre", iPres(i), "iTrl", 1, "iTst", -1, "GInt", GInt("D2L", 3)));
                tc.verifyEqual(full(MLin), full(M), "AbsTol", 1e-13, sprintf("iPre = %d", iPres(i)));
            end
            % On boundary edges, the previous solution exists on the positive side and in the trace space only.
            [~, FLin] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF.empty, "preSol", ws, ...
                "Fwv", LSLF(tc.msh, 1, Fcn.cst(1), d0, d0, "EntIdx", tc.BdEdge, "iPre", 1, "GInt", GInt("D2L", 2)));
            tc.verifyEqual(FLin, F, "AbsTol", 1e-13);
            [~, FLin] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF.empty, "preSol", ws, ...
                "Fwv", LSLF(tc.msh, 1, Fcn.cst(1), 0, d0, "EntIdx", tc.BdEdge, "iPre", 2, "GInt", GInt("D2L", 2)));
            tc.verifyEqual(FLin, F, "AbsTol", 1e-13);
            [~, FLin] = assemble(tc.msh, Uh, Uh, DLF.empty, SLF.empty, "preSol", ws, ...
                "Fwv", LSLF(tc.msh, 1, Fcn.cst(1), d0, d0, "EntIdx", tc.BdEdge, "iPre", -1, "GInt", GInt("D2L", 2)));
            tc.verifyEqual(FLin, zeros(Uh.nGlDoF, 1));
        end
    end
end
%% Local functions.
function DoFVal = evalDoF(fES, fcn)
    % evalDoF: global DoF values of function (nodal DoF in FE space).
    DoFVal = zeros(fES.nGlDoF, 1);
    for iDoF = 1:length(fES.GlDoFs)
        val = fES.GlDoFs(iDoF).eval(fcn, "valType", "num");
        DoFVal(fES.GlDoFs.sub2ind(iDoF)) = val(:);
    end
end
