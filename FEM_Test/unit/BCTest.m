classdef BCTest < matlab.unittest.TestCase
    % BCTest: unit tests of BC.
    % Mesh: mshD2TS([0, 1, 0, 1], 2), 9 nodes (8 on boundary), 16 edges (8 on boundary), 8 elements.
    % 3D mesh: mshD3TS([0, 1, 0, 1, 0, 1], 2), 27 nodes (26 on boundary), 98 edges, 120 faces (48 on boundary), 48 elements.
    properties
        msh; % Mesh.
        BdNode; % Indices of boundary nodes.
        BdEdge; % Indices of boundary edges.
        msh3; % 3D mesh.
    end
    methods (TestClassSetup)
        function setMsh(tc)
            tc.msh = mshD2TS([0, 1, 0, 1], 2);
            tc.BdNode = find(tc.msh.node.type ~= 0);
            tc.BdEdge = find(tc.msh.edge.type ~= 0);
            tc.msh3 = mshD3TS([0, 1, 0, 1, 0, 1], 2);
        end
    end
    methods (Test, TestTags = {'Fast'})
        function P1(tc)
            g = Fcn("D2", "x^2 + y^2");
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0])), BC(g, "node", tc.BdNode));
            tc.verifyEqual(fES.BC.DoFIdx, tc.BdNode);
            tc.verifyEqual(fES.BC.DoFVal, evalNode(g, tc.msh.node.coord(:, tc.BdNode)), "AbsTol", 1e-14);
            tc.verifyEqual(fES.BC.nDoF, 8);
        end
        function P2(tc)
            g = Fcn("D2", "x^3 + y^3");
            fES = FES(tc.msh, FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
                NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]), BC(g, "node", tc.BdNode, "edge", tc.BdEdge));
            tc.verifyEqual(fES.BC.DoFIdx, [tc.BdNode, 9 + tc.BdEdge]);
            tc.verifyEqual(fES.BC.DoFVal, [evalNode(g, tc.msh.node.coord(:, tc.BdNode)), evalNode(g, midPnt(tc.msh, tc.BdEdge))], ...
                "AbsTol", 1e-14);
        end
        function DG1(tc)
            g = Fcn("D2", "x^2 + y^2");
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false)), ...
                BC(g, "node", tc.BdNode));
            % Each DoF on a boundary node is constrained, in every element containing the node.
            GlNode = fES.GlDoFs.EntIdx;
            tc.verifyEqual(fES.BC.DoFIdx, find(ismember(GlNode, tc.BdNode)));
            tc.verifyEqual(fES.BC.DoFVal, evalNode(g, tc.msh.node.coord(:, GlNode(fES.BC.DoFIdx))), "AbsTol", 1e-14);
        end
        function CR(tc)
            g = Fcn("D2", "x^2 + y^2");
            fES = FES(tc.msh, FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])), BC(g, "edge", tc.BdEdge));
            tc.verifyEqual(fES.BC.DoFIdx, tc.BdEdge);
            tc.verifyEqual(fES.BC.DoFVal, evalNode(g, midPnt(tc.msh, tc.BdEdge)), "AbsTol", 1e-14);
        end
        function VP1(tc)
            g = [Fcn("D2", "x^2 + y^2"), Fcn("D2", "x^2 - y^2")];
            fES = FES(tc.msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
                [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
                NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]), BC(g, "node", tc.BdNode));
            tc.verifyEqual(fES.BC.DoFIdx, [tc.BdNode, 9 + tc.BdNode]);
            BdCrd = tc.msh.node.coord(:, tc.BdNode);
            tc.verifyEqual(fES.BC.DoFVal, [evalNode(g(1), BdCrd), evalNode(g(2), BdCrd)], "AbsTol", 1e-14);
        end
        function RT0(tc)
            g = dot(Fcn("D2", "[y^2; x^2]"), MshEnt("D2L").UNV);
            fES = FES(tc.msh, FE("D2T", "[1,0; 0,1; x,y].'", ...
                MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, ...
                "GInt", GInt("D2L", 4))), BC(g, "edge", tc.BdEdge));
            tc.verifyEqual(fES.BC.DoFIdx, tc.BdEdge);
            % Flux through edge w.r.t. normal of global edge, by 3-point Gauss rule (exact for quadratics).
            gPnt = [(1 - sqrt(3/5)) / 2, 1/2, (1 + sqrt(3/5)) / 2];
            gWgt = [5/18, 4/9, 5/18];
            flux = zeros(1, length(tc.BdEdge));
            for iBd = 1:length(tc.BdEdge)
                EgNd1 = tc.msh.node.coord(:, tc.msh.edge.node(1, tc.BdEdge(iBd)));
                EgNd2 = tc.msh.node.coord(:, tc.msh.edge.node(2, tc.BdEdge(iBd)));
                tan = EgNd2 - EgNd1;
                nor = [tan(2); -tan(1)] / norm(tan);
                pnt = EgNd1 + tan * gPnt;
                flux(iBd) = sum((pnt(2, :).^2 * nor(1) + pnt(1, :).^2 * nor(2)) .* gWgt) * norm(tan);
            end
            tc.verifyEqual(fES.BC.DoFVal, flux, "AbsTol", 1e-12);
        end
        function P1_3D(tc)
            msh = tc.msh3;
            g = Fcn("D3", "x^2 + y*z");
            fE = stdFE("P1");
            BdNode = msh.bdEnt(0, 1:6);
            tc.verifyEqual(length(BdNode), 26);
            fES = FES(msh, fE, BC(g, "node", BdNode));
            tc.verifyEqual(fES.BC.DoFIdx, BdNode);
            tc.verifyEqual(fES.BC.DoFVal, evalNode(g, msh.node.coord(:, BdNode)), "AbsTol", 1e-14);
            % Dirichlet condition on faces x = 0, y = 0, z = 0 only.
            DirNode = msh.bdEnt(0, [1, 3, 5]);
            tc.verifyEqual(length(DirNode), 19);
            fES = FES(msh, fE, BC(g, "node", DirNode));
            tc.verifyEqual(fES.BC.DoFIdx, DirNode);
            tc.verifyTrue(all(min(msh.node.coord(:, DirNode), [], 1) == 0));
        end
        function P2_3D(tc)
            msh = tc.msh3;
            g = Fcn("D3", "x^3 + y*z");
            fE = stdFE("P2");
            BdNode = msh.bdEnt(0, 1:6);
            BdEdge = msh.bdEnt(1, 1:6);
            fES = FES(msh, fE, BC(g, "node", BdNode, "edge", BdEdge));
            tc.verifyEqual(fES.BC.DoFIdx, [BdNode, 27 + BdEdge]);
            tc.verifyEqual(fES.BC.DoFVal, [evalNode(g, msh.node.coord(:, BdNode)), evalNode(g, midPnt(msh, BdEdge))], ...
                "AbsTol", 1e-14);
            % Global DoFs on boundary (125 - 27 interior = 98).
            tc.verifyEqual(fES.BC.nDoF, 125 - 27);
        end
        function CR_3D(tc)
            msh = tc.msh3;
            g = Fcn("D3", "x^2 + y*z");
            fE = stdFE("CR");
            BdFace = msh.bdEnt(2, 1:6);
            tc.verifyEqual(length(BdFace), 48);
            fES = FES(msh, fE, BC(g, "face", BdFace));
            tc.verifyEqual(fES.BC.DoFIdx, BdFace);
            FcCen = (msh.node.coord(:, msh.face.node(1, BdFace)) + msh.node.coord(:, msh.face.node(2, BdFace)) + ...
                msh.node.coord(:, msh.face.node(3, BdFace))) / 3;
            tc.verifyEqual(fES.BC.DoFVal, evalNode(g, FcCen), "AbsTol", 1e-14);
            % Face DoFs require boundary faces.
            tc.verifyError(@() FES(msh, fE, BC(g, "node", msh.bdEnt(0, 1:6))), "MATLAB:assertion:failed");
        end
        function SingleDoF(tc)
            g = Fcn("D2", "x^3 + y^3");
            fE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
                NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]);
            % Only single DoFs (one node and one edge midpoint), no entity asserts.
            fES = FES(tc.msh, fE, BC(g, "DoF", [3, 9 + 5]));
            tc.verifyEqual(fES.BC.DoFIdx, [3, 9 + 5]);
            tc.verifyEqual(fES.BC.DoFVal, [evalNode(g, tc.msh.node.coord(:, 3)), evalNode(g, midPnt(tc.msh, 5))], "AbsTol", 1e-14);
            % Together with boundary entities: the interior node is added.
            InNode = find(tc.msh.node.type == 0);
            fES = FES(tc.msh, fE, BC(g, "node", tc.BdNode, "edge", tc.BdEdge, "DoF", InNode));
            tc.verifyEqual(fES.BC.DoFIdx, sort([tc.BdNode, InNode, 9 + tc.BdEdge]));
            tc.verifyEqual(fES.BC.nDoF, 8 + 1 + 8);
        end
        function TP1(tc)
            g = Fcn("D2", "x^2 + y^2");
            fES = FES(tc.msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0)), BC(g, "edge", tc.BdEdge));
            tc.verifyEqual(fES.BC.DoFIdx, [tc.BdEdge, 16 + tc.BdEdge]);
            tc.verifyEqual(fES.BC.DoFVal, [evalNode(g, tc.msh.node.coord(:, tc.msh.edge.node(1, tc.BdEdge))), ...
                evalNode(g, tc.msh.node.coord(:, tc.msh.edge.node(2, tc.BdEdge)))], "AbsTol", 1e-14);
        end
    end
end
%% Local functions.
function val = evalNode(fcn, crd)
    % evalNode: evaluate scalar function at points (by column).
    fun = fcn.getFun;
    val = fun(crd);
end
function crd = midPnt(msh, EgIdx)
    % midPnt: midpoints of edges (by column).
    crd = (msh.node.coord(:, msh.edge.node(1, EgIdx)) + msh.node.coord(:, msh.edge.node(2, EgIdx))) / 2;
end
