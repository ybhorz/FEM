classdef MshEntTest < matlab.unittest.TestCase
    % MshEntTest: unit tests of MshEnt.
    properties
        P2T = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1]; % Vertices of a counter-clockwise triangle (by column).
        P2L = [0, 3; 0, 4]; % Vertices of a line (by column).
        P3T = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3]; % Vertices of a positively oriented tetrahedron.
        P3F = [0.1, 1.2, 0.3; -0.1, 0.2, 1.1; 0.05, 0.1, 0.9]; % Vertices of a triangle in 3D.
    end
    methods (Test, TestTags = {'Fast'})
        function varParm(tc)
            tc.verifyTrue(symEq(MshEnt("D2").var, str2sym("[x; y]")));
            tc.verifyTrue(symEq(MshEnt("D2R").var, str2sym("[l; m]")));
            tc.verifyTrue(symEq(MshEnt("D2R1").var, str2sym("s")));
            tc.verifyTrue(symEq(MshEnt("D2T").var, str2sym("[x; y]")));
            tc.verifyTrue(symEq(MshEnt("D2TR").var, str2sym("[l; m]")));
            tc.verifyTrue(symEq(MshEnt("D2L").var, str2sym("[x; y]")));
            tc.verifyTrue(symEq(MshEnt("D2LR").var, str2sym("s")));
            tc.verifyEmpty(MshEnt("D2").parm);
            tc.verifyEmpty(MshEnt("D2R").parm);
            tc.verifyEmpty(MshEnt("D2R1").parm);
            tc.verifyTrue(symEq(MshEnt("D2T").parm, str2sym("[x1, x2, x3; y1, y2, y3]")));
            tc.verifyTrue(symEq(MshEnt("D2TR").parm, str2sym("[x1, x2, x3; y1, y2, y3]")));
            tc.verifyTrue(symEq(MshEnt("D2L").parm, str2sym("[x1, x2; y1, y2]")));
            tc.verifyTrue(symEq(MshEnt("D2LR").parm, str2sym("[x1, x2; y1, y2]")));
            tc.verifyEqual(MshEnt("D2T").nVar, 2);
            tc.verifyEqual(MshEnt("D2LR").nVar, 1);
            tc.verifyEqual(MshEnt("D2T").sParm, [2, 3]);
            tc.verifyEqual(MshEnt("D2LR").sParm, [2, 2]);
        end
        function dimension(tc)
            types = ["D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
                "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"];
            dims = [2, 2, 1, 2, 2, 2, 1, 3, 3, 2, 1, 3, 3, 3, 2, 3, 1];
            for iType = 1:length(types)
                tc.verifyEqual(MshEnt(types(iType)).dim, dims(iType), types(iType));
            end
            tc.verifyEmpty(MshEnt("VOID").dim);
        end
        function dualType(tc)
            types = ["D2", "D2R", "D2T", "D2TR", "D2L", "D2LR", "D2R1", "VOID", ...
                "D3", "D3R", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR", "D3R2", "D3R1"];
            dlTypes = ["D2R", "D2", "D2TR", "D2T", "D2LR", "D2L", "VOID", "VOID", ...
                "D3R", "D3", "D3TR", "D3T", "D3FR", "D3F", "D3LR", "D3L", "VOID", "VOID"];
            for iType = 1:length(types)
                tc.verifyEqual(MshEnt(types(iType)).dual.type, dlTypes(iType), types(iType));
            end
        end
        function ismemberTable(tc)
            % Row: `mshEnt1`, column: `mshEnt2`, entry: ismember(mshEnt1, mshEnt2).
            types = ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR"];
            tab = logical([1, 1, 1, 1, 1, 1, 1, 1; ... % VOID
                           0, 1, 0, 0, 1, 0, 1, 0; ... % D2
                           0, 0, 1, 0, 0, 1, 0, 0; ... % D2R
                           0, 0, 0, 1, 0, 0, 0, 1; ... % D2R1
                           0, 0, 0, 0, 1, 0, 0, 0; ... % D2T
                           0, 0, 0, 0, 0, 1, 0, 0; ... % D2TR
                           0, 0, 0, 0, 1, 0, 1, 0; ... % D2L
                           0, 0, 0, 0, 0, 0, 0, 1]);   % D2LR
            for i = 1:length(types)
                for j = 1:length(types)
                    tc.verifyEqual(ismember(MshEnt(types(i)), MshEnt(types(j))), tab(i, j), ...
                        "ismember(" + types(i) + ", " + types(j) + ")");
                end
            end
        end
        function varParm3D(tc)
            tc.verifyTrue(symEq(MshEnt("D3").var, str2sym("[x; y; z]")));
            tc.verifyTrue(symEq(MshEnt("D3R").var, str2sym("[l; m; n]")));
            tc.verifyTrue(symEq(MshEnt("D3R2").var, str2sym("[s; t]")));
            tc.verifyTrue(symEq(MshEnt("D3R1").var, str2sym("s")));
            tc.verifyTrue(symEq(MshEnt("D3T").var, str2sym("[x; y; z]")));
            tc.verifyTrue(symEq(MshEnt("D3TR").var, str2sym("[l; m; n]")));
            tc.verifyTrue(symEq(MshEnt("D3F").var, str2sym("[x; y; z]")));
            tc.verifyTrue(symEq(MshEnt("D3FR").var, str2sym("[s; t]")));
            tc.verifyTrue(symEq(MshEnt("D3L").var, str2sym("[x; y; z]")));
            tc.verifyTrue(symEq(MshEnt("D3LR").var, str2sym("s")));
            tc.verifyEmpty(MshEnt("D3").parm);
            tc.verifyEmpty(MshEnt("D3R2").parm);
            tc.verifyTrue(symEq(MshEnt("D3T").parm, str2sym("[x1, x2, x3, x4; y1, y2, y3, y4; z1, z2, z3, z4]")));
            tc.verifyTrue(symEq(MshEnt("D3TR").parm, MshEnt("D3T").parm));
            tc.verifyTrue(symEq(MshEnt("D3F").parm, str2sym("[x1, x2, x3; y1, y2, y3; z1, z2, z3]")));
            tc.verifyTrue(symEq(MshEnt("D3FR").parm, MshEnt("D3F").parm));
            tc.verifyTrue(symEq(MshEnt("D3L").parm, str2sym("[x1, x2; y1, y2; z1, z2]")));
            tc.verifyTrue(symEq(MshEnt("D3LR").parm, MshEnt("D3L").parm));
            tc.verifyEqual(MshEnt("D3T").sParm, [3, 4]);
            tc.verifyEqual(MshEnt("D3F").sParm, [3, 3]);
            tc.verifyEqual(MshEnt("D3L").sParm, [3, 2]);
        end
        function ismemberTable3D(tc)
            % Row: `mshEnt1`, column: `mshEnt2`, entry: ismember(mshEnt1, mshEnt2).
            types = ["VOID", "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"];
            tab = logical([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; ... % VOID
                           0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0; ... % D3
                           0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0; ... % D3R
                           0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0; ... % D3R2
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1; ... % D3R1
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; ... % D3T
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0; ... % D3TR
                           0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0; ... % D3F
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0; ... % D3FR
                           0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0; ... % D3L
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);   % D3LR
            for i = 1:length(types)
                for j = 1:length(types)
                    tc.verifyEqual(ismember(MshEnt(types(i)), MshEnt(types(j))), tab(i, j), ...
                        "ismember(" + types(i) + ", " + types(j) + ")");
                end
            end
            % 2D and 3D types are not compatible.
            types2 = ["D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR"];
            for i = 1:length(types2)
                for j = 2:length(types)
                    tc.verifyFalse(ismember(MshEnt(types2(i)), MshEnt(types(j))), types2(i) + ", " + types(j));
                    tc.verifyFalse(ismember(MshEnt(types(j)), MshEnt(types2(i))), types(j) + ", " + types2(i));
                end
            end
        end
        function typeInfo(tc)
            % Consistency of data of all mesh entity types.
            types = MshEnt.typeLst;
            tc.verifyEqual(length(unique(types)), length(types));
            tc.verifyEqual(types(1), "VOID");
            for iType = 1:length(types)
                entType = types(iType);
                info = MshEnt.getInfo(entType);
                % Type is member of itself.
                tc.verifyTrue(ismember(MshEnt(entType), MshEnt(entType)), entType);
                % Dual is an involution.
                if info.dual ~= "VOID"
                    tc.verifyEqual(MshEnt(info.dual).dual.type, entType, entType);
                end
                % Type after clearing parameter has no parameter, same variable, and is member of the type.
                tc.verifyEmpty(MshEnt(info.free).parm, entType);
                tc.verifyTrue(symEq(MshEnt(info.free).var, info.var), entType);
                tc.verifyTrue(ismember(MshEnt(info.free), MshEnt(entType)), entType);
                % Dimension equals number of variables.
                if entType ~= "VOID"
                    tc.verifyEqual(info.dim, length(info.var), entType);
                end
            end
            % Literal lists of validators in MshEnt and Fcn accept all types in `typeLst`.
            for iType = 1:length(types)
                tc.verifyEqual(MshEnt(types(iType)).type, types(iType));
                tc.verifyEqual(Fcn(types(iType), "1").domn, types(iType));
            end
            tc.verifyError(@() MshEnt("D9"), "MATLAB:validators:mustBeMember");
            tc.verifyError(@() Fcn("D9", "x"), "MATLAB:validators:mustBeMember");
        end
        function topology(tc)
            mshEnt = MshEnt("D2T");
            tc.verifyTrue(symEq(mshEnt.node.coord, mshEnt.parm));
            tc.verifyEqual(mshEnt.elem.node, [1; 2; 3]);
            tc.verifyEqual(mshEnt.edge.node, [1, 2, 3; 2, 3, 1]);
            tc.verifyEqual(mshEnt.nNode, 3);
            tc.verifyEqual(mshEnt.nEdge, 3);
            tc.verifyEqual(mshEnt.msh.type, "D2T");
            tc.verifyTrue(symEq(MshEnt("D2TR").node.coord, sym([0, 1, 0; 0, 0, 1])));
            tc.verifyTrue(symEq(MshEnt("D2L").node.coord, MshEnt("D2L").parm));
            tc.verifyTrue(symEq(MshEnt("D2LR").node.coord, sym([0, 1])));
            tc.verifyEqual(MshEnt("D2LR").msh.type, "D1");
            tc.verifyEqual(MshEnt("D2LR").elem.node, [1; 2]);
        end
        function topology3D(tc)
            mshEnt = MshEnt("D3T");
            tc.verifyTrue(symEq(mshEnt.node.coord, mshEnt.parm));
            tc.verifyEqual(mshEnt.elem.node, [1; 2; 3; 4]);
            tc.verifyEqual(mshEnt.edge.node, [1, 1, 1, 2, 2, 3; 2, 3, 4, 3, 4, 4]);
            tc.verifyEqual(mshEnt.face.node, [2, 1, 1, 1; 3, 4, 2, 3; 4, 3, 4, 2]);
            tc.verifyEqual([mshEnt.nNode, mshEnt.nEdge, mshEnt.nFace], [4, 6, 4]);
            % Face i is opposite to vertex i.
            for iFace = 1:4
                tc.verifyFalse(ismember(iFace, mshEnt.face.node(:, iFace)));
            end
            % Each edge belongs to exactly two faces.
            FcEgNode = [mshEnt.face.node([1, 2], :), mshEnt.face.node([2, 3], :), mshEnt.face.node([3, 1], :)];
            [isEg, iEg] = ismember(sort(FcEgNode, 1)', mshEnt.edge.node', "rows");
            tc.verifyTrue(all(isEg));
            tc.verifyEqual(accumarray(iEg, 1)', 2 * ones(1, 6));
            tc.verifyTrue(symEq(MshEnt("D3TR").node.coord, sym([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1])));
            tc.verifyEqual(MshEnt("D3TR").face.node, mshEnt.face.node);
            % Face and line.
            tc.verifyEqual([MshEnt("D3F").nNode, MshEnt("D3F").nEdge], [3, 3]);
            tc.verifyEqual(MshEnt("D3F").face.node, [1; 2; 3]);
            tc.verifyTrue(symEq(MshEnt("D3FR").node.coord, sym([0, 1, 0; 0, 0, 1])));
            tc.verifyEqual(MshEnt("D3FR").msh.type, "D2T");
            tc.verifyEqual(MshEnt("D3FR").msh.nElem, 1);
            tc.verifyEqual(MshEnt("D3L").edge.node, [1; 2]);
            tc.verifyTrue(symEq(MshEnt("D3LR").node.coord, sym([0, 1])));
            tc.verifyEqual(MshEnt("D3LR").msh.type, "D1");
        end
        function geometryD3T(tc)
            P = tc.P3T;
            mshEnt = MshEnt("D3T");
            UNV = mshEnt.UNV; area = mshEnt.area; UTV = mshEnt.UTV; len = mshEnt.len;
            tc.verifyEqual([length(UNV), length(area), length(UTV), length(len)], [4, 4, 6, 6]);
            cen = mean(P, 2);
            sumNor = zeros(3, 1);
            for iFace = 1:4
                FcNd = P(:, mshEnt.face.node(:, iFace));
                nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
                valUNV = numSym(UNV(iFace).fun, "D3T", P);
                valArea = numSym(area(iFace).fun, "D3T", P);
                tc.verifyEqual(UNV(iFace).domn, "D3T");
                tc.verifyEqual(valUNV, nor / norm(nor), "AbsTol", 1e-12);
                tc.verifyEqual(valArea, norm(nor) / 2, "AbsTol", 1e-12);
                % Unit normal vector is outward.
                tc.verifyGreaterThan(dot(valUNV, mean(FcNd, 2) - cen), 0);
                sumNor = sumNor + valArea * valUNV;
            end
            % Closed surface: sum of area-weighted normals vanishes.
            tc.verifyEqual(sumNor, zeros(3, 1), "AbsTol", 1e-12);
            for iEdge = 1:6
                tan = P(:, mshEnt.edge.node(2, iEdge)) - P(:, mshEnt.edge.node(1, iEdge));
                tc.verifyEqual(numSym(UTV(iEdge).fun, "D3T", P), tan / norm(tan), "AbsTol", 1e-12);
                tc.verifyEqual(numSym(len(iEdge).fun, "D3T", P), norm(tan), "AbsTol", 1e-12);
            end
        end
        function geometryD3F(tc)
            P = tc.P3F;
            mshEnt = MshEnt("D3F");
            nor = cross(P(:, 2) - P(:, 1), P(:, 3) - P(:, 1));
            tc.verifyEqual(numSym(mshEnt.UNV.fun, "D3F", P), nor / norm(nor), "AbsTol", 1e-12);
            tc.verifyEqual(numSym(mshEnt.area.fun, "D3F", P), norm(nor) / 2, "AbsTol", 1e-12);
            tc.verifyEqual(length(mshEnt.UTV), 3);
            tan = P(:, 1) - P(:, 3);
            tc.verifyEqual(numSym(mshEnt.len(3).fun, "D3F", P), norm(tan), "AbsTol", 1e-12);
            % Counter-clockwise triangle in plane z = 0: normal is e_z, area equals that of D2T.
            P2 = tc.P2T;
            areaD2T = ((P2(1, 2) - P2(1, 1)) * (P2(2, 3) - P2(2, 1)) - (P2(1, 3) - P2(1, 1)) * (P2(2, 2) - P2(2, 1))) / 2;
            tc.verifyEqual(numSym(mshEnt.UNV.fun, "D3F", [P2; 0, 0, 0]), [0; 0; 1], "AbsTol", 1e-12);
            tc.verifyEqual(numSym(mshEnt.area.fun, "D3F", [P2; 0, 0, 0]), areaD2T, "AbsTol", 1e-12);
        end
        function geometryD3L(tc)
            P = [0, 1; 0, 2; 0, 2];
            mshEnt = MshEnt("D3L");
            tc.verifyEqual(numSym(mshEnt.UTV.fun, "D3L", P), [1; 2; 2] / 3, "AbsTol", 1e-12);
            tc.verifyEqual(numSym(mshEnt.len.fun, "D3L", P), 3, "AbsTol", 1e-12);
            tc.verifyEmpty(mshEnt.UNV);
            tc.verifyEmpty(mshEnt.area);
        end
        function geometryD2T(tc)
            P = tc.P2T;
            mshEnt = MshEnt("D2T");
            UNV = mshEnt.UNV; UTV = mshEnt.UTV; len = mshEnt.len;
            cen = mean(P, 2);
            for iEdge = 1:3
                EgNd1 = P(:, mshEnt.edge.node(1, iEdge));
                EgNd2 = P(:, mshEnt.edge.node(2, iEdge));
                tan = EgNd2 - EgNd1;
                tc.verifyEqual(UNV(iEdge).domn, "D2T");
                tc.verifyEqual(numSym(UNV(iEdge).fun, "D2T", P), [tan(2); -tan(1)] / norm(tan), "AbsTol", 1e-12);
                tc.verifyEqual(numSym(UTV(iEdge).fun, "D2T", P), tan / norm(tan), "AbsTol", 1e-12);
                tc.verifyEqual(numSym(len(iEdge).fun, "D2T", P), norm(tan), "AbsTol", 1e-12);
                % Unit normal vector is outward.
                tc.verifyGreaterThan(dot(numSym(UNV(iEdge).fun, "D2T", P), (EgNd1 + EgNd2) / 2 - cen), 0);
            end
        end
        function geometryD2L(tc)
            P = tc.P2L;
            mshEnt = MshEnt("D2L");
            tc.verifyEqual(numSym(mshEnt.UNV.fun, "D2L", P), [4; -3] / 5, "AbsTol", 1e-12);
            tc.verifyEqual(numSym(mshEnt.UTV.fun, "D2L", P), [3; 4] / 5, "AbsTol", 1e-12);
            tc.verifyEqual(numSym(mshEnt.len.fun, "D2L", P), 5, "AbsTol", 1e-12);
        end
    end
end
