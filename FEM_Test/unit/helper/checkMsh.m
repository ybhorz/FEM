function [flag, msg] = checkMsh(msh, options)
    % checkMsh: check consistency of mesh data.
    % flag = true if all checks pass, otherwise `msg` describes the first failed check.

    % Check                                                           | D2T | D3T
    % --------------------------------------------------------------- | --- | ---
    % elements are positively oriented (counter-clockwise in 2D)      | Y   | Y
    % elem.edge matches edge.node (with orientation sign)             | Y   | Y
    % edge.elem matches elem.edge                                     | Y   |
    % elem.face, elem.facePerm match face.node (orientation sign)    |     | Y
    % face.elem matches elem.face                                     |     | Y
    % face.edge matches edge.node (with orientation sign)             |     | Y
    % no duplicated edge (face)                                       | Y   | Y
    % facet has one (boundary) or two (interior) elements             | Y   | Y
    % interior facet: one positive and one negative element           | Y   | Y
    % boundary facet: positive element (option `bdOrien`)             | Y   | Y
    % boundary facet: normal is outward (geometric check)             |     | Y
    % boundary facet <=> real(facet.type) > 0 (option `bdType`)       | Y   | Y
    % Euler characteristic V - E + F (- T) = 1 (option `euler`)       | Y   | Y
    % total area or volume (option `vol`)                             | Y   | Y

    arguments (Input)
        msh Msh; % Mesh.
        options.bdOrien logical = true; % Whether check boundary facet is positive in boundary element.
        options.bdType logical = true; % Whether check boundary facet is labeled by positive type.
        options.euler logical = true; % Whether check Euler characteristic (simply connected domain).
        options.vol = []; % Expected total area (volume) of domain.
        options.tol = 1e-12; % Tolerance of area (volume).
    end
    arguments (Output)
        flag (1, 1) logical; % Whether all checks pass.
        msg (1, 1) string; % Message of the first failed check.
    end
    flag = true; msg = "";
    switch msh.type
        case "D2T"
            ElNode = msh.elem.node; ElEdge = msh.elem.edge;
            EgNode = msh.edge.node; EgElem = msh.edge.elem;
            nElem = msh.nElem; nEdge = msh.nEdge;
            % Size.
            if ~isequal(size(ElNode), [3, nElem]) || ~isequal(size(ElEdge), [3, nElem]) || ...
                    ~isequal(size(EgNode), [2, nEdge]) || ~isequal(size(EgElem), [2, nEdge])
                [flag, msg] = setFail("Invalid size of elem.node, elem.edge, edge.node or edge.elem.");
                return;
            end
            % Orientation of element.
            x = msh.node.coord(1, :); y = msh.node.coord(2, :);
            area = ((x(ElNode(2, :)) - x(ElNode(1, :))) .* (y(ElNode(3, :)) - y(ElNode(1, :))) ...
                - (x(ElNode(3, :)) - x(ElNode(1, :))) .* (y(ElNode(2, :)) - y(ElNode(1, :)))) / 2;
            if any(area <= 0)
                [flag, msg] = setFail("Elements are not counter-clockwise: " + idxStr(find(area <= 0)));
                return;
            end
            % Element - edge.
            if any(ElEdge == 0, "all") || any(abs(ElEdge) > nEdge, "all")
                [flag, msg] = setFail("Invalid index in elem.edge.");
                return;
            end
            EgIdx = abs(ElEdge); sgn = sign(ElEdge);
            ElEgNd1 = ElNode; ElEgNd2 = ElNode([2, 3, 1], :);
            EgNd1 = reshape(EgNode(1, EgIdx), size(EgIdx)); EgNd2 = reshape(EgNode(2, EgIdx), size(EgIdx));
            isMatch = (sgn > 0 & EgNd1 == ElEgNd1 & EgNd2 == ElEgNd2) | (sgn < 0 & EgNd1 == ElEgNd2 & EgNd2 == ElEgNd1);
            if ~all(isMatch, "all")
                [~, iElem] = find(~isMatch);
                [flag, msg] = setFail("elem.edge does not match edge.node in elements: " + idxStr(unique(iElem)));
                return;
            end
            % Duplicated edge.
            if size(unique(sort(EgNode, 1)', "rows"), 1) < nEdge
                [flag, msg] = setFail("Duplicated edges.");
                return;
            end
            % Edge - element.
            ElIdx = repmat(1:nElem, 3, 1);
            pairElEg = sortrows([EgIdx(:), sgn(:) .* ElIdx(:)]);
            [iEg, ~, val] = find(EgElem');
            pairEgEl = sortrows([iEg(:), val(:)]);
            if ~isequal(pairElEg, pairEgEl)
                [flag, msg] = setFail("edge.elem does not match elem.edge.");
                return;
            end
            % Connected elements: boundary edge has one, interior edge has two (one positive, one negative).
            % Position of zero entry in edge.elem is arbitrary.
            nCnElem = sum(EgElem ~= 0, 1);
            if any(nCnElem == 0)
                [flag, msg] = setFail("Edges without connected element: " + idxStr(find(nCnElem == 0)));
                return;
            end
            isBd = nCnElem == 1;
            isInvalid = ~isBd & sign(EgElem(1, :)) .* sign(EgElem(2, :)) ~= -1;
            if any(isInvalid)
                [flag, msg] = setFail("Interior edges without one positive and one negative element: " + idxStr(find(isInvalid)));
                return;
            end
            isInvalid = isBd & sum(EgElem, 1) < 0;
            if options.bdOrien && any(isInvalid)
                [flag, msg] = setFail("Boundary edges are not counter-clockwise in boundary element: " + idxStr(find(isInvalid)));
                return;
            end
            if options.bdType && ~isempty(msh.edge.type)
                isInvalid = isBd ~= (real(msh.edge.type) > 0);
                if any(isInvalid)
                    [flag, msg] = setFail("Boundary edges and positive edge types do not match: " + idxStr(find(isInvalid)));
                    return;
                end
            end
            % Euler characteristic.
            if options.euler && msh.nNode - nEdge + nElem ~= 1
                [flag, msg] = setFail(sprintf("Euler characteristic V - E + F = %d.", msh.nNode - nEdge + nElem));
                return;
            end
            % Total area.
            if ~isempty(options.vol) && abs(sum(area) - options.vol) > options.tol
                [flag, msg] = setFail(sprintf("Total area %.15g is not %.15g.", sum(area), options.vol));
                return;
            end
        case "D3T"
            ElNode = msh.elem.node; ElEdge = msh.elem.edge; ElFace = msh.elem.face; ElFcPerm = msh.elem.facePerm;
            EgNode = msh.edge.node;
            FcNode = msh.face.node; FcElem = msh.face.elem; FcEdge = msh.face.edge;
            nElem = msh.nElem; nEdge = msh.nEdge; nFace = msh.nFace;
            % Size.
            if ~isequal(size(ElNode), [4, nElem]) || ~isequal(size(ElEdge), [6, nElem]) || ~isequal(size(ElFace), [4, nElem]) || ...
                    ~isequal(size(ElFcPerm), [4, nElem]) || ~isequal(size(EgNode), [2, nEdge]) || ...
                    ~isequal(size(FcNode), [3, nFace]) || ~isequal(size(FcElem), [2, nFace]) || ~isequal(size(FcEdge), [3, nFace])
                [flag, msg] = setFail("Invalid size of elem.node/edge/face/facePerm, edge.node or face.node/elem/edge.");
                return;
            end
            % Orientation of element.
            X = msh.node.coord;
            v1 = X(:, ElNode(1, :));
            vol = dot(X(:, ElNode(2, :)) - v1, cross(X(:, ElNode(3, :)) - v1, X(:, ElNode(4, :)) - v1)) / 6;
            if any(vol <= 0)
                [flag, msg] = setFail("Elements are not positively oriented: " + idxStr(find(vol <= 0)));
                return;
            end
            % Element - edge.
            if any(ElEdge == 0, "all") || any(abs(ElEdge) > nEdge, "all")
                [flag, msg] = setFail("Invalid index in elem.edge.");
                return;
            end
            LcEgNode = MshEnt("D3T").edge.node;
            EgIdx = abs(ElEdge); sgn = sign(ElEdge);
            ElEgNd1 = ElNode(LcEgNode(1, :), :); ElEgNd2 = ElNode(LcEgNode(2, :), :);
            EgNd1 = reshape(EgNode(1, EgIdx), size(EgIdx)); EgNd2 = reshape(EgNode(2, EgIdx), size(EgIdx));
            isMatch = (sgn > 0 & EgNd1 == ElEgNd1 & EgNd2 == ElEgNd2) | (sgn < 0 & EgNd1 == ElEgNd2 & EgNd2 == ElEgNd1);
            if ~all(isMatch, "all")
                [~, iElem] = find(~isMatch);
                [flag, msg] = setFail("elem.edge does not match edge.node in elements: " + idxStr(unique(iElem)));
                return;
            end
            % Duplicated edge and face.
            if size(unique(sort(EgNode, 1)', "rows"), 1) < nEdge || size(unique(sort(FcNode, 1)', "rows"), 1) < nFace
                [flag, msg] = setFail("Duplicated edges or faces.");
                return;
            end
            % Element - face: local face nodes L = G(perm(c, :)) for global face nodes G, sign by parity of c.
            if any(ElFace == 0, "all") || any(abs(ElFace) > nFace, "all") || any(~ismember(ElFcPerm, 1:6), "all")
                [flag, msg] = setFail("Invalid index in elem.face or elem.facePerm.");
                return;
            end
            LcFcNode = MshEnt("D3T").face.node;
            perm = [1, 2, 3; 2, 3, 1; 3, 1, 2; 1, 3, 2; 3, 2, 1; 2, 1, 3];
            isMatch = true(4, nElem);
            for iLcFc = 1:4
                LcNd = ElNode(LcFcNode(:, iLcFc), :);
                GlNd = FcNode(:, abs(ElFace(iLcFc, :)));
                GlNdPerm = GlNd(sub2ind(size(GlNd), perm(ElFcPerm(iLcFc, :), :)', repmat(1:nElem, 3, 1)));
                isMatch(iLcFc, :) = all(LcNd == GlNdPerm, 1) & sign(ElFace(iLcFc, :)) == 1 - 2 * (ElFcPerm(iLcFc, :) > 3);
            end
            if ~all(isMatch, "all")
                [~, iElem] = find(~isMatch);
                [flag, msg] = setFail("elem.face or elem.facePerm does not match face.node in elements: " + idxStr(unique(iElem)));
                return;
            end
            % Face - element.
            ElIdx = repmat(1:nElem, 4, 1);
            pairElFc = sortrows([abs(ElFace(:)), sign(ElFace(:)) .* ElIdx(:)]);
            [iFc, ~, val] = find(FcElem');
            pairFcEl = sortrows([iFc(:), val(:)]);
            if ~isequal(pairElFc, pairFcEl)
                [flag, msg] = setFail("face.elem does not match elem.face.");
                return;
            end
            nCnElem = sum(FcElem ~= 0, 1);
            if any(nCnElem == 0)
                [flag, msg] = setFail("Faces without connected element: " + idxStr(find(nCnElem == 0)));
                return;
            end
            isBd = nCnElem == 1;
            isInvalid = ~isBd & sign(FcElem(1, :)) .* sign(FcElem(2, :)) ~= -1;
            if any(isInvalid)
                [flag, msg] = setFail("Interior faces without one positive and one negative element: " + idxStr(find(isInvalid)));
                return;
            end
            isInvalid = isBd & sum(FcElem, 1) < 0;
            if options.bdOrien && any(isInvalid)
                [flag, msg] = setFail("Boundary faces are not positive in boundary element: " + idxStr(find(isInvalid)));
                return;
            end
            % Normal of boundary face is outward of domain.
            BdIdx = find(isBd);
            BdElem = abs(sum(FcElem(:, BdIdx), 1));
            BdFcNd1 = X(:, FcNode(1, BdIdx)); BdFcNd2 = X(:, FcNode(2, BdIdx)); BdFcNd3 = X(:, FcNode(3, BdIdx));
            ElCen = (X(:, ElNode(1, BdElem)) + X(:, ElNode(2, BdElem)) + X(:, ElNode(3, BdElem)) + X(:, ElNode(4, BdElem))) / 4;
            isInvalid = dot(cross(BdFcNd2 - BdFcNd1, BdFcNd3 - BdFcNd1), (BdFcNd1 + BdFcNd2 + BdFcNd3) / 3 - ElCen) <= 0;
            if options.bdOrien && any(isInvalid)
                [flag, msg] = setFail("Normals of boundary faces are not outward: " + idxStr(BdIdx(isInvalid)));
                return;
            end
            % Face - edge.
            if any(FcEdge == 0, "all") || any(abs(FcEdge) > nEdge, "all")
                [flag, msg] = setFail("Invalid index in face.edge.");
                return;
            end
            FcEgNd1 = FcNode([1, 2, 3], :); FcEgNd2 = FcNode([2, 3, 1], :);
            EgIdx = abs(FcEdge); sgn = sign(FcEdge);
            EgNd1 = reshape(EgNode(1, EgIdx), size(EgIdx)); EgNd2 = reshape(EgNode(2, EgIdx), size(EgIdx));
            isMatch = (sgn > 0 & EgNd1 == FcEgNd1 & EgNd2 == FcEgNd2) | (sgn < 0 & EgNd1 == FcEgNd2 & EgNd2 == FcEgNd1);
            if ~all(isMatch, "all")
                [~, iFc] = find(~isMatch);
                [flag, msg] = setFail("face.edge does not match edge.node in faces: " + idxStr(unique(iFc)));
                return;
            end
            % Type of boundary face.
            if options.bdType && ~isempty(msh.face.type)
                isInvalid = isBd ~= (real(msh.face.type) > 0);
                if any(isInvalid)
                    [flag, msg] = setFail("Boundary faces and positive face types do not match: " + idxStr(find(isInvalid)));
                    return;
                end
            end
            % Euler characteristic.
            if options.euler && msh.nNode - nEdge + nFace - nElem ~= 1
                [flag, msg] = setFail(sprintf("Euler characteristic V - E + F - T = %d.", msh.nNode - nEdge + nFace - nElem));
                return;
            end
            % Total volume.
            if ~isempty(options.vol) && abs(sum(vol) - options.vol) > options.tol
                [flag, msg] = setFail(sprintf("Total volume %.15g is not %.15g.", sum(vol), options.vol));
                return;
            end
        otherwise
            [flag, msg] = setFail("Unsupported mesh type: " + msh.type);
    end
end
%% Local functions.
function [flag, msg] = setFail(msg)
    % setFail: return failed flag with message.
    flag = false;
    msg = "checkMsh: " + msg;
end
function str = idxStr(idx)
    % idxStr: string of (at most 10) indices.
    idx = idx(:)';
    str = strjoin(string(idx(1:min(10, end))), ", ");
    if length(idx) > 10
        str = str + ", ...";
    end
end
