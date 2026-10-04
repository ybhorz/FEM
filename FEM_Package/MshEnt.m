classdef MshEnt
    % MshEnt: mesh entity.

    % Type   | Dim | Geometry                | Mesh       | Variable | Parameter
    % ------ | --- | ----------------------  | ---------- | -------- | -------------
    % D2     | 2   |                         |            | x; y     |
    % D2R    | 2   |                         |            | l; m     |
    % D2R1   | 1   |                         |            | s        |
    % D2T    | 2   | triangle [ x1, x2, x3;  | 2D element | x; y     | [ x1, x2, x3;
    %        |     |            y1, y2, y3 ] |            |          |   y1, y2, y3 ]
    % D2TR   | 2   | triangle [ 0, 1, 0;     | 2D element | l; m     | [ x1, x2, x3;
    %        |     |            0, 0, 1 ]    |            |          |   y1, y2, y3 ]
    % D2L    | 2   | line [ x1, x2;          | 2D edge    | x; y     | [ x1, x2;
    %        |     |        y1, y2 ]         |            |          |   y1, y2 ]
    % D2LR   | 1   | line [ 0, 1 ]           | 1D element | s        | [ x1, x2;
    %        |     |                         |            |          |   y1, y2 ]
    % D3     | 3   |                         |            | x; y; z  |
    % D3R    | 3   |                         |            | l; m; n  |
    % D3R2   | 2   |                         |            | s; t     |
    % D3R1   | 1   |                         |            | s        |
    % D3T    | 3   | tetrahedron             | 3D element | x; y; z  | [ x1, x2, x3, x4;
    %        |     | [ x1, x2, x3, x4;       |            |          |   y1, y2, y3, y4;
    %        |     |   y1, y2, y3, y4;       |            |          |   z1, z2, z3, z4 ]
    %        |     |   z1, z2, z3, z4 ]      |            |          |
    % D3TR   | 3   | tetrahedron             | 3D element | l; m; n  | [ x1, x2, x3, x4;
    %        |     | [ 0, 1, 0, 0;           |            |          |   y1, y2, y3, y4;
    %        |     |   0, 0, 1, 0;           |            |          |   z1, z2, z3, z4 ]
    %        |     |   0, 0, 0, 1 ]          |            |          |
    % D3F    | 3   | triangle [ x1, x2, x3;  | 3D face    | x; y; z  | [ x1, x2, x3;
    %        |     |            y1, y2, y3;  |            |          |   y1, y2, y3;
    %        |     |            z1, z2, z3 ] |            |          |   z1, z2, z3 ]
    % D3FR   | 2   | triangle [ 0, 1, 0;     | 2D element | s; t     | [ x1, x2, x3;
    %        |     |            0, 0, 1 ]    |            |          |   y1, y2, y3;
    %        |     |                         |            |          |   z1, z2, z3 ]
    % D3L    | 3   | line [ x1, x2;          | 3D edge    | x; y; z  | [ x1, x2;
    %        |     |        y1, y2;          |            |          |   y1, y2;
    %        |     |        z1, z2 ]         |            |          |   z1, z2 ]
    % D3LR   | 1   | line [ 0, 1 ]           | 1D element | s        | [ x1, x2;
    %        |     |                         |            |          |   y1, y2;
    %        |     |                         |            |          |   z1, z2 ]
    %
    % Local numbering of tetrahedron (D3T, D3TR):
    % - edge: [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4].
    % - face: i-th face is opposite to i-th vertex, [2, 3, 4], [1, 4, 3], [1, 2, 4], [1, 3, 2],
    %   nodes are counter-clockwise viewed from outside (normal is outward) for positively oriented tetrahedron.

    properties (Constant)
        typeLst = ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
            "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"]; % List of mesh entity types.
        % Remark: property and argument validators require literal lists, so the list is repeated there
        % (MshEnt.type and Fcn.domn); MshEntTest checks that they accept all types in `typeLst`.
    end
    properties
        type {mustBeMember(type, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
            "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"])} = "VOID"; % Type of mesh entity.
        % D: dimension.
        % T: triangle.
        % T: tetrahedron (3D).
        % F: face.
        % L: line.
        % R: reference.
    end
    properties (Dependent)
        dim; % Dimension.
        % Mesh.
        msh Msh; % Mesh.
        node Node; % Node.
        elem Elem; % Element.
        edge Edge; % Edge.
        face Face; % Face.
        nNode; % Number of nodes.
        nEdge; % Number of edges.
        nFace; % Number of faces.
        % Function.
        var (:, 1) sym; % Standard symbolic variable.
        parm (:, :) sym; % Standard symbolic parameter.
        nVar; % Number of variables.
        sParm; % Size of parameter.
        % Geometry.
        UNV Fcn; % Unit normal vector on facet (edge in 2D, face in 3D).
        UTV Fcn; % Unit tangent vector on edge.
        len Fcn; % Length of edge.
        area Fcn; % Area of face.
    end
    methods
        %% Constructor.
        function mshEnt = MshEnt(type)
            arguments
                type {mustBeMember(type, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
                    "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"])} = "VOID";
            end
            mshEnt.type = type;
        end
        %% Get functions.
        function dim = get.dim(mshEnt)
            dim = MshEnt.getInfo(mshEnt.type).dim;
        end
        function msh = get.msh(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    msh = Msh("D2T", mshEnt.node, mshEnt.elem, mshEnt.edge);
                case {"D2LR", "D3LR"}
                    msh = Msh("D1", mshEnt.node, mshEnt.elem);
                case {"D3T", "D3TR"}
                    msh = Msh("D3T", mshEnt.node, mshEnt.elem, mshEnt.edge, mshEnt.face);
                case "D3FR"
                    msh = Msh("D2T", mshEnt.node, mshEnt.elem, mshEnt.edge);
                otherwise
                    msh = Msh.empty;
            end
        end
        function node = get.node(mshEnt)
            switch mshEnt.type
                case "D2T"
                    node = Node(str2sym("[x1, x2, x3; y1, y2, y3]"));
                case "D2TR"
                    node = Node(sym([0, 1, 0; 0, 0, 1]));
                case "D2L"
                    node = Node(str2sym("[x1, x2; y1, y2]"));
                case "D2LR"
                    node = Node(sym([0, 1]));
                case {"D3T", "D3F", "D3L"}
                    node = Node(mshEnt.parm);
                case "D3TR"
                    node = Node(sym([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]));
                case "D3FR"
                    node = Node(sym([0, 1, 0; 0, 0, 1]));
                case "D3LR"
                    node = Node(sym([0, 1]));
                otherwise
                    node = Node.empty;
            end
        end
        function elem = get.elem(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    elem = Elem([1; 2; 3]);
                case "D2LR"
                    elem = Elem([1; 2]);
                case {"D3T", "D3TR"}
                    elem = Elem([1; 2; 3; 4]);
                case "D3FR"
                    elem = Elem([1; 2; 3]);
                case "D3LR"
                    elem = Elem([1; 2]);
                otherwise
                    elem = Elem.empty;
            end
        end
        function edge = get.edge(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    edge = Edge([1, 2, 3; 2, 3, 1]);
                case "D2L"
                    edge = Edge([1; 2]);
                case {"D3T", "D3TR"}
                    edge = Edge([1, 1, 1, 2, 2, 3; 2, 3, 4, 3, 4, 4]);
                case {"D3F", "D3FR"}
                    edge = Edge([1, 2, 3; 2, 3, 1]);
                case "D3L"
                    edge = Edge([1; 2]);
                otherwise
                    edge = Edge.empty;
            end
        end
        function face = get.face(mshEnt)
            switch mshEnt.type
                case {"D3T", "D3TR"}
                    face = Face([2, 1, 1, 1; 3, 4, 2, 3; 4, 3, 4, 2]);
                case "D3F"
                    face = Face([1; 2; 3]);
                otherwise
                    face = Face.empty;
            end
        end
        function nNode = get.nNode(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR", "D2L", "D2LR", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"}
                    nNode = mshEnt.node.nNode;
                otherwise
                    nNode = [];
            end
        end
        function nEdge = get.nEdge(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR", "D3T", "D3TR", "D3F", "D3FR"}
                    nEdge = mshEnt.edge.nEdge;
                otherwise
                    nEdge = [];
            end
        end
        function nFace = get.nFace(mshEnt)
            switch mshEnt.type
                case {"D3T", "D3TR"}
                    nFace = mshEnt.face.nFace;
                otherwise
                    nFace = [];
            end
        end
        function var = get.var(mshEnt)
            var = MshEnt.getInfo(mshEnt.type).var;
        end
        function parm = get.parm(mshEnt)
            parm = MshEnt.getInfo(mshEnt.type).parm;
        end
        function nVar = get.nVar(mshEnt)
            nVar = length(mshEnt.var);
        end
        function sParm = get.sParm(mshEnt)
            sParm = size(mshEnt.parm);
        end
        function UNV = get.UNV(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2L"}
                    UNV(1:mshEnt.edge.nEdge) = Fcn.cst(0);
                    for iEdge = 1:mshEnt.edge.nEdge
                        EgNd1 = mshEnt.node.coord(:, mshEnt.edge.node(1, iEdge));
                        EgNd2 = mshEnt.node.coord(:, mshEnt.edge.node(2, iEdge));
                        tan = EgNd2 - EgNd1;
                        norm = sqrt(tan(1)^2 + tan(2)^2);
                        UNV(iEdge) = Fcn(mshEnt.type, [tan(2); -tan(1)] / norm);
                    end
                case {"D3T", "D3F"}
                    UNV(1:mshEnt.face.nFace) = Fcn.cst(0);
                    for iFace = 1:mshEnt.face.nFace
                        nor = crsFace(mshEnt, iFace);
                        norm = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
                        UNV(iFace) = Fcn(mshEnt.type, nor / norm);
                    end
                otherwise
                    UNV = Fcn.empty;
            end
        end
        function UTV = get.UTV(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2L"}
                    UTV(1:mshEnt.edge.nEdge) = Fcn.cst(0);
                    for iEdge = 1:mshEnt.edge.nEdge
                        EgNd1 = mshEnt.node.coord(:, mshEnt.edge.node(1, iEdge));
                        EgNd2 = mshEnt.node.coord(:, mshEnt.edge.node(2, iEdge));
                        tan = EgNd2 - EgNd1;
                        norm = sqrt(tan(1)^2 + tan(2)^2);
                        UTV(iEdge) = Fcn(mshEnt.type, tan / norm);
                    end
                case {"D3T", "D3F", "D3L"}
                    UTV(1:mshEnt.edge.nEdge) = Fcn.cst(0);
                    for iEdge = 1:mshEnt.edge.nEdge
                        EgNd1 = mshEnt.node.coord(:, mshEnt.edge.node(1, iEdge));
                        EgNd2 = mshEnt.node.coord(:, mshEnt.edge.node(2, iEdge));
                        tan = EgNd2 - EgNd1;
                        norm = sqrt(tan(1)^2 + tan(2)^2 + tan(3)^2);
                        UTV(iEdge) = Fcn(mshEnt.type, tan / norm);
                    end
                otherwise
                    UTV = Fcn.empty;
            end
        end
        function len = get.len(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2L"}
                    len(1:mshEnt.edge.nEdge) = Fcn.cst(0);
                    for iEdge = 1:mshEnt.edge.nEdge
                        EgNd1 = mshEnt.node.coord(:, mshEnt.edge.node(1, iEdge));
                        EgNd2 = mshEnt.node.coord(:, mshEnt.edge.node(2, iEdge));
                        tan = EgNd2 - EgNd1;
                        norm = sqrt(tan(1)^2 + tan(2)^2);
                        len(iEdge) = Fcn(mshEnt.type, norm);
                    end
                case {"D3T", "D3F", "D3L"}
                    len(1:mshEnt.edge.nEdge) = Fcn.cst(0);
                    for iEdge = 1:mshEnt.edge.nEdge
                        EgNd1 = mshEnt.node.coord(:, mshEnt.edge.node(1, iEdge));
                        EgNd2 = mshEnt.node.coord(:, mshEnt.edge.node(2, iEdge));
                        tan = EgNd2 - EgNd1;
                        norm = sqrt(tan(1)^2 + tan(2)^2 + tan(3)^2);
                        len(iEdge) = Fcn(mshEnt.type, norm);
                    end
                otherwise
                    len = Fcn.empty;
            end
        end
        function area = get.area(mshEnt)
            switch mshEnt.type
                case {"D3T", "D3F"}
                    area(1:mshEnt.face.nFace) = Fcn.cst(0);
                    for iFace = 1:mshEnt.face.nFace
                        nor = crsFace(mshEnt, iFace);
                        norm = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
                        area(iFace) = Fcn(mshEnt.type, norm / 2);
                    end
                otherwise
                    area = Fcn.empty;
            end
        end
        %% Public functions.
        function mshEnt = dual(mshEnt)
            % MshEnt.dual: dual of mesh entity (original <-> reference).
            arguments (Input)
                mshEnt MshEnt;
            end
            arguments (Output)
                mshEnt MshEnt;
            end
            mshEnt.type = MshEnt.getInfo(mshEnt.type).dual;
        end
        function flag = ismember(mshEnt1, mshEnt2)
            % MshEnt.ismember: check if `mshEnt1` is member of `mshEnt2`.

            %      | D2 | D2R | D2R1| D2T | D2TR | D2L | D2LR
            %------|----|-----|-----|-----|------|-----|------
            %  D2  | Y  |  N  |  N  |  Y  |  N   |  Y  |  N
            %  D2R | N  |  Y  |  N  |  N  |  Y   |  N  |  N
            % D2R1 | N  |  N  |  Y  |  N  |  N   |  N  |  Y
            %  D2T | N  |  N  |  N  |  Y  |  N   |  N  |  N
            % D2TR | N  |  N  |  N  |  N  |  Y   |  N  |  N
            %  D2L | N  |  N  |  N  |  Y  |  N   |  Y  |  N
            % D2LR | N  |  N  |  N  |  N  |  N   |  N  |  Y
            % For all types (including 3D), see column `super` of the table in `genInfo`.

            arguments
                mshEnt1 MshEnt;
                mshEnt2 MshEnt;
            end
            flag = ismember(mshEnt2.type, MshEnt.getInfo(mshEnt1.type).super);
        end
    end
    %% Static functions.
    methods (Static)
        function info = getInfo(type)
            % MshEnt.getInfo: get data of mesh entity type.
            % Data are generated once and cached.
            arguments
                type (1, 1) string;
            end
            persistent infoLst;
            if isempty(infoLst)
                infoLst = genInfo();
            end
            info = infoLst.(type);
        end
    end
end
%% Local functions.
function infoLst = genInfo()
    % genInfo: generate data of all mesh entity types.
    % - dim: dimension.
    % - var: standard symbolic variable.
    % - parm: standard symbolic parameter.
    % - dual: dual type (original <-> reference).
    % - free: type after clearing parameter.
    % - super: types of which this type is member.
    %
    % Type | dim | var       | parm                                 | dual | free | super
    % ---- | --- | --------- | ------------------------------------ | ---- | ---- | -----------------
    % VOID |     |           |                                      | VOID | VOID | all types
    % D2   | 2   | [x; y]    |                                      | D2R  | D2   | D2, D2T, D2L
    % D2R  | 2   | [l; m]    |                                      | D2   | D2R  | D2R, D2TR
    % D2R1 | 1   | s         |                                      | VOID | D2R1 | D2R1, D2LR
    % D2T  | 2   | [x; y]    | [x1, x2, x3; y1, y2, y3]             | D2TR | D2   | D2T
    % D2TR | 2   | [l; m]    | [x1, x2, x3; y1, y2, y3]             | D2T  | D2R  | D2TR
    % D2L  | 2   | [x; y]    | [x1, x2; y1, y2]                     | D2LR | D2   | D2L, D2T
    % D2LR | 1   | s         | [x1, x2; y1, y2]                     | D2L  | D2R1 | D2LR
    % D3   | 3   | [x; y; z] |                                      | D3R  | D3   | D3, D3T, D3F, D3L
    % D3R  | 3   | [l; m; n] |                                      | D3   | D3R  | D3R, D3TR
    % D3R2 | 2   | [s; t]    |                                      | VOID | D3R2 | D3R2, D3FR
    % D3R1 | 1   | s         |                                      | VOID | D3R1 | D3R1, D3LR
    % D3T  | 3   | [x; y; z] | [x1, .., x4; y1, .., y4; z1, .., z4] | D3TR | D3   | D3T
    % D3TR | 3   | [l; m; n] | [x1, .., x4; y1, .., y4; z1, .., z4] | D3T  | D3R  | D3TR
    % D3F  | 3   | [x; y; z] | [x1, x2, x3; y1, y2, y3; z1, z2, z3] | D3FR | D3   | D3F, D3T
    % D3FR | 2   | [s; t]    | [x1, x2, x3; y1, y2, y3; z1, z2, z3] | D3F  | D3R2 | D3FR
    % D3L  | 3   | [x; y; z] | [x1, x2; y1, y2; z1, z2]             | D3LR | D3   | D3L, D3F, D3T
    % D3LR | 1   | s         | [x1, x2; y1, y2; z1, z2]             | D3L  | D3R1 | D3LR
    infoLst.VOID = setInfo([], "", "", "VOID", "VOID", MshEnt.typeLst);
    infoLst.D2 = setInfo(2, "[x; y]", "", "D2R", "D2", ["D2", "D2T", "D2L"]);
    infoLst.D2R = setInfo(2, "[l; m]", "", "D2", "D2R", ["D2R", "D2TR"]);
    infoLst.D2R1 = setInfo(1, "s", "", "VOID", "D2R1", ["D2R1", "D2LR"]);
    infoLst.D2T = setInfo(2, "[x; y]", "[x1, x2, x3; y1, y2, y3]", "D2TR", "D2", "D2T");
    infoLst.D2TR = setInfo(2, "[l; m]", "[x1, x2, x3; y1, y2, y3]", "D2T", "D2R", "D2TR");
    infoLst.D2L = setInfo(2, "[x; y]", "[x1, x2; y1, y2]", "D2LR", "D2", ["D2L", "D2T"]);
    infoLst.D2LR = setInfo(1, "s", "[x1, x2; y1, y2]", "D2L", "D2R1", "D2LR");
    infoLst.D3 = setInfo(3, "[x; y; z]", "", "D3R", "D3", ["D3", "D3T", "D3F", "D3L"]);
    infoLst.D3R = setInfo(3, "[l; m; n]", "", "D3", "D3R", ["D3R", "D3TR"]);
    infoLst.D3R2 = setInfo(2, "[s; t]", "", "VOID", "D3R2", ["D3R2", "D3FR"]);
    infoLst.D3R1 = setInfo(1, "s", "", "VOID", "D3R1", ["D3R1", "D3LR"]);
    infoLst.D3T = setInfo(3, "[x; y; z]", "[x1, x2, x3, x4; y1, y2, y3, y4; z1, z2, z3, z4]", "D3TR", "D3", "D3T");
    infoLst.D3TR = setInfo(3, "[l; m; n]", "[x1, x2, x3, x4; y1, y2, y3, y4; z1, z2, z3, z4]", "D3T", "D3R", "D3TR");
    infoLst.D3F = setInfo(3, "[x; y; z]", "[x1, x2, x3; y1, y2, y3; z1, z2, z3]", "D3FR", "D3", ["D3F", "D3T"]);
    infoLst.D3FR = setInfo(2, "[s; t]", "[x1, x2, x3; y1, y2, y3; z1, z2, z3]", "D3F", "D3R2", "D3FR");
    infoLst.D3L = setInfo(3, "[x; y; z]", "[x1, x2; y1, y2; z1, z2]", "D3LR", "D3", ["D3L", "D3F", "D3T"]);
    infoLst.D3LR = setInfo(1, "s", "[x1, x2; y1, y2; z1, z2]", "D3L", "D3R1", "D3LR");
end
function nor = crsFace(mshEnt, iFace)
    % crsFace: cross product (v2 - v1) x (v3 - v1) of vertices of i-th face.
    FcNd1 = mshEnt.node.coord(:, mshEnt.face.node(1, iFace));
    FcNd2 = mshEnt.node.coord(:, mshEnt.face.node(2, iFace));
    FcNd3 = mshEnt.node.coord(:, mshEnt.face.node(3, iFace));
    a = FcNd2 - FcNd1;
    b = FcNd3 - FcNd1;
    nor = [a(2) * b(3) - a(3) * b(2); a(3) * b(1) - a(1) * b(3); a(1) * b(2) - a(2) * b(1)];
end
function info = setInfo(dim, var, parm, dual, free, super)
    % setInfo: set data of mesh entity type.
    info.dim = dim;
    if var == ""
        info.var = sym([]);
    else
        info.var = str2sym(var);
    end
    if parm == ""
        info.parm = sym([]);
    else
        info.parm = str2sym(parm);
    end
    info.dual = dual;
    info.free = free;
    info.super = super;
end
