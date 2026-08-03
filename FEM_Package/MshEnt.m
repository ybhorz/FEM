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

    properties
        type {mustBeMember(type, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR"])} = "VOID"; % Type of mesh entity.
        % D: dimension.
        % T: triangle.
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
        nNode; % Number of nodes.
        nEdge; % Number of edges.
        % Function.
        var (:, 1) sym; % Standard symbolic variable.
        parm (:, :) sym; % Standard symbolic parameter.
        nVar; % Number of variables.
        sParm; % Size of parameter.
        % Geometry.
        UNV Fcn; % Unit normal vector on edge.
        UTV Fcn; % Unit tangent vector on edge.
        len Fcn; % Length of edge.
    end
    methods
        %% Constructor.
        function mshEnt = MshEnt(type)
            arguments
                type {mustBeMember(type, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR"])} = "VOID";
            end
            mshEnt.type = type;
        end
        %% Get functions.
        function dim = get.dim(mshEnt)
            switch mshEnt.type
                case {"D2", "D2R", "D2T", "D2TR", "D2L"}
                    dim = 2;
                case {"D2R1", "D2LR"}
                    dim = 1;
                otherwise
                    dim = [];
            end
        end
        function msh = get.msh(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    msh = Msh("D2T", mshEnt.node, mshEnt.elem, mshEnt.edge);
                case "D2LR"
                    msh = Msh("D1", mshEnt.node, mshEnt.elem);
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
                otherwise
                    edge = Edge.empty;
            end
        end
        function nNode = get.nNode(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR", "D2L", "D2LR"}
                    nNode = mshEnt.node.nNode;
                otherwise
                    nNode = [];
            end
        end
        function nEdge = get.nEdge(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    nEdge = mshEnt.edge.nEdge;
                otherwise
                    nEdge = [];
            end
        end
        function var = get.var(mshEnt)
            switch mshEnt.type
                case {"D2", "D2T", "D2L"}
                    var = str2sym("[x; y]");
                case {"D2R", "D2TR"}
                    var = str2sym("[l; m]");
                case {"D2R1", "D2LR"}
                    var = str2sym("s");
                otherwise
                    var = sym([]);
            end
        end
        function parm = get.parm(mshEnt)
            switch mshEnt.type
                case {"D2T", "D2TR"}
                    parm = str2sym("[x1, x2, x3; y1, y2, y3]");
                case {"D2L", "D2LR"}
                    parm = str2sym("[x1, x2; y1, y2]");
                otherwise
                    parm = sym([]);
            end
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
                otherwise
                    len = Fcn.empty;
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
            switch mshEnt.type
                case "D2"
                    mshEnt.type = "D2R";
                case "D2R"
                    mshEnt.type = "D2";
                case "D2T"
                    mshEnt.type = "D2TR";
                case "D2TR"
                    mshEnt.type = "D2T";
                case "D2L"
                    mshEnt.type = "D2LR";
                case "D2LR"
                    mshEnt.type = "D2L";
                otherwise
                    mshEnt.type = "VOID";
            end
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

            arguments
                mshEnt1 MshEnt;
                mshEnt2 MshEnt;
            end
            switch mshEnt1.type
                case "VOID"
                    flag = true;
                case "D2"
                    if ismember(mshEnt2.type, ["D2", "D2T", "D2L"])
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2R"
                    if ismember(mshEnt2.type, ["D2R", "D2TR"])
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2R1"
                    if ismember(mshEnt2.type, ["D2R1", "D2LR"])
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2T"
                    if ismember(mshEnt2.type, "D2T")
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2TR"
                    if ismember(mshEnt2.type, "D2TR")
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2L"
                    if ismember(mshEnt2.type, ["D2L", "D2T"])
                        flag = true;
                    else
                        flag = false;
                    end
                case "D2LR"
                    if ismember(mshEnt2.type, "D2LR")
                        flag = true;
                    else
                        flag = false;
                    end
            end
        end
    end
end
