function plotFEF(FEFcn)
    % plotFEF: plot Finite Element function.
    arguments
        FEFcn FEF; % Finite Element function.
    end
    assert(ismember(FEFcn.msh.type, "D2T"));
    assert(ismember(FEFcn.elem, ["D2T", "D2LR"]));
    msh = FEFcn.msh;
    fun = FEFcn.getFun;
    switch FEFcn.elem
        case "D2T"
            V = zeros(msh.nElem * 3, 3, FEFcn.nFun);
            F = zeros(msh.nElem, 3);
            for iElem = 1:msh.nElem
                ElParm = FEFcn.ElParm(:, :, iElem);
                ElCoef = FEFcn.ElCoef(:, iElem);
                ElVert = msh.node.coord(:, msh.elem.node(:, iElem));
                val = zeros(FEFcn.nFun, msh.elem.nNode);
                for iVert = 1:msh.elem.nNode
                    temp = fun(ElVert(:, iVert), ElParm, ElCoef);
                    val(:, iVert) = temp(:);
                end
                idx = (iElem - 1) * 3 + (1:3);
                for iFun = 1:FEFcn.nFun
                    V(idx, :, iFun) = [ElVert', val(iFun, :)'];
                end
                F(iElem, :) = idx;
            end
            for iFun = 1:FEFcn.nFun
                figure;
                patch('Faces', F, 'Vertices', V(:, :, iFun), 'FaceVertexCData', V(:, 3, iFun), 'FaceColor', 'interp');
                view(3);
            end
        case "D2LR"
            EgNode1 = msh.edge.node(1, :);
            EgNode2 = msh.edge.node(2, :);
            X = [msh.node.coord(1, EgNode1); msh.node.coord(1, EgNode2)];
            Y = [msh.node.coord(2, EgNode1); msh.node.coord(2, EgNode2)];
            Z = zeros(2, msh.nEdge, FEFcn.nFun);
            for iEdge = 1:msh.nEdge
                ElParm = FEFcn.ElParm(:, :, iEdge);
                ElCoef = FEFcn.ElCoef(:, iEdge);
                EgVert = [0, 1];
                val = zeros(FEFcn.nFun, 2);
                for iVert = 1:2
                    temp = fun(EgVert(iVert), ElParm, ElCoef);
                    val(:, iVert) = temp(:);
                end
                Z(:, iEdge, :) = reshape(val(:, 1:2).', 2, 1, FEFcn.nFun);
            end
            xLine = [X; nan(1, msh.nEdge)];
            yLine = [Y; nan(1, msh.nEdge)];
            for iFun = 1:FEFcn.nFun
                zLine = [Z(:, :, iFun); nan(1, msh.nEdge)];
                figure;
                plot3(xLine(:), yLine(:), zLine(:), '-');
                view(3);
            end
    end
end
