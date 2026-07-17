function plotFEF(FEFcn)
    % plotFEF: plot Finite Element function.
    arguments
        FEFcn FEF; % Finite Element function.
    end
    assert(ismember(FEFcn.msh.type, "D2T"));
    assert(ismember(FEFcn.elem, "D2T"));
    msh = FEFcn.msh;
    fun = FEFcn.getFun;
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
end
