function plotFcn(msh, fcn)
    % plotFcn: plot function.
    arguments
        msh Msh; % Mesh.
        fcn Fcn; % Function.
    end
    assert(ismember(msh.type, "D2T"));
    assert(ismember(fcn.domn, ["D2", "D2T"]));
    fun = fcn.getFun;
    V = zeros(msh.nElem * 3, 3, fcn.nFun);
    F = zeros(msh.nElem, 3);
    for iElem = 1:msh.nElem
        ElVert = msh.node.coord(:, msh.elem.node(:, iElem));
        val = zeros(fcn.nFun, msh.elem.nNode);
        switch fcn.domn
            case "D2"
                for iVert = 1:msh.elem.nNode
                    temp = fun(ElVert(:, iVert));
                    val(:, iVert) = temp(:);
                end
            case "D2T"
                ElParm = msh.node.coord(:, msh.elem.node(:, iElem));
                for iVert = 1:msh.elem.nNode
                    temp = fun(ElVert(:, iVert), ElParm);
                    val(:, iVert) = temp(:);
                end
        end
        idx = (iElem - 1) * 3 + (1:3);
        for iFun = 1:fcn.nFun
            V(idx, :, iFun) = [ElVert', val(iFun, :)'];
        end
        F(iElem, :) = idx;
    end
    for iFun = 1:fcn.nFun
        figure;
        patch('Faces', F, 'Vertices', V(:, :, iFun), 'FaceVertexCData', V(:, 3, iFun), 'FaceColor', 'interp');
        view(3);
    end
end
