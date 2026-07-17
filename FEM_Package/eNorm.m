function eNorm = eNorm(msh, reals, apprs, norms)
    % eNorm: compute error norm between real and approximate solutions.
    arguments
        msh Msh; % Mesh.
        reals (1, :) Fcn; % Real solutions.
        apprs (1, :) FEF; % Approximate solutions.
        norms (1, :) Norm; % Norms.
    end
    assert(ismember(msh.type, "D2T"));
    assert(isequal(msh, apprs.getMsh));
    assert(isequal(msh, norms.getMsh));
    eNorms = zeros(1, length(norms));
    for iNorm = 1:length(norms)
        norm = norms(iNorm);
        real = reals(norm.iFcn);
        appr = apprs(norm.iFcn);
        switch norm.EntDim
            case 2
                assert(ismember(real.domn, ["VOID", "D2", "D2T"]));
                assert(ismember(appr.elem, "D2T"));
                errFcn = dif(real - appr, norm.ord);
                intFcn = norm.form(norm.coef, errFcn, norm.pow);
                intFun = intFcn.getFun;
                for iElem = norm.EntIdx
                    ElParm = appr.ElParm(:, :, iElem);
                    ElCoef = appr.ElCoef(:, iElem);
                    intVal = norm.GInt.eval(@(x) intFun(x, ElParm, ElCoef), ElParm);
                    eNorms(iNorm) = eNorms(iNorm) + intVal;
                end
                eNorms(iNorm) = eNorms(iNorm) ^ (1 / norm.pow);
            case 1
                assert(ismember(appr.elem, ["D2LR", "D2T"]));
                switch appr.elem
                    case "D2LR"
                        assert(ismember(real.domn, ["VOID", "D2", "D2L"]));
                        assert(all(norm.ord == 0));
                        assert(isequal(norm.fcnOpr, "none"));
                        errFcn = dif(real.tfm("D2LR") - appr, norm.ord);
                        intFcn = norm.form(norm.coef.tfm("D2LR"), errFcn, norm.pow) .* Tfm("D2L").JNorm;
                        intFun = intFcn.getFun;
                        for iElem = norm.EntIdx
                            ElParm = appr.ElParm(:, :, iElem);
                            ElCoef = appr.ElCoef(:, iElem);
                            intVal = norm.GInt.eval(@(x) intFun(x, ElParm, ElCoef));
                            eNorms(iNorm) = eNorms(iNorm) + intVal;
                        end
                        eNorms(iNorm) = eNorms(iNorm) ^ (1 / norm.pow);
                    case "D2T"
                        assert(ismember(real.domn, ["VOID", "D2"]));
                        assert(all(norm.ord == 0, "all"));
                        assert(isequal(norm.fcnOpr, "jump"));
                        EgPmSym = MshEnt("D2LR").parm;
                        El1PmSym = sym('El1Parm', appr.sElParm);
                        El2PmSym = sym('El2Parm', appr.sElParm);
                        El1CfSym = sym('El1Coef', [appr.nElCoef, 1]);
                        El2CfSym = sym('El2Coef', [appr.nElCoef, 1]);
                        errFcn = appr.dif(norm.ord).subParm(El1PmSym).subCoef(El1CfSym) - appr.dif(norm.ord).subParm(El2PmSym).subCoef(El2CfSym);
                        intFcn = norm.form(norm.coef, errFcn, norm.pow);
                        intFun = intFcn.getFun("parm", {EgPmSym, El1PmSym, El2PmSym}, "coef", {El1CfSym, El2CfSym});
                        for iEdge = norm.EntIdx
                            EgParm = msh.node.coord(:, msh.edge.node(:, iEdge));
                            CnElem = msh.edge.elem(:, iEdge);
                            iElem1 = abs(CnElem(1)); iElem2 = abs(CnElem(2));
                            if iElem1 == 0 || iElem2 == 0
                                continue
                            end
                            El1Parm = appr.ElParm(:, :, iElem1);
                            El2Parm = appr.ElParm(:, :, iElem2);
                            El1Coef = appr.ElCoef(:, iElem1);
                            El2Coef = appr.ElCoef(:, iElem2);
                            intVal = norm.GInt.eval(@(x) intFun(x, EgParm, El1Parm, El2Parm, El1Coef, El2Coef), EgParm);
                            eNorms(iNorm) = eNorms(iNorm) + intVal;
                        end
                        eNorms(iNorm) = eNorms(iNorm) ^ (1 / norm.pow);
                end
        end
    end
    eNorm = sum(eNorms);
end
