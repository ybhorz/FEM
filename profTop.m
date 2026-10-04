function profTop(script, nTop)
    % profTop: profile a script and print the functions with largest total time and self time (diagnostic).
    % Example (in FEM_Example): profTop("Poisson_SDG_3D", 25)
    arguments
        script (1, 1) string;
        nTop (1, 1) = 25;
    end
    profile clear;
    profile on;
    tic;
    evalin("base", "run(""" + script + """);");
    t = toc;
    profile off;
    info = profile("info");
    F = info.FunctionTable;
    total = [F.TotalTime]';
    self = total;
    for i = 1:numel(F)
        if ~isempty(F(i).Children)
            self(i) = total(i) - sum([F(i).Children.TotalTime]);
        end
    end
    name = string({F.FunctionName})';
    calls = [F.NumCalls]';
    fprintf("Total %.1f s\n\n%-60s %8s %10s\n", t, "Function (by total time)", "Calls", "Total");
    [~, idx] = sort(total, "descend");
    for i = idx(1:min(nTop, end))'
        fprintf("%-60s %8d %10.2f\n", extractBefore(name(i) + string(blanks(60)), 61), calls(i), total(i));
    end
    fprintf("\n%-60s %8s %10s\n", "Function (by self time)", "Calls", "Self");
    [~, idx] = sort(self, "descend");
    for i = idx(1:min(nTop, end))'
        fprintf("%-60s %8d %10.2f\n", extractBefore(name(i) + string(blanks(60)), 61), calls(i), self(i));
    end
end
