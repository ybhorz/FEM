function makeBaseline2D(options)
    % makeBaseline2D: generate 2D regression baseline `baseline_2D.mat`.
    % Run only on unmodified 2D code (the commit before Phase 0), and commit the generated file.
    % With "update", only the error norms of the given examples are recomputed and replaced in the existing baseline;
    % use it only for intentional changes of an example (recorded in `info.update`).
    % Example: makeBaseline2D("update", "Posisson_HDG");
    arguments
        options.example logical = true; % Whether run examples (slow).
        options.update (1, :) string = string.empty; % Examples whose baseline is updated.
    end
    root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(fullfile(root, "FEM_Package"), fullfile(root, "FEM_Test", "baseline"));
    figVis = get(groot, "DefaultFigureVisible");
    set(groot, "DefaultFigureVisible", "off");
    cleanup = onCleanup(@() set(groot, "DefaultFigureVisible", figVis));
    tic;
    file = fullfile(root, "FEM_Test", "baseline", "baseline_2D.mat");
    if ~isempty(options.update)
        base = load(file, "data", "info");
        data = base.data;
        info = base.info;
        new = calcRegression2D("small", false, "example", true, "exName", options.update);
        for iEx = 1:length(options.update)
            name = options.update(iEx);
            fprintf("%s: %s -> %s\n", name, mat2str(data.example.(name), 7), mat2str(new.example.(name), 7));
            data.example.(name) = new.example.(name);
        end
        [~, commit] = system("git -C """ + root + """ rev-parse HEAD");
        upd = struct("date", string(datetime("now")), "commit", strtrim(string(commit)), "example", options.update);
        if isfield(info, "update")
            info.update(end + 1) = upd;
        else
            info.update = upd;
        end
        save(file, "data", "info");
        fprintf("Baseline updated in %s (%.1f s).\n", file, toc);
        return;
    end
    data = calcRegression2D("small", true, "example", options.example);
    info.date = string(datetime("now"));
    info.version = string(version);
    [status, commit] = system("git -C """ + root + """ rev-parse HEAD");
    if status == 0
        info.commit = strtrim(string(commit));
    else
        info.commit = "unknown";
    end
    info.time = toc;
    save(file, "data", "info");
    fprintf("Baseline saved to %s (commit %s, %.1f s).\n", file, info.commit, info.time);
end
