classdef StyleTest < matlab.unittest.TestCase
    % StyleTest: check coding style of source files.

    % Kind        | Files                                                        | Line ending
    % ----------- | ------------------------------------------------------------ | -----------
    % code        | FEM_Package/*.m, FEM_Test/unit/*.m, FEM_Test/unit/helper/*.m | CRLF
    %             | FEM_Test/baseline/*.m, runTests.m                            |
    % live script | FEM_Test/demo/*.m, FEM_Example/*.m                           | LF
    %
    % All files: ASCII only, no tab.
    % Code files: no trailing whitespace (except legacy lines listed in `legacyTrail`).

    properties
        root; % Root folder of repository.
        legacyTrail = struct("mshD2TS", 1); % Number of legacy lines with trailing whitespace.
    end
    methods (TestClassSetup)
        function setRoot(tc)
            tc.root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
        end
    end
    methods (Test, TestTags = {'Fast', 'Style'})
        function codeFile(tc)
            files = [lstFile(fullfile(tc.root, "FEM_Package")), lstFile(fullfile(tc.root, "FEM_Test", "unit")), ...
                lstFile(fullfile(tc.root, "FEM_Test", "unit", "helper")), lstFile(fullfile(tc.root, "FEM_Test", "baseline")), string(fullfile(tc.root, "runTests.m"))];
            tc.assertNotEmpty(files);
            for iFile = 1:length(files)
                bytes = readByte(files(iFile));
                [~, name] = fileparts(files(iFile));
                tc.verifyTrue(all(bytes < 128), files(iFile) + ": non-ASCII character.");
                tc.verifyFalse(any(bytes == 9), files(iFile) + ": tab character.");
                LFIdx = find(bytes == 10);
                tc.verifyTrue(all(LFIdx > 1) && all(bytes(max(LFIdx - 1, 1)) == 13) && nnz(bytes == 13) == length(LFIdx), ...
                    files(iFile) + ": line ending is not CRLF.");
                nTrail = cntTrail(bytes);
                if isfield(tc.legacyTrail, name)
                    nAllow = tc.legacyTrail.(name);
                else
                    nAllow = 0;
                end
                tc.verifyLessThanOrEqual(nTrail, nAllow, files(iFile) + ": trailing whitespace.");
            end
        end
        function liveScript(tc)
            files = [lstFile(fullfile(tc.root, "FEM_Test", "demo")), lstFile(fullfile(tc.root, "FEM_Example"))];
            tc.assertNotEmpty(files);
            for iFile = 1:length(files)
                bytes = readByte(files(iFile));
                tc.verifyTrue(all(bytes < 128), files(iFile) + ": non-ASCII character.");
                tc.verifyFalse(any(bytes == 9), files(iFile) + ": tab character.");
                tc.verifyFalse(any(bytes == 13), files(iFile) + ": line ending is not LF.");
            end
        end
    end
end
%% Local functions.
function files = lstFile(folder, pattern)
    % lstFile: list files in folder (by row).
    arguments
        folder (1, 1) string;
        pattern (1, 1) string = "*.m";
    end
    lst = dir(fullfile(folder, pattern));
    files = strings(1, length(lst));
    for iLst = 1:length(lst)
        files(iLst) = string(fullfile(lst(iLst).folder, lst(iLst).name));
    end
end
function bytes = readByte(file)
    % readByte: read file as bytes (by row).
    fid = fopen(file, "r");
    assert(fid > 0, "Cannot open file: " + file);
    bytes = fread(fid, Inf, "*uint8")';
    fclose(fid);
end
function nTrail = cntTrail(bytes)
    % cntTrail: count lines with trailing whitespace.
    bytes = bytes(bytes ~= 13);
    LFIdx = find(bytes == 10);
    LFIdx = LFIdx(LFIdx > 1);
    nTrail = nnz(bytes(LFIdx - 1) == 32);
    if ~isempty(bytes) && bytes(end) == 32
        nTrail = nTrail + 1;
    end
end
