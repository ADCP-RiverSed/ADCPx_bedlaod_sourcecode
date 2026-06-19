function DepthRNGout = processDepthsRNG(dataADCP1200, ADCPin, varargin)
%PROCESSDEPTHSRNG Process only depth/RNG variables from ADCP files.
%
% IMPORTANT:
%   StartADCP and EndADCP are treated as samples to remove:
%
%       start  = StartADCP(i) + 1;
%       finish = N - EndADCP(i);
%
% Example:
%   StartADCP = 30 and EndADCP = 30 means use samples 31:N-30.
%
% This function:
%   - does NOT save anything to disk
%   - matches ADCPin.ADCPfilename with dataADCP1200.name
%   - uses StartADCP(i) and EndADCP(i) from ADCPin for each file
%   - extracts WF_12, WS_12, WT_12, Temp1
%   - extracts and despikes BT.RNG1-4
%   - calculates Depth1200, Depth1200std, stdSpatial_1200
%   - creates cell arrays such as RNG11_{i}, RNG22_{i}, etc.
%   - optionally plots raw/despiked beam depth series
%   - keeps old/raw variable names in workspace if AssignToCaller/AssignToBase = true
%
% Required inputs:
%   dataADCP1200 = dir(filePattern1200);
%   ADCPin       = output from readInputTable_PhSamples
%
% Example:
%   dataADCP1200 = dir(filePattern1200);
%
%   DepthRNGout = processDepthsRNG(dataADCP1200, ADCPin, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
% Plot selected rows:
%   DepthRNGout = processDepthsRNG(dataADCP1200, ADCPin, ...
%       'PlotResults', true, ...
%       'PlotRows', [1 3 7], ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true);
%
% Main workspace variables created:
%   WF_12, WS_12, WT_12, Temp1
%   RNG1_1200, RNG2_1200, RNG3_1200, RNG4_1200
%   RNG11, RNG22, RNG33, RNG44
%   RNG11_1200, RNG22_1200, RNG33_1200, RNG44_1200
%   RNG11_, RNG22_, RNG33_, RNG44_
%   RNG11_1200_, RNG22_1200_, RNG33_1200_, RNG44_1200_
%   Depth1200, Depth1200std
%   Depth1200series_
%   stdSpatial_1200, stdSpatialALL_1200
%   DepthRNGall, DepthRNGtable, DepthRNGout

%% Parse inputs

p = inputParser;

addRequired(p, 'dataADCP1200', @(x) isstruct(x));
addRequired(p, 'ADCPin', @(x) isstruct(x));

addParameter(p, 'DespikeT', 2, @(x) isnumeric(x) && isscalar(x));

% Plot options
addParameter(p, 'PlotResults', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotRows', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'PlotRaw', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotDespiked', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotSpatialStd', true, @(x) islogical(x) || isnumeric(x));

% Workspace options
addParameter(p, 'AssignToCaller', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'AssignToBase', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, dataADCP1200, ADCPin, varargin{:});
opt = p.Results;

%% Checks

if isempty(dataADCP1200)
    error('dataADCP1200 is empty. Check filePattern1200 and dir(filePattern1200).');
end

if ~isfield(ADCPin, 'ADCPfilename')
    error('ADCPin.ADCPfilename is missing.');
end

%% Variables from ADCPin
% Only ADCPfilename, StartADCP, and EndADCP are used for depth processing.
% Other variables are copied only to preserve row matching for later comparison.

ADCPfilename = string(ADCPin.ADCPfilename(:));
nFiles = numel(ADCPfilename);

InputRow = (1:nFiles)';

StartADCP = getVectorField(ADCPin, 'StartADCP', zeros(nFiles,1), nFiles);
EndADCP   = getVectorField(ADCPin, 'EndADCP',   zeros(nFiles,1), nFiles);

DistanceFromLeftBank = getVectorField(ADCPin, 'DistanceFromLeftBank', nan(nFiles,1), nFiles);
TrateMean            = getVectorField(ADCPin, 'TrateMean',            nan(nFiles,1), nFiles);

D50mobile  = getVectorField(ADCPin, 'D50mobile',  nan(nFiles,1), nFiles);
D50grab    = getVectorField(ADCPin, 'D50grab',    nan(nFiles,1), nFiles);
D50        = getVectorField(ADCPin, 'D50',        nan(nFiles,1), nFiles);
D90        = getVectorField(ADCPin, 'D90',        nan(nFiles,1), nFiles);
DmeaNUMpsd = getVectorField(ADCPin, 'DmeaNUMpsd', nan(nFiles,1), nFiles);
SSC        = getVectorField(ADCPin, 'SSC',        nan(nFiles,1), nFiles);

%% Initial in-memory table

DepthRNGtable = table( ...
    InputRow, ...
    ADCPfilename, ...
    StartADCP, ...
    EndADCP, ...
    DistanceFromLeftBank, ...
    TrateMean, ...
    D50mobile, ...
    D50grab, ...
    D50, ...
    D90, ...
    DmeaNUMpsd, ...
    SSC, ...
    'VariableNames', { ...
    'InputRow', ...
    'ADCPfilename', ...
    'StartADCP', ...
    'EndADCP', ...
    'DistanceFromLeftBank', ...
    'TrateMean', ...
    'D50mobile', ...
    'D50grab', ...
    'D50', ...
    'D90', ...
    'DmeaNUMpsd', ...
    'SSC'});

%% Preallocate old variables

WF_12 = nan(nFiles,1);
WS_12 = nan(nFiles,1);
WT_12 = nan(nFiles,1);

Temp1 = NaN;

% Mean despiked range/depth per beam, one value per file
RNG1_1200 = nan(nFiles,1);
RNG2_1200 = nan(nFiles,1);
RNG3_1200 = nan(nFiles,1);
RNG4_1200 = nan(nFiles,1);

% Mean and standard deviation per beam, one row per file
Depth1200 = nan(nFiles,4);
Depth1200std = nan(nFiles,4);

% Spatial standard deviation across beams
stdSpatial_1200 = nan(nFiles,1);
stdSpatialALL_1200 = cell(nFiles,1);

% Full beam vectors per file - clear names
RNG11_all = cell(nFiles,1);
RNG22_all = cell(nFiles,1);
RNG33_all = cell(nFiles,1);
RNG44_all = cell(nFiles,1);

RNG11_1200_all = cell(nFiles,1);
RNG22_1200_all = cell(nFiles,1);
RNG33_1200_all = cell(nFiles,1);
RNG44_1200_all = cell(nFiles,1);

% Old-style cell-array names
% Access by:
%   RNG11_{i}, RNG22_{i}, RNG33_{i}, RNG44_{i}
%   RNG11_1200_{i}, RNG22_1200_{i}, RNG33_1200_{i}, RNG44_1200_{i}
RNG11_ = cell(nFiles,1);
RNG22_ = cell(nFiles,1);
RNG33_ = cell(nFiles,1);
RNG44_ = cell(nFiles,1);

RNG11_1200_ = cell(nFiles,1);
RNG22_1200_ = cell(nFiles,1);
RNG33_1200_ = cell(nFiles,1);
RNG44_1200_ = cell(nFiles,1);

% Full cropped four-beam series per file
% Depth1200series_{i} = [RNG11, RNG22, RNG33, RNG44]
Depth1200series_ = cell(nFiles,1);

FileFound = false(nFiles,1);
MatchedFile = strings(nFiles,1);
Message = strings(nFiles,1);

startAll = nan(nFiles,1);
finishAll = nan(nFiles,1);

DepthRNGall = repmat(struct( ...
    'InputRow', [], ...
    'ADCPfilename', "", ...
    'matchedFile', "", ...
    'fileFound', false, ...
    'start', [], ...
    'finish', [], ...
    'StartADCP', [], ...
    'EndADCP', [], ...
    'WF_12', [], ...
    'WS_12', [], ...
    'WT_12', [], ...
    'Temp1', [], ...
    'RNG1_1200', [], ...
    'RNG2_1200', [], ...
    'RNG3_1200', [], ...
    'RNG4_1200', [], ...
    'RNG11', [], ...
    'RNG22', [], ...
    'RNG33', [], ...
    'RNG44', [], ...
    'RNG11_1200', [], ...
    'RNG22_1200', [], ...
    'RNG33_1200', [], ...
    'RNG44_1200', [], ...
    'Depth1200', [], ...
    'Depth1200std', [], ...
    'Depth1200series', [], ...
    'stdSpatial_1200', [], ...
    'stdSpatialALL_1200', [], ...
    'message', ""), nFiles, 1);

%% Process files row by row

for i = 1:nFiles

    thisADCPfilename = ADCPfilename(i);

    matchedFile = matchADCPfile(dataADCP1200, thisADCPfilename);

    DepthRNGall(i).InputRow = InputRow(i);
    DepthRNGall(i).ADCPfilename = thisADCPfilename;
    DepthRNGall(i).StartADCP = StartADCP(i);
    DepthRNGall(i).EndADCP = EndADCP(i);

    if matchedFile == ""

        FileFound(i) = false;
        Message(i) = "File not found";

        DepthRNGall(i).fileFound = false;
        DepthRNGall(i).message = Message(i);

        if opt.Verbose
            warning('No matching file found for row %d: %s', i, thisADCPfilename);
        end

        continue
    end

    FileFound(i) = true;
    MatchedFile(i) = matchedFile;

    try
        S = load(matchedFile);

        if ~isfield(S, 'BT')
            error('BT structure is missing in file: %s', matchedFile);
        end

        if ~isfield(S, 'Data')
            error('Data structure is missing in file: %s', matchedFile);
        end

        BT = S.BT;
        Data = S.Data;

        N = getBTlength(BT);

        startRemove = round(StartADCP(i));
        endRemove   = round(EndADCP(i));

        if isnan(startRemove) || startRemove < 0
            startRemove = 0;
        end

        if isnan(endRemove) || endRemove < 0
            endRemove = 0;
        end

        start = startRemove + 1;
        finish = N - endRemove;

        if finish < start
            error('After trimming StartADCP/EndADCP, no samples remain. Row %d.', i);
        end

        startAll(i) = start;
        finishAll(i) = finish;

        if opt.Verbose
            fprintf('\nProcessing row %d/%d\n', i, nFiles);
            fprintf('ADCPfilename: %s\n', thisADCPfilename);
            fprintf('Matched file: %s\n', matchedFile);
            fprintf('N samples: %d | using samples %d:%d\n', N, start, finish);
        end

        %% Crop Data fields using the same start/finish

        WF = getDataFieldOptional(Data, 'WF');
        WS = getDataFieldOptional(Data, 'WS');
        WT = getDataFieldOptional(Data, 'WT');

        WF = cropVectorSafe(WF, start, finish);
        WS = cropVectorSafe(WS, start, finish);
        WT = cropVectorSafe(WT, start, finish);

        WF_12(i) = safeMean(WF) ./ 100;
        WS_12(i) = safeMean(WS) ./ 100;

        try
            WT_12(i) = safeMean(WT);
            Temp1 = WT_12(i);
        catch
            Temp1 = Temp1;
        end

        %% Crop BT.RNG1-4 using same start/finish

        RNG11 = getBTsegment(BT, 'RNG1', start, finish);
        RNG22 = getBTsegment(BT, 'RNG2', start, finish);
        RNG33 = getBTsegment(BT, 'RNG3', start, finish);
        RNG44 = getBTsegment(BT, 'RNG4', start, finish);

        RNG11(RNG11 == 0) = NaN;
        RNG22(RNG22 == 0) = NaN;
        RNG33(RNG33 == 0) = NaN;
        RNG44(RNG44 == 0) = NaN;

        %% Despike cropped beam ranges

        RNG11_1200 = despikeBT(RNG11, opt.DespikeT);
        RNG22_1200 = despikeBT(RNG22, opt.DespikeT);
        RNG33_1200 = despikeBT(RNG33, opt.DespikeT);
        RNG44_1200 = despikeBT(RNG44, opt.DespikeT);

        %% Mean despiked range/depth per beam, one value per file

        RNG1_1200(i) = safeMean(RNG11_1200);
        RNG2_1200(i) = safeMean(RNG22_1200);
        RNG3_1200(i) = safeMean(RNG33_1200);
        RNG4_1200(i) = safeMean(RNG44_1200);

        %% Old depth logic

        R = [RNG11(:), RNG22(:), RNG33(:), RNG44(:)];

        Depth1200(i,:) = nanmeanLocal(R, 1);
        Depth1200std(i,:) = nanstdLocal(R, 0, 1);

        stdSpatialALL_1200{i} = nanstdLocal(R, 0, 2);
        stdSpatial_1200(i) = safeMean(stdSpatialALL_1200{i});

        %% Store full cropped series

        RNG11_all{i} = RNG11;
        RNG22_all{i} = RNG22;
        RNG33_all{i} = RNG33;
        RNG44_all{i} = RNG44;

        RNG11_1200_all{i} = RNG11_1200;
        RNG22_1200_all{i} = RNG22_1200;
        RNG33_1200_all{i} = RNG33_1200;
        RNG44_1200_all{i} = RNG44_1200;

        % Same values with old-style cell-array names
        RNG11_{i} = RNG11;
        RNG22_{i} = RNG22;
        RNG33_{i} = RNG33;
        RNG44_{i} = RNG44;

        RNG11_1200_{i} = RNG11_1200;
        RNG22_1200_{i} = RNG22_1200;
        RNG33_1200_{i} = RNG33_1200;
        RNG44_1200_{i} = RNG44_1200;

        Depth1200series_{i} = R;

        %% Optional plotting

        if opt.PlotResults && shouldPlotThisRow(i, opt.PlotRows)

            plotDepthRNGdiagnostics( ...
                i, ...
                thisADCPfilename, ...
                start, ...
                finish, ...
                RNG11, RNG22, RNG33, RNG44, ...
                RNG11_1200, RNG22_1200, RNG33_1200, RNG44_1200, ...
                stdSpatialALL_1200{i}, ...
                opt.PlotRaw, ...
                opt.PlotDespiked, ...
                opt.PlotSpatialStd);
        end

        %% Store row structure

        DepthRNGall(i).matchedFile = matchedFile;
        DepthRNGall(i).fileFound = true;
        DepthRNGall(i).start = start;
        DepthRNGall(i).finish = finish;

        DepthRNGall(i).WF_12 = WF_12(i);
        DepthRNGall(i).WS_12 = WS_12(i);
        DepthRNGall(i).WT_12 = WT_12(i);
        DepthRNGall(i).Temp1 = Temp1;

        DepthRNGall(i).RNG1_1200 = RNG1_1200(i);
        DepthRNGall(i).RNG2_1200 = RNG2_1200(i);
        DepthRNGall(i).RNG3_1200 = RNG3_1200(i);
        DepthRNGall(i).RNG4_1200 = RNG4_1200(i);

        DepthRNGall(i).RNG11 = RNG11;
        DepthRNGall(i).RNG22 = RNG22;
        DepthRNGall(i).RNG33 = RNG33;
        DepthRNGall(i).RNG44 = RNG44;

        DepthRNGall(i).RNG11_1200 = RNG11_1200;
        DepthRNGall(i).RNG22_1200 = RNG22_1200;
        DepthRNGall(i).RNG33_1200 = RNG33_1200;
        DepthRNGall(i).RNG44_1200 = RNG44_1200;

        DepthRNGall(i).Depth1200 = Depth1200(i,:);
        DepthRNGall(i).Depth1200std = Depth1200std(i,:);
        DepthRNGall(i).Depth1200series = Depth1200series_{i};

        DepthRNGall(i).stdSpatial_1200 = stdSpatial_1200(i);
        DepthRNGall(i).stdSpatialALL_1200 = stdSpatialALL_1200{i};

        DepthRNGall(i).message = "OK";
        Message(i) = "OK";

    catch ME

        Message(i) = "Error: " + string(ME.message);

        DepthRNGall(i).matchedFile = matchedFile;
        DepthRNGall(i).fileFound = true;
        DepthRNGall(i).message = Message(i);

        warning('Processing failed for row %d, file %s:\n%s', ...
            i, thisADCPfilename, ME.message);
    end
end

%% Add outputs to table

DepthRNGtable.FileFound = FileFound;
DepthRNGtable.MatchedFile = MatchedFile;
DepthRNGtable.start = startAll;
DepthRNGtable.finish = finishAll;

DepthRNGtable.WF_12 = WF_12;
DepthRNGtable.WS_12 = WS_12;
DepthRNGtable.WT_12 = WT_12;

DepthRNGtable.RNG1_1200 = RNG1_1200;
DepthRNGtable.RNG2_1200 = RNG2_1200;
DepthRNGtable.RNG3_1200 = RNG3_1200;
DepthRNGtable.RNG4_1200 = RNG4_1200;

DepthRNGtable.Depth1200_B1 = Depth1200(:,1);
DepthRNGtable.Depth1200_B2 = Depth1200(:,2);
DepthRNGtable.Depth1200_B3 = Depth1200(:,3);
DepthRNGtable.Depth1200_B4 = Depth1200(:,4);

DepthRNGtable.Depth1200std_B1 = Depth1200std(:,1);
DepthRNGtable.Depth1200std_B2 = Depth1200std(:,2);
DepthRNGtable.Depth1200std_B3 = Depth1200std(:,3);
DepthRNGtable.Depth1200std_B4 = Depth1200std(:,4);

DepthRNGtable.stdSpatial_1200 = stdSpatial_1200;
DepthRNGtable.Message = Message;

%% Last successful row variables

lastOK = find(FileFound & Message == "OK", 1, 'last');

if isempty(lastOK)
    RNG11 = [];
    RNG22 = [];
    RNG33 = [];
    RNG44 = [];

    RNG11_1200 = [];
    RNG22_1200 = [];
    RNG33_1200 = [];
    RNG44_1200 = [];

    start = [];
    finish = [];
else
    RNG11 = RNG11_all{lastOK};
    RNG22 = RNG22_all{lastOK};
    RNG33 = RNG33_all{lastOK};
    RNG44 = RNG44_all{lastOK};

    RNG11_1200 = RNG11_1200_all{lastOK};
    RNG22_1200 = RNG22_1200_all{lastOK};
    RNG33_1200 = RNG33_1200_all{lastOK};
    RNG44_1200 = RNG44_1200_all{lastOK};

    start = startAll(lastOK);
    finish = finishAll(lastOK);
end

%% Output structure

DepthRNGout = struct();

DepthRNGout.InputRow = InputRow;
DepthRNGout.ADCPfilename = ADCPfilename;
DepthRNGout.StartADCP = StartADCP;
DepthRNGout.EndADCP = EndADCP;

DepthRNGout.DistanceFromLeftBank = DistanceFromLeftBank;
DepthRNGout.TrateMean = TrateMean;
DepthRNGout.D50mobile = D50mobile;
DepthRNGout.D50grab = D50grab;
DepthRNGout.D50 = D50;
DepthRNGout.D90 = D90;
DepthRNGout.DmeaNUMpsd = DmeaNUMpsd;
DepthRNGout.SSC = SSC;

DepthRNGout.WF_12 = WF_12;
DepthRNGout.WS_12 = WS_12;
DepthRNGout.WT_12 = WT_12;
DepthRNGout.Temp1 = Temp1;

DepthRNGout.RNG1_1200 = RNG1_1200;
DepthRNGout.RNG2_1200 = RNG2_1200;
DepthRNGout.RNG3_1200 = RNG3_1200;
DepthRNGout.RNG4_1200 = RNG4_1200;

DepthRNGout.RNG11 = RNG11;
DepthRNGout.RNG22 = RNG22;
DepthRNGout.RNG33 = RNG33;
DepthRNGout.RNG44 = RNG44;

DepthRNGout.RNG11_1200 = RNG11_1200;
DepthRNGout.RNG22_1200 = RNG22_1200;
DepthRNGout.RNG33_1200 = RNG33_1200;
DepthRNGout.RNG44_1200 = RNG44_1200;

DepthRNGout.RNG11_all = RNG11_all;
DepthRNGout.RNG22_all = RNG22_all;
DepthRNGout.RNG33_all = RNG33_all;
DepthRNGout.RNG44_all = RNG44_all;

DepthRNGout.RNG11_1200_all = RNG11_1200_all;
DepthRNGout.RNG22_1200_all = RNG22_1200_all;
DepthRNGout.RNG33_1200_all = RNG33_1200_all;
DepthRNGout.RNG44_1200_all = RNG44_1200_all;

% Old-style cell-array names
DepthRNGout.RNG11_ = RNG11_;
DepthRNGout.RNG22_ = RNG22_;
DepthRNGout.RNG33_ = RNG33_;
DepthRNGout.RNG44_ = RNG44_;

DepthRNGout.RNG11_1200_ = RNG11_1200_;
DepthRNGout.RNG22_1200_ = RNG22_1200_;
DepthRNGout.RNG33_1200_ = RNG33_1200_;
DepthRNGout.RNG44_1200_ = RNG44_1200_;

DepthRNGout.Depth1200series_ = Depth1200series_;

DepthRNGout.Depth1200 = Depth1200;
DepthRNGout.Depth1200std = Depth1200std;
DepthRNGout.stdSpatial_1200 = stdSpatial_1200;
DepthRNGout.stdSpatialALL_1200 = stdSpatialALL_1200;

DepthRNGout.DepthRNGall = DepthRNGall;
DepthRNGout.DepthRNGtable = DepthRNGtable;

DepthRNGout.start = start;
DepthRNGout.finish = finish;
DepthRNGout.startAll = startAll;
DepthRNGout.finishAll = finishAll;

%% Assign to workspace

if opt.AssignToCaller
    assignStructToWorkspace(DepthRNGout, 'caller');
end

if opt.AssignToBase
    assignStructToWorkspace(DepthRNGout, 'base');
end

%% Display

if opt.Verbose
    fprintf('\n=== Depth/RNG-only processing summary ===\n');
    fprintf('Input rows: %d\n', nFiles);
    fprintf('Matched files: %d/%d\n', sum(FileFound), nFiles);
    fprintf('Successfully processed: %d/%d\n', sum(Message == "OK"), nFiles);
    fprintf('StartADCP/EndADCP were used as trim counts.\n');
    fprintf('No files were saved.\n');
    fprintf('=========================================\n');

    disp(DepthRNGtable);
end

end

%% Helper functions

function x = getVectorField(S, fieldName, defaultValue, n)

if isfield(S, fieldName)
    x = S.(fieldName);
else
    x = defaultValue;
end

if isstring(x) || iscell(x)
    x = str2double(string(x));
end

x = double(x);
x = x(:);

if numel(x) == 1 && n > 1
    x = repmat(x, n, 1);
elseif numel(x) < n
    x = [x; nan(n - numel(x), 1)];
elseif numel(x) > n
    x = x(1:n);
end

end


function matchedFile = matchADCPfile(dataADCP1200, thisADCPfilename)

matchedFile = "";

targetName = cleanFileNameFromTable(thisADCPfilename);

if targetName == ""
    return
end

fileNames = string({dataADCP1200.name}');
cleanNames = strings(numel(fileNames),1);

for i = 1:numel(fileNames)
    cleanNames(i) = cleanFileNameFromTable(fileNames(i));
end

[~, targetBase, ~] = fileparts(char(targetName));

fileBases = strings(numel(fileNames),1);

for i = 1:numel(fileNames)
    [~, b, ~] = fileparts(char(cleanNames(i)));
    fileBases(i) = string(b);
end

idx = find(strcmpi(cleanNames, targetName), 1, 'first');

if isempty(idx)
    idx = find(strcmpi(fileBases, string(targetBase)), 1, 'first');
end

if isempty(idx)
    idx = find(startsWith(lower(fileBases), lower(string(targetBase))), 1, 'first');
end

if isempty(idx)
    idx = find(contains(lower(fileBases), lower(string(targetBase))), 1, 'first');
end

if isempty(idx)
    return
end

if isfield(dataADCP1200, 'folder') && ~isempty(dataADCP1200(idx).folder)
    matchedFile = string(fullfile(dataADCP1200(idx).folder, dataADCP1200(idx).name));
else
    matchedFile = string(dataADCP1200(idx).name);
end

end


function s = cleanFileNameFromTable(s)

s = string(s);
s = strtrim(s);

if ismissing(s) || s == "" || lower(s) == "nan"
    s = "";
    return
end

s = erase(s, "'");
s = erase(s, '"');
s = strtrim(s);

end


function v = getDataFieldOptional(Data, fieldName)

if isfield(Data, fieldName)
    v = double(Data.(fieldName));
else
    v = NaN;
end

end


function n = getBTlength(BT)

lens = [];

if isfield(BT, 'RNG1')
    lens(end+1) = numel(BT.RNG1); %#ok<AGROW>
end

if isfield(BT, 'RNG2')
    lens(end+1) = numel(BT.RNG2); %#ok<AGROW>
end

if isfield(BT, 'RNG3')
    lens(end+1) = numel(BT.RNG3); %#ok<AGROW>
end

if isfield(BT, 'RNG4')
    lens(end+1) = numel(BT.RNG4); %#ok<AGROW>
end

if isempty(lens)
    error('BT.RNG1-RNG4 are missing.');
end

n = min(lens);

end


function x = getBTsegment(BT, fieldName, start, finish)

if ~isfield(BT, fieldName)
    error('BT.%s is missing.', fieldName);
end

x = double(BT.(fieldName));
x = x(:);

finish = min(finish, numel(x));

if start > numel(x)
    warning('start exceeds length of BT.%s. Returning NaN.', fieldName);
    x = NaN;
    return
end

x = x(start:finish);

end


function x = cropVectorSafe(x, start, finish)

x = double(x);
x = x(:);

if isempty(x) || all(isnan(x))
    return
end

% If variable is scalar, keep it unchanged.
% This avoids turning scalar Data.WF / Data.WS / Data.WT into NaN.
if numel(x) == 1
    return
end

finish = min(finish, numel(x));

if start > numel(x)
    x = NaN;
else
    x = x(start:finish);
end

end


function BT = despikeBT(BTv,T)
%DESPIKEBT Simple function for despiking BT ADCP data.
%
% Based on the original script logic:
%   T1 = mean(abs(A)) - T * std(abs(A))
%   To = mean(abs(A)) + T * std(abs(A))
%   values outside [T1, To] become NaN.

A = double(BTv);

if isempty(A)
    BT = A;
    return
end

valid = ~isnan(A);

if sum(valid) < 3
    BT = A;
    return
end

T1 = mean(abs(A(valid))) - T .* std(abs(A(valid)));
To = mean(abs(A(valid))) + T .* std(abs(A(valid)));

peaks = (A > To | A < T1);
A(peaks) = NaN;

BT = A;

end


function y = safeMean(x)

if isempty(x)
    y = NaN;
    return
end

x = double(x(:));
x = x(~isnan(x));

if isempty(x)
    y = NaN;
else
    y = mean(x);
end

end


function y = nanmeanLocal(A, dim)

if nargin < 2
    dim = 1;
end

A = double(A);

valid = ~isnan(A);
A(~valid) = 0;

n = sum(valid, dim);
y = sum(A, dim) ./ n;
y(n == 0) = NaN;

end


function y = nanstdLocal(A, flag, dim)

if nargin < 2 || isempty(flag)
    flag = 0;
end

if nargin < 3
    dim = 1;
end

A = double(A);

mu = nanmeanLocal(A, dim);

if dim == 1
    muRep = repmat(mu, size(A,1), 1);
else
    muRep = repmat(mu, 1, size(A,2));
end

valid = ~isnan(A);
diff2 = (A - muRep).^2;
diff2(~valid) = 0;

n = sum(valid, dim);

if flag == 0
    denom = n - 1;
else
    denom = n;
end

denom(denom <= 0) = NaN;

y = sqrt(sum(diff2, dim) ./ denom);

end


function tf = shouldPlotThisRow(i, PlotRows)
%SHOULDPLOTTHISROW Decide if current row should be plotted.

if isempty(PlotRows)
    tf = true;
else
    tf = any(i == PlotRows);
end

end


function plotDepthRNGdiagnostics( ...
    i, ...
    thisADCPfilename, ...
    start, ...
    finish, ...
    RNG11, RNG22, RNG33, RNG44, ...
    RNG11_1200, RNG22_1200, RNG33_1200, RNG44_1200, ...
    stdSpatialALL_1200_i, ...
    PlotRaw, ...
    PlotDespiked, ...
    PlotSpatialStd)
%PLOTDEPTHRNGDIAGNOSTICS Plot raw/despiked BT RNG depth series.
%
% This function only plots. It does not save anything.

nPanels = double(logical(PlotRaw)) + ...
          double(logical(PlotDespiked)) + ...
          double(logical(PlotSpatialStd));

if nPanels == 0
    return
end

figure('Name', sprintf('Depth/RNG diagnostics row %d', i), 'Color', 'w');

panel = 0;

x = (start:finish)';

if numel(x) ~= numel(RNG11)
    x = (1:numel(RNG11))';
end

if PlotRaw

    panel = panel + 1;
    subplot(nPanels, 1, panel);

    plot(x, RNG11, 'DisplayName', 'RNG11 raw'); hold on
    plot(x, RNG22, 'DisplayName', 'RNG22 raw');
    plot(x, RNG33, 'DisplayName', 'RNG33 raw');
    plot(x, RNG44, 'DisplayName', 'RNG44 raw');

    grid on
    xlabel('Ensemble/sample index');
    ylabel('RNG/depth');
    title(sprintf('Raw BT RNG/depths — row %d: %s', ...
        i, cleanTitleString(thisADCPfilename)), ...
        'Interpreter', 'none');

    legend('Location', 'best');
end

if PlotDespiked

    panel = panel + 1;
    subplot(nPanels, 1, panel);

    plot(x, RNG11_1200, 'DisplayName', 'RNG11\_1200'); hold on
    plot(x, RNG22_1200, 'DisplayName', 'RNG22\_1200');
    plot(x, RNG33_1200, 'DisplayName', 'RNG33\_1200');
    plot(x, RNG44_1200, 'DisplayName', 'RNG44\_1200');

    grid on
    xlabel('Ensemble/sample index');
    ylabel('Despiked RNG/depth');
    title('Despiked BT RNG/depths', 'Interpreter', 'none');

    legend('Location', 'best');
end

if PlotSpatialStd

    panel = panel + 1;
    subplot(nPanels, 1, panel);

    plot(x, stdSpatialALL_1200_i, 'DisplayName', 'stdSpatialALL\_1200');

    grid on
    xlabel('Ensemble/sample index');
    ylabel('Std across beams');
    title('Spatial standard deviation across four beams', 'Interpreter', 'none');

    legend('Location', 'best');
end

end


function s = cleanTitleString(s)
%CLEANTITLESTRING Clean string for figure title.

s = string(s);
s = erase(s, "'");
s = erase(s, '"');
s = char(s);

end


function assignStructToWorkspace(S, workspaceName)

f = fieldnames(S);

for i = 1:numel(f)
    assignin(workspaceName, f{i}, S.(f{i}));
end

end