function WaterProfileOut = processWaterProfile(dataADCP1200, ADCPin, varargin)
%PROCESSWATERPROFILE1200 Process water velocity profiles from RDI/RiverPro files.
%
% This function:
%   - matches ADCPin.ADCPfilename with dataADCP1200.name
%   - uses StartADCP(i) and EndADCP(i) as trim counts:
%
%       start  = StartADCP(i) + 1;
%       finish = N - EndADCP(i);
%
%   - crops water-profile velocities:
%
%       V1x_ = v.v1(:,start:finish);
%       V1y_ = v.v2(:,start:finish);
%       V1w_ = v.v3(:,start:finish);
%       V1e_ = v.v4(:,start:finish);
%
%   - converts WF and WS from cm to m:
%
%       Wf_1200(i) = mean(Data.WF)/100;
%       Ws_1200(i) = mean(Data.WS)/100;
%
%   - calculates physical cell-centre coordinate:
%
%       H_1200 = distanceSonTek + Wf_1200(i) + 2/4*Ws_1200(i) ...
%                + [0:nCells-1].*Ws_1200(i);
%
%       z_1200 = Depth_1200(i) - H_1200;
%
%   - keeps only z_1200 > 0
%   - extracts the lower profile using LowerFraction
%   - stores variable-length outputs in cells
%   - assigns all outputs to workspace if requested
%
% Recommended call:
%
%   WaterProfileOut = processWaterProfile1200(dataADCP1200, ADCPin, ...
%       'DepthRNGout', DepthRNGout, ...
%       'distanceSonTek', 0, ...
%       'LowerFraction', 0.50, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
% Important output variables:
%
%   V1x_, V1y_, V1w_, V1e_                     % last successful file
%   V1x__{i}, V1y__{i}, V1w__{i}, V1e__{i}     % all files
%
%   V1x_1, V1y_1, V1w_1, V1e_1
%   V1x_1_{i}, V1y_1_{i}, V1w_1_{i}, V1e_1_{i}
%
%   V1x_1200, V1y_1200
%   Vmod_1200_1
%   Vmod_1200_measured
%   Vmod_1200
%   Vmod1_1200
%   VmodLower_1200
%
%   H_1200, z_1200, z1_1200
%   H_1200_{i}, z_1200_{i}, z1_1200_{i}
%
%   VlastCell_1200
%   Vlast2Cells_1200
%   VlastCell_1200_{i}
%   Vlast2Cells_1200_{i}
%
% No files are saved by this function.

%% Parse inputs

p = inputParser;

addRequired(p, 'dataADCP1200', @(x) isstruct(x));
addRequired(p, 'ADCPin', @(x) isstruct(x));

% Inputs from previous depth processing
addParameter(p, 'DepthRNGout', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'Depth_1200', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Depth1200', [], @(x) isempty(x) || isnumeric(x));

% Geometry
addParameter(p, 'distanceSonTek', 0, @(x) isnumeric(x) && isscalar(x));

% Lower profile selection
% LowerFraction = fraction of the valid profile closest to the bed.
% Example:
%   LowerFraction = 0.50 means use approximately the lower half.
%   LowerFraction = 0.30 means use approximately the lower 30%.
addParameter(p, 'LowerFraction', 0.50, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);

% Backward-compatible alias.
% If logLawpartition is provided, it is used as the start fraction:
%   z1_1200 = z_1200(round(logLawpartition*length(z_1200)):end)
% Otherwise, start fraction is calculated as:
%   lowerStartFraction = 1 - LowerFraction
addParameter(p, 'logLawpartition', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && x <= 1));

% Extrapolation
addParameter(p, 'ApplyTopExtrapolation', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'NTopCellsExtrap', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'TopExtraOffset', 0.00, @(x) isnumeric(x) && isscalar(x));

% Velocity filtering
addParameter(p, 'ApplyVelocityFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'MaxAbsVelocity', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'DespikeVelocity', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'DespikeT', 3, @(x) isnumeric(x) && isscalar(x));

% Plot options
addParameter(p, 'PlotResults', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotRows', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'PlotMeasuredProfile', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotExtrapolatedProfile', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotLastCellSeries', true, @(x) islogical(x) || isnumeric(x));

% Workspace/output
addParameter(p, 'AssignToCaller', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'AssignToBase', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, dataADCP1200, ADCPin, varargin{:});
opt = p.Results;

% Decide lower-profile start fraction
if isempty(opt.logLawpartition)
    lowerStartFraction = 1 - opt.LowerFraction;
else
    lowerStartFraction = opt.logLawpartition;
end

lowerStartFraction = max(0, min(1, lowerStartFraction));

%% Checks

if isempty(dataADCP1200)
    error('dataADCP1200 is empty. Check filePattern1200 and dir(filePattern1200).');
end

if ~isfield(ADCPin, 'ADCPfilename')
    error('ADCPin.ADCPfilename is missing.');
end

%% Variables from ADCPin

ADCPfilename = string(ADCPin.ADCPfilename(:));
nFiles = numel(ADCPfilename);

InputRow = (1:nFiles)';

StartADCP = getVectorField(ADCPin, 'StartADCP', zeros(nFiles,1), nFiles);
EndADCP   = getVectorField(ADCPin, 'EndADCP',   zeros(nFiles,1), nFiles);

DistanceFromLeftBank = getVectorField(ADCPin, 'DistanceFromLeftBank', nan(nFiles,1), nFiles);
TrateMean            = getVectorField(ADCPin, 'TrateMean',            nan(nFiles,1), nFiles);

D50        = getVectorField(ADCPin, 'D50',        nan(nFiles,1), nFiles);
D90        = getVectorField(ADCPin, 'D90',        nan(nFiles,1), nFiles);
DmeaNUMpsd = getVectorField(ADCPin, 'DmeaNUMpsd', nan(nFiles,1), nFiles);
SSC        = getVectorField(ADCPin, 'SSC',        nan(nFiles,1), nFiles);

% Try to get Depth_1200 from previous processing or user input.
Depth_1200 = getDepth1200Input(opt, ADCPin, nFiles);

%% Initial output table

WaterProfileTable = table( ...
    InputRow, ...
    ADCPfilename, ...
    StartADCP, ...
    EndADCP, ...
    DistanceFromLeftBank, ...
    TrateMean, ...
    D50, ...
    D90, ...
    DmeaNUMpsd, ...
    SSC, ...
    Depth_1200, ...
    'VariableNames', { ...
    'InputRow', ...
    'ADCPfilename', ...
    'StartADCP', ...
    'EndADCP', ...
    'DistanceFromLeftBank', ...
    'TrateMean', ...
    'D50', ...
    'D90', ...
    'DmeaNUMpsd', ...
    'SSC', ...
    'Depth_1200'});

%% Preallocate file-level variables

FileFound = false(nFiles,1);
MatchedFile = strings(nFiles,1);
Message = strings(nFiles,1);

startAll = nan(nFiles,1);
finishAll = nan(nFiles,1);

WF_12 = nan(nFiles,1);
WS_12 = nan(nFiles,1);
WT_12 = nan(nFiles,1);

% Old-style aliases in meters
Wf_1200 = nan(nFiles,1);
Ws_1200 = nan(nFiles,1);
Temp1 = NaN;

Z_upValue_all = nan(nFiles,1);

Vmod_1200stdMOD = nan(nFiles,1);
VmodAll_1200 = nan(nFiles,1);
Vmod1low_1200stdMOD = nan(nFiles,1);

VlastCellMean_1200 = nan(nFiles,1);
Vlast2CellsMean_1200 = nan(nFiles,1);

VlastCellMeanX_1200 = nan(nFiles,1);
VlastCellMeanY_1200 = nan(nFiles,1);

Vlast2CellsMeanX_1200 = nan(nFiles,1);
Vlast2CellsMeanY_1200 = nan(nFiles,1);

nCellsUsed_1200 = nan(nFiles,1);
nEnsUsed_1200 = nan(nFiles,1);

LowerFraction_used = repmat(opt.LowerFraction, nFiles, 1);
lowerStartFraction_used = repmat(lowerStartFraction, nFiles, 1);

%% Preallocate cell outputs

V1x__ = cell(nFiles,1);
V1y__ = cell(nFiles,1);
V1w__ = cell(nFiles,1);
V1e__ = cell(nFiles,1);

V1x_1_ = cell(nFiles,1);
V1y_1_ = cell(nFiles,1);
V1w_1_ = cell(nFiles,1);
V1e_1_ = cell(nFiles,1);

V1x_1200_ = cell(nFiles,1);
V1y_1200_ = cell(nFiles,1);
V1w_1200_ = cell(nFiles,1);
V1e_1200_ = cell(nFiles,1);

Vmod_1200_ = cell(nFiles,1);
Vmod_1200_1_ = cell(nFiles,1);
Vmod_1200_measured_ = cell(nFiles,1);
Vmod1_1200_ = cell(nFiles,1);
VmodLower_1200_ = cell(nFiles,1);
Vmod_1200std_ = cell(nFiles,1);

VlastCell_1200_ = cell(nFiles,1);
Vlast2Cells_1200_ = cell(nFiles,1);

VlastCellX_1200_ = cell(nFiles,1);
VlastCellY_1200_ = cell(nFiles,1);

Vlast2CellsX_1200_ = cell(nFiles,1);
Vlast2CellsY_1200_ = cell(nFiles,1);

H_1200_ = cell(nFiles,1);
z_1200_ = cell(nFiles,1);
z1_1200_ = cell(nFiles,1);

V1200extrap_ = cell(nFiles,1);
xfit_ = cell(nFiles,1);
coeffs_ = cell(nFiles,1);

WaterProfileAll = repmat(struct( ...
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
    'Wf_1200', [], ...
    'Ws_1200', [], ...
    'Temp1', [], ...
    'Depth_1200', [], ...
    'LowerFraction', [], ...
    'lowerStartFraction', [], ...
    'H_1200', [], ...
    'z_1200', [], ...
    'z1_1200', [], ...
    'V1x_', [], ...
    'V1y_', [], ...
    'V1w_', [], ...
    'V1e_', [], ...
    'V1x_1', [], ...
    'V1y_1', [], ...
    'V1w_1', [], ...
    'V1e_1', [], ...
    'V1x_1200', [], ...
    'V1y_1200', [], ...
    'V1w_1200', [], ...
    'V1e_1200', [], ...
    'Vmod_1200_1', [], ...
    'Vmod_1200_measured', [], ...
    'Vmod_1200', [], ...
    'Vmod1_1200', [], ...
    'VmodLower_1200', [], ...
    'Z_upValue', [], ...
    'coeffs', [], ...
    'xfit', [], ...
    'V1200extrap', [], ...
    'Vmod_1200std', [], ...
    'Vmod_1200stdMOD', [], ...
    'VmodAll_1200', [], ...
    'Vmod1low_1200stdMOD', [], ...
    'V1x_1std', [], ...
    'V1y_1std', [], ...
    'VlastCell_1200', [], ...
    'Vlast2Cells_1200', [], ...
    'VlastCellMean_1200', [], ...
    'Vlast2CellsMean_1200', [], ...
    'message', ""), nFiles, 1);

%% Process files

for i = 1:nFiles

    thisADCPfilename = ADCPfilename(i);
    matchedFile = matchADCPfile(dataADCP1200, thisADCPfilename);

    WaterProfileAll(i).InputRow = InputRow(i);
    WaterProfileAll(i).ADCPfilename = thisADCPfilename;
    WaterProfileAll(i).StartADCP = StartADCP(i);
    WaterProfileAll(i).EndADCP = EndADCP(i);

    if matchedFile == ""

        FileFound(i) = false;
        Message(i) = "File not found";

        WaterProfileAll(i).fileFound = false;
        WaterProfileAll(i).message = Message(i);

        if opt.Verbose
            warning('No matching file found for row %d: %s', i, thisADCPfilename);
        end

        continue
    end

    FileFound(i) = true;
    MatchedFile(i) = matchedFile;

    try

        S = load(matchedFile);

        if ~isfield(S, 'v')
            error('Velocity structure v is missing in file: %s', matchedFile);
        end

        if ~isfield(S, 'Data')
            error('Data structure is missing in file: %s', matchedFile);
        end

        v = S.v;
        Data = S.Data;

        checkVelocityFields(v);

        N = getVelocityLength(v);

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
            fprintf('\nProcessing water profile row %d/%d\n', i, nFiles);
            fprintf('ADCPfilename: %s\n', thisADCPfilename);
            fprintf('Matched file: %s\n', matchedFile);
            fprintf('N ensembles: %d | using ensembles %d:%d\n', N, start, finish);
        end

        %% Crop velocity matrices

        V1x_ = double(v.v1(:, start:finish));
        V1y_ = double(v.v2(:, start:finish));
        V1w_ = double(v.v3(:, start:finish));
        V1e_ = double(v.v4(:, start:finish));

        if opt.ApplyVelocityFilter
            V1x_ = filterVelocityMatrix(V1x_, opt.MaxAbsVelocity);
            V1y_ = filterVelocityMatrix(V1y_, opt.MaxAbsVelocity);
            V1w_ = filterVelocityMatrix(V1w_, opt.MaxAbsVelocity);
            V1e_ = filterVelocityMatrix(V1e_, opt.MaxAbsVelocity);
        end

        if opt.DespikeVelocity
            V1x_ = despikeVelocityMatrix(V1x_, opt.DespikeT);
            V1y_ = despikeVelocityMatrix(V1y_, opt.DespikeT);
            V1w_ = despikeVelocityMatrix(V1w_, opt.DespikeT);
            V1e_ = despikeVelocityMatrix(V1e_, opt.DespikeT);
        end

        %% Extract/crop instrument geometry fields

        WF = getDataFieldOptional(Data, 'WF');
        WS = getDataFieldOptional(Data, 'WS');
        WT = getDataFieldOptional(Data, 'WT');

        WFcrop = cropVectorSafe(WF, start, finish);
        WScrop = cropVectorSafe(WS, start, finish);
        WTcrop = cropVectorSafe(WT, start, finish);

        % Data.WF and Data.WS are in cm, so divide by 100 to get meters.
        WF_12(i) = safeMean(WFcrop) ./ 100;
        WS_12(i) = safeMean(WScrop) ./ 100;

        Wf_1200(i) = WF_12(i);
        Ws_1200(i) = WS_12(i);

        try
            WT_12(i) = safeMean(WTcrop);
            Temp1 = WT_12(i);
        catch
            Temp1 = Temp1;
        end

        if isnan(Wf_1200(i))
            Wf_1200(i) = getFirstExistingScalar(Data, {'Wf_1200','WF','BlankingDistance','blankingDistance'}) ./ 100;
        end

        if isnan(Ws_1200(i))
            Ws_1200(i) = getFirstExistingScalar(Data, {'Ws_1200','WS','CellSize','cellSize','BinSize','binSize'}) ./ 100;
        end

        if isnan(Wf_1200(i))
            error('Could not determine Wf_1200 / blanking distance for row %d.', i);
        end

        if isnan(Ws_1200(i)) || Ws_1200(i) <= 0
            error('Could not determine Ws_1200 / cell size for row %d.', i);
        end

        %% Mean water profiles before z filtering

        V1x_1200 = nanmeanLocal(V1x_, 2);
        V1y_1200 = nanmeanLocal(V1y_, 2);
        V1w_1200 = nanmeanLocal(V1w_, 2);
        V1e_1200 = nanmeanLocal(V1e_, 2);

        Vmod_1200_1 = sqrt(V1x_1200.^2 + V1y_1200.^2);
        Vmod_1200_measured = nanmeanLocal(sqrt(V1x_.^2 + V1y_.^2), 2);

        %% Determine Depth_1200(i)

        if isnan(Depth_1200(i))
            Depth_1200(i) = estimateDepth1200FromBT(S, start, finish);
        end

        if isnan(Depth_1200(i)) || Depth_1200(i) <= 0
            error('Could not determine valid Depth_1200 for row %d.', i);
        end

        %% Physical cell-center coordinates

        nCells = size(V1x_1200,1);

        H_1200 = opt.distanceSonTek + Wf_1200(i) + 2/4 .* Ws_1200(i) + ...
            (0:nCells-1)' .* Ws_1200(i);

        H_1200 = H_1200(:);

        % z_1200 is distance above bed, not a cell index.
        z_1200 = Depth_1200(i) - H_1200;

        % Keep only valid cells inside the water depth and with valid velocity.
        validZ = z_1200 > 0;
        validV = ~isnan(Vmod_1200_measured);

        validCells = validZ & validV;

        if ~any(validCells)
            error('No valid water-profile cells remain after z_1200 > 0 filtering. Row %d.', i);
        end

        H_1200 = H_1200(validCells);
        z_1200 = z_1200(validCells);

        V1x_ = V1x_(validCells,:);
        V1y_ = V1y_(validCells,:);
        V1w_ = V1w_(validCells,:);
        V1e_ = V1e_(validCells,:);

        V1x_1200 = V1x_1200(validCells);
        V1y_1200 = V1y_1200(validCells);
        V1w_1200 = V1w_1200(validCells);
        V1e_1200 = V1e_1200(validCells);

        Vmod_1200_1 = Vmod_1200_1(validCells);
        Vmod_1200_measured = Vmod_1200_measured(validCells);

        %% Lower/log-law portion

        idxLowerStart = max(1, round(lowerStartFraction .* length(z_1200)));

        z1_1200 = z_1200(idxLowerStart:end);
        Vmod1_1200 = Vmod_1200_measured(idxLowerStart:end);

        VmodLower_1200 = Vmod1_1200;

        %% Upper extrapolation

        % Distance from instrument/reference to first valid cell centre.
        Z_upValue = opt.distanceSonTek + Wf_1200(i) + 2/4 .* Ws_1200(i) + opt.TopExtraOffset;

        V1200extrap = [];
        xfit = [];
        coeffs = [NaN NaN];

        Vmod_1200 = Vmod_1200_measured;

        if opt.ApplyTopExtrapolation

            nTop = min(opt.NTopCellsExtrap, numel(z_1200));

            validTop = ~isnan(z_1200(1:nTop)) & ~isnan(Vmod_1200_measured(1:nTop));

            if sum(validTop) >= 2

                zTop = z_1200(1:nTop);
                vTop = Vmod_1200_measured(1:nTop);

                zTop = zTop(validTop);
                vTop = vTop(validTop);

                coeffs = polyfit(zTop(:)', vTop(:)', 1);

                xfit = linspace(z_1200(1), z_1200(1) + Z_upValue, 2);

                V1200extrap = polyval(coeffs, xfit);

                % Keep old logic: prepend extrapolated points.
                Vmod_1200 = [V1200extrap(:); Vmod_1200_measured(:)];

            else
                if opt.Verbose
                    warning('Row %d: not enough valid top cells for extrapolation. Keeping measured Vmod_1200 only.', i);
                end
            end
        end

        %% Variables for shear-stress step

        V1x_1 = V1x_(1:length(z_1200),:);
        V1y_1 = V1y_(1:length(z_1200),:);
        V1w_1 = V1w_(1:length(z_1200),:);
        V1e_1 = V1e_(1:length(z_1200),:);

        Vmod_1200std = nanstdLocal(sqrt(V1x_'.^2 + V1y_'.^2), 0, 1);
        Vmod_1200std = Vmod_1200std(:);

        Vmod_1200stdMOD(i) = safeMean(Vmod_1200std);
        VmodAll_1200(i) = safeMean(Vmod_1200);

        if isempty(Vmod_1200std)
            Vmod1low_1200stdMOD(i) = NaN;
        else
            idxLowerStdStart = max(1, round(lowerStartFraction .* length(Vmod_1200std)));
            Vmod1low_1200stdMOD(i) = safeMean(Vmod_1200std(idxLowerStdStart:end));
        end

        V1x_1std = nanstdLocal(V1x_1', 0, 1);
        V1y_1std = nanstdLocal(V1y_1', 0, 1);

        V1x_1std = V1x_1std(:);
        V1y_1std = V1y_1std(:);

        %% Last-cell / last-two-cell velocity series for BT-bias comparison

        VmagSeries = sqrt(V1x_1.^2 + V1y_1.^2);

        lastCellIdx = size(VmagSeries,1);

        if lastCellIdx >= 1
            VlastCell_1200 = VmagSeries(lastCellIdx, :)';
            VlastCellX_1200 = V1x_1(lastCellIdx, :)';
            VlastCellY_1200 = V1y_1(lastCellIdx, :)';
        else
            VlastCell_1200 = [];
            VlastCellX_1200 = [];
            VlastCellY_1200 = [];
        end

        if lastCellIdx >= 2
            Vlast2Cells_1200 = nanmeanLocal(VmagSeries(lastCellIdx-1:lastCellIdx, :), 1)';
            Vlast2CellsX_1200 = nanmeanLocal(V1x_1(lastCellIdx-1:lastCellIdx, :), 1)';
            Vlast2CellsY_1200 = nanmeanLocal(V1y_1(lastCellIdx-1:lastCellIdx, :), 1)';
        else
            Vlast2Cells_1200 = VlastCell_1200;
            Vlast2CellsX_1200 = VlastCellX_1200;
            Vlast2CellsY_1200 = VlastCellY_1200;
        end

        VlastCellMean_1200(i) = safeMean(VlastCell_1200);
        Vlast2CellsMean_1200(i) = safeMean(Vlast2Cells_1200);

        VlastCellMeanX_1200(i) = safeMean(VlastCellX_1200);
        VlastCellMeanY_1200(i) = safeMean(VlastCellY_1200);

        Vlast2CellsMeanX_1200(i) = safeMean(Vlast2CellsX_1200);
        Vlast2CellsMeanY_1200(i) = safeMean(Vlast2CellsY_1200);

        nCellsUsed_1200(i) = length(z_1200);
        nEnsUsed_1200(i) = size(V1x_1,2);
        Z_upValue_all(i) = Z_upValue;

        %% Store cells

        V1x__{i} = V1x_;
        V1y__{i} = V1y_;
        V1w__{i} = V1w_;
        V1e__{i} = V1e_;

        V1x_1_{i} = V1x_1;
        V1y_1_{i} = V1y_1;
        V1w_1_{i} = V1w_1;
        V1e_1_{i} = V1e_1;

        V1x_1200_{i} = V1x_1200;
        V1y_1200_{i} = V1y_1200;
        V1w_1200_{i} = V1w_1200;
        V1e_1200_{i} = V1e_1200;

        Vmod_1200_{i} = Vmod_1200;
        Vmod_1200_1_{i} = Vmod_1200_1;
        Vmod_1200_measured_{i} = Vmod_1200_measured;

        Vmod1_1200_{i} = Vmod1_1200;
        VmodLower_1200_{i} = VmodLower_1200;

        Vmod_1200std_{i} = Vmod_1200std;

        VlastCell_1200_{i} = VlastCell_1200;
        Vlast2Cells_1200_{i} = Vlast2Cells_1200;

        VlastCellX_1200_{i} = VlastCellX_1200;
        VlastCellY_1200_{i} = VlastCellY_1200;

        Vlast2CellsX_1200_{i} = Vlast2CellsX_1200;
        Vlast2CellsY_1200_{i} = Vlast2CellsY_1200;

        H_1200_{i} = H_1200;
        z_1200_{i} = z_1200;
        z1_1200_{i} = z1_1200;

        V1200extrap_{i} = V1200extrap;
        xfit_{i} = xfit;
        coeffs_{i} = coeffs;

        %% Optional plots

        if opt.PlotResults && shouldPlotThisRow(i, opt.PlotRows)

            plotWaterProfileDiagnostics( ...
                i, ...
                thisADCPfilename, ...
                z_1200, ...
                z1_1200, ...
                Vmod_1200_measured, ...
                Vmod_1200, ...
                Vmod1_1200, ...
                VlastCell_1200, ...
                Vlast2Cells_1200, ...
                opt.PlotMeasuredProfile, ...
                opt.PlotExtrapolatedProfile, ...
                opt.PlotLastCellSeries);
        end

        %% Store structure row

        WaterProfileAll(i).matchedFile = matchedFile;
        WaterProfileAll(i).fileFound = true;
        WaterProfileAll(i).start = start;
        WaterProfileAll(i).finish = finish;

        WaterProfileAll(i).WF_12 = WF_12(i);
        WaterProfileAll(i).WS_12 = WS_12(i);
        WaterProfileAll(i).WT_12 = WT_12(i);
        WaterProfileAll(i).Wf_1200 = Wf_1200(i);
        WaterProfileAll(i).Ws_1200 = Ws_1200(i);
        WaterProfileAll(i).Temp1 = Temp1;
        WaterProfileAll(i).Depth_1200 = Depth_1200(i);

        WaterProfileAll(i).LowerFraction = opt.LowerFraction;
        WaterProfileAll(i).lowerStartFraction = lowerStartFraction;

        WaterProfileAll(i).H_1200 = H_1200;
        WaterProfileAll(i).z_1200 = z_1200;
        WaterProfileAll(i).z1_1200 = z1_1200;

        WaterProfileAll(i).V1x_ = V1x_;
        WaterProfileAll(i).V1y_ = V1y_;
        WaterProfileAll(i).V1w_ = V1w_;
        WaterProfileAll(i).V1e_ = V1e_;

        WaterProfileAll(i).V1x_1 = V1x_1;
        WaterProfileAll(i).V1y_1 = V1y_1;
        WaterProfileAll(i).V1w_1 = V1w_1;
        WaterProfileAll(i).V1e_1 = V1e_1;

        WaterProfileAll(i).V1x_1200 = V1x_1200;
        WaterProfileAll(i).V1y_1200 = V1y_1200;
        WaterProfileAll(i).V1w_1200 = V1w_1200;
        WaterProfileAll(i).V1e_1200 = V1e_1200;

        WaterProfileAll(i).Vmod_1200_1 = Vmod_1200_1;
        WaterProfileAll(i).Vmod_1200_measured = Vmod_1200_measured;
        WaterProfileAll(i).Vmod_1200 = Vmod_1200;
        WaterProfileAll(i).Vmod1_1200 = Vmod1_1200;
        WaterProfileAll(i).VmodLower_1200 = VmodLower_1200;

        WaterProfileAll(i).Z_upValue = Z_upValue;
        WaterProfileAll(i).coeffs = coeffs;
        WaterProfileAll(i).xfit = xfit;
        WaterProfileAll(i).V1200extrap = V1200extrap;

        WaterProfileAll(i).Vmod_1200std = Vmod_1200std;
        WaterProfileAll(i).Vmod_1200stdMOD = Vmod_1200stdMOD(i);
        WaterProfileAll(i).VmodAll_1200 = VmodAll_1200(i);
        WaterProfileAll(i).Vmod1low_1200stdMOD = Vmod1low_1200stdMOD(i);

        WaterProfileAll(i).V1x_1std = V1x_1std;
        WaterProfileAll(i).V1y_1std = V1y_1std;

        WaterProfileAll(i).VlastCell_1200 = VlastCell_1200;
        WaterProfileAll(i).Vlast2Cells_1200 = Vlast2Cells_1200;
        WaterProfileAll(i).VlastCellMean_1200 = VlastCellMean_1200(i);
        WaterProfileAll(i).Vlast2CellsMean_1200 = Vlast2CellsMean_1200(i);

        WaterProfileAll(i).message = "OK";
        Message(i) = "OK";

    catch ME

        Message(i) = "Error: " + string(ME.message);

        WaterProfileAll(i).matchedFile = matchedFile;
        WaterProfileAll(i).fileFound = true;
        WaterProfileAll(i).message = Message(i);

        warning('Water-profile processing failed for row %d, file %s:\n%s', ...
            i, thisADCPfilename, ME.message);
    end
end

%% Add outputs to table

WaterProfileTable.FileFound = FileFound;
WaterProfileTable.MatchedFile = MatchedFile;
WaterProfileTable.start = startAll;
WaterProfileTable.finish = finishAll;

WaterProfileTable.WF_12 = WF_12;
WaterProfileTable.WS_12 = WS_12;
WaterProfileTable.WT_12 = WT_12;

WaterProfileTable.Wf_1200 = Wf_1200;
WaterProfileTable.Ws_1200 = Ws_1200;

WaterProfileTable.Depth_1200 = Depth_1200;

WaterProfileTable.LowerFraction = LowerFraction_used;
WaterProfileTable.lowerStartFraction = lowerStartFraction_used;

WaterProfileTable.nCellsUsed_1200 = nCellsUsed_1200;
WaterProfileTable.nEnsUsed_1200 = nEnsUsed_1200;

WaterProfileTable.Z_upValue = Z_upValue_all;

WaterProfileTable.VmodAll_1200 = VmodAll_1200;
WaterProfileTable.Vmod_1200stdMOD = Vmod_1200stdMOD;
WaterProfileTable.Vmod1low_1200stdMOD = Vmod1low_1200stdMOD;

WaterProfileTable.VlastCellMean_1200 = VlastCellMean_1200;
WaterProfileTable.Vlast2CellsMean_1200 = Vlast2CellsMean_1200;

WaterProfileTable.VlastCellMeanX_1200 = VlastCellMeanX_1200;
WaterProfileTable.VlastCellMeanY_1200 = VlastCellMeanY_1200;

WaterProfileTable.Vlast2CellsMeanX_1200 = Vlast2CellsMeanX_1200;
WaterProfileTable.Vlast2CellsMeanY_1200 = Vlast2CellsMeanY_1200;

WaterProfileTable.Message = Message;

%% Last successful row variables

lastOK = find(FileFound & Message == "OK", 1, 'last');

if isempty(lastOK)

    V1x_ = [];
    V1y_ = [];
    V1w_ = [];
    V1e_ = [];

    V1x_1200 = [];
    V1y_1200 = [];
    V1w_1200 = [];
    V1e_1200 = [];

    Vmod_1200_1 = [];
    Vmod_1200_measured = [];
    Vmod_1200 = [];
    Vmod1_1200 = [];
    VmodLower_1200 = [];

    Z_upValue = [];
    coeffs = [];
    xfit = [];
    V1200extrap = [];

    V1x_1 = [];
    V1y_1 = [];
    V1w_1 = [];
    V1e_1 = [];

    Vmod_1200std = [];
    V1x_1std = [];
    V1y_1std = [];

    H_1200 = [];
    z_1200 = [];
    z1_1200 = [];

    VlastCell_1200 = [];
    Vlast2Cells_1200 = [];

    VlastCellX_1200 = [];
    VlastCellY_1200 = [];
    Vlast2CellsX_1200 = [];
    Vlast2CellsY_1200 = [];

    start = [];
    finish = [];

else

    V1x_ = V1x__{lastOK};
    V1y_ = V1y__{lastOK};
    V1w_ = V1w__{lastOK};
    V1e_ = V1e__{lastOK};

    V1x_1200 = V1x_1200_{lastOK};
    V1y_1200 = V1y_1200_{lastOK};
    V1w_1200 = V1w_1200_{lastOK};
    V1e_1200 = V1e_1200_{lastOK};

    Vmod_1200_1 = Vmod_1200_1_{lastOK};
    Vmod_1200_measured = Vmod_1200_measured_{lastOK};
    Vmod_1200 = Vmod_1200_{lastOK};
    Vmod1_1200 = Vmod1_1200_{lastOK};
    VmodLower_1200 = VmodLower_1200_{lastOK};

    Z_upValue = Z_upValue_all(lastOK);
    coeffs = coeffs_{lastOK};
    xfit = xfit_{lastOK};
    V1200extrap = V1200extrap_{lastOK};

    V1x_1 = V1x_1_{lastOK};
    V1y_1 = V1y_1_{lastOK};
    V1w_1 = V1w_1_{lastOK};
    V1e_1 = V1e_1_{lastOK};

    Vmod_1200std = Vmod_1200std_{lastOK};
    V1x_1std = nanstdLocal(V1x_1', 0, 1);
    V1y_1std = nanstdLocal(V1y_1', 0, 1);

    V1x_1std = V1x_1std(:);
    V1y_1std = V1y_1std(:);

    H_1200 = H_1200_{lastOK};
    z_1200 = z_1200_{lastOK};
    z1_1200 = z1_1200_{lastOK};

    VlastCell_1200 = VlastCell_1200_{lastOK};
    Vlast2Cells_1200 = Vlast2Cells_1200_{lastOK};

    VlastCellX_1200 = VlastCellX_1200_{lastOK};
    VlastCellY_1200 = VlastCellY_1200_{lastOK};
    Vlast2CellsX_1200 = Vlast2CellsX_1200_{lastOK};
    Vlast2CellsY_1200 = Vlast2CellsY_1200_{lastOK};

    start = startAll(lastOK);
    finish = finishAll(lastOK);
end

%% Build output structure

WaterProfileOut = struct();

WaterProfileOut.InputRow = InputRow;
WaterProfileOut.ADCPfilename = ADCPfilename;
WaterProfileOut.StartADCP = StartADCP;
WaterProfileOut.EndADCP = EndADCP;

WaterProfileOut.DistanceFromLeftBank = DistanceFromLeftBank;
WaterProfileOut.TrateMean = TrateMean;
WaterProfileOut.D50 = D50;
WaterProfileOut.D90 = D90;
WaterProfileOut.DmeaNUMpsd = DmeaNUMpsd;
WaterProfileOut.SSC = SSC;

WaterProfileOut.Depth_1200 = Depth_1200;

WaterProfileOut.LowerFraction = opt.LowerFraction;
WaterProfileOut.lowerStartFraction = lowerStartFraction;

WaterProfileOut.WF_12 = WF_12;
WaterProfileOut.WS_12 = WS_12;
WaterProfileOut.WT_12 = WT_12;
WaterProfileOut.Wf_1200 = Wf_1200;
WaterProfileOut.Ws_1200 = Ws_1200;
WaterProfileOut.Temp1 = Temp1;

WaterProfileOut.V1x_ = V1x_;
WaterProfileOut.V1y_ = V1y_;
WaterProfileOut.V1w_ = V1w_;
WaterProfileOut.V1e_ = V1e_;

WaterProfileOut.V1x_1200 = V1x_1200;
WaterProfileOut.V1y_1200 = V1y_1200;
WaterProfileOut.V1w_1200 = V1w_1200;
WaterProfileOut.V1e_1200 = V1e_1200;

WaterProfileOut.Vmod_1200_1 = Vmod_1200_1;
WaterProfileOut.Vmod_1200_measured = Vmod_1200_measured;
WaterProfileOut.Vmod_1200 = Vmod_1200;
WaterProfileOut.Vmod1_1200 = Vmod1_1200;
WaterProfileOut.VmodLower_1200 = VmodLower_1200;

WaterProfileOut.Z_upValue = Z_upValue;
WaterProfileOut.coeffs = coeffs;
WaterProfileOut.xfit = xfit;
WaterProfileOut.V1200extrap = V1200extrap;

WaterProfileOut.V1x_1 = V1x_1;
WaterProfileOut.V1y_1 = V1y_1;
WaterProfileOut.V1w_1 = V1w_1;
WaterProfileOut.V1e_1 = V1e_1;

WaterProfileOut.Vmod_1200std = Vmod_1200std;
WaterProfileOut.Vmod_1200stdMOD = Vmod_1200stdMOD;
WaterProfileOut.VmodAll_1200 = VmodAll_1200;
WaterProfileOut.Vmod1low_1200stdMOD = Vmod1low_1200stdMOD;

WaterProfileOut.V1x_1std = V1x_1std;
WaterProfileOut.V1y_1std = V1y_1std;

WaterProfileOut.H_1200 = H_1200;
WaterProfileOut.z_1200 = z_1200;
WaterProfileOut.z1_1200 = z1_1200;

WaterProfileOut.VlastCell_1200 = VlastCell_1200;
WaterProfileOut.Vlast2Cells_1200 = Vlast2Cells_1200;

WaterProfileOut.VlastCellX_1200 = VlastCellX_1200;
WaterProfileOut.VlastCellY_1200 = VlastCellY_1200;
WaterProfileOut.Vlast2CellsX_1200 = Vlast2CellsX_1200;
WaterProfileOut.Vlast2CellsY_1200 = Vlast2CellsY_1200;

WaterProfileOut.VlastCellMean_1200 = VlastCellMean_1200;
WaterProfileOut.Vlast2CellsMean_1200 = Vlast2CellsMean_1200;

WaterProfileOut.VlastCellMeanX_1200 = VlastCellMeanX_1200;
WaterProfileOut.VlastCellMeanY_1200 = VlastCellMeanY_1200;

WaterProfileOut.Vlast2CellsMeanX_1200 = Vlast2CellsMeanX_1200;
WaterProfileOut.Vlast2CellsMeanY_1200 = Vlast2CellsMeanY_1200;

WaterProfileOut.V1x__ = V1x__;
WaterProfileOut.V1y__ = V1y__;
WaterProfileOut.V1w__ = V1w__;
WaterProfileOut.V1e__ = V1e__;

WaterProfileOut.V1x_1_ = V1x_1_;
WaterProfileOut.V1y_1_ = V1y_1_;
WaterProfileOut.V1w_1_ = V1w_1_;
WaterProfileOut.V1e_1_ = V1e_1_;

WaterProfileOut.V1x_1200_ = V1x_1200_;
WaterProfileOut.V1y_1200_ = V1y_1200_;
WaterProfileOut.V1w_1200_ = V1w_1200_;
WaterProfileOut.V1e_1200_ = V1e_1200_;

WaterProfileOut.Vmod_1200_ = Vmod_1200_;
WaterProfileOut.Vmod_1200_1_ = Vmod_1200_1_;
WaterProfileOut.Vmod_1200_measured_ = Vmod_1200_measured_;
WaterProfileOut.Vmod1_1200_ = Vmod1_1200_;
WaterProfileOut.VmodLower_1200_ = VmodLower_1200_;
WaterProfileOut.Vmod_1200std_ = Vmod_1200std_;

WaterProfileOut.VlastCell_1200_ = VlastCell_1200_;
WaterProfileOut.Vlast2Cells_1200_ = Vlast2Cells_1200_;

WaterProfileOut.VlastCellX_1200_ = VlastCellX_1200_;
WaterProfileOut.VlastCellY_1200_ = VlastCellY_1200_;
WaterProfileOut.Vlast2CellsX_1200_ = Vlast2CellsX_1200_;
WaterProfileOut.Vlast2CellsY_1200_ = Vlast2CellsY_1200_;

WaterProfileOut.H_1200_ = H_1200_;
WaterProfileOut.z_1200_ = z_1200_;
WaterProfileOut.z1_1200_ = z1_1200_;

WaterProfileOut.V1200extrap_ = V1200extrap_;
WaterProfileOut.xfit_ = xfit_;
WaterProfileOut.coeffs_ = coeffs_;

WaterProfileOut.WaterProfileAll = WaterProfileAll;
WaterProfileOut.WaterProfileTable = WaterProfileTable;

WaterProfileOut.FileFound = FileFound;
WaterProfileOut.MatchedFile = MatchedFile;
WaterProfileOut.Message = Message;

WaterProfileOut.start = start;
WaterProfileOut.finish = finish;
WaterProfileOut.startAll = startAll;
WaterProfileOut.finishAll = finishAll;

%% Assign to workspace

if opt.AssignToCaller
    assignStructToWorkspace(WaterProfileOut, 'caller');
end

if opt.AssignToBase
    assignStructToWorkspace(WaterProfileOut, 'base');
end

%% Display summary

if opt.Verbose
    fprintf('\n=== Water-profile processing summary ===\n');
    fprintf('Input rows: %d\n', nFiles);
    fprintf('Matched files: %d/%d\n', sum(FileFound), nFiles);
    fprintf('Successfully processed: %d/%d\n', sum(Message == "OK"), nFiles);
    fprintf('WF/WS were converted from cm to m using /100.\n');
    fprintf('z_1200 is physical distance above bed, calculated from Depth_1200 - H_1200.\n');
    fprintf('LowerFraction = %.2f, lowerStartFraction = %.2f.\n', opt.LowerFraction, lowerStartFraction);
    fprintf('StartADCP/EndADCP were used as trim counts.\n');
    fprintf('No files were saved.\n');
    fprintf('========================================\n');

    disp(WaterProfileTable);
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


function Depth_1200 = getDepth1200Input(opt, ADCPin, nFiles)

Depth_1200 = nan(nFiles,1);

if ~isempty(opt.Depth_1200)
    Depth_1200 = prepareDepthVector(opt.Depth_1200, nFiles);
    return
end

if ~isempty(opt.Depth1200)
    Depth_1200 = depthFromDepth1200Matrix(opt.Depth1200, nFiles);
    return
end

if ~isempty(opt.DepthRNGout)

    D = opt.DepthRNGout;

    if isfield(D, 'Depth_1200')
        Depth_1200 = prepareDepthVector(D.Depth_1200, nFiles);
        return
    end

    if isfield(D, 'Depth1200')
        Depth_1200 = depthFromDepth1200Matrix(D.Depth1200, nFiles);
        return
    end
end

if isfield(ADCPin, 'Depth_1200')
    Depth_1200 = prepareDepthVector(ADCPin.Depth_1200, nFiles);
    return
end

if isfield(ADCPin, 'Depth1200')
    Depth_1200 = depthFromDepth1200Matrix(ADCPin.Depth1200, nFiles);
    return
end

end


function d = prepareDepthVector(x, nFiles)

x = double(x);

if isvector(x)
    d = x(:);
else
    d = nanmeanLocal(x, 2);
    d = d(:);
end

if numel(d) == 1 && nFiles > 1
    d = repmat(d, nFiles, 1);
elseif numel(d) < nFiles
    d = [d; nan(nFiles - numel(d),1)];
elseif numel(d) > nFiles
    d = d(1:nFiles);
end

end


function d = depthFromDepth1200Matrix(Depth1200, nFiles)

Depth1200 = double(Depth1200);

if isvector(Depth1200)
    d = prepareDepthVector(Depth1200, nFiles);
else
    d = nanmeanLocal(Depth1200, 2);
    d = prepareDepthVector(d, nFiles);
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


function checkVelocityFields(v)

requiredFields = {'v1','v2','v3','v4'};

for k = 1:numel(requiredFields)
    if ~isfield(v, requiredFields{k})
        error('Velocity field v.%s is missing.', requiredFields{k});
    end
end

end


function n = getVelocityLength(v)

n1 = size(v.v1,2);
n2 = size(v.v2,2);
n3 = size(v.v3,2);
n4 = size(v.v4,2);

n = min([n1 n2 n3 n4]);

end


function v = getDataFieldOptional(Data, fieldName)

if isfield(Data, fieldName)
    v = double(Data.(fieldName));
else
    v = NaN;
end

end


function x = cropVectorSafe(x, start, finish)

x = double(x);
x = x(:);

if isempty(x) || all(isnan(x))
    return
end

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


function V = filterVelocityMatrix(V, maxAbsVelocity)

V = double(V);

V(abs(V) > maxAbsVelocity) = NaN;
V(V == -32768) = NaN;
V(V == 32767) = NaN;
V(V == 9999) = NaN;
V(V == -9999) = NaN;

end


function Vout = despikeVelocityMatrix(V, T)

Vout = V;

for r = 1:size(V,1)
    Vout(r,:) = despikeVector(V(r,:), T);
end

end


function y = despikeVector(x, T)

x = double(x);

valid = ~isnan(x);

if sum(valid) < 3
    y = x;
    return
end

T1 = mean(abs(x(valid))) - T .* std(abs(x(valid)));
To = mean(abs(x(valid))) + T .* std(abs(x(valid)));

peaks = (x > To | x < T1);
x(peaks) = NaN;

y = x;

end


function x = getFirstExistingScalar(S, possibleFields)

x = NaN;

for k = 1:numel(possibleFields)
    if isfield(S, possibleFields{k})
        val = double(S.(possibleFields{k}));
        if ~isempty(val)
            x = val(1);
            return
        end
    end
end

end


function Depth_1200_i = estimateDepth1200FromBT(S, start, finish)

Depth_1200_i = NaN;

if ~isfield(S, 'BT')
    return
end

BT = S.BT;

R = [];

if isfield(BT, 'RNG1')
    r = cropVectorSafe(BT.RNG1, start, finish);
    r(r == 0) = NaN;
    R = [R, r(:)]; %#ok<AGROW>
end

if isfield(BT, 'RNG2')
    r = cropVectorSafe(BT.RNG2, start, finish);
    r(r == 0) = NaN;
    R = [R, r(:)]; %#ok<AGROW>
end

if isfield(BT, 'RNG3')
    r = cropVectorSafe(BT.RNG3, start, finish);
    r(r == 0) = NaN;
    R = [R, r(:)]; %#ok<AGROW>
end

if isfield(BT, 'RNG4')
    r = cropVectorSafe(BT.RNG4, start, finish);
    r(r == 0) = NaN;
    R = [R, r(:)]; %#ok<AGROW>
end

if ~isempty(R)
    Depth_1200_i = safeMean(R(:));
end

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

if isempty(PlotRows)
    tf = true;
else
    tf = any(i == PlotRows);
end

end


function plotWaterProfileDiagnostics( ...
    i, ...
    thisADCPfilename, ...
    z_1200, ...
    z1_1200, ...
    Vmod_1200_measured, ...
    Vmod_1200, ...
    Vmod1_1200, ...
    VlastCell_1200, ...
    Vlast2Cells_1200, ...
    PlotMeasuredProfile, ...
    PlotExtrapolatedProfile, ...
    PlotLastCellSeries)

nPanels = double(logical(PlotMeasuredProfile)) + ...
          double(logical(PlotExtrapolatedProfile)) + ...
          double(logical(PlotLastCellSeries));

if nPanels == 0
    return
end

figure('Name', sprintf('Water profile diagnostics row %d', i), 'Color', 'w');

panel = 0;

if PlotMeasuredProfile

    panel = panel + 1;
    subplot(nPanels,1,panel);

    plot(Vmod_1200_measured, z_1200, '-o', 'DisplayName', 'Measured Vmod'); hold on
    plot(Vmod1_1200, z1_1200, 's', 'DisplayName', 'Lower part');

    grid on
    xlabel('Velocity magnitude');
    ylabel('z above bed [m]');
    title(sprintf('Measured water velocity profile — row %d: %s', ...
        i, cleanTitleString(thisADCPfilename)), ...
        'Interpreter', 'none');
    legend('Location','best');
end

if PlotExtrapolatedProfile

    panel = panel + 1;
    subplot(nPanels,1,panel);

    plot(Vmod_1200, 1:numel(Vmod_1200), '-o', 'DisplayName', 'Vmod\_1200 incl. extrap');

    grid on
    xlabel('Velocity magnitude');
    ylabel('Profile index');
    title('Velocity profile with optional upper extrapolation', 'Interpreter','none');
    legend('Location','best');
end

if PlotLastCellSeries

    panel = panel + 1;
    subplot(nPanels,1,panel);

    plot(VlastCell_1200, 'DisplayName', 'Last cell Vmod'); hold on
    plot(Vlast2Cells_1200, 'DisplayName', 'Mean last two cells Vmod');

    grid on
    xlabel('Trimmed ensemble index');
    ylabel('Velocity magnitude');
    title('Near-bed / last-cell velocity series', 'Interpreter','none');
    legend('Location','best');
end

end


function s = cleanTitleString(s)

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