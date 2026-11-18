function BTbedloadOut = processBTbedloadVelocity(dataADCP1200, ADCPin, varargin)
%PROCESSBTBEDLOADVELOCITY1200 Filter apparent bedload velocity from RDI BT data.
%
% This function:
%   - matches ADCPin.ADCPfilename with dataADCP1200.name
%   - uses StartADCP(i) and EndADCP(i) as trim counts
%   - extracts BT velocities from BT structure
%   - uses WaterProfileOut from processWaterProfile1200 for water-profile data
%   - keeps old variable names as much as possible
%   - stores every filtering step as cell arrays in the workspace
%   - exports one explicit final rejection mask per file:
%         BTfilterMask_1200{i}  true  = rejected/invalid ensemble
%         BTvalidMask_1200{i}   true  = retained final ensemble
%     The mask length is identical to the trimmed BT time series and is
%     intended for directly aligning backscatter and other processing.
%   - allows user to switch each filter on/off
%   - optionally plots apparent velocity va from all filtering steps
%     versus TrateMean and MaxTrate, with correlation coefficients.
%
% Important:
%   This function does not recalculate the water profile from v.v1/v.v2.
%   For water-bias and direction filters, provide:
%
%       'WaterProfileOut', WaterProfileOut
%
% New optional plot:
%
%   'PlotVaTrate', true
%
% This plots va from each filtering step against:
%   - TrateMean, marker 'o'
%   - MaxTrate, marker '^'
%
% and writes correlation coefficients inside each subplot.
%
% Recommended call:
%
%   BTbedloadOut = processBTbedloadVelocity1200(dataADCP1200, ADCPin, ...
%       'WaterProfileOut', WaterProfileOut, ...
%       'CoordinateSystem', 'instrument', ...
%       'WaterBiasReference', 'last2cells', ...
%       'BTerrorMultiplier', 4, ...
%       'PlotVaTrate', true, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
% No files are saved by this function.

%% Parse inputs

p = inputParser;

addRequired(p, 'dataADCP1200', @(x) isstruct(x));
addRequired(p, 'ADCPin', @(x) isstruct(x));

% Previous water-profile output
addParameter(p, 'WaterProfileOut', [], @(x) isempty(x) || isstruct(x));

% Coordinate choice: 'instrument', 'earth', or 'both'
addParameter(p, 'CoordinateSystem', 'instrument', @(x) ischar(x) || isstring(x));

% Old control parameters
addParameter(p, 'filtOUTwin', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'filtOUTstd', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'movVelTresh', 0.2, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'RDICutOFFtrim', 0.009, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Ws_min', 0.005, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ExtraFilter', 'no', @(x) ischar(x) || isstring(x));

% Scaling
addParameter(p, 'DummyValue', 1, @(x) isnumeric(x) && isscalar(x));

% Filter switches
addParameter(p, 'ApplySTDDespike', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyBoatMotionFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyMovingMedianFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyExtraPercentileFilter', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyErrorVelocityFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyGPSFilter', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyAttitudeFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyWaterBiasFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyDirectionFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ApplyInstrumentDirectionFilter', true, @(x) islogical(x) || isnumeric(x));

% Error velocity filter
addParameter(p, 'BTerrorMultiplier', 4, @(x) isnumeric(x) && isscalar(x));

% Direction / dot-product filter
addParameter(p, 'DirectionCosThreshold', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'DirectionStdT', 4, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'DirectionEndCutOffset', 2, @(x) isnumeric(x) && isscalar(x));

% Instrument-direction filter
addParameter(p, 'InstrumentEndCutOffset', 3, @(x) isnumeric(x) && isscalar(x));

% Water-bias filter
% WaterBiasReference:
%   'lastcell'   uses VlastCellX_1200_{i}, VlastCellY_1200_{i}
%   'last2cells' uses Vlast2CellsX_1200_{i}, Vlast2CellsY_1200_{i}
%
% WaterBiasComparison:
%   'component' compares x with x and y with y
%   'magnitude' compares sqrt(BTx^2+BTy^2) with sqrt(Vx^2+Vy^2)
addParameter(p, 'WaterBiasReference', 'last2cells', @(x) ischar(x) || isstring(x));
addParameter(p, 'WaterBiasComparison', 'component', @(x) ischar(x) || isstring(x));
addParameter(p, 'WaterBiasUseBothComponents', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'WaterBiasMultiplier', 1.0, @(x) isnumeric(x) && isscalar(x));

% GPS filter
addParameter(p, 'GPSdata', [], @(x) true);
addParameter(p, 'TresholdGPS', 0.2, @(x) isnumeric(x) && isscalar(x));

% Per-file diagnostic plots
addParameter(p, 'PlotResults', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotRows', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'PlotRawFiltered', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotFinalBT', true, @(x) islogical(x) || isnumeric(x));

% Optional va vs transport-rate plots
addParameter(p, 'PlotVaTrate', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotVaTrateSteps', [], @(x) isempty(x) || isnumeric(x) || iscellstr(x) || isstring(x));
addParameter(p, 'PlotVaTrateUseLogX', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotVaTrateUseLogY', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotVaTrateMinN', 3, @(x) isnumeric(x) && isscalar(x));

% Workspace/output
addParameter(p, 'AssignToCaller', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'AssignToBase', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, dataADCP1200, ADCPin, varargin{:});
opt = p.Results;

CoordinateSystem = lower(string(opt.CoordinateSystem));

if ~ismember(CoordinateSystem, ["instrument","earth","both"])
    error('CoordinateSystem must be ''instrument'', ''earth'', or ''both''.');
end

WaterBiasReference = lower(string(opt.WaterBiasReference));
WaterBiasComparison = lower(string(opt.WaterBiasComparison));

if ~ismember(WaterBiasReference, ["lastcell","last2cells"])
    error('WaterBiasReference must be ''lastcell'' or ''last2cells''.');
end

if ~ismember(WaterBiasComparison, ["component","magnitude"])
    error('WaterBiasComparison must be ''component'' or ''magnitude''.');
end

% Old alias for extra percentile filter
if contains(lower(string(opt.ExtraFilter)), "yes")
    opt.ApplyExtraPercentileFilter = true;
end

%% Checks

if isempty(dataADCP1200)
    error('dataADCP1200 is empty. Check filePattern1200 and dir(filePattern1200).');
end

if ~isfield(ADCPin, 'ADCPfilename')
    error('ADCPin.ADCPfilename is missing.');
end

if (opt.ApplyWaterBiasFilter || opt.ApplyDirectionFilter || opt.ApplyInstrumentDirectionFilter) && isempty(opt.WaterProfileOut)
    error(['WaterProfileOut is required when ApplyWaterBiasFilter, ', ...
           'ApplyDirectionFilter, or ApplyInstrumentDirectionFilter is true. ', ...
           'Run processWaterProfile1200 first and pass ''WaterProfileOut'', WaterProfileOut.']);
end

%% Variables from ADCPin

ADCPfilename = string(ADCPin.ADCPfilename(:));
nFiles = numel(ADCPfilename);

InputRow = (1:nFiles)';

StartADCP = getVectorField(ADCPin, 'StartADCP', zeros(nFiles,1), nFiles);
EndADCP   = getVectorField(ADCPin, 'EndADCP',   zeros(nFiles,1), nFiles);

DistanceFromLeftBank = getVectorField(ADCPin, 'DistanceFromLeftBank', nan(nFiles,1), nFiles);

TrateMean = getVectorFieldAny(ADCPin, ...
    {'TrateMean','MeanTrate','MTrateAll','MeanTrateAll'}, ...
    nan(nFiles,1), nFiles);

MaxTrate = getVectorFieldAny(ADCPin, ...
    {'MaxTrate','TrateMax','MTrateMax','MaxTrateAll'}, ...
    nan(nFiles,1), nFiles);

D50        = getVectorField(ADCPin, 'D50',        nan(nFiles,1), nFiles);
D90        = getVectorField(ADCPin, 'D90',        nan(nFiles,1), nFiles);
DmeaNUMpsd = getVectorField(ADCPin, 'DmeaNUMpsd', nan(nFiles,1), nFiles);
SSC        = getVectorField(ADCPin, 'SSC',        nan(nFiles,1), nFiles);

%% Main table

BTbedloadTable = table( ...
    InputRow, ...
    ADCPfilename, ...
    StartADCP, ...
    EndADCP, ...
    DistanceFromLeftBank, ...
    TrateMean, ...
    MaxTrate, ...
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
    'MaxTrate', ...
    'D50', ...
    'D90', ...
    'DmeaNUMpsd', ...
    'SSC'});

%% Status

FileFound = false(nFiles,1);
MatchedFile = strings(nFiles,1);
Message = strings(nFiles,1);

startAll = nan(nFiles,1);
finishAll = nan(nFiles,1);

%% Old scalar outputs

NANsbeforeX1200 = nan(nFiles,1);
NANsbeforeY1200 = nan(nFiles,1);
NANsbeforeX1200perc = nan(nFiles,1);
NANsbeforeY1200perc = nan(nFiles,1);

BTmod1200en = nan(nFiles,1);
BTmod1200org = nan(nFiles,1);
BTxSTD_1200 = nan(nFiles,1);
BTySTD_1200 = nan(nFiles,1);
BTstd_1200 = nan(nFiles,1);

BTmod1200a = nan(nFiles,1);
BTstd1200a = nan(nFiles,1);

BTerror1200m = nan(nFiles,1);

BTstd1200aW = nan(nFiles,1);
BTmod1200aW = nan(nFiles,1);

FilteredDataN_1200 = nan(nFiles,1);
FilteredDataE_1200 = nan(nFiles,1);

BTx_1200f = nan(nFiles,1);
BTy_1200f = nan(nFiles,1);
BTmodf_1200 = nan(nFiles,1);
BT_stdf_1200 = nan(nFiles,1);

FilteredData1_1200 = nan(nFiles,1);
FilteredData2_1200 = nan(nFiles,1);

BT_1_1200mean = nan(nFiles,1);
BT_2_1200mean = nan(nFiles,1);
BT_1_1200std = nan(nFiles,1);
BT_2_1200std = nan(nFiles,1);

%% Mean apparent velocity va for each filtering step

va_raw_1200 = nan(nFiles,1);
va_std_1200 = nan(nFiles,1);
va_boat_1200 = nan(nFiles,1);
va_median_1200 = nan(nFiles,1);
va_extra_1200 = nan(nFiles,1);
va_error_1200 = nan(nFiles,1);
va_gps_1200 = nan(nFiles,1);
va_attitude_1200 = nan(nFiles,1);
va_waterbias_1200 = nan(nFiles,1);
va_direction_1200 = nan(nFiles,1);
va_final_1200 = nan(nFiles,1);

%% Main step-by-step time-series cell outputs

BTx_raw_1200 = cell(nFiles,1);
BTy_raw_1200 = cell(nFiles,1);

BTx_std_1200 = cell(nFiles,1);
BTy_std_1200 = cell(nFiles,1);

BTx_boat_1200 = cell(nFiles,1);
BTy_boat_1200 = cell(nFiles,1);

BTx_median_1200 = cell(nFiles,1);
BTy_median_1200 = cell(nFiles,1);

BTx_extra_1200 = cell(nFiles,1);
BTy_extra_1200 = cell(nFiles,1);

BTx_error_1200 = cell(nFiles,1);
BTy_error_1200 = cell(nFiles,1);

BTx_gps_1200 = cell(nFiles,1);
BTy_gps_1200 = cell(nFiles,1);

BTx_attitude_1200 = cell(nFiles,1);
BTy_attitude_1200 = cell(nFiles,1);

BTx_waterbias_1200 = cell(nFiles,1);
BTy_waterbias_1200 = cell(nFiles,1);

BTx_direction_1200 = cell(nFiles,1);
BTy_direction_1200 = cell(nFiles,1);

BTx_final_1200 = cell(nFiles,1);
BTy_final_1200 = cell(nFiles,1);

%% Magnitude for every step

BTmod_raw_1200 = cell(nFiles,1);
BTmod_std_1200 = cell(nFiles,1);
BTmod_boat_1200 = cell(nFiles,1);
BTmod_median_1200 = cell(nFiles,1);
BTmod_extra_1200 = cell(nFiles,1);
BTmod_error_1200 = cell(nFiles,1);
BTmod_gps_1200 = cell(nFiles,1);
BTmod_attitude_1200 = cell(nFiles,1);
BTmod_waterbias_1200 = cell(nFiles,1);
BTmod_direction_1200 = cell(nFiles,1);
BTmod_final_1200 = cell(nFiles,1);

%% Old cell names preserved

BTx_1200 = cell(nFiles,1);
BTy_1200 = cell(nFiles,1);

BTx_1200a = cell(nFiles,1);
BTy_1200a = cell(nFiles,1);

BTmod1200_ = cell(nFiles,1);
BTmodea_1200_ = cell(nFiles,1);

BTmod1200aW_ = cell(nFiles,1);
BTmodf_1200ALL = cell(nFiles,1);

%% Earth E-N outputs, if available

BTE_raw_1200 = cell(nFiles,1);
BTN_raw_1200 = cell(nFiles,1);

BTE_std_1200 = cell(nFiles,1);
BTN_std_1200 = cell(nFiles,1);

BTE_boat_1200 = cell(nFiles,1);
BTN_boat_1200 = cell(nFiles,1);

BTE_median_1200 = cell(nFiles,1);
BTN_median_1200 = cell(nFiles,1);

BTE_extra_1200 = cell(nFiles,1);
BTN_extra_1200 = cell(nFiles,1);

BTE_gps_1200 = cell(nFiles,1);
BTN_gps_1200 = cell(nFiles,1);

% Old E-N names
BTE_1200 = cell(nFiles,1);
BTN_1200 = cell(nFiles,1);

BTE_1200a = cell(nFiles,1);
BTN_1200a = cell(nFiles,1);

%% Water-profile reused outputs

V1x_1_BT_1200 = cell(nFiles,1);
V1y_1_BT_1200 = cell(nFiles,1);

VwaterBiasX_1200 = cell(nFiles,1);
VwaterBiasY_1200 = cell(nFiles,1);
VwaterBiasMod_1200 = cell(nFiles,1);

%% Other cell outputs

% Explicit per-file masks for downstream processing. These are aligned to
% the trimmed interval startAll(i):finishAll(i), not to the raw full file.
% true in BTfilterMask_1200 means that the final BT velocity is invalid and
% the same ensemble should be set to NaN in comparable processing.
BTfilterMask_1200 = cell(nFiles,1);
BTvalidMask_1200 = cell(nFiles,1);
BTfilterMaskLength_1200 = nan(nFiles,1);
NfilteredFinal_1200 = nan(nFiles,1);
NvalidFinal_1200 = nan(nFiles,1);

peaksBTxy_1200 = cell(nFiles,1);
peaksBT_EN_S = cell(nFiles,1);

BTerror1200all = cell(nFiles,1);
idxBTer_1200 = cell(nFiles,1);

GPSidxOUT_1200 = cell(nFiles,1);
AttitudeidxOUT_1200 = cell(nFiles,1);

WaterBiasALL = cell(nFiles,1);

ErrorVel1200valid = cell(nFiles,1);

alphaV_1200 = cell(nFiles,1);
alphaBTf_1200 = cell(nFiles,1);
alphaBT_1200 = cell(nFiles,1);
VmodALL_1200 = cell(nFiles,1);

BT_1_1200_ALL = cell(nFiles,1);
BT_2_1200_ALL = cell(nFiles,1);

alphaV12_1200 = cell(nFiles,1);
alphaBT12 = cell(nFiles,1);
Vmod12ins_1200 = cell(nFiles,1);

%% Structure per row

BTbedloadAll = repmat(struct( ...
    'InputRow', [], ...
    'ADCPfilename', "", ...
    'matchedFile', "", ...
    'fileFound', false, ...
    'start', [], ...
    'finish', [], ...
    'CoordinateSystem', "", ...
    'BTx_raw', [], ...
    'BTy_raw', [], ...
    'BTx_std', [], ...
    'BTy_std', [], ...
    'BTx_boat', [], ...
    'BTy_boat', [], ...
    'BTx_median', [], ...
    'BTy_median', [], ...
    'BTx_extra', [], ...
    'BTy_extra', [], ...
    'BTx_error', [], ...
    'BTy_error', [], ...
    'BTx_gps', [], ...
    'BTy_gps', [], ...
    'BTx_attitude', [], ...
    'BTy_attitude', [], ...
    'BTx_waterbias', [], ...
    'BTy_waterbias', [], ...
    'BTx_direction', [], ...
    'BTy_direction', [], ...
    'BTx_final', [], ...
    'BTy_final', [], ...
    'BTmod_final', [], ...
    'BTE_raw', [], ...
    'BTN_raw', [], ...
    'BTE_std', [], ...
    'BTN_std', [], ...
    'BTE_boat', [], ...
    'BTN_boat', [], ...
    'BTE_median', [], ...
    'BTN_median', [], ...
    'BT_x_1200', [], ...
    'BT_y_1200', [], ...
    'BT_1_1200', [], ...
    'BT_2_1200', [], ...
    'V1x_1', [], ...
    'V1y_1', [], ...
    'VwaterBiasX', [], ...
    'VwaterBiasY', [], ...
    'WaterBias', [], ...
    'message', ""), nFiles, 1);

%% Process files

for i = 1:nFiles

    thisADCPfilename = ADCPfilename(i);
    matchedFile = matchADCPfile(dataADCP1200, thisADCPfilename);

    BTbedloadAll(i).InputRow = InputRow(i);
    BTbedloadAll(i).ADCPfilename = thisADCPfilename;
    BTbedloadAll(i).CoordinateSystem = CoordinateSystem;

    if matchedFile == ""

        FileFound(i) = false;
        Message(i) = "File not found";

        BTbedloadAll(i).fileFound = false;
        BTbedloadAll(i).message = Message(i);

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

        BT = S.BT;

        N = getBTvelocityLength(BT);

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
            fprintf('\nProcessing BT apparent bedload velocity row %d/%d\n', i, nFiles);
            fprintf('ADCPfilename: %s\n', thisADCPfilename);
            fprintf('Matched file: %s\n', matchedFile);
            fprintf('N ensembles: %d | using ensembles %d:%d\n', N, start, finish);
            fprintf('CoordinateSystem: %s\n', CoordinateSystem);
        end

        %% Extract raw BT components

        [BTx_raw, BTy_raw] = getBTmainComponents(BT, start, finish, CoordinateSystem);
        [BTE_raw, BTN_raw] = getBTearthComponents(BT, start, finish);

        BTx_raw_1200{i} = BTx_raw;
        BTy_raw_1200{i} = BTy_raw;
        BTmod_raw_1200{i} = calcMod(BTx_raw, BTy_raw);
        va_raw_1200(i) = safeMean(BTmod_raw_1200{i}) .* opt.DummyValue;

        BTE_raw_1200{i} = BTE_raw;
        BTN_raw_1200{i} = BTN_raw;

        %% Get water profile from WaterProfileOut only

        [V1x_1, V1y_1, V1w_1, V1e_1, z_1200, ...
            VbiasX, VbiasY, VbiasMod] = getWaterProfileInputsForBT( ...
            i, opt.WaterProfileOut, WaterBiasReference);

        V1x_1_BT_1200{i} = V1x_1;
        V1y_1_BT_1200{i} = V1y_1;

        VwaterBiasX_1200{i} = VbiasX;
        VwaterBiasY_1200{i} = VbiasY;
        VwaterBiasMod_1200{i} = VbiasMod;

        ErrorVel1200valid{i} = V1e_1;

        %% Step 0: initial std despike

        BTx_std = BTx_raw;
        BTy_std = BTy_raw;

        BTE_std = BTE_raw;
        BTN_std = BTN_raw;

        if opt.ApplySTDDespike
            BTx_std = despikeBT(BTx_std, opt.filtOUTstd);
            BTy_std = despikeBT(BTy_std, opt.filtOUTstd);

            if ~isempty(BTE_std)
                BTE_std = despikeBT(BTE_std, opt.filtOUTstd);
            end

            if ~isempty(BTN_std)
                BTN_std = despikeBT(BTN_std, opt.filtOUTstd);
            end
        end

        BTx_std_1200{i} = BTx_std;
        BTy_std_1200{i} = BTy_std;
        BTmod_std_1200{i} = calcMod(BTx_std, BTy_std);
        va_std_1200(i) = safeMean(BTmod_std_1200{i}) .* opt.DummyValue;

        BTE_std_1200{i} = BTE_std;
        BTN_std_1200{i} = BTN_std;

        % Old names after std despike
        BTx_1200{i} = BTx_std;
        BTy_1200{i} = BTy_std;

        BTE_1200{i} = BTE_std;
        BTN_1200{i} = BTN_std;

        %% Count NaNs before main filtering

        NANsbeforeX1200(i) = sum(isnan(BTx_std));
        NANsbeforeY1200(i) = sum(isnan(BTy_std));

        NANsbeforeX1200perc(i) = NANsbeforeX1200(i) ./ numel(BTx_std);
        NANsbeforeY1200perc(i) = NANsbeforeY1200(i) ./ numel(BTy_std);

        %% Step 1: boat-motion threshold

        BTx_boat = BTx_std;
        BTy_boat = BTy_std;

        BTE_boat = BTE_std;
        BTN_boat = BTN_std;

        if opt.ApplyBoatMotionFilter

            [BTx_boat, BTy_boat, peaksBTxy_1200{i}] = ...
                despikeBTmoveBoattresh(BTx_boat, BTy_boat, opt.movVelTresh);

            if ~isempty(BTE_boat) && ~isempty(BTN_boat)
                [BTE_boat, BTN_boat, peaksBT_EN_S{i}] = ...
                    despikeBTmoveBoattresh(BTE_boat, BTN_boat, opt.movVelTresh);
            else
                peaksBT_EN_S{i} = [];
            end

        else
            peaksBTxy_1200{i} = false(size(BTx_boat));
            peaksBT_EN_S{i} = false(size(BTE_boat));
        end

        BTx_boat_1200{i} = BTx_boat;
        BTy_boat_1200{i} = BTy_boat;
        BTmod_boat_1200{i} = calcMod(BTx_boat, BTy_boat);
        va_boat_1200(i) = safeMean(BTmod_boat_1200{i}) .* opt.DummyValue;

        BTE_boat_1200{i} = BTE_boat;
        BTN_boat_1200{i} = BTN_boat;

        BTE_1200{i} = BTE_boat;
        BTN_1200{i} = BTN_boat;

        %% Values before main filtering

        if ~isempty(BTE_boat) && ~isempty(BTN_boat)
            BTmod1200en(i) = safeMean(calcMod(BTE_boat, BTN_boat)) .* opt.DummyValue;
        end

        BTmod1200org(i) = safeMean(calcMod(BTx_boat, BTy_boat)) .* opt.DummyValue;
        BTmod1200_{i} = calcMod(BTx_boat, BTy_boat);

        BTxSTD_1200(i) = nanstdLocal(BTx_boat, 0, 1);
        BTySTD_1200(i) = nanstdLocal(BTy_boat, 0, 1);
        BTstd_1200(i) = nanstdLocal(calcMod(BTx_boat, BTy_boat), 0, 1) .* opt.DummyValue;

        %% Step 2: moving median filtering

        BTx_median = BTx_boat;
        BTy_median = BTy_boat;

        BTE_median = BTE_boat;
        BTN_median = BTN_boat;

        if opt.ApplyMovingMedianFilter
            BTx_median = filloutliersSafe(BTx_median, NaN, 'movmedian', opt.filtOUTwin);
            BTy_median = filloutliersSafe(BTy_median, NaN, 'movmedian', opt.filtOUTwin);

            if ~isempty(BTE_median)
                BTE_median = filloutliersSafe(BTE_median, NaN, 'movmedian', opt.filtOUTwin);
            end

            if ~isempty(BTN_median)
                BTN_median = filloutliersSafe(BTN_median, NaN, 'movmedian', opt.filtOUTwin);
            end
        end

        BTx_median_1200{i} = BTx_median;
        BTy_median_1200{i} = BTy_median;
        BTmod_median_1200{i} = calcMod(BTx_median, BTy_median);
        va_median_1200(i) = safeMean(BTmod_median_1200{i}) .* opt.DummyValue;

        BTE_median_1200{i} = BTE_median;
        BTN_median_1200{i} = BTN_median;

        BTE_1200a{i} = BTE_median;
        BTN_1200a{i} = BTN_median;

        %% Step 2b: optional percentile filter

        BTx_extra = BTx_median;
        BTy_extra = BTy_median;

        BTE_extra = BTE_median;
        BTN_extra = BTN_median;

        if opt.ApplyExtraPercentileFilter
            BTx_extra = filloutliersSafe(BTx_extra, NaN, 'percentiles', [5 95]);
            BTy_extra = filloutliersSafe(BTy_extra, NaN, 'percentiles', [5 95]);

            if ~isempty(BTE_extra)
                BTE_extra = filloutliersSafe(BTE_extra, NaN, 'percentiles', [5 95]);
            end

            if ~isempty(BTN_extra)
                BTN_extra = filloutliersSafe(BTN_extra, NaN, 'percentiles', [5 95]);
            end
        end

        BTx_extra_1200{i} = BTx_extra;
        BTy_extra_1200{i} = BTy_extra;
        BTmod_extra_1200{i} = calcMod(BTx_extra, BTy_extra);
        va_extra_1200(i) = safeMean(BTmod_extra_1200{i}) .* opt.DummyValue;

        BTE_extra_1200{i} = BTE_extra;
        BTN_extra_1200{i} = BTN_extra;

        BTx_1200a{i} = BTx_extra;
        BTy_1200a{i} = BTy_extra;

        BTE_1200a{i} = BTE_extra;
        BTN_1200a{i} = BTN_extra;

        BTmodea_1200 = calcMod(BTx_extra, BTy_extra);
        BTmodea_1200_{i} = BTmodea_1200;

        BTmod1200a(i) = safeMean(BTmodea_1200) .* opt.DummyValue;
        BTstd1200a(i) = nanstdLocal(BTmodea_1200, 0, 1) .* opt.DummyValue;

        %% Step 3: error-velocity filter

        BTx_error = BTx_extra;
        BTy_error = BTy_extra;

        BTerror = getBTerrorVelocity(BT, start, finish);
        BTerror1200m(i) = safeMean(BTerror);
        BTerror1200all{i} = BTerror;

        idxBTer = false(size(BTx_error));

        if opt.ApplyErrorVelocityFilter && ~isempty(BTerror)

            idxBTer = abs(BTerror) > opt.BTerrorMultiplier .* abs(safeMean(BTmodea_1200));
            idxBTer = idxBTer(:);

            nCommon = min([numel(idxBTer), numel(BTx_error), numel(BTy_error)]);
            idxBTer = idxBTer(1:nCommon);

            BTx_error(logical(idxBTer)) = NaN;
            BTy_error(logical(idxBTer)) = NaN;
        end

        idxBTer_1200{i} = idxBTer;

        BTx_error_1200{i} = BTx_error;
        BTy_error_1200{i} = BTy_error;
        BTmod_error_1200{i} = calcMod(BTx_error, BTy_error);
        va_error_1200(i) = safeMean(BTmod_error_1200{i}) .* opt.DummyValue;

        %% Step 4: GPS moving elimination

        BTx_gps = BTx_error;
        BTy_gps = BTy_error;

        idxOUTgps = false(size(BTx_gps));

        if opt.ApplyGPSFilter && ~isempty(opt.GPSdata)

            idxOUTgps = getGPSoutlierIndex(opt.GPSdata, start, finish, opt.TresholdGPS);

            nCommon = min([numel(idxOUTgps), numel(BTx_gps), numel(BTy_gps)]);
            idxOUTgps = idxOUTgps(1:nCommon);

            BTx_gps(logical(idxOUTgps)) = NaN;
            BTy_gps(logical(idxOUTgps)) = NaN;
        end

        GPSidxOUT_1200{i} = idxOUTgps;

        BTx_gps_1200{i} = BTx_gps;
        BTy_gps_1200{i} = BTy_gps;
        BTmod_gps_1200{i} = calcMod(BTx_gps, BTy_gps);
        va_gps_1200(i) = safeMean(BTmod_gps_1200{i}) .* opt.DummyValue;

        %% Step 4b: heading/pitch/roll filter

        BTx_attitude = BTx_gps;
        BTy_attitude = BTy_gps;

        idxOUTattitude = false(size(BTx_attitude));

        if opt.ApplyAttitudeFilter

            [EH, Pitch, Roll] = getAttitudeSeries(S, start, finish);

            idxOUTattitude = getAttitudeOutlierIndex(EH, Pitch, Roll);

            nCommon = min([numel(idxOUTattitude), numel(BTx_attitude), numel(BTy_attitude)]);
            idxOUTattitude = idxOUTattitude(1:nCommon);

            BTx_attitude(logical(idxOUTattitude)) = NaN;
            BTy_attitude(logical(idxOUTattitude)) = NaN;
        end

        AttitudeidxOUT_1200{i} = idxOUTattitude;

        BTx_attitude_1200{i} = BTx_attitude;
        BTy_attitude_1200{i} = BTy_attitude;
        BTmod_attitude_1200{i} = calcMod(BTx_attitude, BTy_attitude);
        va_attitude_1200(i) = safeMean(BTmod_attitude_1200{i}) .* opt.DummyValue;

        %% Step 5: water-bias filter using processWaterProfile1200 outputs

        BTx_waterbias = BTx_attitude;
        BTy_waterbias = BTy_attitude;

        WaterBias = repmat({'NO bias'}, numel(BTx_waterbias), 1);

        if opt.ApplyWaterBiasFilter

            [BTx_waterbias, BTy_waterbias, WaterBias] = applyWaterBiasFilterFromReference( ...
                BTx_waterbias, ...
                BTy_waterbias, ...
                VbiasX, ...
                VbiasY, ...
                WaterBiasComparison, ...
                opt.WaterBiasUseBothComponents, ...
                opt.WaterBiasMultiplier);
        else
            WaterBias(:) = {'not applied'};
        end

        WaterBiasALL{i} = WaterBias;

        BTx_waterbias_1200{i} = BTx_waterbias;
        BTy_waterbias_1200{i} = BTy_waterbias;
        BTmod_waterbias_1200{i} = calcMod(BTx_waterbias, BTy_waterbias);
        va_waterbias_1200(i) = safeMean(BTmod_waterbias_1200{i}) .* opt.DummyValue;

        BTmod1200aW_{i} = BTmod_waterbias_1200{i};
        BTstd1200aW(i) = nanstdLocal(BTmod1200aW_{i}, 0, 1) .* opt.DummyValue;
        BTmod1200aW(i) = safeMean(BTmod1200aW_{i}) .* opt.DummyValue;

        %% Step 6: main direction filter using WaterProfileOut V1x_1/V1y_1

        BT_x_1200 = BTx_waterbias;
        BT_y_1200 = BTy_waterbias;

        if opt.ApplyDirectionFilter

            directionEndCut = max(2, size(V1x_1,1) - opt.DirectionEndCutOffset);

            [BT_x_1200, BT_y_1200, ...
                FilteredDataN_1200(i), FilteredDataE_1200(i), ...
                alphaV_1200{i}, alphaBTf_1200{i}, alphaBT_1200{i}, VmodALL_1200{i}] = ...
                filterBT2org( ...
                V1x_1, ...
                V1y_1, ...
                BTx_waterbias, ...
                BTy_waterbias, ...
                opt.DirectionCosThreshold, ...
                opt.DirectionStdT, ...
                directionEndCut, ...
                NANsbeforeX1200(i), ...
                NANsbeforeY1200(i), ...
                opt.RDICutOFFtrim);
        end

        % filterBT2org can shorten the output when the supplied water-profile
        % series is shorter than the trimmed BT series. Restore the original
        % trimmed length so every downstream mask retains exact ensemble
        % alignment with start:finish.
        nTrimmedBT = numel(BTx_waterbias);
        BT_x_1200 = padVectorToLength(BT_x_1200, nTrimmedBT);
        BT_y_1200 = padVectorToLength(BT_y_1200, nTrimmedBT);

        BTx_direction_1200{i} = BT_x_1200;
        BTy_direction_1200{i} = BT_y_1200;
        BTmod_direction_1200{i} = calcMod(BT_x_1200, BT_y_1200);
        va_direction_1200(i) = safeMean(BTmod_direction_1200{i}) .* opt.DummyValue;

        BTx_final_1200{i} = BT_x_1200;
        BTy_final_1200{i} = BT_y_1200;
        BTmod_final_1200{i} = BTmod_direction_1200{i};
        va_final_1200(i) = safeMean(BTmod_final_1200{i}) .* opt.DummyValue;

        BTx_1200f(i) = safeMean(BT_x_1200);
        BTy_1200f(i) = safeMean(BT_y_1200);

        BTmodff_1200 = BTmod_final_1200{i} .* opt.DummyValue;
        BTmodf_1200(i) = safeMean(BTmodff_1200);
        BT_stdf_1200(i) = nanstdLocal(BTmodff_1200, 0, 1);
        BTmodf_1200ALL{i} = BTmodff_1200;

        % Exact final per-file mask. It includes raw missing values and every
        % selected BT filter because all rejected ensembles are represented
        % by NaN in the final x/y components.
        BTfilterMask_1200{i} = ...
            ~isfinite(BTx_final_1200{i}) | ~isfinite(BTy_final_1200{i});
        BTfilterMask_1200{i} = logical(BTfilterMask_1200{i}(:));
        BTvalidMask_1200{i} = ~BTfilterMask_1200{i};
        BTfilterMaskLength_1200(i) = numel(BTfilterMask_1200{i});
        NfilteredFinal_1200(i) = sum(BTfilterMask_1200{i});
        NvalidFinal_1200(i) = sum(BTvalidMask_1200{i});

        %% Step 6b: optional instrument-direction/homogeneity filter

        BT_1_1200 = [];
        BT_2_1200 = [];

        if opt.ApplyInstrumentDirectionFilter && isfield(BT, 'BTV1') && isfield(BT, 'BTV2')

            BTraw1 = cropVectorSafe(BT.BTV1, start, finish);
            BTraw2 = cropVectorSafe(BT.BTV2, start, finish);

            instrumentEndCut = max(2, size(V1x_1,1) - opt.InstrumentEndCutOffset);

            [BT_1_1200, BT_2_1200, ...
                FilteredData1_1200(i), FilteredData2_1200(i), ...
                alphaV12_1200{i}, alphaBTf_1200{i}, alphaBT12{i}, Vmod12ins_1200{i}] = ...
                filterBT2org( ...
                V1x_1, ...
                V1y_1, ...
                BTraw1, ...
                BTraw2, ...
                opt.DirectionCosThreshold, ...
                opt.DirectionStdT, ...
                instrumentEndCut, ...
                NANsbeforeX1200(i), ...
                NANsbeforeY1200(i), ...
                opt.RDICutOFFtrim);

            BT_1_1200mean(i) = safeMean(BT_1_1200);
            BT_2_1200mean(i) = safeMean(BT_2_1200) .* opt.DummyValue;

            BT_1_1200std(i) = nanstdLocal(BT_1_1200, 0, 1);
            BT_2_1200std(i) = nanstdLocal(BT_2_1200, 0, 1) .* opt.DummyValue;

            BT_1_1200_ALL{i} = BT_1_1200;
            BT_2_1200_ALL{i} = BT_2_1200;
        end

        %% Optional per-file plots

        if opt.PlotResults && shouldPlotThisRow(i, opt.PlotRows)
            plotBTfilterDiagnostics( ...
                i, ...
                thisADCPfilename, ...
                BTx_raw, ...
                BTy_raw, ...
                BTx_median, ...
                BTy_median, ...
                BT_x_1200, ...
                BT_y_1200, ...
                opt.PlotRawFiltered, ...
                opt.PlotFinalBT);
        end

        %% Store row structure

        BTbedloadAll(i).matchedFile = matchedFile;
        BTbedloadAll(i).fileFound = true;
        BTbedloadAll(i).start = start;
        BTbedloadAll(i).finish = finish;

        BTbedloadAll(i).BTx_raw = BTx_raw;
        BTbedloadAll(i).BTy_raw = BTy_raw;

        BTbedloadAll(i).BTx_std = BTx_std;
        BTbedloadAll(i).BTy_std = BTy_std;

        BTbedloadAll(i).BTx_boat = BTx_boat;
        BTbedloadAll(i).BTy_boat = BTy_boat;

        BTbedloadAll(i).BTx_median = BTx_median;
        BTbedloadAll(i).BTy_median = BTy_median;

        BTbedloadAll(i).BTx_extra = BTx_extra;
        BTbedloadAll(i).BTy_extra = BTy_extra;

        BTbedloadAll(i).BTx_error = BTx_error;
        BTbedloadAll(i).BTy_error = BTy_error;

        BTbedloadAll(i).BTx_gps = BTx_gps;
        BTbedloadAll(i).BTy_gps = BTy_gps;

        BTbedloadAll(i).BTx_attitude = BTx_attitude;
        BTbedloadAll(i).BTy_attitude = BTy_attitude;

        BTbedloadAll(i).BTx_waterbias = BTx_waterbias;
        BTbedloadAll(i).BTy_waterbias = BTy_waterbias;

        BTbedloadAll(i).BTx_direction = BT_x_1200;
        BTbedloadAll(i).BTy_direction = BT_y_1200;

        BTbedloadAll(i).BTx_final = BT_x_1200;
        BTbedloadAll(i).BTy_final = BT_y_1200;
        BTbedloadAll(i).BTmod_final = BTmod_final_1200{i};

        BTbedloadAll(i).BTE_raw = BTE_raw;
        BTbedloadAll(i).BTN_raw = BTN_raw;
        BTbedloadAll(i).BTE_std = BTE_std;
        BTbedloadAll(i).BTN_std = BTN_std;
        BTbedloadAll(i).BTE_boat = BTE_boat;
        BTbedloadAll(i).BTN_boat = BTN_boat;
        BTbedloadAll(i).BTE_median = BTE_median;
        BTbedloadAll(i).BTN_median = BTN_median;

        BTbedloadAll(i).BT_x_1200 = BT_x_1200;
        BTbedloadAll(i).BT_y_1200 = BT_y_1200;

        BTbedloadAll(i).BT_1_1200 = BT_1_1200;
        BTbedloadAll(i).BT_2_1200 = BT_2_1200;

        BTbedloadAll(i).V1x_1 = V1x_1;
        BTbedloadAll(i).V1y_1 = V1y_1;
        BTbedloadAll(i).VwaterBiasX = VbiasX;
        BTbedloadAll(i).VwaterBiasY = VbiasY;

        BTbedloadAll(i).WaterBias = WaterBias;

        Message(i) = "OK";
        BTbedloadAll(i).message = Message(i);

    catch ME

        Message(i) = "Error: " + string(ME.message);

        BTbedloadAll(i).matchedFile = matchedFile;
        BTbedloadAll(i).fileFound = true;
        BTbedloadAll(i).message = Message(i);

        warning('BT bedload-velocity processing failed for row %d, file %s:\n%s', ...
            i, thisADCPfilename, ME.message);
    end
end

%% Add table outputs

BTbedloadTable.FileFound = FileFound;
BTbedloadTable.MatchedFile = MatchedFile;
BTbedloadTable.start = startAll;
BTbedloadTable.finish = finishAll;

BTbedloadTable.NANsbeforeX1200 = NANsbeforeX1200;
BTbedloadTable.NANsbeforeY1200 = NANsbeforeY1200;
BTbedloadTable.NANsbeforeX1200perc = NANsbeforeX1200perc;
BTbedloadTable.NANsbeforeY1200perc = NANsbeforeY1200perc;

BTbedloadTable.BTmod1200en = BTmod1200en;
BTbedloadTable.BTmod1200org = BTmod1200org;

BTbedloadTable.BTxSTD_1200 = BTxSTD_1200;
BTbedloadTable.BTySTD_1200 = BTySTD_1200;
BTbedloadTable.BTstd_1200 = BTstd_1200;

BTbedloadTable.BTmod1200a = BTmod1200a;
BTbedloadTable.BTstd1200a = BTstd1200a;

BTbedloadTable.BTerror1200m = BTerror1200m;

BTbedloadTable.BTmod1200aW = BTmod1200aW;
BTbedloadTable.BTstd1200aW = BTstd1200aW;

BTbedloadTable.FilteredDataN_1200 = FilteredDataN_1200;
BTbedloadTable.FilteredDataE_1200 = FilteredDataE_1200;

BTbedloadTable.BTx_1200f = BTx_1200f;
BTbedloadTable.BTy_1200f = BTy_1200f;
BTbedloadTable.BTmodf_1200 = BTmodf_1200;
BTbedloadTable.BT_stdf_1200 = BT_stdf_1200;

BTbedloadTable.FilteredData1_1200 = FilteredData1_1200;
BTbedloadTable.FilteredData2_1200 = FilteredData2_1200;

BTbedloadTable.BT_1_1200mean = BT_1_1200mean;
BTbedloadTable.BT_2_1200mean = BT_2_1200mean;
BTbedloadTable.BT_1_1200std = BT_1_1200std;
BTbedloadTable.BT_2_1200std = BT_2_1200std;

BTbedloadTable.va_raw_1200 = va_raw_1200;
BTbedloadTable.va_std_1200 = va_std_1200;
BTbedloadTable.va_boat_1200 = va_boat_1200;
BTbedloadTable.va_median_1200 = va_median_1200;
BTbedloadTable.va_extra_1200 = va_extra_1200;
BTbedloadTable.va_error_1200 = va_error_1200;
BTbedloadTable.va_gps_1200 = va_gps_1200;
BTbedloadTable.va_attitude_1200 = va_attitude_1200;
BTbedloadTable.va_waterbias_1200 = va_waterbias_1200;
BTbedloadTable.va_direction_1200 = va_direction_1200;
BTbedloadTable.va_final_1200 = va_final_1200;

BTbedloadTable.BTfilterMaskLength_1200 = BTfilterMaskLength_1200;
BTbedloadTable.NfilteredFinal_1200 = NfilteredFinal_1200;
BTbedloadTable.NvalidFinal_1200 = NvalidFinal_1200;

BTbedloadTable.Message = Message;

%% Optional va vs transport-rate plots and correlation table

BTvaStepNames = [
    "raw"
    "std"
    "boat"
    "median"
    "extra"
    "error"
    "gps"
    "attitude"
    "waterbias"
    "direction"
    "final"];

BTvaAllSteps_1200 = [ ...
    va_raw_1200, ...
    va_std_1200, ...
    va_boat_1200, ...
    va_median_1200, ...
    va_extra_1200, ...
    va_error_1200, ...
    va_gps_1200, ...
    va_attitude_1200, ...
    va_waterbias_1200, ...
    va_direction_1200, ...
    va_final_1200];

BTvaTrateCorrTable = makeVaTrateCorrTable( ...
    BTvaStepNames, ...
    BTvaAllSteps_1200, ...
    TrateMean, ...
    MaxTrate, ...
    opt.PlotVaTrateMinN);

if opt.PlotVaTrate
    plotVaVsTrateAllSteps( ...
        BTvaStepNames, ...
        BTvaAllSteps_1200, ...
        TrateMean, ...
        MaxTrate, ...
        opt.PlotVaTrateSteps, ...
        opt.PlotVaTrateUseLogX, ...
        opt.PlotVaTrateUseLogY, ...
        opt.PlotVaTrateMinN);
end

%% Last successful row variables

lastOK = find(FileFound & Message == "OK", 1, 'last');

if isempty(lastOK)

    BTx_ = [];
    BTy_ = [];

    BTE_1200_ = [];
    BTN_1200_ = [];

    BT_x_1200 = [];
    BT_y_1200 = [];

    BT_1_1200 = [];
    BT_2_1200 = [];

    BTmodea_1200 = [];
    BTmodff_1200 = [];
    BTerror = [];
    WaterBias = {};

    start = [];
    finish = [];

else

    BTx_ = BTx_final_1200{lastOK};
    BTy_ = BTy_final_1200{lastOK};

    BTE_1200_ = BTE_1200{lastOK};
    BTN_1200_ = BTN_1200{lastOK};

    BT_x_1200 = BTx_direction_1200{lastOK};
    BT_y_1200 = BTy_direction_1200{lastOK};

    BT_1_1200 = BT_1_1200_ALL{lastOK};
    BT_2_1200 = BT_2_1200_ALL{lastOK};

    BTmodea_1200 = BTmodea_1200_{lastOK};
    BTmodff_1200 = BTmodf_1200ALL{lastOK};
    BTerror = BTerror1200all{lastOK};
    WaterBias = WaterBiasALL{lastOK};

    start = startAll(lastOK);
    finish = finishAll(lastOK);
end

%% Build output structure

BTbedloadOut = struct();

BTbedloadOut.InputRow = InputRow;
BTbedloadOut.ADCPfilename = ADCPfilename;
BTbedloadOut.StartADCP = StartADCP;
BTbedloadOut.EndADCP = EndADCP;

BTbedloadOut.DistanceFromLeftBank = DistanceFromLeftBank;
BTbedloadOut.TrateMean = TrateMean;
BTbedloadOut.MaxTrate = MaxTrate;
BTbedloadOut.D50 = D50;
BTbedloadOut.D90 = D90;
BTbedloadOut.DmeaNUMpsd = DmeaNUMpsd;
BTbedloadOut.SSC = SSC;

BTbedloadOut.CoordinateSystem = CoordinateSystem;

BTbedloadOut.filtOUTwin = opt.filtOUTwin;
BTbedloadOut.filtOUTstd = opt.filtOUTstd;
BTbedloadOut.movVelTresh = opt.movVelTresh;
BTbedloadOut.RDICutOFFtrim = opt.RDICutOFFtrim;
BTbedloadOut.Ws_min = opt.Ws_min;
BTbedloadOut.ExtraFilter = opt.ExtraFilter;
BTbedloadOut.DummyValue = opt.DummyValue;
BTbedloadOut.BTerrorMultiplier = opt.BTerrorMultiplier;
BTbedloadOut.WaterBiasReference = WaterBiasReference;
BTbedloadOut.WaterBiasComparison = WaterBiasComparison;

% Last successful row old names
BTbedloadOut.BTx_ = BTx_;
BTbedloadOut.BTy_ = BTy_;

BTbedloadOut.BTE_1200_ = BTE_1200_;
BTbedloadOut.BTN_1200_ = BTN_1200_;

BTbedloadOut.BT_x_1200 = BT_x_1200;
BTbedloadOut.BT_y_1200 = BT_y_1200;

BTbedloadOut.BT_1_1200 = BT_1_1200;
BTbedloadOut.BT_2_1200 = BT_2_1200;

BTbedloadOut.BTmodea_1200 = BTmodea_1200;
BTbedloadOut.BTmodff_1200 = BTmodff_1200;
BTbedloadOut.BTerror = BTerror;
BTbedloadOut.WaterBias = WaterBias;

% Scalars
BTbedloadOut.NANsbeforeX1200 = NANsbeforeX1200;
BTbedloadOut.NANsbeforeY1200 = NANsbeforeY1200;
BTbedloadOut.NANsbeforeX1200perc = NANsbeforeX1200perc;
BTbedloadOut.NANsbeforeY1200perc = NANsbeforeY1200perc;

BTbedloadOut.BTmod1200en = BTmod1200en;
BTbedloadOut.BTmod1200org = BTmod1200org;
BTbedloadOut.BTmod1200_ = BTmod1200_;

BTbedloadOut.BTxSTD_1200 = BTxSTD_1200;
BTbedloadOut.BTySTD_1200 = BTySTD_1200;
BTbedloadOut.BTstd_1200 = BTstd_1200;

BTbedloadOut.BTmod1200a = BTmod1200a;
BTbedloadOut.BTstd1200a = BTstd1200a;

BTbedloadOut.BTerror1200m = BTerror1200m;
BTbedloadOut.BTerror1200all = BTerror1200all;

BTbedloadOut.BTstd1200aW = BTstd1200aW;
BTbedloadOut.BTmod1200aW = BTmod1200aW;
BTbedloadOut.BTmod1200aW_ = BTmod1200aW_;

BTbedloadOut.FilteredDataN_1200 = FilteredDataN_1200;
BTbedloadOut.FilteredDataE_1200 = FilteredDataE_1200;

BTbedloadOut.BTx_1200f = BTx_1200f;
BTbedloadOut.BTy_1200f = BTy_1200f;
BTbedloadOut.BTmodf_1200 = BTmodf_1200;
BTbedloadOut.BT_stdf_1200 = BT_stdf_1200;
BTbedloadOut.BTmodf_1200ALL = BTmodf_1200ALL;

BTbedloadOut.FilteredData1_1200 = FilteredData1_1200;
BTbedloadOut.FilteredData2_1200 = FilteredData2_1200;

BTbedloadOut.BT_1_1200mean = BT_1_1200mean;
BTbedloadOut.BT_2_1200mean = BT_2_1200mean;
BTbedloadOut.BT_1_1200std = BT_1_1200std;
BTbedloadOut.BT_2_1200std = BT_2_1200std;

% va outputs
BTbedloadOut.va_raw_1200 = va_raw_1200;
BTbedloadOut.va_std_1200 = va_std_1200;
BTbedloadOut.va_boat_1200 = va_boat_1200;
BTbedloadOut.va_median_1200 = va_median_1200;
BTbedloadOut.va_extra_1200 = va_extra_1200;
BTbedloadOut.va_error_1200 = va_error_1200;
BTbedloadOut.va_gps_1200 = va_gps_1200;
BTbedloadOut.va_attitude_1200 = va_attitude_1200;
BTbedloadOut.va_waterbias_1200 = va_waterbias_1200;
BTbedloadOut.va_direction_1200 = va_direction_1200;
BTbedloadOut.va_final_1200 = va_final_1200;

BTbedloadOut.BTvaStepNames = BTvaStepNames;
BTbedloadOut.BTvaAllSteps_1200 = BTvaAllSteps_1200;
BTbedloadOut.BTvaTrateCorrTable = BTvaTrateCorrTable;

% Old cell outputs
BTbedloadOut.BTx_1200 = BTx_1200;
BTbedloadOut.BTy_1200 = BTy_1200;

BTbedloadOut.BTx_1200a = BTx_1200a;
BTbedloadOut.BTy_1200a = BTy_1200a;

BTbedloadOut.BTE_1200 = BTE_1200;
BTbedloadOut.BTN_1200 = BTN_1200;

BTbedloadOut.BTE_1200a = BTE_1200a;
BTbedloadOut.BTN_1200a = BTN_1200a;

% Step-by-step time series
BTbedloadOut.BTx_raw_1200 = BTx_raw_1200;
BTbedloadOut.BTy_raw_1200 = BTy_raw_1200;

BTbedloadOut.BTx_std_1200 = BTx_std_1200;
BTbedloadOut.BTy_std_1200 = BTy_std_1200;

BTbedloadOut.BTx_boat_1200 = BTx_boat_1200;
BTbedloadOut.BTy_boat_1200 = BTy_boat_1200;

BTbedloadOut.BTx_median_1200 = BTx_median_1200;
BTbedloadOut.BTy_median_1200 = BTy_median_1200;

BTbedloadOut.BTx_extra_1200 = BTx_extra_1200;
BTbedloadOut.BTy_extra_1200 = BTy_extra_1200;

BTbedloadOut.BTx_error_1200 = BTx_error_1200;
BTbedloadOut.BTy_error_1200 = BTy_error_1200;

BTbedloadOut.BTx_gps_1200 = BTx_gps_1200;
BTbedloadOut.BTy_gps_1200 = BTy_gps_1200;

BTbedloadOut.BTx_attitude_1200 = BTx_attitude_1200;
BTbedloadOut.BTy_attitude_1200 = BTy_attitude_1200;

BTbedloadOut.BTx_waterbias_1200 = BTx_waterbias_1200;
BTbedloadOut.BTy_waterbias_1200 = BTy_waterbias_1200;

BTbedloadOut.BTx_direction_1200 = BTx_direction_1200;
BTbedloadOut.BTy_direction_1200 = BTy_direction_1200;

BTbedloadOut.BTx_final_1200 = BTx_final_1200;
BTbedloadOut.BTy_final_1200 = BTy_final_1200;

% Step-by-step magnitudes
BTbedloadOut.BTmod_raw_1200 = BTmod_raw_1200;
BTbedloadOut.BTmod_std_1200 = BTmod_std_1200;
BTbedloadOut.BTmod_boat_1200 = BTmod_boat_1200;
BTbedloadOut.BTmod_median_1200 = BTmod_median_1200;
BTbedloadOut.BTmod_extra_1200 = BTmod_extra_1200;
BTbedloadOut.BTmod_error_1200 = BTmod_error_1200;
BTbedloadOut.BTmod_gps_1200 = BTmod_gps_1200;
BTbedloadOut.BTmod_attitude_1200 = BTmod_attitude_1200;
BTbedloadOut.BTmod_waterbias_1200 = BTmod_waterbias_1200;
BTbedloadOut.BTmod_direction_1200 = BTmod_direction_1200;
BTbedloadOut.BTmod_final_1200 = BTmod_final_1200;

% E-N outputs
BTbedloadOut.BTE_raw_1200 = BTE_raw_1200;
BTbedloadOut.BTN_raw_1200 = BTN_raw_1200;

BTbedloadOut.BTE_std_1200 = BTE_std_1200;
BTbedloadOut.BTN_std_1200 = BTN_std_1200;

BTbedloadOut.BTE_boat_1200 = BTE_boat_1200;
BTbedloadOut.BTN_boat_1200 = BTN_boat_1200;

BTbedloadOut.BTE_median_1200 = BTE_median_1200;
BTbedloadOut.BTN_median_1200 = BTN_median_1200;

BTbedloadOut.BTE_extra_1200 = BTE_extra_1200;
BTbedloadOut.BTN_extra_1200 = BTN_extra_1200;

BTbedloadOut.BTE_gps_1200 = BTE_gps_1200;
BTbedloadOut.BTN_gps_1200 = BTN_gps_1200;

% Water profile reused outputs
BTbedloadOut.V1x_1_BT_1200 = V1x_1_BT_1200;
BTbedloadOut.V1y_1_BT_1200 = V1y_1_BT_1200;
BTbedloadOut.VwaterBiasX_1200 = VwaterBiasX_1200;
BTbedloadOut.VwaterBiasY_1200 = VwaterBiasY_1200;
BTbedloadOut.VwaterBiasMod_1200 = VwaterBiasMod_1200;

% Other outputs
% Exact aligned masks for reuse by processBackscatter1200 and other
% ensemble-wise processing. true in BTfilterMask_1200 means rejected.
BTbedloadOut.BTfilterMask_1200 = BTfilterMask_1200;
BTbedloadOut.BTvalidMask_1200 = BTvalidMask_1200;
BTbedloadOut.BTfilterMaskLength_1200 = BTfilterMaskLength_1200;
BTbedloadOut.NfilteredFinal_1200 = NfilteredFinal_1200;
BTbedloadOut.NvalidFinal_1200 = NvalidFinal_1200;
BTbedloadOut.BTfilterMaskConvention = "true_is_rejected";
BTbedloadOut.BTfilterMaskStart_1200 = startAll;
BTbedloadOut.BTfilterMaskEnd_1200 = finishAll;

BTbedloadOut.peaksBTxy_1200 = peaksBTxy_1200;
BTbedloadOut.peaksBT_EN_S = peaksBT_EN_S;

BTbedloadOut.idxBTer_1200 = idxBTer_1200;
BTbedloadOut.GPSidxOUT_1200 = GPSidxOUT_1200;
BTbedloadOut.AttitudeidxOUT_1200 = AttitudeidxOUT_1200;

BTbedloadOut.WaterBiasALL = WaterBiasALL;

BTbedloadOut.ErrorVel1200valid = ErrorVel1200valid;

BTbedloadOut.alphaV_1200 = alphaV_1200;
BTbedloadOut.alphaBTf_1200 = alphaBTf_1200;
BTbedloadOut.alphaBT_1200 = alphaBT_1200;
BTbedloadOut.VmodALL_1200 = VmodALL_1200;

BTbedloadOut.BT_1_1200_ALL = BT_1_1200_ALL;
BTbedloadOut.BT_2_1200_ALL = BT_2_1200_ALL;

BTbedloadOut.alphaV12_1200 = alphaV12_1200;
BTbedloadOut.alphaBT12 = alphaBT12;
BTbedloadOut.Vmod12ins_1200 = Vmod12ins_1200;

BTbedloadOut.BTbedloadAll = BTbedloadAll;
BTbedloadOut.BTbedloadTable = BTbedloadTable;

BTbedloadOut.FileFound = FileFound;
BTbedloadOut.MatchedFile = MatchedFile;
BTbedloadOut.Message = Message;

BTbedloadOut.start = start;
BTbedloadOut.finish = finish;
BTbedloadOut.startAll = startAll;
BTbedloadOut.finishAll = finishAll;

%% Assign to workspace

if opt.AssignToCaller
    assignStructToWorkspace(BTbedloadOut, 'caller');
end

if opt.AssignToBase
    assignStructToWorkspace(BTbedloadOut, 'base');
end

%% Display

if opt.Verbose
    fprintf('\n=== BT apparent bedload velocity processing summary ===\n');
    fprintf('Input rows: %d\n', nFiles);
    fprintf('Matched files: %d/%d\n', sum(FileFound), nFiles);
    fprintf('Successfully processed: %d/%d\n', sum(Message == "OK"), nFiles);
    fprintf('CoordinateSystem: %s\n', CoordinateSystem);
    fprintf('BTerrorMultiplier: %.3f\n', opt.BTerrorMultiplier);
    fprintf('WaterBiasReference: %s\n', WaterBiasReference);
    fprintf('WaterBiasComparison: %s\n', WaterBiasComparison);
    fprintf('Water-profile quantities were reused from WaterProfileOut.\n');
    fprintf('Every filter step was saved as workspace cell arrays.\n');
    fprintf('Exact final BT masks exported: BTfilterMask_1200 (true = rejected).\n');
    fprintf('Final retained ensembles: %d; rejected ensembles: %d.\n', ...
        sum(NvalidFinal_1200,'omitnan'),sum(NfilteredFinal_1200,'omitnan'));
    fprintf('va values for every filter step were saved and correlated with TrateMean/MaxTrate.\n');
    fprintf('StartADCP/EndADCP were used as trim counts.\n');
    fprintf('No files were saved.\n');
    fprintf('=====================================================\n');

    disp(BTbedloadTable);
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


function x = getVectorFieldAny(S, fieldNames, defaultValue, n)

x = defaultValue;

for k = 1:numel(fieldNames)
    if isfield(S, fieldNames{k})
        x = S.(fieldNames{k});
        break
    end
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


function n = getBTvelocityLength(BT)

lens = [];

possibleFields = {'BTV1','BTV2','BTV4','btvE','btvN'};

for k = 1:numel(possibleFields)
    if isfield(BT, possibleFields{k})
        lens(end+1) = numel(BT.(possibleFields{k})); %#ok<AGROW>
    end
end

if isempty(lens)
    error('No usable BT velocity fields found. Expected BT.BTV1/BT.BTV2 or BT.btvE/BT.btvN.');
end

n = min(lens);

end


function [BTx_raw, BTy_raw] = getBTmainComponents(BT, start, finish, CoordinateSystem)

BTx_raw = [];
BTy_raw = [];

if CoordinateSystem == "earth"
    if isfield(BT, 'btvE') && isfield(BT, 'btvN')
        BTx_raw = cropVectorSafe(BT.btvE, start, finish);
        BTy_raw = cropVectorSafe(BT.btvN, start, finish);
        return
    else
        warning('Earth BT fields BT.btvE/BT.btvN are missing. Falling back to BT.BTV1/BT.BTV2.');
    end
end

if isfield(BT, 'BTV1') && isfield(BT, 'BTV2')
    BTx_raw = cropVectorSafe(BT.BTV1, start, finish);
    BTy_raw = cropVectorSafe(BT.BTV2, start, finish);
elseif isfield(BT, 'btvE') && isfield(BT, 'btvN')
    warning('BT.BTV1/BT.BTV2 missing. Falling back to BT.btvE/BT.btvN.');
    BTx_raw = cropVectorSafe(BT.btvE, start, finish);
    BTy_raw = cropVectorSafe(BT.btvN, start, finish);
else
    error('Could not find BT velocity components. Need BT.BTV1/BT.BTV2 or BT.btvE/BT.btvN.');
end

end


function [BTE_raw, BTN_raw] = getBTearthComponents(BT, start, finish)

BTE_raw = [];
BTN_raw = [];

if isfield(BT, 'btvE') && isfield(BT, 'btvN')
    BTE_raw = cropVectorSafe(BT.btvE, start, finish);
    BTN_raw = cropVectorSafe(BT.btvN, start, finish);
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


function x = padVectorToLength(x, n)

x = double(x(:));

if numel(x) < n
    x = [x; nan(n-numel(x),1)];
elseif numel(x) > n
    x = x(1:n);
end

end


function m = calcMod(x, y)

if isempty(x) || isempty(y)
    m = [];
    return
end

n = min(numel(x), numel(y));

x = double(x(:));
y = double(y(:));

m = sqrt(x(1:n).^2 + y(1:n).^2);

end


function BT = despikeBT(BTv,T)

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


function [BTx_1,BTy_1,peaks] = despikeBTmoveBoattresh(BTx_,BTy_,T)

BTx_ = double(BTx_(:));
BTy_ = double(BTy_(:));

BTv = sqrt(BTx_.^2 + BTy_.^2);
A = BTv;

T1 = T;
To = -T;

peaks = (A > T1 | A < To);

BTx_1 = BTx_;
BTy_1 = BTy_;

BTx_1(peaks) = NaN;
BTy_1(peaks) = NaN;

end


function BTerror = getBTerrorVelocity(BT, start, finish)

BTerror = [];

if isfield(BT, 'BTV4')
    BTerror = cropVectorSafe(BT.BTV4, start, finish);
elseif isfield(BT, 'btvU')
    BTerror = cropVectorSafe(BT.btvU, start, finish);
end

end


function [V1x_1, V1y_1, V1w_1, V1e_1, z_1200, VbiasX, VbiasY, VbiasMod] = ...
    getWaterProfileInputsForBT(i, WaterProfileOut, WaterBiasReference)

if isempty(WaterProfileOut)
    error('WaterProfileOut is empty.');
end

if ~isfield(WaterProfileOut, 'V1x_1_') || ~isfield(WaterProfileOut, 'V1y_1_')
    error('WaterProfileOut must contain V1x_1_ and V1y_1_ from processWaterProfile1200.');
end

if numel(WaterProfileOut.V1x_1_) < i || isempty(WaterProfileOut.V1x_1_{i})
    error('WaterProfileOut.V1x_1_{%d} is missing or empty.', i);
end

V1x_1 = WaterProfileOut.V1x_1_{i};
V1y_1 = WaterProfileOut.V1y_1_{i};

if isfield(WaterProfileOut, 'V1w_1_') && numel(WaterProfileOut.V1w_1_) >= i
    V1w_1 = WaterProfileOut.V1w_1_{i};
else
    V1w_1 = NaN(size(V1x_1));
end

if isfield(WaterProfileOut, 'V1e_1_') && numel(WaterProfileOut.V1e_1_) >= i
    V1e_1 = WaterProfileOut.V1e_1_{i};
else
    V1e_1 = NaN(size(V1x_1));
end

if isfield(WaterProfileOut, 'z_1200_') && numel(WaterProfileOut.z_1200_) >= i
    z_1200 = WaterProfileOut.z_1200_{i};
else
    z_1200 = (1:size(V1x_1,1))';
end

switch WaterBiasReference

    case "lastcell"

        if isfield(WaterProfileOut, 'VlastCellX_1200_') && ...
                isfield(WaterProfileOut, 'VlastCellY_1200_') && ...
                numel(WaterProfileOut.VlastCellX_1200_) >= i && ...
                ~isempty(WaterProfileOut.VlastCellX_1200_{i})

            VbiasX = WaterProfileOut.VlastCellX_1200_{i};
            VbiasY = WaterProfileOut.VlastCellY_1200_{i};

        else
            VbiasX = V1x_1(end,:)';
            VbiasY = V1y_1(end,:)';
        end

    case "last2cells"

        if isfield(WaterProfileOut, 'Vlast2CellsX_1200_') && ...
                isfield(WaterProfileOut, 'Vlast2CellsY_1200_') && ...
                numel(WaterProfileOut.Vlast2CellsX_1200_) >= i && ...
                ~isempty(WaterProfileOut.Vlast2CellsX_1200_{i})

            VbiasX = WaterProfileOut.Vlast2CellsX_1200_{i};
            VbiasY = WaterProfileOut.Vlast2CellsY_1200_{i};

        else
            if size(V1x_1,1) >= 2
                VbiasX = nanmeanLocal(V1x_1(end-1:end,:), 1)';
                VbiasY = nanmeanLocal(V1y_1(end-1:end,:), 1)';
            else
                VbiasX = V1x_1(end,:)';
                VbiasY = V1y_1(end,:)';
            end
        end

    otherwise
        error('Unknown WaterBiasReference: %s', WaterBiasReference);
end

VbiasX = VbiasX(:);
VbiasY = VbiasY(:);
VbiasMod = sqrt(VbiasX.^2 + VbiasY.^2);

end


function [BTx_, BTy_, WaterBias] = applyWaterBiasFilterFromReference( ...
    BTx_, BTy_, VbiasX, VbiasY, WaterBiasComparison, WaterBiasUseBothComponents, WaterBiasMultiplier)

BTx_ = double(BTx_(:));
BTy_ = double(BTy_(:));

VbiasX = double(VbiasX(:));
VbiasY = double(VbiasY(:));

WaterBias = repmat({'NO bias'}, numel(BTx_), 1);

nEns = min([numel(BTx_), numel(BTy_), numel(VbiasX), numel(VbiasY)]);

WaterBiasComparison = lower(string(WaterBiasComparison));

for bb = 1:nEns

    removeNow = false;

    switch WaterBiasComparison

        case "component"

            IwbX = abs(BTx_(bb)) > WaterBiasMultiplier .* abs(VbiasX(bb));
            IwbY = abs(BTy_(bb)) > WaterBiasMultiplier .* abs(VbiasY(bb));

            if WaterBiasUseBothComponents
                removeNow = IwbX || IwbY;
            else
                removeNow = IwbX;
            end

        case "magnitude"

            BTmod = sqrt(BTx_(bb).^2 + BTy_(bb).^2);
            Vmod = sqrt(VbiasX(bb).^2 + VbiasY(bb).^2);

            removeNow = abs(BTmod) > WaterBiasMultiplier .* abs(Vmod);

        otherwise
            error('Unknown WaterBiasComparison: %s', WaterBiasComparison);
    end

    if removeNow
        WaterBias{bb} = 'possible water bias';
        BTx_(bb) = NaN;
        BTy_(bb) = NaN;
    end
end

end


function [BT_N,BT_E, FilteredDataN, FilteredDataE,alphaV,alphaBTf,alphaBT,Vmod] = ...
    filterBT2org(V_E, V_N, BT_E, BT_N, Fi, T, endCut, Nx, Ny, CutOff)

V_E = double(V_E);
V_N = double(V_N);

BT_E = double(BT_E(:));
BT_N = double(BT_N(:));

nEns = min([size(V_E,2), size(V_N,2), numel(BT_E), numel(BT_N)]);

V_E = V_E(:,1:nEns);
V_N = V_N(:,1:nEns);
BT_E = BT_E(1:nEns);
BT_N = BT_N(1:nEns);

BTmod = sqrt(BT_E.^2 + BT_N.^2);

alphaBT = atan2d(BT_N, BT_E);
alphaBT = wrapTo360local(alphaBT);

V_EE = nan(nEns,1);
V_NN = nan(nEns,1);

endCut = round(endCut);
endCut = max(2, min(endCut, size(V_E,1)));

for ii = 1:nEns

    vE = V_E(:,ii);
    vN = V_N(:,ii);

    vE = vE(2:endCut);
    vN = vN(2:endCut);

    V_EE(ii) = safeMean(vE);
    V_NN(ii) = safeMean(vN);
end

alphaV = atan2d(V_NN, V_EE);
alphaV = wrapTo360local(alphaV);

Vmod = sqrt(V_EE.^2 + V_NN.^2);

COSalpha = (V_EE .* BT_E + BT_N .* V_NN) ./ (Vmod .* BTmod);
COSalpha(isnan(COSalpha)) = -Inf;

BT_N(COSalpha < Fi) = NaN;
BT_E(COSalpha < Fi) = NaN;

T1N = safeMean(BT_N) + T .* nanstdLocal(BT_N, 0, 1);
ToN = safeMean(BT_N) - T .* nanstdLocal(abs(BT_N), 0, 1);

BT_N(BT_N < ToN | BT_N > T1N) = NaN;

T1E = safeMean(BT_E) + T .* nanstdLocal(BT_E, 0, 1);
ToE = safeMean(BT_E) - T .* nanstdLocal(BT_E, 0, 1);

BT_E(BT_E < ToE | BT_E > T1E) = NaN;

alphaBTf = atan2d(BT_N, BT_E);
alphaBTf = wrapTo360local(alphaBTf);

TotalNAN_N = sum(isnan(BT_N)) - Nx;
TotalNAN_E = sum(isnan(BT_E)) - Ny;

FilteredDataN = TotalNAN_N ./ numel(BT_N);
FilteredDataE = TotalNAN_E ./ numel(BT_E);

BTmod = sqrt(BT_E.^2 + BT_N.^2);
BTmod(BTmod < CutOff) = NaN;

end


function idxOUT = getGPSoutlierIndex(GPSdata, start, finish, TresholdGPS)

idxOUT = false(finish - start + 1, 1);

if isempty(GPSdata)
    return
end

if ischar(GPSdata) || isstring(GPSdata)
    return
end

if isnumeric(GPSdata) && all(isnan(GPSdata))
    return
end

idxOUTg = false(size(idxOUT));
idxOUTvtg = false(size(idxOUT));

if isstruct(GPSdata) && isfield(GPSdata, 'GPVfinal')
    GPV = cropVectorSafe(GPSdata.GPVfinal, start, finish);
    idxOUTg = abs(GPV) > TresholdGPS;
end

if isstruct(GPSdata) && isfield(GPSdata, 'VTGV1')
    VTG = cropVectorSafe(GPSdata.VTGV1, start, finish);
    idxOUTvtg = abs(VTG) > TresholdGPS;
elseif isstruct(GPSdata) && isfield(GPSdata, 'VTGv_ext')
    VTG = cropVectorSafe(GPSdata.VTGv_ext, start, finish);
    idxOUTvtg = abs(VTG) > TresholdGPS;
end

nCommon = min([numel(idxOUT), numel(idxOUTg), numel(idxOUTvtg)]);

idxOUT = false(nCommon,1);
idxOUT = logical(idxOUTg(1:nCommon)) | logical(idxOUTvtg(1:nCommon));

end


function [EH, Pitch, Roll] = getAttitudeSeries(S, start, finish)

EH = [];
Pitch = [];
Roll = [];

if isfield(S, 'Data')
    Data = S.Data;

    EH = getFirstExistingVector(Data, {'EH','Heading','heading','HEADING'}, start, finish);
    Pitch = getFirstExistingVector(Data, {'Pitch','PITCH','pitch'}, start, finish);
    Roll = getFirstExistingVector(Data, {'Roll','ROLL','roll'}, start, finish);
end

if isempty(EH) && isfield(S, 'Sensor')
    EH = getFirstExistingVector(S.Sensor, {'EH','Heading','heading','HEADING'}, start, finish);
    Pitch = getFirstExistingVector(S.Sensor, {'Pitch','PITCH','pitch'}, start, finish);
    Roll = getFirstExistingVector(S.Sensor, {'Roll','ROLL','roll'}, start, finish);
end

end


function x = getFirstExistingVector(S, fields, start, finish)

x = [];

for k = 1:numel(fields)
    if isfield(S, fields{k})
        x = cropVectorSafe(S.(fields{k}), start, finish);
        return
    end
end

end


function idx = getAttitudeOutlierIndex(EH, Pitch, Roll)

n = max([numel(EH), numel(Pitch), numel(Roll), 1]);
idx = false(n,1);

if ~isempty(EH)
    idxEH = isoutlierSafe(EH, 'quartiles');
    idx(1:numel(idxEH)) = idx(1:numel(idxEH)) | idxEH(:);
end

if ~isempty(Pitch)
    idxP = isoutlierSafe(Pitch, 'quartiles');
    idx(1:numel(idxP)) = idx(1:numel(idxP)) | idxP(:);
end

if ~isempty(Roll)
    idxR = isoutlierSafe(Roll, 'quartiles');
    idx(1:numel(idxR)) = idx(1:numel(idxR)) | idxR(:);
end

end


function y = filloutliersSafe(x, fillValue, method, param)

x = double(x(:));

try
    y = filloutliers(x, fillValue, method, param);
catch
    y = x;

    switch lower(string(method))

        case "movmedian"

            win = param;
            medx = movmedian(x, win, 'omitnan');
            dev = abs(x - medx);
            madx = movmedian(dev, win, 'omitnan');
            idx = dev > 3 .* madx;
            y(idx) = fillValue;

        case "percentiles"

            xx = x(~isnan(x));

            if isempty(xx)
                return
            end

            prc = prctile(xx, param);
            idx = x < prc(1) | x > prc(2);
            y(idx) = fillValue;
    end
end

end


function idx = isoutlierSafe(x, method)

x = double(x(:));

try
    idx = isoutlier(x, method);
catch
    switch lower(string(method))
        case "quartiles"

            xx = x(~isnan(x));

            if isempty(xx)
                idx = false(size(x));
                return
            end

            q1 = prctile(xx, 25);
            q3 = prctile(xx, 75);
            iqrVal = q3 - q1;
            L = q1 - 1.5 .* iqrVal;
            U = q3 + 1.5 .* iqrVal;
            idx = x < L | x > U;

        otherwise
            idx = false(size(x));
    end
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


function a = wrapTo360local(a)

a = mod(a, 360);
a(a < 0) = a(a < 0) + 360;

end


function tf = shouldPlotThisRow(i, PlotRows)

if isempty(PlotRows)
    tf = true;
else
    tf = any(i == PlotRows);
end

end


function plotBTfilterDiagnostics( ...
    i, ...
    thisADCPfilename, ...
    BTx_raw, ...
    BTy_raw, ...
    BTx_afterMedian, ...
    BTy_afterMedian, ...
    BT_x_1200, ...
    BT_y_1200, ...
    PlotRawFiltered, ...
    PlotFinalBT)

nPanels = double(logical(PlotRawFiltered)) + double(logical(PlotFinalBT));

if nPanels == 0
    return
end

figure('Name', sprintf('BT filtering diagnostics row %d', i), 'Color', 'w');

panel = 0;

if PlotRawFiltered

    panel = panel + 1;
    subplot(nPanels,1,panel);

    plot(BTx_raw, 'DisplayName', 'BTx raw'); hold on
    plot(BTy_raw, 'DisplayName', 'BTy raw');
    plot(BTx_afterMedian, 'DisplayName', 'BTx after median');
    plot(BTy_afterMedian, 'DisplayName', 'BTy after median');

    grid on
    xlabel('Trimmed ensemble index');
    ylabel('BT velocity');
    title(sprintf('BT raw and intermediate filtering — row %d: %s', ...
        i, cleanTitleString(thisADCPfilename)), ...
        'Interpreter', 'none');
    legend('Location','best');
end

if PlotFinalBT

    panel = panel + 1;
    subplot(nPanels,1,panel);

    plot(BT_x_1200, 'DisplayName', 'BT\_x\_1200 final'); hold on
    plot(BT_y_1200, 'DisplayName', 'BT\_y\_1200 final');
    plot(sqrt(BT_x_1200.^2 + BT_y_1200.^2), 'DisplayName', 'BTmod final');

    grid on
    xlabel('Trimmed ensemble index');
    ylabel('BT velocity');
    title('Final BT apparent bedload velocity after selected filters', 'Interpreter','none');
    legend('Location','best');
end

end


function BTvaTrateCorrTable = makeVaTrateCorrTable( ...
    BTvaStepNames, ...
    BTvaAllSteps_1200, ...
    TrateMean, ...
    MaxTrate, ...
    minN)

nSteps = numel(BTvaStepNames);

R_MeanTrate = nan(nSteps,1);
N_MeanTrate = zeros(nSteps,1);
R_MaxTrate = nan(nSteps,1);
N_MaxTrate = zeros(nSteps,1);

for k = 1:nSteps

    va = BTvaAllSteps_1200(:,k);

    [R_MeanTrate(k), N_MeanTrate(k)] = corrNoNaN(TrateMean, va, minN);
    [R_MaxTrate(k),  N_MaxTrate(k)]  = corrNoNaN(MaxTrate,  va, minN);
end

BTvaTrateCorrTable = table( ...
    BTvaStepNames(:), ...
    R_MeanTrate, ...
    N_MeanTrate, ...
    R_MaxTrate, ...
    N_MaxTrate, ...
    'VariableNames', { ...
    'Step', ...
    'R_MeanTrate', ...
    'N_MeanTrate', ...
    'R_MaxTrate', ...
    'N_MaxTrate'});

end


function plotVaVsTrateAllSteps( ...
    BTvaStepNames, ...
    BTvaAllSteps_1200, ...
    TrateMean, ...
    MaxTrate, ...
    PlotVaTrateSteps, ...
    useLogX, ...
    useLogY, ...
    minN)

stepIdx = selectVaPlotSteps(BTvaStepNames, PlotVaTrateSteps);

if isempty(stepIdx)
    warning('PlotVaTrate is true, but no valid plotting steps were selected.');
    return
end

nSteps = numel(stepIdx);
nRows = ceil(sqrt(nSteps));
nCols = ceil(nSteps ./ nRows);

figure('Name', 'BT apparent velocity va vs transport rate', 'Color', 'w');

try
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    useTiled = true;
catch
    useTiled = false;
end

for jj = 1:nSteps

    k = stepIdx(jj);
    stepName = BTvaStepNames(k);
    va = BTvaAllSteps_1200(:,k);

    [rMean, nMean] = corrNoNaN(TrateMean, va, minN);
    [rMax,  nMax]  = corrNoNaN(MaxTrate,  va, minN);

    if useTiled
        nexttile;
    else
        subplot(nRows, nCols, jj);
    end

    h1 = plot(TrateMean, va, 'o', 'DisplayName', 'MeanTrate');
    hold on
    h2 = plot(MaxTrate, va, '^', 'DisplayName', 'MaxTrate');

    grid on
    xlabel('Transport rate');
    ylabel('v_a');
    title(sprintf('v_a %s', stepName), 'Interpreter', 'none');

    if useLogX
        set(gca, 'XScale', 'log');
    end

    if useLogY
        set(gca, 'YScale', 'log');
    end

    legend([h1 h2], 'Location', 'best');

    txt = sprintf('r Mean = %.2f (n=%d)\nr Max = %.2f (n=%d)', ...
        rMean, nMean, rMax, nMax);

    addTextInsideAxes(txt);
end

end


function stepIdx = selectVaPlotSteps(BTvaStepNames, PlotVaTrateSteps)

if isempty(PlotVaTrateSteps)
    stepIdx = 1:numel(BTvaStepNames);
    return
end

if isnumeric(PlotVaTrateSteps)
    stepIdx = PlotVaTrateSteps(:)';
    stepIdx = stepIdx(stepIdx >= 1 & stepIdx <= numel(BTvaStepNames));
    return
end

steps = lower(string(PlotVaTrateSteps(:)));
allSteps = lower(string(BTvaStepNames(:)));

stepIdx = [];

for i = 1:numel(steps)
    idx = find(allSteps == steps(i), 1, 'first');

    if ~isempty(idx)
        stepIdx(end+1) = idx; %#ok<AGROW>
    end
end

stepIdx = unique(stepIdx, 'stable');

end


function [r, n] = corrNoNaN(x, y, minN)

x = double(x(:));
y = double(y(:));

nCommon = min(numel(x), numel(y));

x = x(1:nCommon);
y = y(1:nCommon);

idx = ~isnan(x) & ~isnan(y) & isfinite(x) & isfinite(y);

n = sum(idx);

if n < minN
    r = NaN;
    return
end

C = corrcoef(x(idx), y(idx));
r = C(1,2);

end


function addTextInsideAxes(txt)

xl = xlim;
yl = ylim;

if strcmpi(get(gca, 'XScale'), 'log')
    xpos = 10.^(log10(xl(1)) + 0.05 .* diff(log10(xl)));
else
    xpos = xl(1) + 0.05 .* diff(xl);
end

if strcmpi(get(gca, 'YScale'), 'log')
    ypos = 10.^(log10(yl(2)) - 0.08 .* diff(log10(yl)));
else
    ypos = yl(2) - 0.08 .* diff(yl);
end

text(xpos, ypos, txt, ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'Margin', 4);

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