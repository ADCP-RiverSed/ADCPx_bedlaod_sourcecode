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
%   - extracts the RDI/M9 vertical-beam depth when available
%   - converts vertical-beam depth to metres and rejects values > 10 m
%   - creates cell arrays such as RNG11_{i}, RNG22_{i}, etc.
%   - optionally plots raw/despiked beam depth series
%   - keeps time series only in cell arrays, one cell per measurement
%   - does not expose loose variables from only the last processed file
%
% NaN handling:
%   - NaN samples are omitted from all mean, standard-deviation, and
%     spatial depth calculations.
%   - RNG1_1200 to RNG4_1200 and Depth1200 use the despiked series.
%   - If despiking leaves an entire beam without a valid sample, its raw
%     cropped series is used as a fallback, still omitting NaN values.
%   - A mean remains NaN only when the complete cropped beam contains no
%     valid non-NaN sample.
%
% Required inputs:
%   dataADCP1200 = dir(filePattern1200);
%   ADCPin       = output from readInputTable_PhSamples
%
% Example:
%   DepthRNGout = processDepthsRNG(dataADCP1200, ADCPin, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
% Main per-measurement workspace variables created:
%   WF_12, WS_12, WT_12, Temp1 (vectors, one value per measurement)
%   RNG1_1200, RNG2_1200, RNG3_1200, RNG4_1200
%   RNG11_, RNG22_, RNG33_, RNG44_
%   RNG11_1200_, RNG22_1200_, RNG33_1200_, RNG44_1200_
%   Depth1200, Depth1200std
%   Depth1200series_
%   stdSpatial_1200, stdSpatialALL_1200
%   VBDepth, VBDepthstd, VBDepth_, VBAvailable
%   VBSource, VBInstrument, VBOriginalUnit, VBScaleToMeters
%   DepthRNGall, DepthRNGtable, DepthRNGout

%% Parse inputs

p = inputParser;

addRequired(p, 'dataADCP1200', @(x) isstruct(x));
addRequired(p, 'ADCPin', @(x) isstruct(x));

addParameter(p, 'DespikeT', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'VBMaxDepth', 10, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

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

% Temp1 was previously overwritten in every loop iteration and therefore
% represented only the last measurement. It is now one value per file.
Temp1 = nan(nFiles,1);

RNG1_1200 = nan(nFiles,1);
RNG2_1200 = nan(nFiles,1);
RNG3_1200 = nan(nFiles,1);
RNG4_1200 = nan(nFiles,1);

Depth1200 = nan(nFiles,4);
Depth1200std = nan(nFiles,4);
Depth1200Mean = nan(nFiles,1);

stdSpatial_1200 = nan(nFiles,1);
stdSpatialALL_1200 = cell(nFiles,1);

RNG11_ = cell(nFiles,1);
RNG22_ = cell(nFiles,1);
RNG33_ = cell(nFiles,1);
RNG44_ = cell(nFiles,1);

RNG11_1200_ = cell(nFiles,1);
RNG22_1200_ = cell(nFiles,1);
RNG33_1200_ = cell(nFiles,1);
RNG44_1200_ = cell(nFiles,1);

Depth1200series_ = cell(nFiles,1);

% Vertical-beam outputs. VBDepth_ is the only stored vertical-beam time
% series and contains processed depth in metres, one cell per measurement.
VBDepth_ = cell(nFiles,1);
VBDepth = nan(nFiles,1);
VBDepthstd = nan(nFiles,1);
VBAvailable = false(nFiles,1);
VBSource = strings(nFiles,1);
VBInstrument = strings(nFiles,1);
VBOriginalUnit = strings(nFiles,1);
VBScaleToMeters = nan(nFiles,1);
VBDifferenceFromFourBeam = nan(nFiles,1);
VBRelativeDifference = nan(nFiles,1);
VBMessage = repmat("Vertical beam not available", nFiles, 1);

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
    'Depth1200', [], ...
    'Depth1200std', [], ...
    'Depth1200Mean', [], ...
    'stdSpatial_1200', [], ...
    'VBAvailable', false, ...
    'VBDepth', [], ...
    'VBDepthstd', [], ...
    'VBSource', "", ...
    'VBInstrument', "", ...
    'VBOriginalUnit', "", ...
    'VBScaleToMeters', [], ...
    'VBDifferenceFromFourBeam', [], ...
    'VBRelativeDifference', [], ...
    'VBMessage', "", ...
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
            Temp1(i) = WT_12(i);
        catch
            Temp1(i) = NaN;
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
        % NaN values are explicitly removed before calculating each mean.
        % If a complete despiked beam is NaN, the raw cropped beam is used.

        RNG1_1200(i) = safeMeanWithFallback(RNG11_1200, RNG11);
        RNG2_1200(i) = safeMeanWithFallback(RNG22_1200, RNG22);
        RNG3_1200(i) = safeMeanWithFallback(RNG33_1200, RNG33);
        RNG4_1200(i) = safeMeanWithFallback(RNG44_1200, RNG44);

        %% Final depth statistics from despiked data, omitting NaNs
        % Keep one authoritative mean depth per beam. These values are
        % based on the despiked series and use the raw cropped series only
        % when despiking leaves no valid sample in the complete beam.

        Depth1200(i,:) = [ ...
            RNG1_1200(i), ...
            RNG2_1200(i), ...
            RNG3_1200(i), ...
            RNG4_1200(i)];

        % Calculate temporal standard deviation from the despiked series.
        % If the despiked series cannot produce a standard deviation, use
        % the raw cropped series as fallback. NaNs are omitted in both.
        Depth1200std(i,1) = safeStdWithFallback(RNG11_1200, RNG11, 0);
        Depth1200std(i,2) = safeStdWithFallback(RNG22_1200, RNG22, 0);
        Depth1200std(i,3) = safeStdWithFallback(RNG33_1200, RNG33, 0);
        Depth1200std(i,4) = safeStdWithFallback(RNG44_1200, RNG44, 0);

        % Use despiked beam series for spatial variability. A raw cropped
        % beam is substituted only when its entire despiked series is NaN.
        R1 = seriesWithFallback(RNG11_1200, RNG11);
        R2 = seriesWithFallback(RNG22_1200, RNG22);
        R3 = seriesWithFallback(RNG33_1200, RNG33);
        R4 = seriesWithFallback(RNG44_1200, RNG44);

        R = [R1(:), R2(:), R3(:), R4(:)];

        % Standard deviation across the four beams for every sample.
        % NaNs are omitted. A sample remains NaN only when fewer than two
        % beams contain valid values because sample std uses flag = 0.
        stdSpatialALL_1200{i} = nanstdLocal(R, 0, 2);

        % Mean spatial standard deviation, omitting NaN samples.
        stdSpatial_1200(i) = safeMean(stdSpatialALL_1200{i});

        % Mean of the four slanted beams, used as the independent reference
        % for vertical-beam unit detection.
        Depth1200Mean(i) = safeMean(Depth1200(i,:));

        %% Vertical beam: RDI and SonTek M9

        [VBfull, VBSource(i), VBInstrument(i)] = ...
            findVerticalBeamSeries(S);

        if ~isempty(VBfull)
            VBcropped = cropVectorSafe(VBfull, start, finish);
            VBcropped(VBcropped == 0) = NaN;

            [VBmetres, VBOriginalUnit(i), VBScaleToMeters(i), ...
                unitMessage] = convertVerticalBeamToMetres( ...
                VBcropped, ...
                Depth1200Mean(i), ...
                opt.VBMaxDepth, ...
                VBSource(i));

            if any(isfinite(VBmetres))
                VBdespiked = despikeBT(VBmetres, opt.DespikeT);
                VBDepth_{i} = seriesWithFallback(VBdespiked, VBmetres);

                VBDepth(i) = safeMeanWithFallback(VBdespiked, VBmetres);
                VBDepthstd(i) = safeStdWithFallback( ...
                    VBdespiked, VBmetres, 0);
                VBAvailable(i) = isfinite(VBDepth(i));

                if VBAvailable(i) && isfinite(Depth1200Mean(i))
                    VBDifferenceFromFourBeam(i) = ...
                        VBDepth(i) - Depth1200Mean(i);
                    VBRelativeDifference(i) = ...
                        abs(VBDifferenceFromFourBeam(i)) ./ ...
                        max(abs(Depth1200Mean(i)), eps);
                end

                VBMessage(i) = unitMessage;
            else
                VBDepth_{i} = nan(size(VBcropped));
                VBMessage(i) = unitMessage;
            end
        else
            VBDepth_{i} = [];
        end

        %% Store full cropped series

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
                VBDepth_{i}, ...
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
        DepthRNGall(i).Temp1 = Temp1(i);

        DepthRNGall(i).RNG1_1200 = RNG1_1200(i);
        DepthRNGall(i).RNG2_1200 = RNG2_1200(i);
        DepthRNGall(i).RNG3_1200 = RNG3_1200(i);
        DepthRNGall(i).RNG4_1200 = RNG4_1200(i);

        DepthRNGall(i).Depth1200 = Depth1200(i,:);
        DepthRNGall(i).Depth1200std = Depth1200std(i,:);
        DepthRNGall(i).Depth1200Mean = Depth1200Mean(i);

        DepthRNGall(i).stdSpatial_1200 = stdSpatial_1200(i);

        DepthRNGall(i).VBAvailable = VBAvailable(i);
        DepthRNGall(i).VBDepth = VBDepth(i);
        DepthRNGall(i).VBDepthstd = VBDepthstd(i);
        DepthRNGall(i).VBSource = VBSource(i);
        DepthRNGall(i).VBInstrument = VBInstrument(i);
        DepthRNGall(i).VBOriginalUnit = VBOriginalUnit(i);
        DepthRNGall(i).VBScaleToMeters = VBScaleToMeters(i);
        DepthRNGall(i).VBDifferenceFromFourBeam = ...
            VBDifferenceFromFourBeam(i);
        DepthRNGall(i).VBRelativeDifference = VBRelativeDifference(i);
        DepthRNGall(i).VBMessage = VBMessage(i);

        DepthRNGall(i).message = "OK";
        Message(i) = "OK";

    catch ME

        Message(i) = "Error: " + string(ME.message);
        VBMessage(i) = "Processing error: " + string(ME.message);

        DepthRNGall(i).matchedFile = matchedFile;
        DepthRNGall(i).fileFound = true;
        DepthRNGall(i).VBMessage = VBMessage(i);
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
DepthRNGtable.Depth1200Mean = Depth1200Mean;

DepthRNGtable.Depth1200std_B1 = Depth1200std(:,1);
DepthRNGtable.Depth1200std_B2 = Depth1200std(:,2);
DepthRNGtable.Depth1200std_B3 = Depth1200std(:,3);
DepthRNGtable.Depth1200std_B4 = Depth1200std(:,4);

DepthRNGtable.stdSpatial_1200 = stdSpatial_1200;
DepthRNGtable.VBAvailable = VBAvailable;
DepthRNGtable.VBDepth = VBDepth;
DepthRNGtable.VBDepthstd = VBDepthstd;
DepthRNGtable.VBSource = VBSource;
DepthRNGtable.VBInstrument = VBInstrument;
DepthRNGtable.VBOriginalUnit = VBOriginalUnit;
DepthRNGtable.VBScaleToMeters = VBScaleToMeters;
DepthRNGtable.VBDifferenceFromFourBeam = VBDifferenceFromFourBeam;
DepthRNGtable.VBRelativeDifference = VBRelativeDifference;
DepthRNGtable.VBMessage = VBMessage;
DepthRNGtable.Message = Message;

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
DepthRNGout.Depth1200Mean = Depth1200Mean;
DepthRNGout.stdSpatial_1200 = stdSpatial_1200;
DepthRNGout.stdSpatialALL_1200 = stdSpatialALL_1200;

DepthRNGout.VBDepth_ = VBDepth_;
DepthRNGout.VBDepth = VBDepth;
DepthRNGout.VBDepthstd = VBDepthstd;
DepthRNGout.VBAvailable = VBAvailable;
DepthRNGout.VBSource = VBSource;
DepthRNGout.VBInstrument = VBInstrument;
DepthRNGout.VBOriginalUnit = VBOriginalUnit;
DepthRNGout.VBScaleToMeters = VBScaleToMeters;
DepthRNGout.VBDifferenceFromFourBeam = VBDifferenceFromFourBeam;
DepthRNGout.VBRelativeDifference = VBRelativeDifference;
DepthRNGout.VBMessage = VBMessage;
DepthRNGout.VBMaxDepth = opt.VBMaxDepth;

DepthRNGout.DepthRNGall = DepthRNGall;
DepthRNGout.DepthRNGtable = DepthRNGtable;

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
    fprintf('Vertical beam available: %d/%d\n', sum(VBAvailable), nFiles);
    fprintf('Vertical-beam maximum accepted depth: %.3g m\n', opt.VBMaxDepth);
    fprintf('StartADCP/EndADCP were used as trim counts.\n');
    fprintf('NaN values were omitted from RNG/depth statistics.\n');
    fprintf('Time series are stored only in per-measurement cell arrays.\n');
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
%SAFEMEAN Mean after removing NaN values explicitly.

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


function y = safeMeanWithFallback(primary, fallback)
%SAFEMEANWITHFALLBACK Mean of primary data omitting NaNs.
% If primary contains no valid sample, use fallback data, also omitting NaNs.

y = safeMean(primary);

if isnan(y)
    y = safeMean(fallback);
end

end


function y = safeStd(x, flag)
%SAFESTD Standard deviation after removing NaN values explicitly.

if nargin < 2 || isempty(flag)
    flag = 0;
end

if isempty(x)
    y = NaN;
    return
end

x = double(x(:));
x = x(~isnan(x));

if isempty(x)
    y = NaN;
    return
end

if flag == 0 && numel(x) < 2
    y = NaN;
else
    y = std(x, flag);
end

end


function y = safeStdWithFallback(primary, fallback, flag)
%SAFESTDWITHFALLBACK Standard deviation omitting NaNs, with fallback.
% The fallback is used when the primary series has no valid standard
% deviation, for example when it contains fewer than two valid samples.

if nargin < 3 || isempty(flag)
    flag = 0;
end

y = safeStd(primary, flag);

if isnan(y)
    y = safeStd(fallback, flag);
end

end


function x = seriesWithFallback(primary, fallback)
%SERIESWITHFALLBACK Return primary unless it is entirely invalid.
% NaN locations inside a partly valid despiked series are retained so that
% nanstdLocal can omit them sample by sample across beams.

primary = double(primary(:));
fallback = double(fallback(:));

if any(~isnan(primary))
    x = primary;
else
    x = fallback;
end

end


function [v, source, instrument] = findVerticalBeamSeries(S)
%FINDVERTICALBEAMSERIES Find a vertical-beam depth/range vector.
%
% Priority is given to standardized depth fields. Legacy M9 and RDI aliases
% are retained as fallbacks. The returned vector is still in source units;
% conversion to metres is handled separately.

v = [];
source = "";
instrument = detectInstrument(S);

candidates = { ...
    {'T','VB','StandardDepth'}, "T.VB.StandardDepth"; ...
    {'T','VB','Depth'},         "T.VB.Depth"; ...
    {'T','VB','VBdepth'},       "T.VB.VBdepth"; ...
    {'T','VB','VBDepth'},       "T.VB.VBDepth"; ...
    {'T','VB','Range'},         "T.VB.Range"; ...
    {'T','VB','RNG'},           "T.VB.RNG"; ...
    {'VB','StandardDepth'}, "VB.StandardDepth"; ...
    {'VB','Depth'},         "VB.Depth"; ...
    {'T','BT','VBdepth'},   "T.BT.VBdepth"; ...
    {'T','BT','VBDepth'},   "T.BT.VBDepth"; ...
    {'T','BT','VBRange'},   "T.BT.VBRange"; ...
    {'T','BT','RNG5'},      "T.BT.RNG5"; ...
    {'BT','VBdepth'},       "BT.VBdepth"; ...
    {'BT','VBDepth'},       "BT.VBDepth"; ...
    {'VB','VBdepth'},       "VB.VBdepth"; ...
    {'VB','VBDepth'},       "VB.VBDepth"; ...
    {'VB','Range'},         "VB.Range"; ...
    {'VB','RNG'},           "VB.RNG"; ...
    {'BT','VBRange'},       "BT.VBRange"; ...
    {'BT','RNG5'},          "BT.RNG5"; ...
    {'T','BS_adp','VB_RNGM'}, "T.BS_adp.VB_RNGM"; ...
    {'T','BS_adp','VBdepth'}, "T.BS_adp.VBdepth"; ...
    {'BS_adp','VB_RNGM'},   "BS_adp.VB_RNGM"; ...
    {'BS_adp','VBdepth'},   "BS_adp.VBdepth"; ...
    {'DepthSounder_Range'}, "DepthSounder_Range"; ...
    {'VB_RNGM'},            "VB_RNGM"; ...
    {'VBdepth'},            "VBdepth"};

for k = 1:size(candidates,1)
    [found, candidate] = getNestedNumericField(S, candidates{k,1});

    if found && any(isfinite(candidate(:)))
        v = double(candidate(:));
        source = candidates{k,2};

        if instrument == "Unknown"
            if startsWith(source, "T.", 'IgnoreCase', true) || ...
                    any(contains(source, ["BS_adp", "DepthSounder"], ...
                    'IgnoreCase', true))
                instrument = "SonTek M9";
            elseif startsWith(source, "BT.", 'IgnoreCase', true)
                instrument = "RDI";
            end
        end
        return
    end
end

end


function instrument = detectInstrument(S)

instrument = "Unknown";

if isfield(S, 'M9meta') || ...
        (isfield(S, 'T') && isstruct(S.T) && isfield(S.T, 'M9meta'))
    instrument = "SonTek M9";
    return
end

if isfield(S, 'T') && isstruct(S.T) && ...
        isfield(S.T, 'VB') && isstruct(S.T.VB) && ...
        isfield(S.T.VB, 'SourceInstrument')
    try
        sourceInstrument = string(S.T.VB.SourceInstrument);
        sourceInstrument = sourceInstrument(1);

        if strlength(strtrim(sourceInstrument)) > 0
            instrument = sourceInstrument;
            return
        end
    catch
    end
end

if isfield(S, 'VB') && isstruct(S.VB) && ...
        isfield(S.VB, 'SourceInstrument')
    try
        sourceInstrument = string(S.VB.SourceInstrument);
        sourceInstrument = sourceInstrument(1);

        if strlength(strtrim(sourceInstrument)) > 0
            instrument = sourceInstrument;
            return
        end
    catch
    end
end

if isfield(S, 'BS_adp') && isstruct(S.BS_adp) && ...
        isfield(S.BS_adp, 'VB_RNGM')
    instrument = "SonTek M9";
    return
end

if isfield(S, 'T') && isstruct(S.T) && ...
        isfield(S.T, 'BS_adp') && isstruct(S.T.BS_adp) && ...
        isfield(S.T.BS_adp, 'VB_RNGM')
    instrument = "SonTek M9";
    return
end

if isfield(S, 'Data') && isstruct(S.Data)
    metadataNames = {'InstrumentName','InstrumentType','SystemType'};

    for k = 1:numel(metadataNames)
        if isfield(S.Data, metadataNames{k})
            try
                metadataValue = string(S.Data.(metadataNames{k}));
                metadataValue = metadataValue(1);

                if any(contains(metadataValue, ["M9", "SonTek"], ...
                        'IgnoreCase', true))
                    instrument = "SonTek M9";
                elseif any(contains(metadataValue, ...
                        ["RDI", "RiverPro", "RioGrande", "WorkHorse"], ...
                        'IgnoreCase', true))
                    instrument = "RDI";
                end
                return
            catch
            end
        end
    end
end

end


function [found, value] = getNestedNumericField(S, pathParts)

found = false;
value = [];
current = S;

for k = 1:numel(pathParts)
    fieldName = pathParts{k};

    if ~isstruct(current) || ~isfield(current, fieldName)
        return
    end

    current = current.(fieldName);
end

try
    if isnumeric(current) || islogical(current)
        value = double(current);
    elseif iscell(current) || isstring(current) || ischar(current)
        value = str2double(string(current));
    else
        return
    end
catch
    return
end

value = double(value(:));
found = ~isempty(value);

end


function [metres, originalUnit, scaleToMetres, message] = ...
        convertVerticalBeamToMetres( ...
        rawValues, referenceDepth, maxDepth, source)
%CONVERTVERTICALBEAMTOMETRES Resolve m/dm/cm using the four-beam depth.

rawValues = double(rawValues(:));
rawValues(rawValues == 0) = NaN;

metres = nan(size(rawValues));
originalUnit = "";
scaleToMetres = NaN;
message = "Vertical beam contains no valid numeric data";

finiteRaw = isfinite(rawValues) & rawValues > 0;
if ~any(finiteRaw)
    return
end

candidateScales = [1, 0.1, 0.01];
candidateUnits = ["m", "dm", "cm"];
candidateScores = inf(size(candidateScales));
candidateUsable = false(size(candidateScales));

referenceAvailable = ...
    isfinite(referenceDepth) && referenceDepth > 0 && ...
    referenceDepth <= maxDepth;

for k = 1:numel(candidateScales)
    candidate = rawValues .* candidateScales(k);
    valid = isfinite(candidate) & candidate > 0 & candidate <= maxDepth;

    if ~any(valid)
        continue
    end

    candidateUsable(k) = true;
    candidateMedian = safeMedian(candidate(valid));
    rejectedFraction = sum(finiteRaw & ~valid) ./ sum(finiteRaw);

    if referenceAvailable
        candidateScores(k) = ...
            abs(candidateMedian - referenceDepth) ./ ...
            max(abs(referenceDepth), eps) + rejectedFraction;
    else
        candidateScores(k) = rejectedFraction;
    end
end

if ~any(candidateUsable)
    message = "No m/dm/cm conversion produced depth within 0-" + ...
        string(maxDepth) + " m";
    return
end

if referenceAvailable
    [~, bestIndex] = min(candidateScores);
else
    sourceLikelyMetres = ...
        any(contains(source, ...
        ["StandardDepth", ".Depth", "VBdepth", "VB_RNGM", "RNG5"], ...
        'IgnoreCase', true));

    if sourceLikelyMetres && candidateUsable(1)
        bestIndex = 1;
    else
        usableIndices = find(candidateUsable);

        if numel(usableIndices) ~= 1
            message = "Vertical-beam unit is ambiguous because the " + ...
                "four-beam reference depth is unavailable";
            return
        end

        bestIndex = usableIndices(1);
    end
end

scaleToMetres = candidateScales(bestIndex);
originalUnit = candidateUnits(bestIndex);
metres = rawValues .* scaleToMetres;

% The user-defined river-depth limit is applied to every sample, not only
% to the mean, so impossible values cannot influence despiking/statistics.
invalid = ~isfinite(metres) | metres <= 0 | metres > maxDepth;
metres(invalid) = NaN;

if any(isfinite(metres))
    message = "Vertical beam converted from " + originalUnit + ...
        " to m (scale " + string(scaleToMetres) + ")";
else
    message = "Vertical beam had no values within 0-" + ...
        string(maxDepth) + " m after unit conversion";
end

end


function y = safeMedian(x)

if isempty(x)
    y = NaN;
    return
end

x = double(x(:));
x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = median(x);
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


function plotDepthRNGdiagnostics( ...
    i, ...
    thisADCPfilename, ...
    start, ...
    finish, ...
    RNG11, RNG22, RNG33, RNG44, ...
    RNG11_1200, RNG22_1200, RNG33_1200, RNG44_1200, ...
    VBDepth_i, ...
    stdSpatialALL_1200_i, ...
    PlotRaw, ...
    PlotDespiked, ...
    PlotSpatialStd)

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

    if ~isempty(VBDepth_i) && any(isfinite(VBDepth_i))
        xVB = (start:(start + numel(VBDepth_i) - 1))';
        plot(xVB, VBDepth_i, 'k--', 'LineWidth', 1.25, ...
            'DisplayName', 'Vertical beam [m]');
    end

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
