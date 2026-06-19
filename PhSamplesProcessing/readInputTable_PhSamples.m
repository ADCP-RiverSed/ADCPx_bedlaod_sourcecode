function ADCPin = readInputTable_PhSamples(thisFolder, varargin)
%READINPUTTABLE_PHSAMPLES Robust function to read ADCP input table.
%
% This function reads and extracts all input-table variables.
% Acoustic/physical constants are NOT calculated here.
%
% Usage:
%   ADCPin = readInputTable_PhSamples(thisFolder);
%
% Optional:
%   ADCPin = readInputTable_PhSamples(thisFolder, 'ENDdata', 10);
%   ADCPin = readInputTable_PhSamples(thisFolder, 'AssignToCaller', true);
%   ADCPin = readInputTable_PhSamples(thisFolder, 'AssignToBase', true);
%   ADCPin = readInputTable_PhSamples(thisFolder, 'ExcelFile', fullfile(thisFolder,'MyInputTable.xlsx'));
%
% Main behaviour:
%   - automatically finds Excel file with "InputTable" in name
%   - reads Excel range A:AA
%   - keeps columns only until EndADCP
%   - removes rows where "Name as it is in the folder" is empty
%   - substitutes missing mobile D-values with grab D-values
%   - calculates DnumMean14 and DmeaNUMpsd from grainsize_lognormal(D50,D90)
%   - keeps old variable names for compatibility
%   - optionally assigns variables to caller/base workspace

%% Parse inputs

    p = inputParser;

    addRequired(p, 'thisFolder', @(x) ischar(x) || isstring(x));

    addParameter(p, 'ExcelFile', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Pattern', '*InputTable*.xlsx', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Sheet', 1);
    addParameter(p, 'Range', 'A:AA', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ENDdata', [], @(x) isempty(x) || isnumeric(x));

    addParameter(p, 'AssignToCaller', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'AssignToBase', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

    parse(p, thisFolder, varargin{:});
    opt = p.Results;

    thisFolder = char(thisFolder);

    if ~exist(thisFolder, 'dir')
        error('Folder does not exist:\n%s', thisFolder);
    end

%% Find InputTable Excel file

    if ~isempty(opt.ExcelFile)

        excelFile = char(opt.ExcelFile);

        if ~exist(excelFile, 'file')
            error('Specified Excel file does not exist:\n%s', excelFile);
        end

        [~, excelName, excelExt] = fileparts(excelFile);
        excelDisplayName = [excelName excelExt];

    else

        DIRps = dir(fullfile(thisFolder, char(opt.Pattern)));

        % Ignore temporary Excel lock files
        isTemp = startsWith({DIRps.name}, '~$');
        DIRps = DIRps(~isTemp);

        if isempty(DIRps)
            error('No Excel file matching "%s" was found in:\n%s', ...
                char(opt.Pattern), thisFolder);
        end

        % If several InputTable files exist, use newest one
        if numel(DIRps) > 1

            [~, newestIdx] = max([DIRps.datenum]);

            if opt.Verbose
                fprintf('\nSeveral InputTable Excel files found. Using newest one:\n');
                for i = 1:numel(DIRps)
                    fprintf('  %d: %s\n', i, DIRps(i).name);
                end
            end

            DIRps = DIRps(newestIdx);
        end

        excelFile = fullfile(thisFolder, DIRps(1).name);
        excelDisplayName = DIRps(1).name;
    end

    if opt.Verbose
        fprintf('\nReading ADCP input table:\n%s\n', excelFile);
    end

%% Read Excel file robustly

    opts = detectImportOptions( ...
        excelFile, ...
        'Sheet', opt.Sheet, ...
        'Range', char(opt.Range), ...
        'VariableNamingRule', 'preserve');

    try
        opts = setvaropts(opts, opts.VariableNames, 'WhitespaceRule', 'preserve');
    catch
    end

    DataP = readtable(excelFile, opts);

    % Remove fully empty rows/columns
    DataP = removeEmptyRowsAndColumns(DataP);

    if isempty(DataP) || height(DataP) == 0
        error('Input table is empty after removing empty rows/columns.');
    end

%% Keep only columns until EndADCP

    endIdx = findColumnIndex(DataP, {'EndADCP', 'End ADCP', 'End_ADCP', 'EndADCP '});

    if ~isempty(endIdx)
        DataP = DataP(:, 1:endIdx);
    end

%% Find leading ADCP filename column

    leadColNames = { ...
        'Name as it is in the folder', ...
        'NameAsItIsInTheFolder', ...
        'Number as it is in the folder', ...
        'NumberAsItIsInTheFolder', ...
        'ADCPfilename', ...
        'ADCP file name', ...
        'ADCP filename'};

    leadIdx = findColumnIndex(DataP, leadColNames);

    if isempty(leadIdx)
        fprintf('\nAvailable Excel headers are:\n');
        disp(DataP.Properties.VariableNames');
        error(['Could not find the leading ADCP filename column.\n', ...
               'Expected column name similar to: "Name as it is in the folder".']);
    end

    leadRaw = DataP{:, leadIdx};
    leadStr = string(leadRaw);
    leadStr = strtrim(leadStr);

    validRows = ~(ismissing(leadStr) | leadStr == "" | lower(leadStr) == "nan");

    DataP = DataP(validRows, :);
    leadStr = leadStr(validRows);

    if height(DataP) == 0
        error('No non-empty ADCP rows found using leading column: %s', ...
            DataP.Properties.VariableNames{leadIdx});
    end

%% Limit to ENDdata if requested

    if ~isempty(opt.ENDdata)
        ENDdata = min(opt.ENDdata, height(DataP));
        DataP = DataP(1:ENDdata, :);
    else
        ENDdata = height(DataP);
    end

%% Clean column names and store original columns too

    cleanNames = matlab.lang.makeValidName(DataP.Properties.VariableNames);
    cleanNames = matlab.lang.makeUniqueStrings(cleanNames);

    RawColumns = struct();

    for i = 1:width(DataP)
        RawColumns.(cleanNames{i}) = DataP{:, i};
    end

%% Map table columns to original variable names

    Date_ADCP = getCol(DataP, { ...
        'Date', ...
        'Date ADCP', ...
        'Date_ADCP', ...
        'Date-collection', ...
        'Date collection'}, false, strings(ENDdata,1));

ADCPfilename = getCol(DataP, leadColNames, true, strings(ENDdata,1));
ADCPfilename = string(ADCPfilename);
ADCPfilename = strtrim(ADCPfilename);

% Start/end ADCP indices or flags from input table
StartADCP = getNumCol(DataP, { ...
    'StartADCP', ...
    'Start ADCP', ...
    'Start_ADCP', ...
    'Start'}, false, nan(ENDdata,1));

EndADCP = getNumCol(DataP, { ...
    'EndADCP', ...
    'End ADCP', ...
    'End_ADCP', ...
    'End'}, false, nan(ENDdata,1));


    numName = (1:ENDdata)';
    Position1ADCPs = numName;
    Position2ADCPs = numName;

    DistanceFromLeftBank = getNumCol(DataP, { ...
        'Distance from the reference point', ...
        'DistanceFromTheReferencePoint', ...
        'Distance from left bank', ...
        'DistanceFromLeftBank'}, true, nan(ENDdata,1));

    DistanceToBank = DistanceFromLeftBank;

%% Transport rates

    TrateMax = getNumCol(DataP, { ...
        'TrateMax', ...
        'TrateMax ', ...
        'Trate Max', ...
        'T rate max'}, false, nan(ENDdata,1));

    TrateMin = getNumCol(DataP, { ...
        'TrateMin', ...
        'Trate Min', ...
        'T rate min'}, false, nan(ENDdata,1));

    TrateMean = getNumCol(DataP, { ...
        'TrateMean', ...
        'TrateMean ', ...
        'Trate Mean', ...
        'T rate mean'}, false, nan(ENDdata,1));

    MTrateAll = TrateMean;
    MinTrate = TrateMin;
    MaxTrate = TrateMax;
    MeanTrate = TrateMean;

%% Sand transport

    Tsand = getNumCol(DataP, { ...
        'Tsand', ...
        'Tsand ', ...
        'T sand'}, false, nan(ENDdata,1));

    Tsand_1 = getNumCol(DataP, { ...
        'Tsand.1', ...
        'Tsand_1', ...
        'Tsand 1'}, false, Tsand);

    DistanceFromRefPointSand = getNumCol(DataP, { ...
        'Distance from the reference point for sand', ...
        'DistanceFromTheReferencePointForSand'}, false, nan(ENDdata,1));

    TrateSAND = getNumCol(DataP, { ...
        'T sand (g/s.m)', ...
        'TSand_g_s_m_', ...
        'T sand g/s.m', ...
        'TsandGsm', ...
        'Tsand'}, false, nan(ENDdata,1));

    TloadSSCsand = getNumCol(DataP, { ...
        'T suspended (g/s.m)', ...
        'Tsuspended (g/s.m)', ...
        'TSuspended_g_s_m_', ...
        'T suspended g/s.m', ...
        'TsuspendedGsm'}, false, nan(ENDdata,1));

    try
        FsS = TrateMean ./ TrateSAND;
    catch
        FsS = ones(ENDdata,1);
    end

%% Grain size distributions - mobile

    D10mobile = getNumCol(DataP, {'D10mobile', 'D10 mobile'}, false, nan(ENDdata,1));
    D16mobile = getNumCol(DataP, {'D16mobile', 'D16 mobile'}, false, nan(ENDdata,1));
    D50mobile = getNumCol(DataP, {'D50mobile', 'D50 mobile'}, false, nan(ENDdata,1));
    D84mobile = getNumCol(DataP, {'D84mobile', 'D84 mobile'}, false, nan(ENDdata,1));
    D90mobile = getNumCol(DataP, {'D90mobile', 'D90mobile ', 'D90 mobile'}, false, nan(ENDdata,1));

%% Grain size distributions - grab samples

    D10grab = getNumCol(DataP, {'D10grab', 'D10 grab'}, false, nan(ENDdata,1));
    D16grab = getNumCol(DataP, {'D16grab', 'D16 grab'}, false, nan(ENDdata,1));
    D50grab = getNumCol(DataP, {'D50grab', 'D50 grab'}, false, nan(ENDdata,1));
    D84grab = getNumCol(DataP, {'D84grab', 'D84 grab'}, false, nan(ENDdata,1));
    D90grab = getNumCol(DataP, {'D90grab', 'D90grab ', 'D90 grab'}, false, nan(ENDdata,1));

%% Reassign grain sizes for compatibility with old scripts
% If mobile D-values are missing, use grab D-values.
% If both are missing, use old fallback estimates where possible.

    % Grab variables kept with old names
    D_10grabM = D10grab;
    D_16grabM = D16grab;
    D_50grabM = D50grab;
    D_84grabM = D84grab;
    D_90grabM = D90grab;

    % Mobile variables with grab fallback
    D10 = fillMissingDwithGrab(D10mobile, D10grab);
    D16 = fillMissingDwithGrab(D16mobile, D16grab);
    D50 = fillMissingDwithGrab(D50mobile, D50grab);
    D84 = fillMissingDwithGrab(D84mobile, D84grab);
    D90 = fillMissingDwithGrab(D90mobile, D90grab);

    % Fallbacks from D50 where still missing
    idx = isnan(D10) & ~isnan(D50);
    D10(idx) = 0.1 .* D50(idx);

    idx = isnan(D16) & ~isnan(D50);
    D16(idx) = 0.15 .* D50(idx);

    idx = isnan(D84) & ~isnan(D50);
    D84(idx) = 2.4 .* D50(idx);

    idx = isnan(D90) & ~isnan(D84);
    D90(idx) = D84(idx);

    idx = isnan(D90) & ~isnan(D50);
    D90(idx) = 2.4 .* D50(idx);

    % Lognormal function requires D90 > D50.
    % If D90 <= D50, adjust slightly to avoid crash.
    idxBad = ~isnan(D50) & ~isnan(D90) & D90 <= D50;
    if any(idxBad)
        warning('Some rows had D90 <= D50. Setting D90 = 1.01 * D50 for those rows.');
        D90(idxBad) = 1.01 .* D50(idxBad);
    end

    % D65 estimate, old name preserved
    D65 = D50 .* 1.1;

    % Temporary assignments for compatibility
    D50_05 = D50;
    D50_14 = D50;

    % Lognormal grain-size statistics from D50 and D90
    gs = grainsize_lognormal(D50, D90);

    % Old variables preserved
    DnumMean14 = gs.D_mean_n;
    DmeaNUMpsd = gs.D_mean_n;

    % Additional lognormal variables also kept as raw names
    mu_ln      = gs.mu_ln;
    sigma_ln   = gs.sigma_ln;
    D_mean_n   = gs.D_mean_n;
    D_mean_v   = gs.D_mean_v;
    Dg         = gs.Dg;
    sigma_g    = gs.sigma_g;
    Sauter     = gs.Sauter;
    skewness_n = gs.skewness_n;
    kurt_n     = gs.kurt_n;
    phi_mean   = gs.phi_mean;
    phi_sigma  = gs.phi_sigma;

%% Standard deviation

    Dstd = getNumCol(DataP, {'Dstd', 'Dstd ', 'D std', 'D_std'}, false, nan(ENDdata,1));

    if all(isnan(Dstd))
        Dstd = (D84 - D16) ./ 2;
    end

    Dstd_14 = Dstd;
    Dstd_05 = Dstd;

%% Sediment information

    Dss = getNumCol(DataP, {'Dss', 'D ss', 'D_SS'}, false, nan(ENDdata,1));

    D_SS = Dss;

    try
        D_SS = fillmissing(D_SS, 'linear');
    catch
    end

    SSC = getNumCol(DataP, {'SSC', 'SSC ', 'Suspended sediment concentration'}, false, nan(ENDdata,1));

    try
        SSC = fillmissing(SSC, 'linear');
    catch
    end

    SSCsand = getNumCol(DataP, { ...
        'SSCsand mg/l', ...
        'SSCsandMg_l', ...
        'SSC sand mg/l', ...
        'SSCsand'}, false, zeros(ENDdata,1));

%% Input check table

InputCheckTable = table( ...
    (1:ENDdata)', ...
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

%% Build output structure with old variable names

    ADCPin = struct();

    ADCPin.thisFolder = thisFolder;
    ADCPin.excelFile = excelFile;
    ADCPin.excelDisplayName = excelDisplayName;
    ADCPin.DataP = DataP;
    ADCPin.RawColumns = RawColumns;
    ADCPin.InputCheckTable = InputCheckTable;
    ADCPin.ENDdata = ENDdata;
    ADCPin.Date_ADCP = Date_ADCP;
    ADCPin.ADCPfilename = ADCPfilename;

    ADCPin.StartADCP = StartADCP;
    ADCPin.EndADCP = EndADCP;

    ADCPin.numName = numName;
    ADCPin.Position1ADCPs = Position1ADCPs;
    ADCPin.Position2ADCPs = Position2ADCPs;

    ADCPin.DistanceFromLeftBank = DistanceFromLeftBank;
    ADCPin.DistanceToBank = DistanceToBank;

    ADCPin.TrateMax = TrateMax;
    ADCPin.TrateMin = TrateMin;
    ADCPin.TrateMean = TrateMean;
    ADCPin.MTrateAll = MTrateAll;
    ADCPin.MinTrate = MinTrate;
    ADCPin.MaxTrate = MaxTrate;
    ADCPin.MeanTrate = MeanTrate;

    ADCPin.Tsand = Tsand;
    ADCPin.Tsand_1 = Tsand_1;
    ADCPin.DistanceFromRefPointSand = DistanceFromRefPointSand;
    ADCPin.TrateSAND = TrateSAND;
    ADCPin.TloadSSCsand = TloadSSCsand;
    ADCPin.FsS = FsS;

    ADCPin.D10mobile = D10mobile;
    ADCPin.D16mobile = D16mobile;
    ADCPin.D50mobile = D50mobile;
    ADCPin.D84mobile = D84mobile;
    ADCPin.D90mobile = D90mobile;

    ADCPin.D10grab = D10grab;
    ADCPin.D16grab = D16grab;
    ADCPin.D50grab = D50grab;
    ADCPin.D84grab = D84grab;
    ADCPin.D90grab = D90grab;

    ADCPin.D50 = D50;
    ADCPin.D_50grabM = D_50grabM;
    ADCPin.D90 = D90;
    ADCPin.D_90grabM = D_90grabM;
    ADCPin.D10 = D10;
    ADCPin.D_10grabM = D_10grabM;
    ADCPin.D16 = D16;
    ADCPin.D_16grabM = D_16grabM;
    ADCPin.D84 = D84;
    ADCPin.D_84grabM = D_84grabM;
    ADCPin.D65 = D65;

    ADCPin.D50_05 = D50_05;
    ADCPin.D50_14 = D50_14;
    ADCPin.DnumMean14 = DnumMean14;
    ADCPin.DmeaNUMpsd = DmeaNUMpsd;

    ADCPin.gs = gs;
    ADCPin.mu_ln = mu_ln;
    ADCPin.sigma_ln = sigma_ln;
    ADCPin.D_mean_n = D_mean_n;
    ADCPin.D_mean_v = D_mean_v;
    ADCPin.Dg = Dg;
    ADCPin.sigma_g = sigma_g;
    ADCPin.Sauter = Sauter;
    ADCPin.skewness_n = skewness_n;
    ADCPin.kurt_n = kurt_n;
    ADCPin.phi_mean = phi_mean;
    ADCPin.phi_sigma = phi_sigma;

    ADCPin.Dstd = Dstd;
    ADCPin.Dstd_14 = Dstd_14;
    ADCPin.Dstd_05 = Dstd_05;

    ADCPin.Dss = Dss;
    ADCPin.D_SS = D_SS;
    ADCPin.SSC = SSC;
    ADCPin.SSCsand = SSCsand;

%% Optional assignment to caller/base workspace

    if opt.AssignToCaller
        assignStructToWorkspace(ADCPin, 'caller');
    end

    if opt.AssignToBase
        assignStructToWorkspace(ADCPin, 'base');
    end

%% Display summary

    if opt.Verbose
        fprintf('\n=== Input Table Summary ===\n');
        fprintf('File: %s\n', excelDisplayName);
        fprintf('Rows read: %d\n', ENDdata);
        fprintf('Leading ADCP file column: %s\n', DataP.Properties.VariableNames{leadIdx});
        fprintf('First ADCP file/name: %s\n', string(ADCPfilename(1)));

        try
            fprintf('Date first row: %s\n', string(Date_ADCP(1)));
        catch
        end

        fprintf('Distance range: %.1f - %.1f m from reference/left bank\n', ...
            min(DistanceFromLeftBank, [], 'omitnan'), ...
            max(DistanceFromLeftBank, [], 'omitnan'));

        fprintf('D50 used range: %.2f - %.2f\n', ...
            min(D50, [], 'omitnan'), ...
            max(D50, [], 'omitnan'));

        fprintf('D90 used range: %.2f - %.2f\n', ...
            min(D90, [], 'omitnan'), ...
            max(D90, [], 'omitnan'));

        fprintf('DmeaNUMpsd range: %.2f - %.2f\n', ...
            min(DmeaNUMpsd, [], 'omitnan'), ...
            max(DmeaNUMpsd, [], 'omitnan'));

        fprintf('Mean transport rate range: %.3f - %.3f g/m.s\n', ...
            min(TrateMean, [], 'omitnan'), ...
            max(TrateMean, [], 'omitnan'));

        fprintf('===========================\n');

        disp(InputCheckTable);
    end
end

%% Helper functions

function T = removeEmptyRowsAndColumns(T)

    if isempty(T)
        return
    end

    nRows = height(T);
    nCols = width(T);

    emptyRow = true(nRows,1);

    for i = 1:nRows
        for j = 1:nCols
            val = T{i,j};

            if ~isEmptyValue(val)
                emptyRow(i) = false;
                break
            end
        end
    end

    T(emptyRow,:) = [];

    if isempty(T)
        return
    end

    nRows = height(T);
    nCols = width(T);

    emptyCol = true(1,nCols);

    for j = 1:nCols
        for i = 1:nRows
            val = T{i,j};

            if ~isEmptyValue(val)
                emptyCol(j) = false;
                break
            end
        end
    end

    T(:,emptyCol) = [];
end


function tf = isEmptyValue(val)

    if isempty(val)
        tf = true;
        return
    end

    if iscell(val)
        if isempty(val{1})
            tf = true;
            return
        end

        val = val{1};
    end

    try
        if ismissing(val)
            tf = true;
            return
        end
    catch
    end

    if ischar(val) || isstring(val)
        s = strtrim(string(val));
        tf = ismissing(s) || s == "" || lower(s) == "nan";
        return
    end

    if isnumeric(val)
        tf = all(isnan(val));
        return
    end

    tf = false;
end


function idx = findColumnIndex(T, possibleNames)

    headers = T.Properties.VariableNames;
    idx = [];

    normHeaders = cellfun(@normalizeName, headers, 'UniformOutput', false);

    for k = 1:numel(possibleNames)
        target = normalizeName(possibleNames{k});
        hit = find(strcmp(normHeaders, target), 1, 'first');

        if ~isempty(hit)
            idx = hit;
            return
        end
    end
end


function x = getCol(T, possibleNames, required, defaultValue)

    idx = findColumnIndex(T, possibleNames);

    if isempty(idx)
        if required
            fprintf('\nAvailable Excel headers are:\n');
            disp(T.Properties.VariableNames');
            error('Required column not found. Expected one of: %s', ...
                strjoin(possibleNames, ', '));
        else
            x = defaultValue;
            return
        end
    end

    x = T{:, idx};
end


function x = getNumCol(T, possibleNames, required, defaultValue)

    x = getCol(T, possibleNames, required, defaultValue);

    if isempty(x)
        x = defaultValue;
        return
    end

    if isnumeric(x)
        x = double(x);
        return
    end

    if islogical(x)
        x = double(x);
        return
    end

    if iscell(x)
        x = string(x);
    end

    if ischar(x)
        x = string(x);
    end

    if isstring(x)
        x = strtrim(x);
        x = replace(x, ",", ".");
        x = regexprep(x, '[^\d\.\-\+eE]', '');
        x = str2double(x);
        return
    end

    try
        x = double(x);
    catch
        warning('Could not convert column to numeric. Using default value.');
        x = defaultValue;
    end
end


function s = normalizeName(s)

    s = string(s);
    s = strtrim(s);
    s = lower(s);
    s = regexprep(s, '[^a-z0-9]', '');
    s = char(s);
end


function assignStructToWorkspace(S, workspaceName)

    f = fieldnames(S);

    skipFields = {'RawColumns'};

    for i = 1:numel(f)

        if any(strcmp(f{i}, skipFields))
            continue
        end

        assignin(workspaceName, f{i}, S.(f{i}));
    end
end


function Dout = fillMissingDwithGrab(Dmobile, Dgrab)
%FILLMISSINGDWITHGRAB Fill missing mobile grain-size values using grab values.
%
% Dmobile is used where available.
% Dgrab is used where Dmobile is NaN.

    Dmobile = double(Dmobile);
    Dgrab = double(Dgrab);

    Dout = Dmobile;

    idxMissing = isnan(Dout);

    if numel(Dgrab) == numel(Dout)
        Dout(idxMissing) = Dgrab(idxMissing);
    elseif isscalar(Dgrab)
        Dout(idxMissing) = Dgrab;
    else
        warning('Dgrab size does not match Dmobile size. Missing mobile values not fully replaced.');
    end
end


function gs = grainsize_lognormal(D50, D90)
%GRAINSIZE_LOGNORMAL Lognormal grain-size statistics from D50 and D90.
%
% Inputs:
%   D50 - median grain diameter
%   D90 - 90th percentile grain diameter
%
% Output structure gs:
%   mu_ln
%   sigma_ln
%   D_mean_n
%   D_mean_v
%   Dg
%   sigma_g
%   Sauter
%   skewness_n
%   kurt_n
%   phi_mean
%   phi_sigma
%
% D50 and D90 must have the same unit, normally mm.

    D50 = double(D50);
    D90 = double(D90);

    if ~isequal(size(D50), size(D90))
        error('grainsize_lognormal:badInput', ...
            'D50 and D90 must have the same size.');
    end

    mu_ln      = nan(size(D50));
    sigma_ln   = nan(size(D50));
    D_mean_n   = nan(size(D50));
    D_mean_v   = nan(size(D50));
    Dg         = nan(size(D50));
    sigma_g    = nan(size(D50));
    Sauter     = nan(size(D50));
    skewness_n = nan(size(D50));
    kurt_n     = nan(size(D50));
    phi_mean   = nan(size(D50));
    phi_sigma  = nan(size(D50));

    valid = ~isnan(D50) & ~isnan(D90) & D50 > 0 & D90 > D50;

    if any(~valid)
        warning('grainsize_lognormal:missingOrInvalidRows', ...
            'Some D50/D90 rows are missing or invalid. Output is NaN for those rows.');
    end

    if any(valid)

        z90 = 1.281551566;

        mu_ln(valid) = log(D50(valid));
        sigma_ln(valid) = log(D90(valid) ./ D50(valid)) ./ z90;

        D_mean_n(valid) = exp(mu_ln(valid) + sigma_ln(valid).^2 ./ 2);
        sigma_g(valid) = exp(sigma_ln(valid));

        skewness_n(valid) = (exp(sigma_ln(valid).^2) + 2) .* ...
            sqrt(exp(sigma_ln(valid).^2) - 1);

        kurt_n(valid) = exp(4 .* sigma_ln(valid).^2) + ...
            2 .* exp(3 .* sigma_ln(valid).^2) + ...
            3 .* exp(2 .* sigma_ln(valid).^2) - 6;

        D_mean_v(valid) = exp(mu_ln(valid) + 3 .* sigma_ln(valid).^2 ./ 2);

        Sauter(valid) = exp(mu_ln(valid) + 5 .* sigma_ln(valid).^2 ./ 2);

        Dg(valid) = D50(valid);

        phi_mean(valid) = -log2(D_mean_n(valid));
        phi_sigma(valid) = sigma_ln(valid) ./ log(2);
    end

    gs.mu_ln = mu_ln;
    gs.sigma_ln = sigma_ln;
    gs.D_mean_n = D_mean_n;
    gs.D_mean_v = D_mean_v;
    gs.Dg = Dg;
    gs.sigma_g = sigma_g;
    gs.Sauter = Sauter;
    gs.skewness_n = skewness_n;
    gs.kurt_n = kurt_n;
    gs.phi_mean = phi_mean;
    gs.phi_sigma = phi_sigma;
end