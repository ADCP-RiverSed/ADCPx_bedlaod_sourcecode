function [BT_N,BT_E,FilteredDataN,FilteredDataE,alphaV,alphaBTf,alphaBT,Vmod,stats,TotalNaNDataN,TotalNaNDataE] = ...
    filterBT2WaterProfile(V_E,V_N,BT_E,BT_N,Depth,CellSize,varargin)
% FILTERBT2WATERPROFILE
%
% Direction-based filter for bottom-track (BT) velocities using ADCP water
% velocity profiles.
%
% -------------------------------------------------------------------------
% MAIN IDEA
% -------------------------------------------------------------------------
% For each time sample:
%   1) Smooth the ADCP water-profile velocities V_E and V_N in time
%   2) Limit valid water-profile cells using Depth and CellSize
%   3) Build one reference water vector from the profile using one of:
%        - 'full'     : mean of all valid cells above bottom
%        - 'last'     : last valid cell above bottom
%        - 'last2avg' : mean of the last two valid cells above bottom
%   4) Compute the angle difference between BT and water reference using
%      the dot product
%   5) Remove BT values if the angle difference is larger than the threshold
%
% By default the threshold is 90 degrees, which means BT is removed when
% it points into the opposite half-plane relative to the water vector.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% V_E, V_N
%   Water velocity profile matrices [nCells x N]
%   Rows    = vertical cells
%   Columns = time samples
%
% BT_E, BT_N
%   Bottom-track velocity components, vectors of length N
%
% Depth
%   Water depth for each time sample
%   Can be:
%     - scalar
%     - vector length N
%
% CellSize
%   Cell size for each time sample
%   Can be:
%     - scalar
%     - vector length N
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE INPUTS
% -------------------------------------------------------------------------
% 'AngleThreshold' : angle threshold in degrees, default = 90
%                    BT is removed when angle difference > AngleThreshold
%
% 'Fi'             : cosine threshold for the dot product
%                    If given, it overrides AngleThreshold
%                    Examples:
%                       Fi = 0    -> 90 degrees
%                       Fi = 0.5  -> 60 degrees
%                       Fi = -0.5 -> 120 degrees
%
% 'ProfileMethod'  : 'full', 'last', or 'last2avg'
%                    default = 'full'
%
% 'SmoothWindow'   : moving-average window along time, default = 5
%
% 'StartCell'      : first cell to consider, default = 1
%
% 'EndCell'        : last cell to consider before depth limitation
%                    default = []  -> all cells
%
% 'FirstCellCenter': first cell-center distance from transducer
%                    Can be:
%                      []         -> default = CellSize/2
%                      scalar
%                      vector length N
%
% 'UseAbsDepth'    : true/false, default = false
%                    If true, abs(Depth) is used
%
% -------------------------------------------------------------------------
% DEPTH LIMITATION
% -------------------------------------------------------------------------
% For each time sample, a cell is considered valid if its center lies above
% the bottom:
%
%   z_cell(k,t) = FirstCellCenter(t) + (k-1)*CellSize(t)
%   valid if z_cell(k,t) < Depth(t)
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% BT_N, BT_E
%   Filtered BT components
%
% FilteredDataN, FilteredDataE
%   Percentage of newly filtered data in BT_N and BT_E relative to
%   originally valid samples
%
% alphaV
%   Direction angle of reference water vector [deg, 0-360]
%
% alphaBTf
%   Direction angle of filtered BT vector [deg, 0-360]
%
% alphaBT
%   Direction angle of original BT vector [deg, 0-360]
%
% Vmod
%   Magnitude of reference water vector
%
% stats
%   Diagnostics structure
%
% TotalNaNDataN, TotalNaNDataE
%   Total percentage of NaNs in final BT_N and BT_E relative to all samples
%
% -------------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % Keep original orientation of BT
    % ---------------------------------------------------------------------
    wasRowE = isrow(BT_E);
    wasRowN = isrow(BT_N);

    BT_E_in = BT_E(:);
    BT_N_in = BT_N(:);

    N = numel(BT_E_in);

    if numel(BT_N_in) ~= N
        error('BT_E and BT_N must have the same length.');
    end

    % ---------------------------------------------------------------------
    % Check water-profile matrices
    % ---------------------------------------------------------------------
    if ~ismatrix(V_E) || ~ismatrix(V_N)
        error('V_E and V_N must be 2D matrices of size [nCells x N].');
    end

    [nCellsE,nTimeE] = size(V_E);
    [nCellsN,nTimeN] = size(V_N);

    if nCellsE ~= nCellsN || nTimeE ~= nTimeN
        error('V_E and V_N must have the same size.');
    end

    if nTimeE ~= N
        error('Number of columns in V_E/V_N must equal length(BT_E/BT_N).');
    end

    % ---------------------------------------------------------------------
    % Depth handling: scalar or vector length N
    % ---------------------------------------------------------------------
    if isscalar(Depth)
        Depth = repmat(Depth,N,1);
    else
        Depth = Depth(:);
    end

    if numel(Depth) ~= N
        error('Depth must be scalar or length N.');
    end

    % ---------------------------------------------------------------------
    % CellSize handling: scalar or vector length N
    % ---------------------------------------------------------------------
    if isscalar(CellSize)
        CellSize = repmat(CellSize,N,1);
    else
        CellSize = CellSize(:);
    end

    if numel(CellSize) ~= N
        error('CellSize must be scalar or length N.');
    end

    if any(~isfinite(CellSize) | CellSize <= 0)
        error('All CellSize values must be finite and > 0.');
    end

    % ---------------------------------------------------------------------
    % Parse options
    % ---------------------------------------------------------------------
    p = inputParser;
    addParameter(p,'AngleThreshold',90,@(x)isnumeric(x) && isscalar(x) && x>=0 && x<=180);
    addParameter(p,'Fi',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=-1 && x<=1));
    addParameter(p,'ProfileMethod','full',@(x)ischar(x) || isstring(x));
    addParameter(p,'SmoothWindow',5,@(x)isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p,'StartCell',1,@(x)isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p,'EndCell',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
    addParameter(p,'FirstCellCenter',[],@(x)isempty(x) || isscalar(x) || isvector(x));
    addParameter(p,'UseAbsDepth',false,@(x)islogical(x) || isnumeric(x));
    parse(p,varargin{:});

    AngleThreshold  = p.Results.AngleThreshold;
    Fi              = p.Results.Fi;
    ProfileMethod   = lower(char(p.Results.ProfileMethod));
    SmoothWindow    = max(1,round(p.Results.SmoothWindow));
    StartCell       = round(p.Results.StartCell);
    EndCell         = p.Results.EndCell;
    FirstCellCenter = p.Results.FirstCellCenter;
    UseAbsDepth     = logical(p.Results.UseAbsDepth);

    if isempty(Fi)
        Fi = cosd(AngleThreshold);
    else
        AngleThreshold = acosd(Fi);
    end

    if ~ismember(ProfileMethod,{'full','last','last2avg'})
        error('ProfileMethod must be ''full'', ''last'', or ''last2avg''.');
    end

    if isempty(EndCell)
        EndCell = nCellsE;
    else
        EndCell = round(EndCell);
    end

    if StartCell > EndCell || StartCell < 1 || EndCell > nCellsE
        error('StartCell/EndCell must define a valid range within the number of cells.');
    end

    % ---------------------------------------------------------------------
    % FirstCellCenter handling
    % ---------------------------------------------------------------------
    if isempty(FirstCellCenter)
        FirstCellCenter = CellSize/2;
    elseif isscalar(FirstCellCenter)
        FirstCellCenter = repmat(FirstCellCenter,N,1);
    else
        FirstCellCenter = FirstCellCenter(:);
    end

    if numel(FirstCellCenter) ~= N
        error('FirstCellCenter must be empty, scalar, or length N.');
    end

    if any(~isfinite(FirstCellCenter) | FirstCellCenter < 0)
        error('All FirstCellCenter values must be finite and >= 0.');
    end

    % ---------------------------------------------------------------------
    % Smooth water profiles in time
    % Each row = one cell time series
    % ---------------------------------------------------------------------
    V_E_smooth = movmean(V_E,SmoothWindow,2,'omitnan');
    V_N_smooth = movmean(V_N,SmoothWindow,2,'omitnan');

    % ---------------------------------------------------------------------
    % Preallocate reference vectors
    % ---------------------------------------------------------------------
    VE_ref               = NaN(N,1);
    VN_ref               = NaN(N,1);
    lastValidCell        = NaN(N,1);
    secondLastValidCell  = NaN(N,1);
    validCellCount       = zeros(N,1);
    cellCentersLastValid = NaN(N,1);

    % ---------------------------------------------------------------------
    % Build reference water vector for each time sample
    % ---------------------------------------------------------------------
    for ii = 1:N

        depthNow = Depth(ii);
        if UseAbsDepth
            depthNow = abs(depthNow);
        end

        cellSizeNow    = CellSize(ii);
        firstCenterNow = FirstCellCenter(ii);

        if ~isfinite(depthNow) || isnan(depthNow) || depthNow <= 0
            continue
        end

        % Cell centers at this time step
        cellCentersNow = firstCenterNow + (0:nCellsE-1)' * cellSizeNow;

        % Restrict to requested cell range
        candidateCells = (StartCell:EndCell)';

        % Keep only cells whose center is above bottom
        validCells = candidateCells(cellCentersNow(candidateCells) < depthNow);

        % Require at least one finite value in E or N
        if ~isempty(validCells)
            keepFinite = isfinite(V_E_smooth(validCells,ii)) | isfinite(V_N_smooth(validCells,ii));
            validCells = validCells(keepFinite);
        end

        if isempty(validCells)
            continue
        end

        validCellCount(ii)  = numel(validCells);
        lastValidCell(ii)   = validCells(end);
        cellCentersLastValid(ii) = cellCentersNow(validCells(end));

        if numel(validCells) >= 2
            secondLastValidCell(ii) = validCells(end-1);
        end

        switch ProfileMethod
            case 'full'
                VE_ref(ii) = mean(V_E_smooth(validCells,ii),'omitnan');
                VN_ref(ii) = mean(V_N_smooth(validCells,ii),'omitnan');

            case 'last'
                k = validCells(end);
                VE_ref(ii) = V_E_smooth(k,ii);
                VN_ref(ii) = V_N_smooth(k,ii);

            case 'last2avg'
                if numel(validCells) >= 2
                    k1 = validCells(end-1);
                    k2 = validCells(end);
                    VE_ref(ii) = mean([V_E_smooth(k1,ii), V_E_smooth(k2,ii)],'omitnan');
                    VN_ref(ii) = mean([V_N_smooth(k1,ii), V_N_smooth(k2,ii)],'omitnan');
                else
                    k = validCells(end);
                    VE_ref(ii) = V_E_smooth(k,ii);
                    VN_ref(ii) = V_N_smooth(k,ii);
                end
        end
    end

    % ---------------------------------------------------------------------
    % Magnitudes
    % ---------------------------------------------------------------------
    BTmod = hypot(BT_E_in,BT_N_in);
    Vmod  = hypot(VE_ref,VN_ref);

    % ---------------------------------------------------------------------
    % Angles
    % ---------------------------------------------------------------------
    alphaBT = mod(atan2d(BT_N_in,BT_E_in),360);
    alphaV  = mod(atan2d(VN_ref,VE_ref),360);

    % ---------------------------------------------------------------------
    % Dot product / cosine of angle
    % ---------------------------------------------------------------------
    COSalpha = NaN(N,1);
    validDot = isfinite(BT_E_in) & isfinite(BT_N_in) & ...
               isfinite(VE_ref)  & isfinite(VN_ref)  & ...
               isfinite(BTmod)   & isfinite(Vmod)    & ...
               (BTmod > 0)       & (Vmod > 0);

    COSalpha(validDot) = (VE_ref(validDot).*BT_E_in(validDot) + ...
                          VN_ref(validDot).*BT_N_in(validDot)) ./ ...
                         (Vmod(validDot).*BTmod(validDot));

    % Numerical safety
    COSalpha = max(-1,min(1,COSalpha));

    AngleDiff = NaN(N,1);
    AngleDiff(validDot) = acosd(COSalpha(validDot));

    % ---------------------------------------------------------------------
    % Filter mask:
    % remove BT when angle difference is larger than threshold
    % equivalent to COSalpha < Fi
    % ---------------------------------------------------------------------
    maskAngle = false(N,1);
    maskAngle(validDot) = COSalpha(validDot) < Fi;

    % ---------------------------------------------------------------------
    % Apply filtering
    % ---------------------------------------------------------------------
    BT_E = BT_E_in;
    BT_N = BT_N_in;

    BT_E(maskAngle) = NaN;
    BT_N(maskAngle) = NaN;

    alphaBTf = mod(atan2d(BT_N,BT_E),360);

    % ---------------------------------------------------------------------
    % Percentages
    % ---------------------------------------------------------------------
    % Originally valid samples
    validE0 = isfinite(BT_E_in);
    validN0 = isfinite(BT_N_in);

    % Newly filtered samples
    nFilteredE = nnz(maskAngle & validE0);
    nFilteredN = nnz(maskAngle & validN0);

    if nnz(validN0) > 0
        FilteredDataN = nFilteredN / nnz(validN0) * 100;
    else
        FilteredDataN = NaN;
    end

    if nnz(validE0) > 0
        FilteredDataE = nFilteredE / nnz(validE0) * 100;
    else
        FilteredDataE = NaN;
    end

    % Total final NaN percentages relative to all samples
    TotalNaNDataN = nnz(isnan(BT_N)) / N * 100;
    TotalNaNDataE = nnz(isnan(BT_E)) / N * 100;

    % Original NaN percentages
    OriginalNaNDataN = nnz(isnan(BT_N_in)) / N * 100;
    OriginalNaNDataE = nnz(isnan(BT_E_in)) / N * 100;

    % ---------------------------------------------------------------------
    % Return original orientation
    % ---------------------------------------------------------------------
    if wasRowE
        BT_E     = BT_E.';
        alphaV   = alphaV.';
        alphaBT  = alphaBT.';
        alphaBTf = alphaBTf.';
        Vmod     = Vmod.';
    end

    if wasRowN
        BT_N = BT_N.';
    end

    % ---------------------------------------------------------------------
    % Diagnostics
    % ---------------------------------------------------------------------
    stats = struct();
    stats.VE_ref               = VE_ref;
    stats.VN_ref               = VN_ref;
    stats.V_E_smooth           = V_E_smooth;
    stats.V_N_smooth           = V_N_smooth;
    stats.BTmod                = BTmod;
    stats.Vmod                 = hypot(VE_ref,VN_ref);
    stats.COSalpha             = COSalpha;
    stats.AngleDiff            = AngleDiff;
    stats.maskAngle            = maskAngle;
    stats.lastValidCell        = lastValidCell;
    stats.secondLastValidCell  = secondLastValidCell;
    stats.validCellCount       = validCellCount;
    stats.cellCentersLastValid = cellCentersLastValid;
    stats.ProfileMethod        = ProfileMethod;
    stats.AngleThreshold       = AngleThreshold;
    stats.Fi                   = Fi;
    stats.StartCell            = StartCell;
    stats.EndCell              = EndCell;
    stats.Depth                = Depth;
    stats.CellSize             = CellSize;
    stats.FirstCellCenter      = FirstCellCenter;

    % Added NaN / filtering diagnostics
    stats.nSamples            = N;
    stats.nValidOriginalN     = nnz(validN0);
    stats.nValidOriginalE     = nnz(validE0);
    stats.nFilteredN          = nFilteredN;
    stats.nFilteredE          = nFilteredE;
    stats.FilteredDataN       = FilteredDataN;
    stats.FilteredDataE       = FilteredDataE;
    stats.OriginalNaNDataN    = OriginalNaNDataN;
    stats.OriginalNaNDataE    = OriginalNaNDataE;
    stats.TotalNaNDataN       = TotalNaNDataN;
    stats.TotalNaNDataE       = TotalNaNDataE;
end