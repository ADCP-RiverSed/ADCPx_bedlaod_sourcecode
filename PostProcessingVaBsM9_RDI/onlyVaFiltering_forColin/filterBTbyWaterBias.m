function [BTx_f,BTy_f,maskFiltered,filteredIdx,pctFiltered,WaterBias,stats] = ...
    filterBTbyWaterBias(BTx,BTy,V1x_1,V1y_1,Depth,CellSize,varargin)
% FILTERBTBYWATERPROFILE
%
% Filter bottom-track (BT) velocities using near-bed water-profile velocity.
%
% -------------------------------------------------------------------------
% MAIN IDEA
% -------------------------------------------------------------------------
% For each time sample (column):
%   1) Smooth the velocity profiles V1x_1 and V1y_1 in time
%   2) Find the last valid profile cell above the bed using Depth and CellSize
%   3) Build a reference near-bed water velocity from:
%        - the last valid cell, or
%        - the average of the last two valid cells
%   4) Filter BT values if they are inconsistent with that reference
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% BTx, BTy
%   Bottom-track velocity components, vectors of length N
%
% V1x_1, V1y_1
%   Velocity profile matrices of size [nCells x N]
%   Rows    = profile cells
%   Columns = time samples
%
% Depth
%   Measured water depth for each time sample
%   Can be:
%     - scalar
%     - vector length N
%
% CellSize
%   Cell size. Can be:
%     - scalar
%     - vector length N
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE INPUTS
% -------------------------------------------------------------------------
% 'SmoothWindow'     : moving-average window along time, default = 5
%
% 'ReferenceMethod'  : 'last' or 'last2avg', default = 'last2avg'
%
% 'CompareMode'      : 'signedless' or 'absgreater', default = 'signedless'
%
% 'ApplyTo'          : 'either', 'both', 'xonly', 'yonly'
%                      default = 'either'
%
% 'FirstCellCenter'  : first cell-center distance from transducer
%                      Can be:
%                        []         -> default = CellSize/2
%                        scalar
%                        vector length N
%
% 'UseAbsDepthCheck' : true/false, default = false
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% BTx_f, BTy_f
%   Filtered BT components with flagged samples replaced by NaN
%
% maskFiltered
%   Final logical mask of filtered samples
%
% filteredIdx
%   Indices of filtered samples
%
% pctFiltered
%   Percentage of filtered samples relative to originally valid BT samples
%
% WaterBias
%   Cell array of labels:
%       'possible water bias'
%       'NO bias'
%       'insufficient valid cells'
%
% stats
%   Structure with diagnostics
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
% [BTx_f,BTy_f,maskFiltered,filteredIdx,pctFiltered,WaterBias,stats] = ...
%     filterBTbyWaterProfile(BTx_,BTy_,V1x_1,V1y_1,Depth,CellSize,...
%     'SmoothWindow',7,...
%     'ReferenceMethod','last2avg',...
%     'CompareMode','signedless',...
%     'ApplyTo','either');
%
% figure; hold on
% plot(BTx_,'k')
% plot(BTx_f,'b','LineWidth',1.2)
% plot(filteredIdx,BTx_(filteredIdx),'ro')
% plot(stats.refX,'g','LineWidth',1)
% legend('BTx raw','BTx filtered','Filtered','Reference X')
% title(['Filtered = ' num2str(pctFiltered,'%.2f') '%'])
% grid on

    % -----------------------------
    % Basic checks
    % -----------------------------
    if ~isvector(BTx) || ~isvector(BTy)
        error('BTx and BTy must be vectors.');
    end

    BTx = BTx(:);
    BTy = BTy(:);

    N = numel(BTx);

    if numel(BTy) ~= N
        error('BTx and BTy must have the same length.');
    end

    if ~ismatrix(V1x_1) || ~ismatrix(V1y_1)
        error('V1x_1 and V1y_1 must be 2D matrices of size [nCells x N].');
    end

    [nCells1,nTime1] = size(V1x_1);
    [nCells2,nTime2] = size(V1y_1);

    if nCells1 ~= nCells2 || nTime1 ~= nTime2
        error('V1x_1 and V1y_1 must have the same size.');
    end

    if nTime1 ~= N
        error('Number of columns in V1x_1 / V1y_1 must equal length(BTx).');
    end

    % Depth: scalar or vector length N
    if isscalar(Depth)
        Depth = repmat(Depth,N,1);
    else
        Depth = Depth(:);
    end

    if numel(Depth) ~= N
        error('Depth must be scalar or length N.');
    end

    % CellSize: scalar or vector length N
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

    % -----------------------------
    % Options
    % -----------------------------
    p = inputParser;
    addParameter(p,'SmoothWindow',5,@(x)isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p,'ReferenceMethod','last2avg',@(x)ischar(x) || isstring(x));
    addParameter(p,'CompareMode','signedless',@(x)ischar(x) || isstring(x));
    addParameter(p,'ApplyTo','either',@(x)ischar(x) || isstring(x));
    addParameter(p,'FirstCellCenter',[],@(x)isempty(x) || isscalar(x) || isvector(x));
    addParameter(p,'UseAbsDepthCheck',false,@(x)islogical(x) || isnumeric(x));
    parse(p,varargin{:});

    SmoothWindow     = max(1,round(p.Results.SmoothWindow));
    ReferenceMethod  = lower(char(p.Results.ReferenceMethod));
    CompareMode      = lower(char(p.Results.CompareMode));
    ApplyTo          = lower(char(p.Results.ApplyTo));
    FirstCellCenter  = p.Results.FirstCellCenter;
    UseAbsDepthCheck = logical(p.Results.UseAbsDepthCheck);

    if ~ismember(ReferenceMethod,{'last','last2avg'})
        error('ReferenceMethod must be ''last'' or ''last2avg''.');
    end

    if ~ismember(CompareMode,{'signedless','absgreater'})
        error('CompareMode must be ''signedless'' or ''absgreater''.');
    end

    if ~ismember(ApplyTo,{'either','both','xonly','yonly'})
        error('ApplyTo must be ''either'', ''both'', ''xonly'', or ''yonly''.');
    end

    % FirstCellCenter: [] or scalar or length N
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

    % -----------------------------
    % Smooth profiles in time
    % Each row = one cell time series
    % -----------------------------
    V1x_smooth = movmean(V1x_1,SmoothWindow,2,'omitnan');
    V1y_smooth = movmean(V1y_1,SmoothWindow,2,'omitnan');

    % -----------------------------
    % Preallocate outputs
    % -----------------------------
    BTx_f = BTx;
    BTy_f = BTy;

    lastValidCell       = NaN(N,1);
    secondLastValidCell = NaN(N,1);
    refX                = NaN(N,1);
    refY                = NaN(N,1);

    maskX = false(N,1);
    maskY = false(N,1);
    maskFiltered = false(N,1);

    WaterBias = repmat({'NO bias'},N,1);

    % Optional diagnostics for variable geometry
    cellCentersLastValid = NaN(N,1);
    cellSizeUsed         = CellSize;
    firstCellCenterUsed  = FirstCellCenter;

    % -----------------------------
    % Loop over time
    % -----------------------------
    for bb = 1:N

        depthNow = Depth(bb);
        if UseAbsDepthCheck
            depthNow = abs(depthNow);
        end

        if ~isfinite(depthNow) || isnan(depthNow) || depthNow <= 0
            WaterBias{bb} = 'insufficient valid cells';
            continue
        end

        cellSizeNow = CellSize(bb);
        firstCenterNow = FirstCellCenter(bb);

        % Cell centers for this time step
        cellCentersNow = firstCenterNow + (0:nCells1-1)' * cellSizeNow;

        % Valid cells = cell centers above the bottom
        validCells = find(cellCentersNow < depthNow);

        % Also require at least one finite profile component in that cell
        validCells = validCells( ...
            isfinite(V1x_smooth(validCells,bb)) | isfinite(V1y_smooth(validCells,bb)) );

        if isempty(validCells)
            WaterBias{bb} = 'insufficient valid cells';
            continue
        end

        % Last valid cell
        lastIdx = validCells(end);
        lastValidCell(bb) = lastIdx;
        cellCentersLastValid(bb) = cellCentersNow(lastIdx);

        % Second last valid cell if available
        if numel(validCells) >= 2
            secondIdx = validCells(end-1);
            secondLastValidCell(bb) = secondIdx;
        else
            secondIdx = NaN;
        end

        % Build reference near-bed velocity
        switch ReferenceMethod
            case 'last'
                refX(bb) = V1x_smooth(lastIdx,bb);
                refY(bb) = V1y_smooth(lastIdx,bb);

            case 'last2avg'
                if isfinite(secondIdx)
                    refX(bb) = mean([V1x_smooth(secondIdx,bb), V1x_smooth(lastIdx,bb)],'omitnan');
                    refY(bb) = mean([V1y_smooth(secondIdx,bb), V1y_smooth(lastIdx,bb)],'omitnan');
                else
                    refX(bb) = V1x_smooth(lastIdx,bb);
                    refY(bb) = V1y_smooth(lastIdx,bb);
                end
        end

        % If no valid reference, skip
        if ~isfinite(refX(bb)) && ~isfinite(refY(bb))
            WaterBias{bb} = 'insufficient valid cells';
            continue
        end

        % Compare BT against reference
        switch CompareMode
            case 'signedless'
                if isfinite(refX(bb)) && isfinite(BTx(bb))
                    maskX(bb) = BTx(bb) < refX(bb);
                end
                if isfinite(refY(bb)) && isfinite(BTy(bb))
                    maskY(bb) = BTy(bb) < refY(bb);
                end

            case 'absgreater'
                if isfinite(refX(bb)) && isfinite(BTx(bb))
                    maskX(bb) = abs(BTx(bb)) > abs(refX(bb));
                end
                if isfinite(refY(bb)) && isfinite(BTy(bb))
                    maskY(bb) = abs(BTy(bb)) > abs(refY(bb));
                end
        end

        % Combine X/Y masks
        switch ApplyTo
            case 'either'
                maskFiltered(bb) = maskX(bb) | maskY(bb);

            case 'both'
                maskFiltered(bb) = maskX(bb) & maskY(bb);

            case 'xonly'
                maskFiltered(bb) = maskX(bb);

            case 'yonly'
                maskFiltered(bb) = maskY(bb);
        end

        % Label
        if maskFiltered(bb)
            WaterBias{bb} = 'possible water bias';
            BTx_f(bb) = NaN;
            BTy_f(bb) = NaN;
        else
            WaterBias{bb} = 'NO bias';
        end
    end

    % Do not count already invalid BT as filtered
    validBT0 = isfinite(BTx) & isfinite(BTy);
    maskFiltered(~validBT0) = false;

    % Ensure outputs are NaN where final mask is true
    BTx_f(maskFiltered) = NaN;
    BTy_f(maskFiltered) = NaN;

    filteredIdx = find(maskFiltered);
    nFiltered   = nnz(maskFiltered);
    nValidBT    = nnz(validBT0);

    if nValidBT > 0
        pctFiltered = 100 * nFiltered / nValidBT;
    else
        pctFiltered = NaN;
    end

    % Diagnostics
    stats = struct();
    stats.V1x_smooth           = V1x_smooth;
    stats.V1y_smooth           = V1y_smooth;
    stats.lastValidCell        = lastValidCell;
    stats.secondLastValidCell  = secondLastValidCell;
    stats.lastValidCellCenter  = cellCentersLastValid;
    stats.refX                 = refX;
    stats.refY                 = refY;
    stats.maskX                = maskX;
    stats.maskY                = maskY;
    stats.maskFiltered         = maskFiltered;
    stats.filteredIdx          = filteredIdx;
    stats.nFiltered            = nFiltered;
    stats.pctFiltered          = pctFiltered;
    stats.ReferenceMethod      = ReferenceMethod;
    stats.CompareMode          = CompareMode;
    stats.ApplyTo              = ApplyTo;
    stats.SmoothWindow         = SmoothWindow;
    stats.Depth                = Depth;
    stats.CellSize             = cellSizeUsed;
    stats.FirstCellCenter      = firstCellCenterUsed;
end