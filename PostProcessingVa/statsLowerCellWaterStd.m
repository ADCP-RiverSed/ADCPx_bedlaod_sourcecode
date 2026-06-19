function stats = statsLowerCellWaterStd(V1x_1,V1y_1,DepthBeam1,DepthBeam2,DepthBeam3,DepthBeam4,CellSize,varargin)
% STATSLOWERCELLWATERSTD
% Standard deviation of water velocity in the lower cells.
%
% Inputs:
%   V1x_1, V1y_1   : water velocity profiles [nCells x N]
%   DepthBeam1..4  : beam depths [1xN or Nx1]
%   CellSize       : scalar or vector [1xN or Nx1], in meters
%
% Optional name-value inputs:
%   'DepthMethod'      : 'median' (default), 'min', 'max', 'mean'
%   'FirstCellCenter'  : [] (default = CellSize/2), scalar, or vector length N
%   'nLowerCells'      : number of deepest valid cells to use, default = 2
%   'UseAbsDepth'      : true/false, default = false
%
% Output:
%   stats : structure with lower-cell means/stds and per-time extracted series

    p = inputParser;
    addParameter(p,'DepthMethod','median',@(x)ischar(x) || isstring(x));
    addParameter(p,'FirstCellCenter',[],@(x)isempty(x) || isscalar(x) || isvector(x));
    addParameter(p,'nLowerCells',2,@(x)isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p,'UseAbsDepth',false,@(x)islogical(x) || isnumeric(x));
    parse(p,varargin{:});

    DepthMethod     = lower(char(p.Results.DepthMethod));
    FirstCellCenter = p.Results.FirstCellCenter;
    nLowerCells     = round(p.Results.nLowerCells);
    UseAbsDepth     = logical(p.Results.UseAbsDepth);

    [nCells,N] = size(V1x_1);

    if ~isequal(size(V1x_1),size(V1y_1))
        error('V1x_1 and V1y_1 must have the same size.');
    end

    D1 = DepthBeam1(:);
    D2 = DepthBeam2(:);
    D3 = DepthBeam3(:);
    D4 = DepthBeam4(:);

    if ~(numel(D1)==N && numel(D2)==N && numel(D3)==N && numel(D4)==N)
        error('DepthBeam1..4 must all have length N.');
    end

    if isscalar(CellSize)
        CellSize = repmat(CellSize,N,1);
    else
        CellSize = CellSize(:);
    end

    if numel(CellSize) ~= N
        error('CellSize must be scalar or length N.');
    end

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

    % Representative depth from the 4 beams
    Dall = [D1 D2 D3 D4];
    switch DepthMethod
        case 'median'
            DepthRef = median(Dall,2,'omitnan');
        case 'mean'
            DepthRef = mean(Dall,2,'omitnan');
        case 'min'
            DepthRef = min(Dall,[],2,'omitnan');
        case 'max'
            DepthRef = max(Dall,[],2,'omitnan');
        otherwise
            error('DepthMethod must be ''median'', ''mean'', ''min'', or ''max''.');
    end

    if UseAbsDepth
        DepthRef = abs(DepthRef);
    end

    % Outputs per time
    lowerX_mean = NaN(N,1);
    lowerY_mean = NaN(N,1);
    lowerMag_mean = NaN(N,1);

    lowerX_std = NaN(N,1);
    lowerY_std = NaN(N,1);
    lowerMag_std = NaN(N,1);

    lastValidCell = NaN(N,1);
    nValidLowerCells = zeros(N,1);

    for t = 1:N

        depthNow = DepthRef(t);
        if ~isfinite(depthNow) || isnan(depthNow) || depthNow <= 0
            continue
        end

        cellCenters = FirstCellCenter(t) + (0:nCells-1)' * CellSize(t);

        validCells = find(cellCenters < depthNow);

        if isempty(validCells)
            continue
        end

        lastValidCell(t) = validCells(end);

        idx = validCells(max(1,end-nLowerCells+1):end);

        vx = V1x_1(idx,t);
        vy = V1y_1(idx,t);
        vm = hypot(vx,vy);

        lowerX_mean(t) = mean(vx,'omitnan');
        lowerY_mean(t) = mean(vy,'omitnan');
        lowerMag_mean(t) = mean(vm,'omitnan');

        lowerX_std(t) = std(vx,'omitnan');
        lowerY_std(t) = std(vy,'omitnan');
        lowerMag_std(t) = std(vm,'omitnan');

        nValidLowerCells(t) = numel(idx);
    end

    stats = struct();
    stats.DepthRef = DepthRef;
    stats.lastValidCell = lastValidCell;
    stats.nValidLowerCells = nValidLowerCells;

    stats.lowerX_mean = lowerX_mean;
    stats.lowerY_mean = lowerY_mean;
    stats.lowerMag_mean = lowerMag_mean;

    stats.lowerX_std = lowerX_std;
    stats.lowerY_std = lowerY_std;
    stats.lowerMag_std = lowerMag_std;

    stats.lowerX_std_mean = mean(lowerX_std,'omitnan');
    stats.lowerY_std_mean = mean(lowerY_std,'omitnan');
    stats.lowerMag_std_mean = mean(lowerMag_std,'omitnan');

    stats.lowerX_std_global = std(lowerX_mean,'omitnan');
    stats.lowerY_std_global = std(lowerY_mean,'omitnan');
    stats.lowerMag_std_global = std(lowerMag_mean,'omitnan');
end