function [BT,nFiltered,pctFiltered,filteredIdx,spikeMask,info] = despikeBT_01(BTv,varargin)
% DESPIKEBT  Despike BT ADCP data using absolute-amplitude thresholds
%
% This version does NOT use signed residuals such as (A - baseline)
% for the main filtering logic. Instead, it works on:
%
%     x = abs(A)
%
% so points are filtered based on their amplitude.
%
% METHODS
% -------------------------------------------------------------------------
% 'globalstd'
%   Global amplitude filter:
%       remove if abs(A) > Threshold * std(A)
%
% 'globalmad'
%   Robust global amplitude filter:
%       x      = abs(A)
%       medX   = median(x)
%       MADX   = median(abs(x - medX))
%       scale  = 1.4826 * MADX
%       remove if x > medX + Threshold * scale
%
% 'localmean'
%   Local amplitude filter using moving mean and moving std on abs(A):
%       x         = abs(A)
%       localMean = movmean(x,Window)
%       localStd  = movstd(x,Window)
%       remove if x > localMean + Threshold * localStd
%
% 'localmad'
%   Robust local amplitude filter using moving median and moving MAD on abs(A):
%       x         = abs(A)
%       localMed  = movmedian(x,Window)
%       localMAD  = movmedian(abs(x - localMed),Window)
%       scale     = 1.4826 * localMAD
%       remove if x > localMed + Threshold * scale
%
% ABOUT 1.4826
% -------------------------------------------------------------------------
% 1.4826 converts MAD to a sigma-like scale estimate:
%
%     sigma_robust ~= 1.4826 * MAD
%
% for normally distributed data.
%
% INPUT
% -------------------------------------------------------------------------
% BTv : input vector
%
% Name-value options:
%   'Method'      : 'globalstd', 'globalmad', 'localmean', 'localmad'
%                   default = 'globalstd'
%   'Threshold'   : threshold multiplier, default = 3
%   'Window'      : moving window size for local methods, default = 11
%   'Iterations'  : number of despiking passes, default = 1
%
% OUTPUT
% -------------------------------------------------------------------------
% BT           : filtered signal, with filtered values replaced by NaN
% nFiltered    : number of filtered values
% pctFiltered  : percent of filtered values relative to original non-NaN data
% filteredIdx  : indices of filtered values
% spikeMask    : logical mask of filtered values
% info         : structure with details
%
% EXAMPLES
% -------------------------------------------------------------------------
%
% % Apply global std and plot
% [BT_std,n_std,pct_std,idx_std,mask_std,info_std] = despikeBT(BTv,...
%     'Method','globalstd','Threshold',3,'Iterations',1);
%
% figure; hold on
% plot(BTv,'k')
% plot(BT_std,'b','LineWidth',1.2)
% plot(idx_std,BTv(idx_std),'ro','MarkerSize',6)
% legend('Raw','Filtered','Removed spikes')
% title(['globalstd | removed = ' num2str(n_std) ...
%        ' (' num2str(pct_std,'%.2f') '%)'])
% xlabel('Sample'); ylabel('BT')
% grid on
%
%
% % Apply local MAD and plot
% [BT_lmad,n_lmad,pct_lmad,idx_lmad,mask_lmad,info_lmad] = despikeBT(BTv,...
%     'Method','localmad','Window',21,'Threshold',3.5,'Iterations',1);
%
% figure; hold on
% plot(BTv,'k')
% plot(BT_lmad,'b','LineWidth',1.2)
% plot(idx_lmad,BTv(idx_lmad),'ro','MarkerSize',6)
% legend('Raw','Filtered','Removed spikes')
% title(['localmad | removed = ' num2str(n_lmad) ...
%        ' (' num2str(pct_lmad,'%.2f') '%)'])
% xlabel('Sample'); ylabel('BT')
% grid on
%
%
% % Apply local mean and plot
% [BT_lmean,n_lmean,pct_lmean,idx_lmean,mask_lmean,info_lmean] = despikeBT(BTv,...
%     'Method','localmean','Window',21,'Threshold',3,'Iterations',1);
%
% figure; hold on
% plot(BTv,'k')
% plot(BT_lmean,'b','LineWidth',1.2)
% plot(idx_lmean,BTv(idx_lmean),'ro','MarkerSize',6)
% legend('Raw','Filtered','Removed spikes')
% title(['localmean | removed = ' num2str(n_lmean) ...
%        ' (' num2str(pct_lmean,'%.2f') '%)'])
% xlabel('Sample'); ylabel('BT')
% grid on

    if ~isvector(BTv)
        error('BTv must be a vector.');
    end

    wasRow = isrow(BTv);
    A0 = BTv(:);
    A  = A0;

    p = inputParser;
    addParameter(p,'Method','globalstd',@(x)ischar(x) || isstring(x));
    addParameter(p,'Threshold',3,@(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p,'Window',11,@(x)isnumeric(x) && isscalar(x) && x>=3);
    addParameter(p,'Iterations',1,@(x)isnumeric(x) && isscalar(x) && x>=1);
    parse(p,varargin{:});

    method = lower(char(p.Results.Method));
    T      = p.Results.Threshold;
    win    = round(p.Results.Window);
    nIter  = round(p.Results.Iterations);

    if mod(win,2)==0
        win = win + 1;
    end

    originalNaN = isnan(A);
    nValid0     = nnz(~originalNaN);

    spikeMask = false(size(A));
    lastCenter = NaN(size(A));
    lastScale  = NaN(size(A));
    lastThresh = NaN(size(A));

    for k = 1:nIter
        x = abs(A);

        switch method

            case 'globalstd'
                center = zeros(size(A));          % no center used
                scale  = nanstd_local(A);         % std of raw series
                thresh = T * scale;
                currentMask = x > thresh;

            case 'globalmad'
                medX   = nanmedian_local(x);
                madX   = nanmedian_local(abs(x - medX));

                % 1.4826 converts MAD to a sigma-like robust scale
                scale  = 1.4826 * madX;

                center = medX * ones(size(A));
                thresh = center + T * scale;
                currentMask = x > thresh;

            case 'localmean'
                center = movmean(x,win,'omitnan');
                scale  = movstd(x,win,'omitnan');
                scale  = enforce_min_scale(scale);
                thresh = center + T .* scale;
                currentMask = x > thresh;

            case 'localmad'
                center = movmedian(x,win,'omitnan');
                localMAD = movmedian(abs(x - center),win,'omitnan');

                % Convert local MAD to a local sigma-like robust scale
                scale  = 1.4826 * localMAD;
                scale  = enforce_min_scale(scale);
                thresh = center + T .* scale;
                currentMask = x > thresh;

            otherwise
                error('Unknown Method. Use globalstd, globalmad, localmean, or localmad.');
        end

        if isscalar(scale)
            scale = scale * ones(size(A));
        end
        if isscalar(thresh)
            thresh = thresh * ones(size(A));
        end

        currentMask(originalNaN | isnan(A)) = false;
        currentMask(~isfinite(thresh)) = false;

        if ~any(currentMask)
            lastCenter = center;
            lastScale  = scale;
            lastThresh = thresh;
            break
        end

        spikeMask = spikeMask | currentMask;
        A(currentMask) = NaN;

        lastCenter = center;
        lastScale  = scale;
        lastThresh = thresh;
    end

    BT = A;
    nFiltered = nnz(spikeMask);
    filteredIdx = find(spikeMask);

    if nValid0 > 0
        pctFiltered = 100 * nFiltered / nValid0;
    else
        pctFiltered = NaN;
    end

    if wasRow
        BT         = BT.';
        spikeMask  = spikeMask.';
        lastCenter = lastCenter.';
        lastScale  = lastScale.';
        lastThresh = lastThresh.';
    end

    info = struct();
    info.method        = method;
    info.threshold     = T;
    info.window        = win;
    info.iterations    = nIter;
    info.nFiltered     = nFiltered;
    info.pctFiltered   = pctFiltered;
    info.filteredIdx   = filteredIdx;
    info.spikeMask     = spikeMask;
    info.lastCenter    = lastCenter;
    info.lastScale     = lastScale;
    info.lastThreshold = lastThresh;
    info.originalNaN   = originalNaN;
    info.nOriginalValid = nValid0;
end

% -------------------------------------------------------------------------
% Helper functions
% -------------------------------------------------------------------------
function s = nanstd_local(x)
    x = x(~isnan(x));
    if numel(x) < 2
        s = NaN;
    else
        s = std(x);
    end
end

function m = nanmedian_local(x)
    x = x(~isnan(x));
    if isempty(x)
        m = NaN;
    else
        m = median(x);
    end
end

function s = enforce_min_scale(s)
    medS = median(s(~isnan(s) & isfinite(s)));
    if isempty(medS) || medS <= 0
        medS = eps;
    end
    minS = max(0.05*medS, eps);
    s(~isfinite(s) | s < minS) = minS;
end