function stats = statsBTrawStd(BTx,BTy)
% STATSBTRAWSTD
% Statistics of raw bottom-track velocity in x, y, and magnitude.
%
% Inputs:
%   BTx, BTy : bottom-track velocity components [1xN or Nx1]
%
% Output:
%   stats : structure with mean, std, median, min, max, and BT magnitude

    BTx = BTx(:);
    BTy = BTy(:);

    BTmag = hypot(BTx,BTy);

    stats = struct();
    stats.BTx_mean   = mean(BTx,'omitnan');
    stats.BTx_std    = std(BTx,'omitnan');
    stats.BTx_median = median(BTx,'omitnan');
    stats.BTx_min    = min(BTx,[],'omitnan');
    stats.BTx_max    = max(BTx,[],'omitnan');

    stats.BTy_mean   = mean(BTy,'omitnan');
    stats.BTy_std    = std(BTy,'omitnan');
    stats.BTy_median = median(BTy,'omitnan');
    stats.BTy_min    = min(BTy,[],'omitnan');
    stats.BTy_max    = max(BTy,[],'omitnan');

    stats.BTmag      = BTmag;
    stats.BTmag_mean = mean(BTmag,'omitnan');
    stats.BTmag_std  = std(BTmag,'omitnan');
    stats.BTmag_median = median(BTmag,'omitnan');
    stats.BTmag_min  = min(BTmag,[],'omitnan');
    stats.BTmag_max  = max(BTmag,[],'omitnan');

    stats.nValid = nnz(isfinite(BTx) & isfinite(BTy));
end