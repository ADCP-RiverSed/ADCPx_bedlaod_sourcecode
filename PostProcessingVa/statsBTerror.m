function stats = statsBTerror(BTerr)
% STATSBTERROR
% Statistics of bottom-track error velocity.
%
% Input:
%   BTerr : bottom-track error velocity [1xN or Nx1]
%
% Output:
%   stats : structure with mean, std, median, min, max, RMS

    BTerr = BTerr(:);

    stats = struct();
    stats.mean   = mean(BTerr,'omitnan');
    stats.std    = std(BTerr,'omitnan');
    stats.median = median(BTerr,'omitnan');
    stats.min    = min(BTerr,[],'omitnan');
    stats.max    = max(BTerr,[],'omitnan');
    stats.rms    = sqrt(mean(BTerr.^2,'omitnan'));
    stats.nValid = nnz(isfinite(BTerr));
end