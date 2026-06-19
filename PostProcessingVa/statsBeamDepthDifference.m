function stats = statsBeamDepthDifference(DepthBeam1,DepthBeam2,DepthBeam3,DepthBeam4)
% STATSBEAMDEPTHDIFFERENCE
% Statistics of beam depth differences.
%
% Inputs:
%   DepthBeam1..DepthBeam4 : beam ranges/depths [1xN or Nx1]
%
% Output:
%   stats : structure with pairwise beam differences and temporal changes

    D1 = DepthBeam1(:);
    D2 = DepthBeam2(:);
    D3 = DepthBeam3(:);
    D4 = DepthBeam4(:);

    if ~(numel(D1)==numel(D2) && numel(D1)==numel(D3) && numel(D1)==numel(D4))
        error('All beam depth vectors must have the same length.');
    end

    % Pairwise differences
    d12 = D1 - D2;
    d13 = D1 - D3;
    d14 = D1 - D4;
    d23 = D2 - D3;
    d24 = D2 - D4;
    d34 = D3 - D4;

    Dall = [D1 D2 D3 D4];
    beamSpread = max(Dall,[],2,'omitnan') - min(Dall,[],2,'omitnan');

    % Temporal change per beam
    dD1_time = [NaN; diff(D1)];
    dD2_time = [NaN; diff(D2)];
    dD3_time = [NaN; diff(D3)];
    dD4_time = [NaN; diff(D4)];

    stats = struct();

    stats.d12 = d12; stats.d13 = d13; stats.d14 = d14;
    stats.d23 = d23; stats.d24 = d24; stats.d34 = d34;

    stats.meanAbs_d12 = mean(abs(d12),'omitnan');
    stats.meanAbs_d13 = mean(abs(d13),'omitnan');
    stats.meanAbs_d14 = mean(abs(d14),'omitnan');
    stats.meanAbs_d23 = mean(abs(d23),'omitnan');
    stats.meanAbs_d24 = mean(abs(d24),'omitnan');
    stats.meanAbs_d34 = mean(abs(d34),'omitnan');

    stats.std_d12 = std(d12,'omitnan');
    stats.std_d13 = std(d13,'omitnan');
    stats.std_d14 = std(d14,'omitnan');
    stats.std_d23 = std(d23,'omitnan');
    stats.std_d24 = std(d24,'omitnan');
    stats.std_d34 = std(d34,'omitnan');

    stats.beamSpread = beamSpread;
    stats.beamSpread_mean = mean(beamSpread,'omitnan');
    stats.beamSpread_std  = std(beamSpread,'omitnan');
    stats.beamSpread_max  = max(beamSpread,[],'omitnan');

    stats.dD1_time = dD1_time;
    stats.dD2_time = dD2_time;
    stats.dD3_time = dD3_time;
    stats.dD4_time = dD4_time;

    stats.dD1_time_std = std(dD1_time,'omitnan');
    stats.dD2_time_std = std(dD2_time,'omitnan');
    stats.dD3_time_std = std(dD3_time,'omitnan');
    stats.dD4_time_std = std(dD4_time,'omitnan');
end