% Testing and applying the functions under a basic scenario, as suggested in the
% monograph, for moderate to intensive sediment transport conditions with
% minimal external disturbances.
%
% Field applications should be tuned accordingly.
% The functions contain many options and can be adapted as needed.

clear all; clc

maindirADD = 'C:\Users\slavenc\OneDrive - NTNU';
addpath(genpath([maindirADD '\ADCPdataColin\functiosforprocesingNIWA']));

% 
BTx  = BT.BTV1; % bedlaod velocity in x 
BTy  = BT.BTV2;
BTerr= BT.BTV4 % error BT velocity
% btvE = BT.btvE;
% btvN = BT.btvN;
V1x_1    = v.v1; % water velocity in x 
V1y_1    = v.v2;
DepthBeam1    = BT.RNG1;
DepthBeam2    = BT.RNG2;
DepthBeam3    = BT.RNG3;
DepthBeam4    = BT.RNG4;

CellSize = Data.WS./100;

%% STEP 0 - Some statistics on the raw data
% a) Raw BT statistics
statsBT = statsBTrawStd(BTx,BTy);

% b) BT error statistics
statsErr = statsBTerror(BTerr);

% c) Beam depth difference statistics
statsDepth = statsBeamDepthDifference(DepthBeam1,DepthBeam2,DepthBeam3,DepthBeam4);

% d) Water velocity statistics in lower cells (in this case the last three,
% not always teh last - putting nLowerCells=1 is vaILBLE OR VALID...  MIN 3
%  IS RECCOMENDED. 
statsLower = statsLowerCellWaterStd(V1x_1,V1y_1,...
    DepthBeam1,DepthBeam2,DepthBeam3,DepthBeam4,mean(CellSize),...
    'DepthMethod','mean',...
    'nLowerCells',3);

%% STEP 1 - Threshold filtering
% This is applied to BTx and BTy (or btvE and btvN). Open teh function for
% details> 'Threshold' is Thresholdxstd
% there is a global and local option, global uses the entire time series
% for std and local you apply on a window. 

[BTx_1,nFilteredX,pctFilteredX,filteredIdX,spikeMaskX,infoX] = despikeBT_01(BTx,...
    'Method','globalstd',...
    'Threshold',3,...
    'Iterations',1);

[BTy_1,nFilteredY,pctFilteredY,filteredIdY,spikeMaskY,infoY] = despikeBT_01(BTy,...
    'Method','globalstd',...
    'Threshold',3,...
    'Iterations',1);

figure; hold on
plot(BTx,'k')
plot(BTx_1,'b','LineWidth',1.2)
plot(find(spikeMaskX),BTx(spikeMaskX),'ro','MarkerSize',6)
legend('Raw','Filtered','Removed spikes')
title(['globalstd | removed = ' num2str(nFilteredX)])
xlabel('Sample')
ylabel('BT velocity')
grid on

% Apply local MAD filtering and plot.
% In some cases this fits better. Usually I get similar results, but if I
% suspect changes in bedload intensity, e.g. dune passage, this is the
% preferred option. Check the function options.
%
% [BT_localmad,n2,p2,idx2,mask2,info2] = despikeBT_01(BTx,...
%     'Method','localmad','Window',60,'Threshold',3,'Iterations',1);
%
% figure; hold on
% plot(BTx,'k')
% plot(BT_localmad,'b','LineWidth',1.2)
% plot(find(mask2),BTx(mask2),'ro','MarkerSize',6)
% legend('Raw','Filtered','Removed spikes')
% title(['localmad | removed = ' num2str(n2)])
% xlabel('Sample')
% ylabel('BT velocity')
% grid on

%% STEP 2 - GPS filter
% The function recognises either UTM coordinates (m) or lon/lat (degrees).
% Tested mainly with M9. The RDI setup in Germany is fixed on the boat.
%
% Typical GPS threshold is 0.2-0.3 m/s. Very rarely should BT velocity be
% larger than 0.3 m/s. The maximum I measured was slightly below 0.2 m/s.

GPSlong = GPS.Longitude;
GPSlat  = GPS.Latitude;

% BTx_1 = BottomTrack.BT_Vel(:,1);
% BTy_1 = BottomTrack.BT_Vel(:,2);

VTG    = [];
deltaT = diff(System.Time);

[BTx_2,BTy_2,BTmag_f,maskFiltered,filteredIdx,pctFiltered,stats] = ...
    filterBTbyGPSVTG(GPSlong,GPSlat,deltaT,BTx_1,BTy_1,VTG,...
    'CoordType','deg', ...
    'GPSThreshold',0.25, ...
    'VTGThreshold',0.25);

figure
subplot(2,1,1)
plot(stats.gpsStepDistance,'k')
ylabel('Distance per step [m]')
title('GPS step distance from lon/lat coordinates')
grid on

subplot(2,1,2)
plot(stats.gpsSpeed,'b')
ylabel('GPS speed')
xlabel('Sample')
title('GPS-derived speed')
grid on

figure; hold on
plot(stats.BTmag_raw,'k')
plot(BTmag_f,'b','LineWidth',1.2)
plot(stats.gpsSpeed,'g')
plot(filteredIdx,stats.BTmag_raw(filteredIdx),'ro')
legend('Raw BT magnitude','Filtered BT magnitude','GPS speed','Filtered')
title(['Filtered = ' num2str(pctFiltered,'%.2f') '%'])
grid on

%% STEP 3 - Upper water-bias filter
% Deletes all BT values with magnitude larger than the water velocity in
% the last valid cell, or optionally the average of the last two valid cells.
%
% Important:
% - v and BT must be in the same coordinate system
% - Depth and CellSize must be in meters
% - The reference is taken from the last valid cell or the average of the
%   last two valid cells above the bed

BTx_2    = BTx_1;
BTy_2    = BTx_1;
V1x_1    = v.v1;
V1y_1    = v.v2;
Depth    = BT.RNG;
CellSize = Data.WS./100;

% In this example I did not have GPS values, so the BT data are read again.

[BTx_3,BTy_3,maskFiltered,filteredIdx,pctFiltered,WaterBias,stats] = ...
    filterBTbyWaterBias(BTx_2,BTy_2,V1x_1,V1y_1,Depth,CellSize,...
    'SmoothWindow',5,...
    'ReferenceMethod','last',...
    'CompareMode','absgreater');

%% STEP 4 - Water-direction filter
% Removes BT values when the angle between the BT vector and the
% corresponding water-velocity reference vector is larger than 90 degrees.

[BT_Nf,BT_Ef,pctFilteredN,pctFilteredE,alphaV,alphaBTf,alphaBT,Vmod,stats,totalNaNPctN,totalNaNPctE] = ...
    filterBT2WaterProfile(V1x_1,V1y_1,BTx_3,BTy_3,Depth,CellSize,...
    'ProfileMethod','full',...
    'AngleThreshold',90,...
    'SmoothWindow',5);

BT_E_raw = BTx(:);
BT_N_raw = BTy(:);
BTmag_raw = hypot(BT_E_raw,BT_N_raw);
BTmag_f   = hypot(BT_Ef(:),BT_Nf(:));
idxFilt   = find(stats.maskAngle);

figure
subplot(3,1,1); hold on
plot(BT_E_raw,'k'); plot(BT_Ef,'b'); plot(idxFilt,BT_E_raw(idxFilt),'ro')
title(['BT_E | filtered = ' num2str(pctFilteredE,'%.2f') '%']); grid on

subplot(3,1,2); hold on
plot(BT_N_raw,'k'); plot(BT_Nf,'b'); plot(idxFilt,BT_N_raw(idxFilt),'ro')
title(['BT_N | filtered = ' num2str(pctFilteredN,'%.2f') '%']); grid on

subplot(3,1,3); hold on
plot(stats.AngleDiff,'k'); yline(stats.AngleThreshold,'r--'); plot(idxFilt,stats.AngleDiff(idxFilt),'ro')
title('Angle difference and threshold'); grid on
xlabel('Sample')

% note that 40-55% filtering is normal and usually a good dataset. Above
% 75% is probabkly wak transport or very noisy data. above 85% I dont count
% as valid.... the error could be anything: compass interference, GPS, etc.
% These data should be marked and the metadata should be checked, not
% immideiately removed. Somtimes it is only the compass an issue. 

