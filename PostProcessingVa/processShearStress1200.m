function ShearStressOut = processShearStress1200(dataADCP1200, ADCPin, WaterProfileOut, varargin)
%PROCESSSHEARSTRESS1200 Calculate shear velocity, shear stress, critical
% shear and excess shear for the 1200-kHz ADCP workflow.
%
% PURPOSE
% -------
% This function contains ONLY shear-related calculations. It does not
% calculate bedload concentration (Cb), active-layer thickness or bedload
% transport. Those calculations belong in processBedloadTrate.
%
% GRAIN-SIZE AND UNIT CONVENTION
% ------------------------------
% D50mobile is the representative mobile-grain D50 and is interpreted as
% millimetres [mm]. It is calculated row-by-row as the mean of D50_14 and
% D50_05, with fallback to D50 and finally D_50grabM when needed. It is
% converted exactly once:
%
%   D50critical_m = D50mobile ./ 1000;        % [m]
%
% D50mobile is used consistently for all Shields stresses, critical-shear
% calculations and excess-shear comparisons. D_50grabM is retained as an
% input/output variable but does not control the critical threshold.
%
% D65 CONVENTION (Manning-WC formula)
% -----------------------------------
% shearMannWC uses D65 in METRES. A measured D65 (read in mm) is used when
% present; otherwise D65mobile_mm is estimated from the grain-size
% distribution by log-linear interpolation of D against the standard-normal
% quantile through D50, D84 and D90 (all mm), then divided by 1000 once.
% D65source_1200 records whether D65 was measured or estimated.
%
% SHEAR-VELOCITY METHODS AND PRESERVED VARIABLE NAMES
% ---------------------------------------------------
% 1) Log-law fit to the near-bed velocity profile
%      UsheaAll_1200       shear velocity [m/s]
%      Tau_All1200         Shields stress [-]
%      tau_All1200_Pa      physical shear stress [Pa]
%      The fit excludes near-bed cells inside the bed side-lobe zone
%      (z <= (1-cos(BeamAngleDeg))*Depth) and, when enough cells remain,
%      is restricted to the log layer (z <= LogLawMaxZFraction*Depth).
%      nLogLawCells_1200 records the number of cells used per row. The raw
%      profiles z1_1200/Vmod1_1200 are never modified.
%
% 2) Reynolds-stress estimates from velocity fluctuations
%      U_shearRE1_1200, U_shearRE2_1200                [m/s]
%      Tau_AllRE1_1200, Tau_AllRE2_1200                [-]
%      tau_AllRE1_1200_Pa, tau_AllRE2_1200_Pa          [Pa]
%
% 3) Hydraulic/Chezy estimate
%      u_shearH1200 = sqrt(g)*U/C                         [m/s]
%      Tau_All1200_H                                     [-]
%      tau_All1200_H_Pa                                 [Pa]
%
% 4) Depth-average Van Rijn estimates
%      [Ushear1_1200, Ushear2_1200, Chezzy_1200] =
%          shearVanRijn(U,Depth,a,D0,D90)
%      aVanRijn may be one scalar or one value per input row.
%      Ushear1 (grain) needs only ks (grain roughness, m). Ushear2 (form)
%      additionally needs D0=HEightDUNE1200 (m). A missing D0 leaves
%      Ushear2_1200 NaN but no longer invalidates Ushear1_1200.
%      D0 is derived from DepthRNGout.stdSpatial_1200 (the mean spatial
%      depth difference among the four simultaneous beam depths); fallbacks
%      use stdSpatialALL_1200 or reconstruct it from RNG11..RNG44. Grain
%      size is never substituted for D0.
%      Tau_All1200_1, Tau_All1200_2                     [-]
%      tau_All1200_1_Pa, tau_All1200_2_Pa              [Pa]
%
% 5) Manning/WC estimate
%      tau_Mann1200       physical shear stress [Pa]
%      uShear_Mann1200    shear velocity [m/s]
%      tau_Mann_star1200  Shields stress [-]
%
% CRITICAL-SHEAR METHODS AND VARIABLES
% ------------------------------------
% CriticalShieldsMethod may be:
%   'vanrijn'  : Van Rijn piecewise fit based on Dstar.
%   'brownlie' : Brownlie (1981) fit based on particle Reynolds number Rep.
%   'soulsby'  : Soulsby & Whitehouse (1997) fit based on Dstar.
%   'both'     : calculates both Brownlie and Soulsby; Brownlie is used as
%                the primary ThetaCr_1200 value, while both are exported.
%   'constant' or 'fixed' : uses ThetaCrConstant.
%
% Dimensionless quantities and explicit fits:
%   Rep   = sqrt(R*g*D50)*D50/ni
%   Dstar = D50*(R*g/ni^2)^(1/3), with Rep = Dstar^(3/2)
%   Brownlie:
%     theta_c = 0.22*Rep^(-0.6) + 0.06*exp(-17.77*Rep^(-0.6))
%   Soulsby-Whitehouse:
%     theta_c = 0.30/(1+1.2*Dstar) +
%               0.055*(1-exp(-0.020*Dstar))
%   Dimensional critical stress and velocity:
%     Tucr_1200 = ThetaCr_1200*(ros-ro)*g*D50
%     ucr_1200  = sqrt(Tucr_1200/ro)
%
%   Rep_1200                  particle Reynolds number [-]
%   Dstar_1200                dimensionless particle diameter [-]
%   ThetaCrVanRijn_1200       Van Rijn critical Shields value [-]
%   ThetaCrBrownlie_1200      Brownlie critical Shields value [-]
%   ThetaCrSoulsby_1200       Soulsby-Whitehouse critical Shields value [-]
%   ThetaCr_1200              selected/primary critical Shields value [-]
%   Tucr_1200                 selected critical physical stress [Pa]
%   Tau_Crn                   unchanged alias of Tucr_1200 [Pa]
%   ucr_1200                  selected critical shear velocity [m/s]
%
% Method-specific Tucr aliases are retained for processBedloadTrate:
%   TucrLog_1200, TucrRE1_1200, TucrRE2_1200, TucrH_1200,
%   Tucr1_1200, Tucr2_1200, TucrMann_1200.
%
% OPTIONAL GRAIN/SKIN-FRICTION PARTITION
% --------------------------------------
% When ApplySkinFrictionPartition=true, the function uses depth-averaged
% velocity and hydraulic radius to estimate grain stress:
%   C_grain = 18*log10(12*Rh/(3*D90))
%   tau_skin = ro*g*(U/C_grain)^2
%   theta_skin = tau_skin/((ros-ro)*g*D50)
% By default, U=VmodAll_1200, Rh=Depth_1200 and D90=2.5*D50. These can be
% overridden using USkin, RhSkin and D90Skin (D90Skin is in metres).
% Additional outputs are C_grain_1200, tau_skin_1200, theta_skin_1200,
% uShear_skin1200, mobility_ratio_skin_1200, moving_skin_1200 and
% TvSkin_1200. Existing variable names are not changed.
%
% EXCESS-SHEAR VARIABLES
% ----------------------
% For every method, excess shear is calculated only from comparable,
% dimensionless quantities:
%
%   Tv = (Theta_method - ThetaCr_1200)./ThetaCr_1200
%
% Preserved/output names:
%   Tv_1200, TvRE1_1200, TvRE2_1200, TvH_1200,
%   Tv1_1200, Tv2_1200, TvMann_1200.
%
% Note: grain/profile-based stresses (log-law, Reynolds, Ushear1, skin)
% represent the stress acting on the grains. On large bedform-dominated
% rivers it is physically normal for these to give theta'/thetaCr < 1 while
% bedload transport still occurs, because (a) transport is driven by the
% total stress those methods do not see and (b) motion is intermittent
% below the continuous-motion threshold. Such rows are not errors; report
% them as partial transport rather than discarding them.
%
% The equivalent velocity ratio and mobility flags are also returned:
%   Uratio... = uStar./ucr_1200
%   Mobile... = Tv > 0
%
% PLOTS
% -----
%   'PlotTrateShear',true
%       TrateMean and MaxTrate versus every shear-velocity method.
%
%   'PlotCriticalComparison',true
%       uStar versus ucr_1200 with a 1:1 line. Every subplot has identical
%       x- and y-axis limits.
%
%   'PlotExcessShear',true
%       TrateMean and MaxTrate versus excess shear. Every subplot has the
%       same x-axis limits and the same y-axis limits.
%
% RECOMMENDED CALL
% ----------------
% ShearStressOut = processShearStress1200( ...
%     dataADCP1200, PhSamplesin, WaterProfileOut, ...
%     'ni', ni, ...
%     'ks', ks, ...
%     'SO', SO, ...
%     'CZ', [], ...
%     'DepthRNGout', DepthRNGout, ...
%     'aVanRijn', 2, ...
%     'BeamAngleDeg', 20, ...
%     'CriticalShieldsMethod', 'vanrijn', ...
%     'ApplySkinFrictionPartition', true, ...
%     'PlotTrateShear', true, ...
%     'PlotCriticalComparison', true, ...
%     'PlotExcessShear', true, ...
%     'AssignToCaller', true, ...
%     'AssignToBase', true, ...
%     'Verbose', true);
%
% No files are saved.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'dataADCP1200',@(x) isempty(x) || isstruct(x));
addRequired(p,'ADCPin',@(x) isstruct(x));
addRequired(p,'WaterProfileOut',@(x) isstruct(x));

addParameter(p,'ro',1000,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'ros',2650,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'g',9.81,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'R',1.65,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'kappa',0.4,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);

% IMPORTANT: explicit kinematic-viscosity input [m^2/s].
addParameter(p,'ni',[],@(x) isempty(x) || isnumeric(x));

addParameter(p,'SO',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'SOIsPermille',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'CZ',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'distanceSonTek',0,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'ks',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'HEightDUNE1200',[],@(x) isempty(x) || isnumeric(x));
% Output structure from processDepthsRNG. Used to derive D0/HEightDUNE1200
% from stdSpatial_1200: the mean spatial depth difference among the four
% simultaneous beam depths. Temporal beam standard deviations are not used.
addParameter(p,'DepthRNGout',[],@(x) isempty(x) || isstruct(x));

addParameter(p,'x0',[0.05 0.001],@(x) isnumeric(x) && ~isempty(x));
addParameter(p,'REcoef',0.19,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'PlotFlag','noplot',@(x) ischar(x) || isstring(x));
addParameter(p,'aVanRijn',2,@(x) isnumeric(x) && ~isempty(x));

% Log-law fit conditioning. BeamAngleDeg sets the bed side-lobe exclusion
% zone (1-cos(BeamAngleDeg)) of the depth; LogLawMaxZFraction caps the fit
% to the log layer when at least LogLawMinCells cells remain after trimming.
addParameter(p,'BeamAngleDeg',20,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x<45);
addParameter(p,'LogLawMaxZFraction',0.2,@(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
addParameter(p,'LogLawMinCells',3,@(x) isnumeric(x) && isscalar(x) && x>=2);

addParameter(p,'CriticalShieldsMethod','vanrijn',@(x) ischar(x) || isstring(x));
addParameter(p,'ThetaCrConstant',0.047,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);

% Optional grain/skin-friction partition. D90Skin is in metres.
addParameter(p,'ApplySkinFrictionPartition',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'D90Skin',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'USkin',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'RhSkin',[],@(x) isempty(x) || isnumeric(x));

addParameter(p,'ApplyLogLaw',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'ApplyReynolds',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'ApplyChezy',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'ApplyDepthAverageVanRijn',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'ApplyManningWC',true,@(x) islogical(x) || isnumeric(x));


addParameter(p,'PlotTrateShear',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'PlotTrateShearMethods',[],@(x) isempty(x) || isnumeric(x) || iscellstr(x) || isstring(x));
addParameter(p,'PlotTrateShearUseLogX',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'PlotTrateShearUseLogY',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'PlotTrateShearMinN',3,@(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p,'PlotCriticalComparison',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'PlotExcessShear',false,@(x) islogical(x) || isnumeric(x));

addParameter(p,'AssignToCaller',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToBase',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'Verbose',true,@(x) islogical(x) || isnumeric(x));

parse(p,dataADCP1200,ADCPin,WaterProfileOut,varargin{:});
opt = p.Results;

if opt.ros <= opt.ro
    error('ros must be greater than ro.');
end
validCriticalMethods = ["vanrijn","van rijn","brownlie","soulsby", ...
    "soulsby-whitehouse","soulsby & whitehouse","both","constant","fixed"];
selectedCriticalMethod = lower(strtrim(string(opt.CriticalShieldsMethod)));
if ~any(selectedCriticalMethod == validCriticalMethods)
    error(['CriticalShieldsMethod must be vanrijn, brownlie, soulsby, ', ...
        'both, constant or fixed.']);
end
% critical_shear_river defines submerged specific gravity from densities.
R_used = (opt.ros-opt.ro)./opt.ro;
if isfinite(opt.R) && abs(opt.R-R_used) > 1e-10
    warning('processShearStress1200:RConsistency', ...
        ['Input R=%.6g differs from (ros-ro)/ro=%.6g. The density-derived ', ...
         'value is used for all Shields and critical-shear calculations.'], ...
        opt.R,R_used);
end

%% Checks and common inputs
if ~isfield(ADCPin,'ADCPfilename')
    error('ADCPin.ADCPfilename is missing.');
end
requiredProfileFields = {'z1_1200_','Vmod1_1200_','V1x_1_','V1y_1_','V1w_1_'};
for k = 1:numel(requiredProfileFields)
    if ~isfield(WaterProfileOut,requiredProfileFields{k})
        error('WaterProfileOut.%s is missing.',requiredProfileFields{k});
    end
end

ADCPfilename = string(ADCPin.ADCPfilename(:));
nFiles = numel(ADCPfilename);
InputRow = (1:nFiles)';

StartADCP = getVectorField(ADCPin,'StartADCP',zeros(nFiles,1),nFiles);
EndADCP = getVectorField(ADCPin,'EndADCP',zeros(nFiles,1),nFiles);
DistanceFromLeftBank = getVectorField(ADCPin,'DistanceFromLeftBank',nan(nFiles,1),nFiles);

TrateMean = getVectorFieldAny(ADCPin,{'TrateMean','MeanTrate','MTrateAll','MeanTrateAll'},nan(nFiles,1),nFiles);
MaxTrate = getVectorFieldAny(ADCPin,{'MaxTrate','TrateMax','MTrateMax','MaxTrateAll'},nan(nFiles,1),nFiles);

D_90grabM = getVectorFieldAny(ADCPin,{'D_90grabM','D90grabM','D90_grabM','D90grab'},nan(nFiles,1),nFiles);
D_50grabM = getVectorFieldAny(ADCPin,{'D_50grabM','D50grabM','D50_grabM','D50grab'},nan(nFiles,1),nFiles);
D50_14 = getVectorFieldAny(ADCPin,{'D50_14','D50_14mm','D50mobile14'},nan(nFiles,1),nFiles);
D50_05 = getVectorFieldAny(ADCPin,{'D50_05','D50_05mm','D50mobile05'},nan(nFiles,1),nFiles);
D50 = getVectorFieldAny(ADCPin,{'D50','D_50'},nan(nFiles,1),nFiles);
D65 = getVectorFieldAny(ADCPin,{'D65','D_65','D65_14'},nan(nFiles,1),nFiles);
D84 = getVectorFieldAny(ADCPin,{'D84','D_84','D84_14'},nan(nFiles,1),nFiles);   % [mm]
FsS = getVectorFieldAny(ADCPin,{'FsS','Fs','sandFraction','SandFraction'},nan(nFiles,1),nFiles);
D_g_std_14 = getVectorFieldAny(ADCPin,{'D_g_std_14','Dstd_14','D_g_std','Dstd'},nan(nFiles,1),nFiles);

% Representative mobile D50 [mm]. This is the D50 used for every critical
% Shields, dimensional critical-stress and excess-shear calculation.
D50mobile = nan(nFiles,1);
for iD50 = 1:nFiles
    mobileCandidates = [D50_14(iD50),D50_05(iD50)];
    mobileCandidates = mobileCandidates(isfinite(mobileCandidates) & mobileCandidates>0);
    if ~isempty(mobileCandidates)
        D50mobile(iD50) = mean(mobileCandidates);
    end
    if ~isfinite(D50mobile(iD50)) || D50mobile(iD50)<=0
        D50mobile(iD50) = D50(iD50);
    end
    if ~isfinite(D50mobile(iD50)) || D50mobile(iD50)<=0
        D50mobile(iD50) = D_50grabM(iD50);
    end
end

D50critical_mm = D50mobile;
D50critical_m = D50critical_mm./1000;

if all(~isfinite(D50critical_mm) | D50critical_mm<=0)
    error(['No valid D50mobile values were found. Provide D50_14/D50_05, ', ...
        'D50, or D_50grabM in millimetres.']);
end

% Representative D65 [mm] for shearMannWC. Prefer a measured D65; otherwise
% estimate it from the grain-size distribution (log-linear in the standard-
% normal quantile through D50, D84 and D90). Converted to metres once at the
% call site, matching the D65m [m] convention in shearMannWC.
D65mobile_mm = nan(nFiles,1);
D65source_1200 = strings(nFiles,1);
D65source_1200(:) = "missing";
for iD65 = 1:nFiles
    if isfinite(D65(iD65)) && D65(iD65)>0
        D65mobile_mm(iD65) = D65(iD65);
        D65source_1200(iD65) = "input_D65";
    else
        [D65mobile_mm(iD65),D65source_1200(iD65)] = estimateD65FromPSD( ...
            D50critical_mm(iD65), D84(iD65), D_90grabM(iD65));
    end
end

SO_1200 = prepareOptionalVector(opt.SO,nFiles);
if all(isnan(SO_1200))
    SO_1200 = getVectorFieldAny(ADCPin,{'SO','S0','Slope','slope','WaterSlope','waterSlope'},nan(nFiles,1),nFiles);
end
CZ_input = prepareOptionalVector(opt.CZ,nFiles);

% Slope-magnitude sanity check (diagnostic only; does not alter any result).
% After applying SOIsPermille, the effective gradient S0 should fall in a
% physically plausible river range. A median far outside ~[1e-5, 1e-1]
% usually signals a unit mismatch between permille and m/m.
if logical(opt.SOIsPermille), S0check = SO_1200./1000; else, S0check = SO_1200; end
S0valid = S0check(isfinite(S0check) & S0check>0);
if ~isempty(S0valid)
    medS0 = median(S0valid);
    if medS0 < 1e-5 || medS0 > 1e-1
        warning('processShearStress1200:slopeMagnitude', ...
            ['Median effective slope S0=%.3g after SOIsPermille=%d is outside ', ...
             'the typical river range [1e-5, 1e-1]. Check whether SO is in ', ...
             'permille or m/m.'], medS0, logical(opt.SOIsPermille));
    end
end

Depth_1200 = getVectorFieldAny(WaterProfileOut,{'Depth_1200','Depth1200'},nan(nFiles,1),nFiles);
VmodAll_1200 = getVectorFieldAny(WaterProfileOut,{'VmodAll_1200','VmodAll1200'},nan(nFiles,1),nFiles);
WT_12 = getVectorFieldAny(WaterProfileOut,{'WT_12','WT12','Temp1'},nan(nFiles,1),nFiles);

if all(isnan(VmodAll_1200)) && isfield(WaterProfileOut,'Vmod_1200_')
    for i = 1:min(nFiles,numel(WaterProfileOut.Vmod_1200_))
        VmodAll_1200(i) = safeMean(WaterProfileOut.Vmod_1200_{i});
    end
end

% Skin-friction inputs. Missing values are filled row-by-row from the
% standard workflow inputs, following critical_shear_river conventions.
D90Skin_1200 = prepareOptionalVector(opt.D90Skin,nFiles); % [m]
USkin_1200 = prepareOptionalVector(opt.USkin,nFiles);     % [m/s]
RhSkin_1200 = prepareOptionalVector(opt.RhSkin,nFiles);   % [m]
for iSkin = 1:nFiles
    if ~isfinite(D90Skin_1200(iSkin)) || D90Skin_1200(iSkin)<=0
        D90Skin_1200(iSkin) = 2.5.*D50critical_m(iSkin);
    end
    if ~isfinite(USkin_1200(iSkin))
        USkin_1200(iSkin) = VmodAll_1200(iSkin);
    end
    if ~isfinite(RhSkin_1200(iSkin)) || RhSkin_1200(iSkin)<=0
        RhSkin_1200(iSkin) = Depth_1200(iSkin);
    end
end

ks = prepareOptionalVector(opt.ks,nFiles);
if all(isnan(ks))
    ks = getVectorFieldAny(ADCPin,{'ks','Ks','k_s','roughnessHeight'},nan(nFiles,1),nFiles);
end
% D_90grabM is in millimetres; shearVanRijn requires grain roughness in m.
missingKs = ~isfinite(ks) | ks<=0;
ks(missingKs) = D_90grabM(missingKs)./1000;

% D0 / dune-height handling. Preserve the variable name HEightDUNE1200.
% Priority: explicit input -> DepthRNGout spatial beam-depth difference ->
% ADCPin dune-height field -> NaN. Grain size is never used as dune height.
HEightDUNE1200 = prepareOptionalVector(opt.HEightDUNE1200,nFiles);
HEightDUNEsource_1200 = strings(nFiles,1);
HEightDUNEsource_1200(:) = "missing";
validD0 = isfinite(HEightDUNE1200) & HEightDUNE1200>0;
HEightDUNEsource_1200(validD0) = "input_HEightDUNE1200";

[D0FromDepthRNG_1200,D0BeamStd_1200,D0DepthRNGSource_1200] = ...
    deriveDuneHeightFromDepthRNG(opt.DepthRNGout,nFiles);
fillD0 = ~validD0 & isfinite(D0FromDepthRNG_1200) & D0FromDepthRNG_1200>0;
HEightDUNE1200(fillD0) = D0FromDepthRNG_1200(fillD0);
HEightDUNEsource_1200(fillD0) = D0DepthRNGSource_1200(fillD0);

D0FromADCPin = getVectorFieldAny(ADCPin, ...
    {'HEightDUNE1200','HeightDUNE1200','DuneHeight1200','D0'}, ...
    nan(nFiles,1),nFiles);
fillD0 = (~isfinite(HEightDUNE1200) | HEightDUNE1200<=0) & ...
    isfinite(D0FromADCPin) & D0FromADCPin>0;
HEightDUNE1200(fillD0) = D0FromADCPin(fillD0);
HEightDUNEsource_1200(fillD0) = "ADCPin_D0";

aVanRijn_1200 = prepareVector(opt.aVanRijn,nFiles);

ni = prepareOptionalVector(opt.ni,nFiles);
niInput = ni;
niSource_1200 = strings(nFiles,1);
niSource_1200(:) = "missing";
niSource_1200(isfinite(ni) & ni>0) = "input_ni";
Temp1 = WT_12;

if isempty(dataADCP1200), dataADCP1200 = struct([]); end
FileFound = false(nFiles,1);
MatchedFile = strings(nFiles,1);
Message = strings(nFiles,1);

%% Preallocate preserved and new shear variables
D90shear = D_90grabM;
D50bottom = D_50grabM./1000;
missingBottom = ~isfinite(D50bottom) | D50bottom<=0;
D50bottom(missingBottom) = D50critical_m(missingBottom);

Rep_1200 = nan(nFiles,1);
Dstar_1200 = nan(nFiles,1);
ThetaCrVanRijn_1200 = nan(nFiles,1);
ThetaCrBrownlie_1200 = nan(nFiles,1);
ThetaCrSoulsby_1200 = nan(nFiles,1);
ThetaCr_1200 = nan(nFiles,1);
Tucr_1200 = nan(nFiles,1);
Tau_Crn = nan(nFiles,1);
ucr_1200 = nan(nFiles,1);

C_grain_1200 = nan(nFiles,1);
tau_skin_1200 = nan(nFiles,1);
theta_skin_1200 = nan(nFiles,1);
uShear_skin1200 = nan(nFiles,1);
mobility_ratio_skin_1200 = nan(nFiles,1);
moving_skin_1200 = false(nFiles,1);
TvSkin_1200 = nan(nFiles,1);
UratioSkin_1200 = nan(nFiles,1);
MobileSkin_1200 = false(nFiles,1);

TucrLog_1200 = nan(nFiles,1);
TucrRE1_1200 = nan(nFiles,1);
TucrRE2_1200 = nan(nFiles,1);
TucrH_1200 = nan(nFiles,1);
Tucr1_1200 = nan(nFiles,1);
Tucr2_1200 = nan(nFiles,1);
TucrMann_1200 = nan(nFiles,1);

x_1200 = cell(nFiles,1);
resnorm_1200 = cell(nFiles,1);
residual_1200 = cell(nFiles,1);
exitflag__1200 = cell(nFiles,1);
output_1200 = cell(nFiles,1);
fun_1200 = cell(nFiles,1);
nLogLawCells_1200 = nan(nFiles,1);   % cells used in each log-law fit (audit)

UsheaAll_1200 = nan(nFiles,1);
Tau_All1200 = nan(nFiles,1);
tau_All1200_Pa = nan(nFiles,1);

U_shearRE1_1200 = nan(nFiles,1);
Tau_AllRE1_1200 = nan(nFiles,1);
tau_AllRE1_1200_Pa = nan(nFiles,1);
U_shearRE2_1200 = nan(nFiles,1);
Tau_AllRE2_1200 = nan(nFiles,1);
tau_AllRE2_1200_Pa = nan(nFiles,1);

uw_mean1200 = cell(nFiles,1);
vw_mean1200 = cell(nFiles,1);
TKE_prof1200 = cell(nFiles,1);
tau_turb_prof1200 = cell(nFiles,1);

CZ_1200 = nan(nFiles,1);
CZsource_1200 = strings(nFiles,1);
u_shearH1200 = nan(nFiles,1);
Tau_All1200_H = nan(nFiles,1);
tau_All1200_H_Pa = nan(nFiles,1);

Ushear1_1200 = nan(nFiles,1);
Ushear2_1200 = nan(nFiles,1);
Chezzy_1200 = nan(nFiles,1);
Tau_All1200_1 = nan(nFiles,1);
Tau_All1200_2 = nan(nFiles,1);
tau_All1200_1_Pa = nan(nFiles,1);
tau_All1200_2_Pa = nan(nFiles,1);

tau_Mann1200 = nan(nFiles,1);
tau_Mann_star1200 = nan(nFiles,1);
tau_rm_star1200 = nan(nFiles,1);
tau_rm_starGaumann1200 = nan(nFiles,1);
tau_rm1200 = nan(nFiles,1);
uShear_Mann1200 = nan(nFiles,1);

Tv_1200 = nan(nFiles,1);
TvRE1_1200 = nan(nFiles,1);
TvRE2_1200 = nan(nFiles,1);
TvH_1200 = nan(nFiles,1);
Tv1_1200 = nan(nFiles,1);
Tv2_1200 = nan(nFiles,1);
TvMann_1200 = nan(nFiles,1);

Uratio_1200 = nan(nFiles,1);
UratioRE1_1200 = nan(nFiles,1);
UratioRE2_1200 = nan(nFiles,1);
UratioH_1200 = nan(nFiles,1);
Uratio1_1200 = nan(nFiles,1);
Uratio2_1200 = nan(nFiles,1);
UratioMann_1200 = nan(nFiles,1);

Mobile_1200 = false(nFiles,1);
MobileRE1_1200 = false(nFiles,1);
MobileRE2_1200 = false(nFiles,1);
MobileH_1200 = false(nFiles,1);
Mobile1_1200 = false(nFiles,1);
Mobile2_1200 = false(nFiles,1);
MobileMann_1200 = false(nFiles,1);

Fr_1200 = nan(nFiles,1);
Tausk_1200 = nan(nFiles,1);

z1_1200_all = cell(nFiles,1);
Vmod1_1200_all = cell(nFiles,1);
V1x_1_all = cell(nFiles,1);
V1y_1_all = cell(nFiles,1);
V1w_1_all = cell(nFiles,1);

%% Per-row structure
ShearStressAll = repmat(struct( ...
    'InputRow',[],'ADCPfilename',"",'matchedFile',"",'fileFound',false, ...
    'D90shear',[],'D50bottom',[],'D50mobile',[],'D50critical_mm',[], ...
    'D50critical_m',[],'ks',[],'HEightDUNE1200',[], ...
    'HEightDUNEsource_1200',"",'D0FromDepthRNG_1200',[], ...
    'D0BeamStd_1200',[],'ni',[],'Temp1',[], ...
    'Rep_1200',[],'Dstar_1200',[],'ThetaCrVanRijn_1200',[], ...
    'ThetaCrBrownlie_1200',[],'ThetaCrSoulsby_1200',[], ...
    'ThetaCr_1200',[],'Tucr_1200',[],'ucr_1200',[], ...
    'D90Skin_1200',[],'USkin_1200',[],'RhSkin_1200',[], ...
    'C_grain_1200',[],'tau_skin_1200',[],'theta_skin_1200',[], ...
    'uShear_skin1200',[],'mobility_ratio_skin_1200',[], ...
    'moving_skin_1200',false,'TvSkin_1200',[],'UratioSkin_1200',[], ...
    'z1_1200',[],'Vmod1_1200',[],'V1x_1',[],'V1y_1',[],'V1w_1',[], ...
    'x_1200',[],'nLogLawCells_1200',[],'UsheaAll_1200',[],'Tau_All1200',[],'tau_All1200_Pa',[], ...
    'U_shearRE1_1200',[],'Tau_AllRE1_1200',[],'tau_AllRE1_1200_Pa',[], ...
    'U_shearRE2_1200',[],'Tau_AllRE2_1200',[],'tau_AllRE2_1200_Pa',[], ...
    'TKE_prof1200',[],'tau_turb_prof1200',[],'uw_mean1200',[],'vw_mean1200',[], ...
    'CZ_1200',[],'CZsource_1200',"",'u_shearH1200',[], ...
    'Tau_All1200_H',[],'tau_All1200_H_Pa',[], ...
    'Ushear1_1200',[],'Ushear2_1200',[],'Chezzy_1200',[], ...
    'Tau_All1200_1',[],'Tau_All1200_2',[], ...
    'tau_All1200_1_Pa',[],'tau_All1200_2_Pa',[], ...
    'tau_Mann1200',[],'tau_Mann_star1200',[],'tau_rm_star1200',[], ...
    'tau_rm_starGaumann1200',[],'tau_rm1200',[],'uShear_Mann1200',[], ...
    'Tv_1200',[],'TvRE1_1200',[],'TvRE2_1200',[],'TvH_1200',[], ...
    'Tv1_1200',[],'Tv2_1200',[],'TvMann_1200',[], ...
    'Fr_1200',[],'Tausk_1200',[],'message',""),nFiles,1);

%% Main loop
for i = 1:nFiles
    ShearStressAll(i).InputRow = InputRow(i);
    ShearStressAll(i).ADCPfilename = ADCPfilename(i);

    try
        if ~isfinite(D50critical_mm(i)) || D50critical_mm(i)<=0
            error('Invalid D50mobile at row %d. Expected a positive value in mm.',i);
        end

        matchedFile = "";
        if ~isempty(dataADCP1200)
            matchedFile = matchADCPfile(dataADCP1200,ADCPfilename(i));
        end
        MatchedFile(i) = matchedFile;
        FileFound(i) = matchedFile~="";

        z1_1200 = getCellField(WaterProfileOut,'z1_1200_',i);
        Vmod1_1200 = getCellField(WaterProfileOut,'Vmod1_1200_',i);
        V1x_1 = getCellField(WaterProfileOut,'V1x_1_',i);
        V1y_1 = getCellField(WaterProfileOut,'V1y_1_',i);
        V1w_1 = getCellField(WaterProfileOut,'V1w_1_',i);

        if isempty(z1_1200) || isempty(Vmod1_1200) || isempty(V1x_1) || isempty(V1y_1) || isempty(V1w_1)
            error('Required water-profile data are missing for row %d.',i);
        end

        z1_1200 = double(z1_1200(:));
        Vmod1_1200 = double(Vmod1_1200(:));
        V1x_1 = double(V1x_1);
        V1y_1 = double(V1y_1);
        V1w_1 = double(V1w_1);

        z1_1200_all{i} = z1_1200;
        Vmod1_1200_all{i} = Vmod1_1200;
        V1x_1_all{i} = V1x_1;
        V1y_1_all{i} = V1y_1;
        V1w_1_all{i} = V1w_1;

        if ~isfinite(ks(i)) || ks(i)<=0
            if isfinite(D90shear(i)) && D90shear(i)>0
                ks(i)=D90shear(i)./1000;
            end
        end
        % Do not replace missing dune height with grain size. A missing D0
        % intentionally leaves Ushear2_1200 as NaN.

        Temp1(i) = WT_12(i);
        if ~isfinite(Temp1(i)) && matchedFile~=""
            try
                STemp = load(matchedFile,'Data');
                if isfield(STemp,'Data') && isfield(STemp.Data,'WT')
                    Temp1(i)=safeMean(STemp.Data.WT);
                end
            catch
            end
        end
        if ~isfinite(ni(i)) || ni(i)<=0
            ni(i)=callWaterViscosity(Temp1(i));
            niSource_1200(i)="calculated_from_temperature";
        end
        if ~isfinite(ni(i)) || ni(i)<=0
            error('Invalid ni at row %d. Expected kinematic viscosity in m^2/s.',i);
        end

        [Rep_1200(i),Dstar_1200(i),ThetaCrVanRijn_1200(i), ...
            ThetaCrBrownlie_1200(i),ThetaCrSoulsby_1200(i), ...
            ThetaCr_1200(i),Tucr_1200(i),ucr_1200(i)] = ...
            calculateCriticalShear(D50critical_m(i),ni(i),R_used,opt.g,opt.ro, ...
                opt.CriticalShieldsMethod,opt.ThetaCrConstant);
        Tau_Crn(i)=Tucr_1200(i);

        if logical(opt.ApplySkinFrictionPartition)
            [C_grain_1200(i),tau_skin_1200(i),theta_skin_1200(i), ...
                uShear_skin1200(i),mobility_ratio_skin_1200(i), ...
                moving_skin_1200(i)] = calculateSkinFriction( ...
                D50critical_m(i),D90Skin_1200(i),USkin_1200(i), ...
                RhSkin_1200(i),opt.ro,R_used,opt.g,ThetaCr_1200(i));
            TvSkin_1200(i)=excessShear(theta_skin_1200(i),ThetaCr_1200(i));
            UratioSkin_1200(i)=safeRatio(uShear_skin1200(i),ucr_1200(i));
            MobileSkin_1200(i)=moving_skin_1200(i);
        end
        TucrLog_1200(i)=Tucr_1200(i);
        TucrRE1_1200(i)=Tucr_1200(i);
        TucrRE2_1200(i)=Tucr_1200(i);
        TucrH_1200(i)=Tucr_1200(i);
        Tucr1_1200(i)=Tucr_1200(i);
        Tucr2_1200(i)=Tucr_1200(i);
        TucrMann_1200(i)=Tucr_1200(i);

        if logical(opt.ApplyLogLaw)
            % Side-lobe / near-bed contamination guard for the log-law fit
            % only. z1_1200 is height above the bed, so cells nearest the
            % bed (small z) within (1-cos(beamAngle)) of the depth are
            % contaminated by the bed side-lobe return; when enough cells
            % remain, the fit is restricted to the log layer
            % (z <= LogLawMaxZFraction*Depth). Raw profiles are preserved.
            zFit=z1_1200; vFit=Vmod1_1200;
            keepFit=isfinite(zFit)&isfinite(vFit);
            if isfinite(Depth_1200(i)) && Depth_1200(i)>0
                zSideLobe=(1-cosd(opt.BeamAngleDeg)).*Depth_1200(i);
                keepFit=keepFit & zFit>zSideLobe;
                inLog=keepFit & zFit<=opt.LogLawMaxZFraction.*Depth_1200(i);
                if sum(inLog)>=opt.LogLawMinCells
                    keepFit=inLog;
                end
            end
            if sum(keepFit)>=opt.LogLawMinCells
                zFit=zFit(keepFit); vFit=vFit(keepFit);
            end
            nLogLawCells_1200(i)=numel(zFit);
            [xTmp,resTmp,residTmp,exitTmp,outTmp,funTmp] = ...
                callCalculateShearVel(zFit,vFit',opt.kappa,ks(i),opt.x0,opt.PlotFlag);
            x_1200{i}=xTmp; resnorm_1200{i}=resTmp; residual_1200{i}=residTmp;
            exitflag__1200{i}=exitTmp; output_1200{i}=outTmp; fun_1200{i}=funTmp;
            if ~isempty(xTmp), UsheaAll_1200(i)=xTmp(1); end
            [Tau_All1200(i),tau_All1200_Pa(i)] = ...
                stressFromUstar(UsheaAll_1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
        end

        if logical(opt.ApplyReynolds)
            [uw_mean,vw_mean,x1,x2,~,~,TKE1200,tau_turb1200] = ...
                callCalculateREShearVelXYcord(V1x_1,V1y_1,V1w_1,z1_1200, ...
                    Depth_1200(i),opt.x0,opt.REcoef,opt.PlotFlag);
            if ~isempty(x1), U_shearRE1_1200(i)=x1(1); end
            if ~isempty(x2), U_shearRE2_1200(i)=x2(1); end
            [Tau_AllRE1_1200(i),tau_AllRE1_1200_Pa(i)] = ...
                stressFromUstar(U_shearRE1_1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
            [Tau_AllRE2_1200(i),tau_AllRE2_1200_Pa(i)] = ...
                stressFromUstar(U_shearRE2_1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
            TKE_prof1200{i}=nanmeanLocal(TKE1200,2);
            tau_turb_prof1200{i}=nanmeanLocal(tau_turb1200,2);
            uw_mean1200{i}=nanmeanLocal(uw_mean,2);
            vw_mean1200{i}=nanmeanLocal(vw_mean,2);
        end

        if logical(opt.ApplyChezy)
            if isfinite(CZ_input(i)) && CZ_input(i)>0
                CZ_1200(i)=CZ_input(i); CZsource_1200(i)="input_CZ";
            else
                [CZ_1200(i),CZsource_1200(i)] = calculateChezyFromSlope( ...
                    VmodAll_1200(i),Depth_1200(i),SO_1200(i),logical(opt.SOIsPermille));
            end
            if isfinite(CZ_1200(i)) && CZ_1200(i)>0
                % Correct Chezy conversion: u* = sqrt(g)*U/C.
                u_shearH1200(i)=sqrt(opt.g).*VmodAll_1200(i)./CZ_1200(i);
                [Tau_All1200_H(i),tau_All1200_H_Pa(i)] = ...
                    stressFromUstar(u_shearH1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
            end
        end

        if logical(opt.ApplyDepthAverageVanRijn)
            [Ushear1_1200(i),Ushear2_1200(i),Chezzy_1200(i)] = ...
                shearVanRijn(VmodAll_1200(i),Depth_1200(i),aVanRijn_1200(i), ...
                    HEightDUNE1200(i),ks(i));
            [Tau_All1200_1(i),tau_All1200_1_Pa(i)] = ...
                stressFromUstar(Ushear1_1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
            [Tau_All1200_2(i),tau_All1200_2_Pa(i)] = ...
                stressFromUstar(Ushear2_1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
        end

        if logical(opt.ApplyManningWC)
            % D65 in metres for shearMannWC. Prefer measured/estimated D65
            % (mm, computed once before the loop); the final safeguard keeps
            % D50critical_m, which is already in metres.
            if isfinite(D65mobile_mm(i)) && D65mobile_mm(i)>0
                D65_i = D65mobile_mm(i)./1000;     % mm -> m
            else
                D65_i = D50critical_m(i);          % already in metres
            end
            [tau_Mann1200(i),~,tau_rm_star1200(i),tau_rm_starGaumann1200(i), ...
                tau_rm1200(i),uShear_Mann1200(i)] = shearMannWC( ...
                    VmodAll_1200(i),SO_1200(i),D65_i,D50critical_m(i), ...
                    FsS(i),D_g_std_14(i),logical(opt.SOIsPermille),opt.ro);
            [tau_Mann_star1200(i),~] = stressFromUstar( ...
                uShear_Mann1200(i),D50critical_m(i),R_used,opt.g,opt.ro);
        end

        Tv_1200(i)=excessShear(Tau_All1200(i),ThetaCr_1200(i));
        TvRE1_1200(i)=excessShear(Tau_AllRE1_1200(i),ThetaCr_1200(i));
        TvRE2_1200(i)=excessShear(Tau_AllRE2_1200(i),ThetaCr_1200(i));
        TvH_1200(i)=excessShear(Tau_All1200_H(i),ThetaCr_1200(i));
        Tv1_1200(i)=excessShear(Tau_All1200_1(i),ThetaCr_1200(i));
        Tv2_1200(i)=excessShear(Tau_All1200_2(i),ThetaCr_1200(i));
        TvMann_1200(i)=excessShear(tau_Mann_star1200(i),ThetaCr_1200(i));

        Uratio_1200(i)=safeRatio(UsheaAll_1200(i),ucr_1200(i));
        UratioRE1_1200(i)=safeRatio(U_shearRE1_1200(i),ucr_1200(i));
        UratioRE2_1200(i)=safeRatio(U_shearRE2_1200(i),ucr_1200(i));
        UratioH_1200(i)=safeRatio(u_shearH1200(i),ucr_1200(i));
        Uratio1_1200(i)=safeRatio(Ushear1_1200(i),ucr_1200(i));
        Uratio2_1200(i)=safeRatio(Ushear2_1200(i),ucr_1200(i));
        UratioMann_1200(i)=safeRatio(uShear_Mann1200(i),ucr_1200(i));

        Mobile_1200(i)=isfinite(Tv_1200(i)) && Tv_1200(i)>0;
        MobileRE1_1200(i)=isfinite(TvRE1_1200(i)) && TvRE1_1200(i)>0;
        MobileRE2_1200(i)=isfinite(TvRE2_1200(i)) && TvRE2_1200(i)>0;
        MobileH_1200(i)=isfinite(TvH_1200(i)) && TvH_1200(i)>0;
        Mobile1_1200(i)=isfinite(Tv1_1200(i)) && Tv1_1200(i)>0;
        Mobile2_1200(i)=isfinite(Tv2_1200(i)) && Tv2_1200(i)>0;
        MobileMann_1200(i)=isfinite(TvMann_1200(i)) && TvMann_1200(i)>0;

        denomFr=opt.g.*Depth_1200(i)+opt.distanceSonTek;
        if isfinite(VmodAll_1200(i)) && isfinite(denomFr) && denomFr>0
            Fr_1200(i)=VmodAll_1200(i)./sqrt(denomFr);
        end
        if isfinite(Tau_All1200_H(i)) && isfinite(Fr_1200(i)) && ...
                Tau_All1200_H(i)>=0 && Fr_1200(i)>=0
            Tausk_1200(i)=0.05+0.7.*(Tau_All1200_H(i).*Fr_1200(i).^0.7).^0.8;
        end

        s=ShearStressAll(i);
        s.matchedFile=matchedFile; s.fileFound=FileFound(i);
        s.D90shear=D90shear(i); s.D50bottom=D50bottom(i); s.D50mobile=D50mobile(i);
        s.D50critical_mm=D50critical_mm(i); s.D50critical_m=D50critical_m(i);
        s.ks=ks(i); s.HEightDUNE1200=HEightDUNE1200(i);
        s.HEightDUNEsource_1200=HEightDUNEsource_1200(i);
        s.D0FromDepthRNG_1200=D0FromDepthRNG_1200(i);
        s.D0BeamStd_1200=D0BeamStd_1200(i,:);
        s.ni=ni(i); s.Temp1=Temp1(i);
        s.Rep_1200=Rep_1200(i); s.Dstar_1200=Dstar_1200(i);
        s.ThetaCrVanRijn_1200=ThetaCrVanRijn_1200(i);
        s.ThetaCrBrownlie_1200=ThetaCrBrownlie_1200(i);
        s.ThetaCrSoulsby_1200=ThetaCrSoulsby_1200(i);
        s.ThetaCr_1200=ThetaCr_1200(i); s.Tucr_1200=Tucr_1200(i); s.ucr_1200=ucr_1200(i);
        s.D90Skin_1200=D90Skin_1200(i); s.USkin_1200=USkin_1200(i); s.RhSkin_1200=RhSkin_1200(i);
        s.C_grain_1200=C_grain_1200(i); s.tau_skin_1200=tau_skin_1200(i);
        s.theta_skin_1200=theta_skin_1200(i); s.uShear_skin1200=uShear_skin1200(i);
        s.mobility_ratio_skin_1200=mobility_ratio_skin_1200(i);
        s.moving_skin_1200=moving_skin_1200(i); s.TvSkin_1200=TvSkin_1200(i);
        s.UratioSkin_1200=UratioSkin_1200(i);
        s.z1_1200=z1_1200; s.Vmod1_1200=Vmod1_1200;
        s.V1x_1=V1x_1; s.V1y_1=V1y_1; s.V1w_1=V1w_1;
        s.x_1200=x_1200{i}; s.nLogLawCells_1200=nLogLawCells_1200(i);
        s.UsheaAll_1200=UsheaAll_1200(i); s.Tau_All1200=Tau_All1200(i); s.tau_All1200_Pa=tau_All1200_Pa(i);
        s.U_shearRE1_1200=U_shearRE1_1200(i); s.Tau_AllRE1_1200=Tau_AllRE1_1200(i); s.tau_AllRE1_1200_Pa=tau_AllRE1_1200_Pa(i);
        s.U_shearRE2_1200=U_shearRE2_1200(i); s.Tau_AllRE2_1200=Tau_AllRE2_1200(i); s.tau_AllRE2_1200_Pa=tau_AllRE2_1200_Pa(i);
        s.TKE_prof1200=TKE_prof1200{i}; s.tau_turb_prof1200=tau_turb_prof1200{i}; s.uw_mean1200=uw_mean1200{i}; s.vw_mean1200=vw_mean1200{i};
        s.CZ_1200=CZ_1200(i); s.CZsource_1200=CZsource_1200(i);
        s.u_shearH1200=u_shearH1200(i); s.Tau_All1200_H=Tau_All1200_H(i); s.tau_All1200_H_Pa=tau_All1200_H_Pa(i);
        s.Ushear1_1200=Ushear1_1200(i); s.Ushear2_1200=Ushear2_1200(i); s.Chezzy_1200=Chezzy_1200(i);
        s.Tau_All1200_1=Tau_All1200_1(i); s.Tau_All1200_2=Tau_All1200_2(i);
        s.tau_All1200_1_Pa=tau_All1200_1_Pa(i); s.tau_All1200_2_Pa=tau_All1200_2_Pa(i);
        s.tau_Mann1200=tau_Mann1200(i); s.tau_Mann_star1200=tau_Mann_star1200(i);
        s.tau_rm_star1200=tau_rm_star1200(i); s.tau_rm_starGaumann1200=tau_rm_starGaumann1200(i);
        s.tau_rm1200=tau_rm1200(i); s.uShear_Mann1200=uShear_Mann1200(i);
        s.Tv_1200=Tv_1200(i); s.TvRE1_1200=TvRE1_1200(i); s.TvRE2_1200=TvRE2_1200(i);
        s.TvH_1200=TvH_1200(i); s.Tv1_1200=Tv1_1200(i); s.Tv2_1200=Tv2_1200(i); s.TvMann_1200=TvMann_1200(i);
        s.Fr_1200=Fr_1200(i); s.Tausk_1200=Tausk_1200(i);
        Message(i)="OK"; s.message=Message(i); ShearStressAll(i)=s;

    catch ME
        Message(i)="Error: "+string(ME.message);
        ShearStressAll(i).message=Message(i);
        warning('Shear-stress processing failed for row %d, file %s:\n%s',i,ADCPfilename(i),ME.message);
    end
end

% Physical-coherence diagnostic retained from critical_shear_river.
lowRepMask = isfinite(Rep_1200) & Rep_1200 < 1e2;
if any(lowRepMask)
    warning('processShearStress1200:lowRe', ...
        ['Rep < approximately 100 for %d row(s). Critical-Shields fits can ', ...
         'diverge more in the transitional/smooth regime.'],sum(lowRepMask));
end

% Plumbing diagnostic: explain NaN Van Rijn rows. Ushear1 needs a valid ks
% (grain roughness); Ushear2 additionally needs a valid D0/HEightDUNE1200.
nNoKs = sum(~(isfinite(ks)&ks>0));
nNoD0 = sum(~(isfinite(HEightDUNE1200)&HEightDUNE1200>0));
if nNoKs>0
    warning('processShearStress1200:ksMissing', ...
        ['ks (grain roughness) is missing/invalid for %d row(s); Van Rijn ', ...
         'u_*1 (Ushear1_1200) is NaN there. Provide ks or D_90grabM [mm].'],nNoKs);
end
if nNoD0>0
    warning('processShearStress1200:D0Missing', ...
        ['D0/HEightDUNE1200 (dune height) is missing/invalid for %d row(s); ', ...
         'Van Rijn u_*2 (Ushear2_1200) is NaN there by design. Provide ', ...
         'DepthRNGout.stdSpatial_1200 or HEightDUNE1200 [m].'],nNoD0);
end

%% Method matrices, diagnostics and plots
ShearVelocityNames1200=["UsheaAll_1200";"U_shearRE1_1200";"U_shearRE2_1200"; ...
    "u_shearH1200";"Ushear1_1200";"Ushear2_1200";"uShear_Mann1200"];
ShearVelocityLabels1200=["Log-law u_*";"Reynolds u_* 1";"Reynolds u_* 2"; ...
    "Hydraulic/Chezy u_*";"Van Rijn u_* 1";"Van Rijn u_* 2";"Manning-WC u_*"];
ShearVelocityAll1200=[UsheaAll_1200,U_shearRE1_1200,U_shearRE2_1200, ...
    u_shearH1200,Ushear1_1200,Ushear2_1200,uShear_Mann1200];

ExcessShearAll1200=[Tv_1200,TvRE1_1200,TvRE2_1200,TvH_1200,Tv1_1200,Tv2_1200,TvMann_1200];
UratioAll1200=[Uratio_1200,UratioRE1_1200,UratioRE2_1200,UratioH_1200,Uratio1_1200,Uratio2_1200,UratioMann_1200];
MobileAll1200=[Mobile_1200,MobileRE1_1200,MobileRE2_1200,MobileH_1200,Mobile1_1200,Mobile2_1200,MobileMann_1200];

% Add the optional skin-friction method to comparisons when it was computed.
if any(isfinite(uShear_skin1200))
    ShearVelocityNames1200(end+1,1)="uShear_skin1200";
    ShearVelocityLabels1200(end+1,1)="Van Rijn grain/skin stress u_*";
    ShearVelocityAll1200(:,end+1)=uShear_skin1200;
    ExcessShearAll1200(:,end+1)=TvSkin_1200;
    UratioAll1200(:,end+1)=UratioSkin_1200;
    MobileAll1200(:,end+1)=MobileSkin_1200;
end

TrateShearCorrTable=makeCorrelationTable(ShearVelocityNames1200,ShearVelocityLabels1200,ShearVelocityAll1200,TrateMean,MaxTrate,opt.PlotTrateShearMinN,'ShearVelocity');
ExcessShearCorrTable=makeCorrelationTable(ShearVelocityNames1200,ShearVelocityLabels1200,ExcessShearAll1200,TrateMean,MaxTrate,opt.PlotTrateShearMinN,'ExcessShear');
CriticalShearSummaryTable=makeCriticalSummary(ShearVelocityNames1200,ShearVelocityLabels1200,ShearVelocityAll1200,ucr_1200,ExcessShearAll1200);

if logical(opt.PlotTrateShear)
    plotMetricVsTransport(ShearVelocityNames1200,ShearVelocityLabels1200,ShearVelocityAll1200, ...
        TrateMean,MaxTrate,opt.PlotTrateShearMethods,logical(opt.PlotTrateShearUseLogX), ...
        logical(opt.PlotTrateShearUseLogY),opt.PlotTrateShearMinN, ...
        'Transport rate versus shear velocity','Shear velocity, u_* [m/s]');
end
if logical(opt.PlotCriticalComparison)
    plotCriticalComparison(ShearVelocityNames1200,ShearVelocityLabels1200,ShearVelocityAll1200,ucr_1200,ExcessShearAll1200);
end
if logical(opt.PlotExcessShear)
    plotMetricVsTransport(ShearVelocityNames1200,ShearVelocityLabels1200,ExcessShearAll1200, ...
        TrateMean,MaxTrate,[],false,logical(opt.PlotTrateShearUseLogY), ...
        opt.PlotTrateShearMinN,'Transport rate versus excess shear','Excess shear, T [-]');
end

%% Output table
ShearStressTable=table(InputRow,ADCPfilename,DistanceFromLeftBank,Depth_1200,VmodAll_1200, ...
    TrateMean,MaxTrate,SO_1200,CZ_input,CZ_1200,CZsource_1200,ks, ...
    HEightDUNE1200,HEightDUNEsource_1200,D0FromDepthRNG_1200, ...
    D0BeamStd_1200,ni,niSource_1200,Temp1, ...
    D_90grabM,D_50grabM,D50critical_m,D65mobile_mm,D65source_1200,Rep_1200,Dstar_1200, ...
    ThetaCrVanRijn_1200,ThetaCrBrownlie_1200,ThetaCrSoulsby_1200, ...
    ThetaCr_1200,Tucr_1200,ucr_1200,D90Skin_1200,USkin_1200,RhSkin_1200, ...
    C_grain_1200,tau_skin_1200,theta_skin_1200,uShear_skin1200, ...
    mobility_ratio_skin_1200,moving_skin_1200,TvSkin_1200, ...
    UsheaAll_1200,nLogLawCells_1200,Tau_All1200,tau_All1200_Pa,U_shearRE1_1200,Tau_AllRE1_1200, ...
    tau_AllRE1_1200_Pa,U_shearRE2_1200,Tau_AllRE2_1200,tau_AllRE2_1200_Pa, ...
    u_shearH1200,Tau_All1200_H,tau_All1200_H_Pa,Ushear1_1200,Ushear2_1200,Chezzy_1200, ...
    Tau_All1200_1,Tau_All1200_2,tau_All1200_1_Pa,tau_All1200_2_Pa, ...
    tau_Mann1200,tau_Mann_star1200,uShear_Mann1200,Tv_1200,TvRE1_1200,TvRE2_1200, ...
    TvH_1200,Tv1_1200,Tv2_1200,TvMann_1200,FileFound,MatchedFile,Message);

%% Last successful profiles
lastOK=find(Message=="OK",1,'last');
if isempty(lastOK)
    z1_1200=[]; Vmod1_1200=[]; V1x_1=[]; V1y_1=[]; V1w_1=[];
else
    z1_1200=z1_1200_all{lastOK}; Vmod1_1200=Vmod1_1200_all{lastOK};
    V1x_1=V1x_1_all{lastOK}; V1y_1=V1y_1_all{lastOK}; V1w_1=V1w_1_all{lastOK};
end

%% Build output structure
ShearStressOut=struct();
varsToExport={ ...
    'InputRow','ADCPfilename','StartADCP','EndADCP','DistanceFromLeftBank', ...
    'Depth_1200','VmodAll_1200','TrateMean','MaxTrate','SO_1200','CZ_input', ...
    'CZ_1200','CZsource_1200','D_90grabM','D_50grabM','D50_14','D50_05','D50', ...
    'D65','D84','D65mobile_mm','D65source_1200','FsS','D_g_std_14','D90shear','D50bottom','D50mobile', ...
    'R_used', ...
    'D50critical_mm','D50critical_m','ni','niInput','niSource_1200','Temp1','ks', ...
    'HEightDUNE1200','HEightDUNEsource_1200','D0FromDepthRNG_1200', ...
    'D0BeamStd_1200','D0DepthRNGSource_1200','aVanRijn_1200', ...
    'Rep_1200','Dstar_1200', ...
    'ThetaCrVanRijn_1200','ThetaCrBrownlie_1200','ThetaCrSoulsby_1200', ...
    'ThetaCr_1200','Tucr_1200','Tau_Crn','ucr_1200', ...
    'D90Skin_1200','USkin_1200','RhSkin_1200','C_grain_1200', ...
    'tau_skin_1200','theta_skin_1200','uShear_skin1200', ...
    'mobility_ratio_skin_1200','moving_skin_1200','TvSkin_1200', ...
    'UratioSkin_1200','MobileSkin_1200', ...
    'TucrLog_1200','TucrRE1_1200','TucrRE2_1200', ...
    'TucrH_1200','Tucr1_1200','Tucr2_1200','TucrMann_1200', ...
    'z1_1200','Vmod1_1200','V1x_1','V1y_1','V1w_1','z1_1200_all', ...
    'Vmod1_1200_all','V1x_1_all','V1y_1_all','V1w_1_all','x_1200','nLogLawCells_1200', ...
    'resnorm_1200','residual_1200','exitflag__1200','output_1200','fun_1200', ...
    'UsheaAll_1200','Tau_All1200','tau_All1200_Pa', ...
    'U_shearRE1_1200','Tau_AllRE1_1200','tau_AllRE1_1200_Pa', ...
    'U_shearRE2_1200','Tau_AllRE2_1200','tau_AllRE2_1200_Pa', ...
    'uw_mean1200','vw_mean1200','TKE_prof1200','tau_turb_prof1200', ...
    'u_shearH1200','Tau_All1200_H','tau_All1200_H_Pa', ...
    'Ushear1_1200','Ushear2_1200','Chezzy_1200','Tau_All1200_1', ...
    'Tau_All1200_2','tau_All1200_1_Pa','tau_All1200_2_Pa', ...
    'tau_Mann1200','tau_Mann_star1200','tau_rm_star1200', ...
    'tau_rm_starGaumann1200','tau_rm1200','uShear_Mann1200', ...
    'Tv_1200','TvRE1_1200','TvRE2_1200','TvH_1200','Tv1_1200','Tv2_1200','TvMann_1200', ...
    'Uratio_1200','UratioRE1_1200','UratioRE2_1200','UratioH_1200', ...
    'Uratio1_1200','Uratio2_1200','UratioMann_1200', ...
    'Mobile_1200','MobileRE1_1200','MobileRE2_1200','MobileH_1200', ...
    'Mobile1_1200','Mobile2_1200','MobileMann_1200','Fr_1200','Tausk_1200', ...
    'ShearVelocityNames1200','ShearVelocityLabels1200','ShearVelocityAll1200', ...
    'ExcessShearAll1200','UratioAll1200','MobileAll1200','TrateShearCorrTable', ...
    'ExcessShearCorrTable','CriticalShearSummaryTable','FileFound','MatchedFile', ...
    'Message','ShearStressAll','ShearStressTable'};
for k=1:numel(varsToExport)
    name=varsToExport{k}; ShearStressOut.(name)=eval(name);
end
ShearStressOut.SOIsPermille=logical(opt.SOIsPermille);
ShearStressOut.CriticalShieldsMethod=string(opt.CriticalShieldsMethod);
ShearStressOut.ThetaCrConstant=opt.ThetaCrConstant;
ShearStressOut.ApplySkinFrictionPartition=logical(opt.ApplySkinFrictionPartition);
ShearStressOut.DepthRNGoutUsed=~isempty(opt.DepthRNGout);
ShearStressOut.BeamAngleDeg=opt.BeamAngleDeg;
ShearStressOut.LogLawMaxZFraction=opt.LogLawMaxZFraction;
ShearStressOut.LogLawMinCells=opt.LogLawMinCells;

if logical(opt.AssignToCaller), assignStructToWorkspace(ShearStressOut,'caller'); end
if logical(opt.AssignToBase), assignStructToWorkspace(ShearStressOut,'base'); end

if logical(opt.Verbose)
    fprintf('\n=== Shear-stress processing summary ===\n');
    fprintf('Rows successfully processed: %d/%d\n',sum(Message=="OK"),nFiles);
    fprintf('Critical D50 source: D50mobile [mm], divided by 1000 once.\n');
    fprintf('Critical method: %s\n',string(opt.CriticalShieldsMethod));
    fprintf('Skin-friction partition applied: %d\n',logical(opt.ApplySkinFrictionPartition));
    fprintf('Log-law side-lobe guard: beam=%g deg, log-layer cap=%.2f, min cells=%d\n', ...
        opt.BeamAngleDeg,opt.LogLawMaxZFraction,opt.LogLawMinCells);
    fprintf('D65 measured rows: %d; estimated rows: %d\n', ...
        sum(D65source_1200=="input_D65"),sum(startsWith(D65source_1200,"estimated")));
    fprintf('No Cb, active-layer thickness or bedload transport was calculated.\n');
    fprintf('PlotTrateShear=%d, PlotCriticalComparison=%d, PlotExcessShear=%d\n', ...
        logical(opt.PlotTrateShear),logical(opt.PlotCriticalComparison),logical(opt.PlotExcessShear));
    fprintf('=======================================\n');
    disp(ShearStressTable);
end
end

%% Helpers
function [D0,D0BeamStd,source]=deriveDuneHeightFromDepthRNG(DepthRNGout,nFiles)
%DERIVEDUNEHEIGHTFROMDEPTHRNG Derive D0 [m] from spatial beam-depth differences.
%
% Preferred definition for every input row:
%   D0 = DepthRNGout.stdSpatial_1200
%
% stdSpatial_1200 is interpreted as the mean, over ensembles, of the
% simultaneous spatial standard deviation across RNG11/RNG22/RNG33/RNG44.
% It therefore represents the differences among the four beam depths, not
% the temporal standard deviation of an individual beam.
%
% Fallbacks:
%   1) mean of the positive finite values in stdSpatialALL_1200 for each row;
%   2) reconstruct the spatial standard-deviation series from the four raw
%      beam-range series and take its mean.
%
% D0BeamStd is retained only as an audit diagnostic containing the temporal
% standard deviation of each individual beam. It is not used to define D0.
% All depth/range values are assumed to be in metres.

D0=nan(nFiles,1);
D0BeamStd=nan(nFiles,4);
source=strings(nFiles,1);
source(:)="missing";
if isempty(DepthRNGout) || ~isstruct(DepthRNGout), return; end

% 1) Use the already processed spatial depth-difference metric directly.
if isfield(DepthRNGout,'stdSpatial_1200')
    direct=prepareVector(DepthRNGout.stdSpatial_1200,nFiles);
    valid=isfinite(direct) & direct>0;
    D0(valid)=direct(valid);
    source(valid)="DepthRNGout_stdSpatial_1200";
end

% 2) If only the ensemble-wise spatial differences are available, average
%    them row by row. This container may be a cell array, a row-by-ensemble
%    matrix, a column-by-ensemble matrix, or a vector for one input row.
if isfield(DepthRNGout,'stdSpatialALL_1200')
    C=DepthRNGout.stdSpatialALL_1200;
    for i=1:nFiles
        if isfinite(D0(i)) && D0(i)>0, continue; end
        vals=getSeriesFromContainerRow(C,i,nFiles);
        vals=double(vals(:));
        vals=vals(isfinite(vals) & vals>0);
        if ~isempty(vals)
            D0(i)=mean(vals);
            source(i)="DepthRNGout_mean_stdSpatialALL_1200";
        end
    end
end

beamRawNames={ ...
    {'RNG11_','RNG11','RNG11_all','RNG11_1200_all'}, ...
    {'RNG22_','RNG22','RNG22_all','RNG22_1200_all'}, ...
    {'RNG33_','RNG33','RNG33_all','RNG33_1200_all'}, ...
    {'RNG44_','RNG44','RNG44_all','RNG44_1200_all'}};
beamStdNames={ ...
    {'RNG11std_1200','RNG11Std_1200','stdRNG11_1200','STD_RNG11_1200'}, ...
    {'RNG22std_1200','RNG22Std_1200','stdRNG22_1200','STD_RNG22_1200'}, ...
    {'RNG33std_1200','RNG33Std_1200','stdRNG33_1200','STD_RNG33_1200'}, ...
    {'RNG44std_1200','RNG44Std_1200','stdRNG44_1200','STD_RNG44_1200'}};

% 3) Reconstruct the same spatial statistic from simultaneous raw ranges.
for i=1:nFiles
    beamSeries=cell(1,4);
    lengths=zeros(1,4);
    for b=1:4
        beamSeries{b}=getSeriesFromStructRow(DepthRNGout,beamRawNames{b},i,nFiles);
        D0BeamStd(i,b)=safeStd(beamSeries{b});
        if ~isfinite(D0BeamStd(i,b)) || D0BeamStd(i,b)<=0
            D0BeamStd(i,b)=getScalarFromStructRow( ...
                DepthRNGout,beamStdNames{b},i,nFiles);
        end
        lengths(b)=numel(beamSeries{b});
    end

    if isfinite(D0(i)) && D0(i)>0, continue; end
    if any(lengths==0), continue; end

    nCommon=min(lengths);
    beamMatrix=nan(nCommon,4);
    for b=1:4
        x=double(beamSeries{b}(:));
        beamMatrix(:,b)=x(1:nCommon);
    end

    spatialDifference=nan(nCommon,1);
    for j=1:nCommon
        vals=beamMatrix(j,:);
        vals=vals(isfinite(vals));
        if numel(vals)>=2
            spatialDifference(j)=std(vals,0);
        end
    end
    spatialDifference=spatialDifference(isfinite(spatialDifference) & spatialDifference>0);
    if ~isempty(spatialDifference)
        D0(i)=mean(spatialDifference);
        source(i)="DepthRNGout_reconstructed_spatial_std";
    end
end
end

function vals=getSeriesFromContainerRow(C,i,nFiles)
vals=[];
if iscell(C)
    if numel(C)>=i && ~isempty(C{i}), vals=C{i}; end
elseif isnumeric(C)
    C=double(C);
    if isvector(C)
        if nFiles==1
            vals=C(:);
        elseif numel(C)==nFiles
            vals=C(i);
        end
    elseif size(C,1)==nFiles
        vals=C(i,:).';
    elseif size(C,2)==nFiles
        vals=C(:,i);
    end
end
end

function x=getSeriesFromStructRow(S,names,i,nFiles)
x=[];
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    if iscell(C)
        if numel(C)>=i && ~isempty(C{i}), x=double(C{i}(:)); return; end
    elseif isnumeric(C)
        C=double(C);
        if ndims(C)==2 && size(C,1)==nFiles && size(C,2)>1
            x=C(i,:).'; return;
        elseif ndims(C)==2 && size(C,2)==nFiles && size(C,1)>1
            x=C(:,i); return;
        elseif nFiles==1 && numel(C)>1
            x=C(:); return;
        end
    end
end
end

function v=getScalarFromStructRow(S,names,i,nFiles)
v=NaN;
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    if iscell(C)
        if numel(C)>=i, v=safeMean(C{i}); end
    elseif isnumeric(C)
        C=double(C);
        if isvector(C) && numel(C)>=i, v=C(i);
        elseif size(C,1)==nFiles && size(C,2)>=1, v=C(i,1);
        elseif size(C,2)==nFiles && size(C,1)>=1, v=C(1,i);
        end
    end
    if isfinite(v), return; end
end
end

function vals=getFourValuesFromContainer(C,i,nFiles)
vals=[];
if iscell(C)
    if numel(C)>=i, vals=C{i}; end
elseif isnumeric(C)
    C=double(C);
    if size(C,1)==nFiles && size(C,2)>=4, vals=C(i,1:4);
    elseif size(C,2)==nFiles && size(C,1)>=4, vals=C(1:4,i).';
    elseif nFiles==1 && numel(C)>=4, vals=C(1:4);
    end
end
end

function s=safeStd(x)
s=NaN;
if isempty(x), return; end
x=double(x(:)); x=x(isfinite(x));
if numel(x)>=2, s=std(x,0); end
end

function [D65mm,source]=estimateD65FromPSD(D50mm,D84mm,D90mm)
%ESTIMATED65FROMPSD Estimate D65 [mm] from D50, D84 and D90 [mm].
% Grain sizes are taken as log-normal, so log(D) is linear in the standard-
% normal quantile z of the cumulative percentile. D65 is obtained from a
% least-squares fit of log(D) on z through the available points
% D50 (z=0), D84 (z=0.9945) and D90 (z=1.2816), evaluated at z65=0.3853.
% Falls back to a two-point fit, then to D50 if no spread is available.
D65mm=NaN; source="missing";
z50=0; z84=0.9945; z90=1.2816; z65=0.3853;

pts=zeros(0,2);                              % [z, log(D)]
if isfinite(D50mm) && D50mm>0, pts(end+1,:)=[z50,log(D50mm)]; end %#ok<AGROW>
if isfinite(D84mm) && D84mm>0, pts(end+1,:)=[z84,log(D84mm)]; end %#ok<AGROW>
if isfinite(D90mm) && D90mm>0, pts(end+1,:)=[z90,log(D90mm)]; end %#ok<AGROW>

if size(pts,1)>=2 && numel(unique(pts(:,1)))>=2
    c=polyfit(pts(:,1),pts(:,2),1);          % log(D) = c(1)*z + c(2)
    D65mm=exp(polyval(c,z65));
    if size(pts,1)>=3, source="estimated_D50_D84_D90";
    else,             source="estimated_two_point"; end
end

if ~isfinite(D65mm) || D65mm<=0               % no spread information
    if isfinite(D50mm) && D50mm>0
        D65mm=D50mm; source="fallback_D50";
    else
        D65mm=NaN;   source="missing";
    end
end
end

function x=getVectorField(S,fieldName,defaultValue,n)
if isfield(S,fieldName), x=S.(fieldName); else, x=defaultValue; end
x=prepareVector(x,n);
end

function x=getVectorFieldAny(S,fieldNames,defaultValue,n)
x=defaultValue;
for k=1:numel(fieldNames)
    if isfield(S,fieldNames{k}), x=S.(fieldNames{k}); break; end
end
x=prepareVector(x,n);
end

function x=prepareOptionalVector(x,n)
if isempty(x), x=nan(n,1); else, x=prepareVector(x,n); end
end

function x=prepareVector(x,n)
if isstring(x) || iscell(x), x=str2double(string(x)); end
x=double(x); x=x(:);
if isempty(x), x=nan(n,1);
elseif numel(x)==1 && n>1, x=repmat(x,n,1);
elseif numel(x)<n, x=[x;nan(n-numel(x),1)];
elseif numel(x)>n, x=x(1:n); end
end

function c=getCellField(S,fieldName,i)
c=[];
if isfield(S,fieldName)
    C=S.(fieldName);
    if iscell(C) && numel(C)>=i, c=C{i}; elseif isnumeric(C) && i==1, c=C; end
end
end

function matchedFile=matchADCPfile(dataADCP1200,thisADCPfilename)
matchedFile="";
if isempty(dataADCP1200) || ~isfield(dataADCP1200,'name'), return; end
target=cleanFileNameFromTable(thisADCPfilename);
if target=="", return; end
names=string({dataADCP1200.name}'); clean=strings(numel(names),1); bases=strings(numel(names),1);
for k=1:numel(names)
    clean(k)=cleanFileNameFromTable(names(k)); [~,b,~]=fileparts(char(clean(k))); bases(k)=string(b);
end
[~,tb,~]=fileparts(char(target));
idx=find(strcmpi(clean,target),1);
if isempty(idx), idx=find(strcmpi(bases,string(tb)),1); end
if isempty(idx), idx=find(startsWith(lower(bases),lower(string(tb))),1); end
if isempty(idx), idx=find(contains(lower(bases),lower(string(tb))),1); end
if isempty(idx), return; end
if isfield(dataADCP1200,'folder') && ~isempty(dataADCP1200(idx).folder)
    matchedFile=string(fullfile(dataADCP1200(idx).folder,dataADCP1200(idx).name));
else
    matchedFile=string(dataADCP1200(idx).name);
end
end

function s=cleanFileNameFromTable(s)
s=strtrim(string(s));
if ismissing(s) || s=="" || lower(s)=="nan", s=""; return; end
s=strtrim(erase(erase(s,"'"),'"'));
end

function [CZ,source]=calculateChezyFromSlope(U,Depth,SO,SOIsPermille)
CZ=NaN; source="missing";
if ~all(isfinite([U,Depth,SO])) || U<=0 || Depth<=0 || SO<=0, return; end
if SOIsPermille, S0=SO./1000; else, S0=SO; end
if ~isfinite(S0) || S0<=0, return; end
CZ=U./sqrt(Depth.*S0); source="calculated_from_SO_Depth_V";
end

function niValue=callWaterViscosity(T)
if ~isfinite(T), T=15; end
try
    if exist('WaterViscosity','file')==2
        niValue=WaterViscosity(T);
        if isfinite(niValue) && niValue>0, return; end
    end
catch
end
mu=2.414e-5.*10.^(247.8./(T+133.15));
niValue=mu./1000;
end

function [Rep,Dstar,thetaVanRijn,thetaBrownlie,thetaSoulsby, ...
    thetaCr,tauCrPa,uCr]=calculateCriticalShear( ...
    D50m,ni,R,g,ro,method,thetaConst)
%CALCULATECRITICALSHEAR Critical Shields stress from selected formulation.
Rep=NaN; Dstar=NaN; thetaVanRijn=NaN; thetaBrownlie=NaN;
thetaSoulsby=NaN; thetaCr=NaN; tauCrPa=NaN; uCr=NaN;
if ~all(isfinite([D50m,ni,R,g,ro])) || ...
        D50m<=0 || ni<=0 || R<=0 || g<=0 || ro<=0
    return
end

Rep=sqrt(R.*g.*D50m).*D50m./ni;
Dstar=D50m.*(R.*g./ni.^2).^(1/3);

thetaVanRijn=criticalShieldsVanRijn(Dstar);
thetaBrownlie=0.22.*Rep.^(-0.6) + ...
    0.06.*exp(-17.77.*Rep.^(-0.6));
thetaSoulsby=0.30./(1+1.2.*Dstar) + ...
    0.055.*(1-exp(-0.020.*Dstar));

method=lower(strtrim(string(method)));
switch method
    case {"vanrijn","van rijn"}
        thetaCr=thetaVanRijn;
    case "brownlie"
        thetaCr=thetaBrownlie;
    case {"soulsby","soulsby-whitehouse","soulsby & whitehouse"}
        thetaCr=thetaSoulsby;
    case "both"
        % Preserve critical_shear_river convention: Brownlie is primary.
        thetaCr=thetaBrownlie;
    case {"constant","fixed"}
        thetaCr=thetaConst;
    otherwise
        error(['Unknown CriticalShieldsMethod: %s. Use vanrijn, brownlie, ', ...
            'soulsby, both, constant or fixed.'],method);
end

tauCrPa=ro.*R.*g.*D50m.*thetaCr;
uCr=sqrt(tauCrPa./ro);
end

function [Cgrain,tauSkin,thetaSkin,uSkin,mobilityRatio,moving]= ...
    calculateSkinFriction(D50m,D90m,U,Rh,ro,R,g,thetaCr)
%CALCULATESKINFRICTION Van Rijn grain/skin-friction stress partition.
Cgrain=NaN; tauSkin=NaN; thetaSkin=NaN; uSkin=NaN;
mobilityRatio=NaN; moving=false;
if ~all(isfinite([D50m,D90m,U,Rh,ro,R,g,thetaCr])) || ...
        D50m<=0 || D90m<=0 || U<0 || Rh<=0 || ro<=0 || ...
        R<=0 || g<=0 || thetaCr<=0
    return
end
arg=12.*Rh./(3.*D90m);
if arg<=1
    return
end
Cgrain=18.*log10(arg);
if ~isfinite(Cgrain) || Cgrain<=0
    return
end
tauSkin=ro.*g.*(U./Cgrain).^2;
thetaSkin=tauSkin./(ro.*R.*g.*D50m);
uSkin=sqrt(tauSkin./ro);
mobilityRatio=thetaSkin./thetaCr;
moving=isfinite(mobilityRatio) && mobilityRatio>=1;
end

function thetaCr=criticalShieldsVanRijn(Dstar)
if ~isfinite(Dstar) || Dstar<=0, thetaCr=NaN;
elseif Dstar<=4, thetaCr=0.24.*Dstar.^(-1);
elseif Dstar<=10, thetaCr=0.14.*Dstar.^(-0.64);
elseif Dstar<=20, thetaCr=0.04.*Dstar.^(-0.1);
elseif Dstar<=150, thetaCr=0.013.*Dstar.^(0.29);
else, thetaCr=0.055; end
end

function [theta,tauPa]=stressFromUstar(uStar,D50m,R,g,ro)
theta=NaN; tauPa=NaN;
if ~all(isfinite([uStar,D50m,R,g,ro])) || uStar<0 || D50m<=0, return; end
theta=uStar.^2./(R.*g.*D50m);
tauPa=ro.*uStar.^2;
end

function T=excessShear(theta,thetaCr)
if all(isfinite([theta,thetaCr])) && thetaCr>0, T=(theta-thetaCr)./thetaCr; else, T=NaN; end
end

function r=safeRatio(a,b)
if all(isfinite([a,b])) && b>0, r=a./b; else, r=NaN; end
end

function [x,resnorm,residual,exitflag,output,fun]=callCalculateShearVel(z,V,kappa,ks,x0,plotFlag)
try
    if exist('calculateShearVel','file')==2
        [x,resnorm,residual,exitflag,output,fun]=calculateShearVel(z,V,kappa,ks,x0,plotFlag); return;
    end
catch
end
z=double(z(:)); V=double(V(:)); valid=isfinite(z)&isfinite(V)&z>0&V>=0;
if sum(valid)<2 || ~isfinite(ks) || ks<=0
    x=[NaN NaN]; resnorm=NaN; residual=nan(size(V)); exitflag=-1;
    output=struct('message','Fallback log-law fit failed.'); fun=[]; return;
end
X=log(30.*z(valid)./ks); Y=V(valid); validX=isfinite(X); X=X(validX); Y=Y(validX);
if numel(X)<2
    x=[NaN NaN]; resnorm=NaN; residual=nan(size(V)); exitflag=-1;
    output=struct('message','Fallback log-law fit failed.'); fun=[]; return;
end
coef=polyfit(X,Y,1); x=[coef(1).*kappa,coef(2)]; Yfit=polyval(coef,X);
residual=Y-Yfit; resnorm=sum(residual.^2); exitflag=1;
output=struct('message','Fallback log-law linear fit used.'); fun=coef;
end

function [uw_mean,vw_mean,x1,x2,fun,fun2,TKE1200,tau_turb1200]= ...
    callCalculateREShearVelXYcord(V1x_1,V1y_1,V1w_1,z1_1200,DepthValue,x0,coef,plotFlag)
try
    if exist('calculateREShearVelXYcord','file')==2
        [uw_mean,vw_mean,x1,x2,fun,fun2,TKE1200,tau_turb1200]= ...
            calculateREShearVelXYcord(V1x_1,V1y_1,V1w_1,z1_1200,DepthValue,x0,coef,plotFlag); return;
    end
catch
end
mx=nanmeanLocal(V1x_1,2); my=nanmeanLocal(V1y_1,2); mw=nanmeanLocal(V1w_1,2);
ux=V1x_1-repmat(mx,1,size(V1x_1,2)); uy=V1y_1-repmat(my,1,size(V1y_1,2)); uz=V1w_1-repmat(mw,1,size(V1w_1,2));
uw_mean=ux.*uz; vw_mean=uy.*uz; TKE1200=0.5.*(ux.^2+uy.^2+uz.^2); tau_turb1200=-uw_mean;
x1=[sqrt(abs(safeMean(uw_mean(:)))),NaN]; x2=[sqrt(abs(safeMean(vw_mean(:)))),NaN]; fun=[]; fun2=[];
end

function [Ushear1,Ushear2,Chezzy]=shearVanRijn(U,Depth,a,D0,D90)
% Ushear1 and Chezzy depend only on U, Depth, a and D90.
% Ushear2 additionally requires D0. Therefore a missing D0 must not
% invalidate the grain-roughness/skin-friction estimate Ushear1.
Ushear1=NaN; Ushear2=NaN; Chezzy=NaN;

% Grain-roughness branch: independent of D0.
if all(isfinite([U,Depth,a,D90])) && U>=0 && Depth>0 && a>0 && D90>0
    arg1=12.*Depth./(a.*D90);
    if isfinite(arg1) && arg1>1
        Chezzy=18.*log10(arg1);
        d1=5.75.*log10(arg1);
        if isfinite(d1) && d1~=0
            Ushear1=U./d1;
        end
    end
end

% Dune/bedform branch: requires a valid D0.
if all(isfinite([U,Depth,a,D0])) && U>=0 && Depth>0 && a>0 && D0>0
    arg2=Depth./(a.*D0./30);
    if isfinite(arg2) && arg2>1
        d2=2.5.*log(arg2);
        if isfinite(d2) && d2~=0
            Ushear2=U./d2;
        end
    end
end
end

function [tau_Mann,tau_Mann_star,tau_rm_star,tau_rm_starGaumann,tau_rm,uShear_Mann]= ...
    shearMannWC(VmodAll_S,SO,D65m,D50m,FsS,Dstd,SOIsPermille,ro)
% ro [kg/m^3] is used to convert tau_Mann to shear velocity. It defaults to
% 1000 if not supplied, preserving the previous behaviour exactly.
if nargin<8 || ~isfinite(ro) || ro<=0, ro=1000; end
tau_Mann=NaN; tau_Mann_star=NaN; tau_rm_star=NaN; tau_rm_starGaumann=NaN; tau_rm=NaN; uShear_Mann=NaN;
if SOIsPermille, S0=SO./1000; else, S0=SO; end
if all(isfinite([VmodAll_S,S0,D65m])) && VmodAll_S>=0 && S0>0 && D65m>0
    tau_Mann=17.*(S0.*D65m).^0.25.*VmodAll_S.^1.5;
    if tau_Mann>=0, uShear_Mann=sqrt(tau_Mann./ro); end
end
if all(isfinite([uShear_Mann,D50m])) && D50m>0
    tau_Mann_star=uShear_Mann.^2./(1.65.*9.81.*D50m);
end
if isfinite(FsS), tau_rm_star=0.021+0.0105.*exp(-20.*FsS); end
if isfinite(Dstd), tau_rm_starGaumann=0.03+(0.052-0.03)./(1+exp(7.1.*(Dstd-1.66))); end
if all(isfinite([tau_rm_star,D50m])) && D50m>0
    tau_rm=tau_rm_star.*(1.65.*9.81.*D50m);
end
end

function T=makeCorrelationTable(names,labels,X,yMean,yMax,minN,prefix)
n=numel(names); Rmean=nan(n,1); Nmean=zeros(n,1); Rmax=nan(n,1); Nmax=zeros(n,1);
for k=1:n
    [Rmean(k),Nmean(k)]=corrNoNaN(X(:,k),yMean,minN);
    [Rmax(k),Nmax(k)]=corrNoNaN(X(:,k),yMax,minN);
end
prefix = char(string(prefix));
varNames = {[prefix 'Variable'], [prefix 'Method'], ...
    'R_TrateMean','N_TrateMean','R_MaxTrate','N_MaxTrate'};
T=table(names(:),labels(:),Rmean,Nmean,Rmax,Nmax,'VariableNames',varNames);
end

function T=makeCriticalSummary(names,labels,U,ucr,Tv)
n=numel(names); Nvalid=zeros(n,1); Nabv=zeros(n,1); PercentAbove=nan(n,1); MedianRatio=nan(n,1); MedianExcess=nan(n,1);
for k=1:n
    valid=isfinite(U(:,k))&isfinite(ucr)&ucr>0; Nvalid(k)=sum(valid);
    if Nvalid(k)>0
        Nabv(k)=sum(U(valid,k)>ucr(valid)); PercentAbove(k)=100.*Nabv(k)./Nvalid(k);
        MedianRatio(k)=median(U(valid,k)./ucr(valid));
        tv=Tv(valid,k); tv=tv(isfinite(tv)); if ~isempty(tv), MedianExcess(k)=median(tv); end
    end
end
T=table(names(:),labels(:),Nvalid,Nabv,PercentAbove,MedianRatio,MedianExcess, ...
    'VariableNames',{'ShearVelocityVariable','ShearVelocityMethod','NValid','NAboveCritical','PercentAboveCritical','MedianUstarOverUcr','MedianExcessShear'});
end

function plotMetricVsTransport(names,labels,X,yMean,yMax,selected,useLogX,useLogY,minN,figureName,xLabelText)
idx=selectMethods(names,labels,selected); if isempty(idx), warning('No valid plotting methods selected.'); return; end
isExcessPlot=contains(lower(string(xLabelText)),'excess');
allX=X(:,idx); xv=allX(isfinite(allX));
if useLogX, xv=xv(xv>0); end
if isExcessPlot && ~useLogX, xv=[xv;0]; end
if isempty(xv), sharedX=[]; else, sharedX=expandedLimits(xv,useLogX); end
yv=[yMean(isfinite(yMean));yMax(isfinite(yMax))]; if useLogY, yv=yv(yv>0); end
if isempty(yv), sharedY=[]; else, sharedY=expandedLimits(yv,useLogY); end
n=numel(idx); nr=ceil(sqrt(n)); nc=ceil(n/nr); figure('Name',figureName,'Color','w');
try, tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact'); tiled=true; catch, tiled=false; end
for jj=1:n
    k=idx(jj); if tiled, nexttile; else, subplot(nr,nc,jj); end; hold on;
    x=X(:,k); vm=isfinite(x)&isfinite(yMean); vx=isfinite(x)&isfinite(yMax);
    if useLogX, vm=vm&x>0; vx=vx&x>0; end
    if useLogY, vm=vm&yMean>0; vx=vx&yMax>0; end
    h1=plot(x(vm),yMean(vm),'o','LineStyle','none','DisplayName','Mean Trate');
    h2=plot(x(vx),yMax(vx),'^','LineStyle','none','DisplayName','Max Trate');
    if isExcessPlot, xline(0,'--','Critical, T=0','HandleVisibility','off'); end
    grid on; box on; xlabel(xLabelText); ylabel('Transport rate'); title(labels(k),'Interpreter','none');
    if useLogX, set(gca,'XScale','log'); end; if useLogY, set(gca,'YScale','log'); end
    if ~isempty(sharedX), xlim(sharedX); end; if ~isempty(sharedY), ylim(sharedY); end
    legend([h1 h2],'Location','best');
    [r1,n1]=corrNoNaN(x,yMean,minN); [r2,n2]=corrNoNaN(x,yMax,minN);
    addTextInsideAxes(sprintf('r Mean = %.2f (n=%d)\nr Max = %.2f (n=%d)',r1,n1,r2,n2));
end
end

function plotCriticalComparison(names,labels,U,ucr,Tv)
vals=[ucr(:);U(:)]; vals=vals(isfinite(vals)&vals>=0); if isempty(vals), warning('No valid critical-comparison data.'); return; end
lim=expandedLimits([0;vals],false); n=numel(names); nr=ceil(sqrt(n)); nc=ceil(n/nr);
figure('Name','Critical shear velocity comparison','Color','w');
try, tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact'); tiled=true; catch, tiled=false; end
for k=1:n
    if tiled, nexttile; else, subplot(nr,nc,k); end; hold on;
    valid=isfinite(ucr)&isfinite(U(:,k))&ucr>=0&U(:,k)>=0;
    below=valid&U(:,k)<=ucr; above=valid&U(:,k)>ucr;
    h1=plot(ucr(below),U(below,k),'o','LineStyle','none','DisplayName','At/below critical');
    h2=plot(ucr(above),U(above,k),'^','LineStyle','none','DisplayName','Above critical');
    plot(lim,lim,'--','DisplayName','u_* = u_{cr}');
    xlim(lim); ylim(lim); axis square; grid on; box on;
    xlabel('Critical shear velocity, u_{cr} [m/s]'); ylabel('Calculated shear velocity, u_* [m/s]');
    title(labels(k),'Interpreter','none'); legend([h1 h2],'Location','best');
    tv=Tv(valid,k); tv=tv(isfinite(tv)); med=NaN; if ~isempty(tv), med=median(tv); end
    addTextInsideAxes(sprintf('n=%d\nAbove=%d (%.1f%%)\nMedian T=%.2f',sum(valid),sum(above),100*sum(above)/max(sum(valid),1),med));
end
end

function idx=selectMethods(names,labels,selected)
if isempty(selected), idx=1:numel(names); return; end
if isnumeric(selected), idx=unique(round(selected(:)'),'stable'); idx=idx(idx>=1&idx<=numel(names)); return; end
req=lower(strtrim(string(selected(:)))); nv=lower(string(names(:))); lv=lower(string(labels(:))); idx=[];
for k=1:numel(req)
    j=find(nv==req(k),1); if isempty(j), j=find(lv==req(k),1); end
    if ~isempty(j), idx(end+1)=j; end %#ok<AGROW>
end
idx=unique(idx,'stable');
end

function lim=expandedLimits(v,useLog)
v=v(isfinite(v));
if useLog
    v=v(v>0); lo=min(v); hi=max(v); if lo==hi, lo=lo/1.2; hi=hi*1.2; end
    lim=[10^(log10(lo)-0.05*(log10(hi)-log10(lo))),10^(log10(hi)+0.05*(log10(hi)-log10(lo)))];
else
    lo=min(v); hi=max(v); if lo==hi, d=max(abs(lo)*0.1,1e-6); else, d=0.05*(hi-lo); end
    lim=[lo-d,hi+d];
end
end

function [r,n]=corrNoNaN(x,y,minN)
x=double(x(:)); y=double(y(:)); m=min(numel(x),numel(y)); x=x(1:m); y=y(1:m); valid=isfinite(x)&isfinite(y); n=sum(valid);
if n<minN || numel(unique(x(valid)))<2 || numel(unique(y(valid)))<2, r=NaN; return; end
C=corrcoef(x(valid),y(valid)); r=C(1,2);
end

function addTextInsideAxes(txt)
drawnow; xl=xlim; yl=ylim;
if strcmpi(get(gca,'XScale'),'log') && all(xl>0), xpos=10^(log10(xl(1))+0.05*diff(log10(xl))); else, xpos=xl(1)+0.05*diff(xl); end
if strcmpi(get(gca,'YScale'),'log') && all(yl>0), ypos=10^(log10(yl(2))-0.08*diff(log10(yl))); else, ypos=yl(2)-0.08*diff(yl); end
text(xpos,ypos,txt,'VerticalAlignment','top','BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Margin',4);
end

function y=safeMean(x)
if isempty(x), y=NaN; return; end
x=double(x(:)); x=x(isfinite(x)); if isempty(x), y=NaN; else, y=mean(x); end
end

function y=nanmeanLocal(A,dim)
if nargin<2, dim=1; end
A=double(A); valid=isfinite(A); A(~valid)=0; count=sum(valid,dim); y=sum(A,dim)./count; y(count==0)=NaN;
end

function assignStructToWorkspace(S,workspaceName)
f=fieldnames(S); for k=1:numel(f), assignin(workspaceName,f{k},S.(f{k})); end
end