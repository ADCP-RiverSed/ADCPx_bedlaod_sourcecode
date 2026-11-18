function BedloadTransportOut = processBedloadTrate(ADCPin, ShearStressOut, BTbedloadOut, varargin)
%PROCESSBEDLOADTRATE Calculate bedload concentration, active-layer
% thickness and kinematic bedload transport.
%
% This function is intentionally separated from processShearStress1200.
% All bedload concentration, calibrated Cb, settling velocity, empirical
% active-layer thickness, Van Rijn transport and kinematic transport
% calculations are contained in this function. Their existing variable
% names are preserved.
%
% Critical shear is NOT independently recalculated here. Tucr_1200 is
% imported from ShearStressOut as the controlling physical threshold. The
% equivalent Dstar, ThetaCr and ucr representations are also imported; if
% an older ShearStressOut lacks one of those equivalent fields, it is
% reconstructed algebraically from Tucr_1200, D50mobile and ni:
%
%   Dstar_1200       dimensionless particle diameter [-]
%   ThetaCr_1200     critical Shields stress [-]
%   Tucr_1200        critical physical shear stress [Pa]
%   Tau_Crn          alias of Tucr_1200 [Pa]
%   ucr_1200         critical shear velocity [m/s]
%
% Method-specific critical-stress aliases are also accepted:
%
%   TucrLog_1200, TucrRE1_1200, TucrRE2_1200, TucrH_1200,
%   Tucr1_1200, Tucr2_1200 and TucrMann_1200.
%
% Excess shear used for every Cb calculation is evaluated from physical
% stresses with identical units:
%
%   T = (tau_method_Pa - Tucr_method_Pa) ./ Tucr_method_Pa
%
% The function also imports the excess-shear variables calculated by
% processShearStress1200 and stores the difference for diagnostics.
%
% IMPORTANT GRAIN-SIZE CONVENTION
% -------------------------------
% D50mobile is interpreted as millimetres [mm]. The preferred source is
% ShearStressOut.D50mobile, but only when the physical-sample input confirms
% that a mobile D50 is available. If needed, it is reconstructed from the
% mobile-sample fields D50mobile, D50_14, D50_05 or the legacy D50 field.
% D_50grabM is retained for diagnostics but is never substituted for a
% missing mobile D50. The mobile size is converted once:
%
%   D50critical_m = D50mobile ./ 1000;
%
% It is the D50 used for:
%   - Dstar / Dv
%   - active-layer thickness
%   - Cb
%   - Van Rijn transport
%   - kinematic transport
%
% If no positive finite mobile D50 is available, every kinematic transport
% output for that row is set to zero. Concentration, active-layer-thickness
% and critical-shear diagnostics remain NaN where they cannot be calculated.
%
% INSTRUMENT COEFFICIENT
% ----------------------
% Ninstrument is dimensionless and multiplies every kinematic transport
% estimate: q_s = Ninstrument*V_b*C_b*D_active*rho_s. It does not modify
% the empirical Van Rijn QvR transport-capacity formulas.
%
% UNITS
% -----
% D50critical_mm            [mm]
% D50critical_m             [m]
% ni                        [m^2/s]
% shear velocity            [m/s]
% tau and Tucr              [Pa]
% Theta and T               [-]
% Cb                        [-], volumetric concentration
% active-layer thickness    [m]
% BTmodf_1200               [m/s]
% qs and QvR                [kg/(m s)]
%
% Cb validity and final concentration
% ------------------------------------
% A Van Rijn reference concentration is calculated for every available
% applied-shear method:
%
%   Log-law, Reynolds 1, Reynolds 2, Hydraulic/Chezy,
%   Van Rijn 1, Van Rijn 2, Manning and explicit grain/skin friction.
%
% A method value is valid when:
%
%   CbMin <= Cb <= CbMax && isreal(Cb) && isfinite(Cb)
%
% The default inclusive limits are:
%
%   0.005 <= Cb <= 0.55
%
% CbFinal_1200 is the mean or median of the valid method values. If no
% method produces a valid value, the optional fallback CbFallback is used
% (default 0.1 when UseCbFallback=true).
%
% Final active-layer thickness
% ----------------------------
% The seven empirical active-layer-thickness estimates returned by
% Dthickness6 are retained in DDt1200. DthicknessFinal_1200 is the mean or
% median of estimates within:
%
%   DthicknessMinFactor*D50mobile <= Dthickness <=
%   DthicknessMaxFactor*D50mobile
%
% Defaults are 1.5 and 10. If no estimate is valid, the optional fallback
% DthicknessFallbackFactor*D50mobile is used (default factor 2).
%
% Kinematic transport calculations
% --------------------------------
% All original method-specific kinematic calculations are retained exactly:
%
%   qs_method = Ninstrument * BTmodf_1200 * Cb_method * ...
%               (n*D50mobile) * ros
%
% A separate final kinematic result is additionally calculated from the
% bounded/fallback final inputs:
%
%   qsKinematicFinal_1200 = Ninstrument * BTmodf_1200 * ...
%       CbFinal_1200 * DthicknessFinal_1200 * ros
%
% Ninstrument is an instrument-specific dimensionless correction
% coefficient. It may be scalar or one value per input row and defaults to
% 1, which reproduces the uncorrected equations.
%
% Thus, the original n*D50 calculations are not replaced. The final value
% is an additional estimate based on the selected final concentration and
% selected final active-layer thickness.
%
% FINAL TRANSPORT COMPARISON PLOTS
% --------------------------------
% When PlotFinalTrateComparison=true, one figure is created with:
%   1) calculated qsKinematicFinal_1200 versus measured TrateMean;
%   2) calculated qsKinematicFinal_1200 versus measured MaxTrate.
% Both panels use the same x/y limits and include a 1:1 agreement line.
% Pearson r, bias, MAE and RMSE are reported in each panel. Calculated
% transport is in kg/(m s). MeasuredTrateToKgFactor converts the measured
% TrateMean and MaxTrate inputs to kg/(m s): use 1 when they are already in
% kg/(m s), or 0.001 when they are supplied in g/(m s).
%
% Recommended workflow:
%
%   ShearStressOut = processShearStress1200( ...
%       dataADCP1200, PhSamplesin, WaterProfileOut, ...
%       'ni', ni, ...
%       'ks', ks, ...
%       'SO', SO, ...
%       'CZ', [], ...
%       'aVanRijn', 2, ...
%       'CriticalShieldsMethod', 'vanrijn', ...
%       'PlotTrateShear', true, ...
%       'PlotCriticalComparison', true, ...
%       'PlotExcessShear', true, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
%   BedloadTransportOut = processBedloadTrate( ...
%       PhSamplesin, ShearStressOut, BTbedloadOut, ...
%       'ni', ni, ...
%       'n', 5, ...
%       'Ninstrument', 1, ...
%       'CbMin', 0.005, ...
%       'CbMax', 0.55, ...
%       'UseCbFallback', true, ...
%       'CbFallback', 0.1, ...
%       'CbFinalStatistic', 'mean', ...
%       'DthicknessMinFactor', 1.5, ...
%       'DthicknessMaxFactor', 10, ...
%       'UseDthicknessFallback', true, ...
%       'DthicknessFallbackFactor', 2, ...
%       'DthicknessFinalStatistic', 'mean', ...
%       'MeasuredTrateToKgFactor', 1, ... % use 0.001 for measured g/(m s)
%       'PlotFinalTrateComparison', true, ...
%       'PlotFinalTrateUseLogAxes', false, ...
%       'AssignToCaller', true, ...
%       'AssignToBase', true, ...
%       'Verbose', true);
%
% No files are saved. Variable-length or per-position series are exported
% only as cell arrays, with cell i corresponding to input measurement i.
% No temporary value from the final loop iteration is exported.

%% Parse inputs

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'ADCPin', @(x) isstruct(x));
addRequired(p, 'ShearStressOut', @(x) isstruct(x));
addRequired(p, 'BTbedloadOut', @(x) isstruct(x));

% Physical constants
addParameter(p, 'ro', 1000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ros', 2650, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'g', 9.81, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'R', 1.65, @(x) isnumeric(x) && isscalar(x) && x > 0);

% Kinematic viscosity [m^2/s], scalar or vector
addParameter(p, 'ni', [], @(x) isempty(x) || isnumeric(x));

% Active-layer multiplier in Dthickness = n*D50
addParameter(p, 'n', 1, @(x) isnumeric(x) && ~isempty(x));

% Instrument-specific correction coefficient for every kinematic model.
% Scalar or one value per input row/instrument; default 1 means no correction.
addParameter(p, 'Ninstrument', 1, ...
    @(x) isnumeric(x) && ~isempty(x));

% Valid concentration limits: inclusive bounds
addParameter(p, 'CbMin', 0.005, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'CbMax', 0.55, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

% Final Cb aggregation and optional fallback
addParameter(p, 'CbFinalStatistic', 'mean', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'UseCbFallback', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'CbFallback', 0.1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

% Include the custom calibrated concentration in the final aggregation.
% The default is false because it is not an independent shear method.
addParameter(p, 'IncludeCalibratedCbInFinal', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

% Final active-layer-thickness aggregation and optional fallback
addParameter(p, 'DthicknessMinFactor', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'DthicknessMaxFactor', 10, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'DthicknessFinalStatistic', 'mean', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'UseDthicknessFallback', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'DthicknessFallbackFactor', 2, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

% Use the custom calibrated concentration calculation
addParameter(p, 'ApplyNewCb', true, @(x) islogical(x) || isnumeric(x));

% Relative tolerance used only for critical-input consistency diagnostics
addParameter(p, 'CriticalConsistencyTolerance', 1e-6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);

% Final calculated-versus-measured transport comparison. Calculated qs is
% in kg/(m s). Multiply measured TrateMean/MaxTrate by this factor before
% plotting and calculating statistics (1 for kg/(m s), 0.001 for g/(m s)).
addParameter(p, 'MeasuredTrateToKgFactor', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'PlotFinalTrateComparison', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'PlotFinalTrateUseLogAxes', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'PlotFinalTrateMinN', 3, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);

% Workspace/output
addParameter(p, 'AssignToCaller', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'AssignToBase', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, ADCPin, ShearStressOut, BTbedloadOut, varargin{:});
opt = p.Results;

if opt.CbMin >= opt.CbMax
    error('CbMin must be smaller than CbMax.');
end

opt.CbFinalStatistic = lower(string(opt.CbFinalStatistic));
if ~any(opt.CbFinalStatistic == ["mean","median"])
    error('CbFinalStatistic must be ''mean'' or ''median''.');
end

if logical(opt.UseCbFallback) && ...
        (opt.CbFallback < opt.CbMin || opt.CbFallback > opt.CbMax)
    error(['CbFallback must lie within the inclusive range ', ...
        '[CbMin,CbMax] when UseCbFallback is true.']);
end

if opt.DthicknessMinFactor >= opt.DthicknessMaxFactor
    error('DthicknessMinFactor must be smaller than DthicknessMaxFactor.');
end

opt.DthicknessFinalStatistic = ...
    lower(string(opt.DthicknessFinalStatistic));

if ~any(opt.DthicknessFinalStatistic == ["mean","median"])
    error('DthicknessFinalStatistic must be ''mean'' or ''median''.');
end

if logical(opt.UseDthicknessFallback) && ...
        (opt.DthicknessFallbackFactor < opt.DthicknessMinFactor || ...
         opt.DthicknessFallbackFactor > opt.DthicknessMaxFactor)
    error(['DthicknessFallbackFactor must lie within the inclusive ', ...
        'range [DthicknessMinFactor,DthicknessMaxFactor] when ', ...
        'UseDthicknessFallback is true.']);
end

%% Input dimensions and identifiers

if isfield(ADCPin, 'ADCPfilename')
    ADCPfilename = string(ADCPin.ADCPfilename(:));
elseif isfield(ShearStressOut, 'ADCPfilename')
    ADCPfilename = string(ShearStressOut.ADCPfilename(:));
else
    error('ADCPfilename is missing from both ADCPin and ShearStressOut.');
end

nFiles = numel(ADCPfilename);
InputRow = (1:nFiles)';

DistanceFromLeftBank = getVectorFieldAny(ADCPin, ...
    {'DistanceFromLeftBank'}, nan(nFiles,1), nFiles);
TrateMean = getVectorFieldAny(ADCPin, ...
    {'TrateMean','MeanTrate','MTrateAll','MeanTrateAll'}, nan(nFiles,1), nFiles);
MaxTrate = getVectorFieldAny(ADCPin, ...
    {'MaxTrate','TrateMax','MTrateMax','MaxTrateAll'}, nan(nFiles,1), nFiles);

MeanTrate14 = getVectorFieldAny(ADCPin, ...
    {'MeanTrate14','TrateMean14','MTrate14'}, nan(nFiles,1), nFiles);
MeanTrate05 = getVectorFieldAny(ADCPin, ...
    {'MeanTrate05','TrateMean05','MTrate05'}, nan(nFiles,1), nFiles);

%% Mobile D50 inputs

D_50grabM = getVectorFieldAny(ADCPin, ...
    {'D_50grabM','D50grabM','D50_grabM','D50grab'}, nan(nFiles,1), nFiles);

% Other grain-size inputs retained for the original calculations.
D50_14 = getVectorFieldAny(ADCPin, ...
    {'D50_14','D50_14mm','D50mobile14'}, nan(nFiles,1), nFiles);

D50_05 = getVectorFieldAny(ADCPin, ...
    {'D50_05','D50_05mm','D50mobile05'}, nan(nFiles,1), nFiles);

D50 = getVectorFieldAny(ADCPin, ...
    {'D50','D_50'}, nan(nFiles,1), nFiles);

% Mobile D50 reported directly by the physical-sample reader, when present.
% If this field exists, its NaN/zero rows are authoritative: no legacy
% mobile-size field is allowed to replace them.
directMobileFieldPresent = any(isfield(ADCPin, ...
    {'D50mobile','D50_mobile','D50Mob'}));
D50mobileInput = getVectorFieldAny(ADCPin, ...
    {'D50mobile','D50_mobile','D50Mob'}, nan(nFiles,1), nFiles);

% Retain the shear-processing value for exact consistency, but use it only
% when at least one physical mobile-sample field is valid. This prevents a
% grab-size fallback in ShearStressOut from being reintroduced here.
D50mobileShear = getVectorFromStruct(ShearStressOut, ...
    'D50mobile', nan(nFiles,1), nFiles);
D50mobile = nan(nFiles,1);

for iD50 = 1:nFiles
    if directMobileFieldPresent
        mobileCandidates = D50mobileInput(iD50);
    else
        mobileCandidates = [ ...
            D50_14(iD50), ...
            D50_05(iD50), ...
            D50(iD50)];
    end
    mobileCandidates = mobileCandidates(isfinite(mobileCandidates) & mobileCandidates>0);

    if ~isempty(mobileCandidates)
        if isfinite(D50mobileShear(iD50)) && D50mobileShear(iD50)>0
            D50mobile(iD50) = D50mobileShear(iD50);
        else
            D50mobile(iD50) = mean(mobileCandidates);
        end
    end
end

D50critical_mm = D50mobile;
D50critical_m = D50critical_mm ./ 1000;
D50mobileM = D50critical_m;
NoMobileD50_1200 = ~isfinite(D50mobile) | D50mobile <= 0;
KinematicZeroDueToNoMobileD50_1200 = NoMobileD50_1200;

if any(NoMobileD50_1200)
    bad = find(NoMobileD50_1200);
    warning('processBedloadTrate:noMobileD50', ...
        ['No positive finite mobile D50 is available for row(s) %s. ', ...
         'All kinematic transport outputs are set to zero for those rows; ', ...
         'grab D50 is not substituted.'], mat2str(bad(:)'));
end

%% Viscosity input

ni = prepareOptionalVector(opt.ni, nFiles);
if all(isnan(ni)) && isfield(ShearStressOut, 'ni')
    ni = prepareVector(ShearStressOut.ni, nFiles);
end
if all(isnan(ni)) && isfield(ADCPin, 'ni')
    ni = prepareVector(ADCPin.ni, nFiles);
end

invalidNi = (~isfinite(ni) | ni <= 0) & ~NoMobileD50_1200;
if any(invalidNi)
    bad = find(invalidNi);
    error(['ni must contain one positive kinematic-viscosity value [m^2/s] ', ...
        'for every row with a valid mobile D50. Invalid rows: %s'], ...
        mat2str(bad(:)'));
end

%% Active-layer multiplier

n = prepareVector(opt.n, nFiles);
if any(~isfinite(n) | n <= 0)
    error('All n values must be positive and finite.');
end

Ninstrument = prepareVector(opt.Ninstrument, nFiles);
if any(~isfinite(Ninstrument) | Ninstrument <= 0)
    bad = find(~isfinite(Ninstrument) | Ninstrument <= 0);
    error(['Ninstrument must contain one positive finite coefficient for ', ...
        'every input row/instrument. Invalid rows: %s'], mat2str(bad(:)'));
end

DthicknessKin_1200 = n .* D50mobileM;

%% Apparent BT bedload velocity

BTmodf_1200 = getBTmodfVector(BTbedloadOut, nFiles);

%% Shear velocities, stresses and critical inputs from ShearStressOut

% Shear velocities [m/s]
UsheaAll_1200 = getVectorFromStruct(ShearStressOut, ...
    'UsheaAll_1200', nan(nFiles,1), nFiles);

U_shearRE1_1200 = getVectorFromStruct(ShearStressOut, ...
    'U_shearRE1_1200', nan(nFiles,1), nFiles);

U_shearRE2_1200 = getVectorFromStruct(ShearStressOut, ...
    'U_shearRE2_1200', nan(nFiles,1), nFiles);

u_shearH1200 = getVectorFromStruct(ShearStressOut, ...
    'u_shearH1200', nan(nFiles,1), nFiles);

Ushear1_1200 = getVectorFromStruct(ShearStressOut, ...
    'Ushear1_1200', nan(nFiles,1), nFiles);

Ushear2_1200 = getVectorFromStruct(ShearStressOut, ...
    'Ushear2_1200', nan(nFiles,1), nFiles);

uShear_Mann1200 = getVectorFromStruct(ShearStressOut, ...
    'uShear_Mann1200', nan(nFiles,1), nFiles);

% Explicit Van Rijn grain/skin-friction method
uShear_skin1200 = getVectorFromStruct(ShearStressOut, ...
    'uShear_skin1200', nan(nFiles,1), nFiles);

theta_skin_1200 = getVectorFromStruct(ShearStressOut, ...
    'theta_skin_1200', nan(nFiles,1), nFiles);

tau_skin_1200 = getVectorFromStruct(ShearStressOut, ...
    'tau_skin_1200', nan(nFiles,1), nFiles);

% Dimensionless Shields stresses from processShearStress1200
Tau_All1200 = getVectorFromStruct(ShearStressOut, ...
    'Tau_All1200', nan(nFiles,1), nFiles);

Tau_AllRE1_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tau_AllRE1_1200', nan(nFiles,1), nFiles);

Tau_AllRE2_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tau_AllRE2_1200', nan(nFiles,1), nFiles);

Tau_All1200_H = getVectorFromStruct(ShearStressOut, ...
    'Tau_All1200_H', nan(nFiles,1), nFiles);

Tau_All1200_1 = getVectorFromStruct(ShearStressOut, ...
    'Tau_All1200_1', nan(nFiles,1), nFiles);

Tau_All1200_2 = getVectorFromStruct(ShearStressOut, ...
    'Tau_All1200_2', nan(nFiles,1), nFiles);

tau_Mann_star1200 = getVectorFromStruct(ShearStressOut, ...
    'tau_Mann_star1200', nan(nFiles,1), nFiles);

% Physical method stresses [Pa]
tau_All1200_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_All1200_Pa', nan(nFiles,1), nFiles);

tau_AllRE1_1200_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_AllRE1_1200_Pa', nan(nFiles,1), nFiles);

tau_AllRE2_1200_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_AllRE2_1200_Pa', nan(nFiles,1), nFiles);

tau_All1200_H_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_All1200_H_Pa', nan(nFiles,1), nFiles);

tau_All1200_1_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_All1200_1_Pa', nan(nFiles,1), nFiles);

tau_All1200_2_Pa = getVectorFromStruct(ShearStressOut, ...
    'tau_All1200_2_Pa', nan(nFiles,1), nFiles);

tau_Mann1200 = getVectorFromStruct(ShearStressOut, ...
    'tau_Mann1200', nan(nFiles,1), nFiles);

% Common critical quantities calculated by processShearStress1200.
% Tucr_1200 is the controlling input [Pa]. Equivalent representations are
% read from ShearStressOut when present. For compatibility with an older
% ShearStressOut, missing Dstar, ThetaCr or ucr are reconstructed directly
% from the imported Tucr, D50mobile and ni; no new critical-shear criterion
% is selected in this function.
Tucr_1200 = requireVectorFromStruct(ShearStressOut, ...
    'Tucr_1200', nFiles);

Dstar_1200 = getVectorFromStruct(ShearStressOut, ...
    'Dstar_1200', nan(nFiles,1), nFiles);

ThetaCr_1200 = getVectorFromStruct(ShearStressOut, ...
    'ThetaCr_1200', nan(nFiles,1), nFiles);

ucr_1200 = getVectorFromStruct(ShearStressOut, ...
    'ucr_1200', nan(nFiles,1), nFiles);

missingDstar = ~isfinite(Dstar_1200) | Dstar_1200 <= 0;
Dstar_1200(missingDstar) = D50critical_m(missingDstar) .* ...
    (opt.R .* opt.g ./ ni(missingDstar).^2).^(1/3);

missingTheta = ~isfinite(ThetaCr_1200) | ThetaCr_1200 <= 0;
ThetaCr_1200(missingTheta) = Tucr_1200(missingTheta) ./ ...
    (opt.ro .* opt.R .* opt.g .* D50critical_m(missingTheta));

missingUcr = ~isfinite(ucr_1200) | ucr_1200 <= 0;
ucr_1200(missingUcr) = sqrt(Tucr_1200(missingUcr) ./ opt.ro);

Tau_Crn = getVectorFromStruct(ShearStressOut, ...
    'Tau_Crn', Tucr_1200, nFiles);

% Method-specific critical-stress aliases. They are normally identical,
% because critical stress depends on sediment and fluid properties, not on
% the method used to estimate the applied shear.
TucrLog_1200 = getVectorFromStruct(ShearStressOut, ...
    'TucrLog_1200', Tucr_1200, nFiles);

TucrRE1_1200 = getVectorFromStruct(ShearStressOut, ...
    'TucrRE1_1200', Tucr_1200, nFiles);

TucrRE2_1200 = getVectorFromStruct(ShearStressOut, ...
    'TucrRE2_1200', Tucr_1200, nFiles);

TucrH_1200 = getVectorFromStruct(ShearStressOut, ...
    'TucrH_1200', Tucr_1200, nFiles);

Tucr1_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tucr1_1200', Tucr_1200, nFiles);

Tucr2_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tucr2_1200', Tucr_1200, nFiles);

TucrMann_1200 = getVectorFromStruct(ShearStressOut, ...
    'TucrMann_1200', Tucr_1200, nFiles);

% Excess-shear values supplied by processShearStress1200 are retained for
% diagnostics. processBedloadTrate independently evaluates T from the
% physical applied stress and the imported physical critical stress.
TvInput_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tv_1200', nan(nFiles,1), nFiles);

TvRE1Input_1200 = getVectorFromStruct(ShearStressOut, ...
    'TvRE1_1200', nan(nFiles,1), nFiles);

TvRE2Input_1200 = getVectorFromStruct(ShearStressOut, ...
    'TvRE2_1200', nan(nFiles,1), nFiles);

TvHInput_1200 = getVectorFromStruct(ShearStressOut, ...
    'TvH_1200', nan(nFiles,1), nFiles);

Tv1Input_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tv1_1200', nan(nFiles,1), nFiles);

Tv2Input_1200 = getVectorFromStruct(ShearStressOut, ...
    'Tv2_1200', nan(nFiles,1), nFiles);

TvMannInput_1200 = getVectorFromStruct(ShearStressOut, ...
    'TvMann_1200', nan(nFiles,1), nFiles);

TvSkinInput_1200 = getVectorFromStruct(ShearStressOut, ...
    'TvSkin_1200', nan(nFiles,1), nFiles);

%% Original matrix inputs for W* variables

QperDall = getMatrixFieldAny(ADCPin, {'QperDall','QperD_all','QperD'});
PSD_percentage = getMatrixFieldAny(ADCPin, ...
    {'PSD_percentage','PSDpercentage','PSD_perc'});
QperPos_14 = getMatrixFieldAny(ADCPin, {'QperPos_14','QperPos14'});
QperPos_05 = getMatrixFieldAny(ADCPin, {'QperPos_05','QperPos05'});

if isempty(QperDall)
    QperDall = nan(nFiles,1);
end
if size(QperDall,1) < nFiles
    QperDall = padRows(QperDall, nFiles);
elseif size(QperDall,1) > nFiles
    QperDall = QperDall(1:nFiles,:);
end

if isempty(PSD_percentage)
    PSD_percentage = nan(size(QperDall));
end
if size(PSD_percentage,1) < nFiles
    PSD_percentage = padRows(PSD_percentage, nFiles);
elseif size(PSD_percentage,1) > nFiles
    PSD_percentage = PSD_percentage(1:nFiles,:);
end
if size(PSD_percentage,2) < size(QperDall,2)
    PSD_percentage(:,end+1:size(QperDall,2)) = NaN;
elseif size(PSD_percentage,2) > size(QperDall,2)
    PSD_percentage = PSD_percentage(:,1:size(QperDall,2));
end

if isempty(QperPos_14)
    QperPos_14 = nan(1,nFiles);
end
if isempty(QperPos_05)
    QperPos_05 = nan(size(QperPos_14));
end
if size(QperPos_14,2) < nFiles
    QperPos_14(:,end+1:nFiles) = NaN;
end
if size(QperPos_05,2) < nFiles
    QperPos_05(:,end+1:nFiles) = NaN;
end
QperPos_14 = QperPos_14(:,1:nFiles);
QperPos_05 = QperPos_05(:,1:nFiles);

nClasses = size(QperDall,2);
nPos = max(size(QperPos_14,1), size(QperPos_05,1));
if size(QperPos_14,1) < nPos
    QperPos_14(end+1:nPos,:) = NaN;
end
if size(QperPos_05,1) < nPos
    QperPos_05(end+1:nPos,:) = NaN;
end

%% Preallocate concentration, thickness and critical diagnostics

% Dv is the Van Rijn dimensionless particle diameter and is identical to
% Dstar imported from processShearStress1200.
Dv_1200 = Dstar_1200;
DvH_1200 = Dstar_1200;
Dv1_1200 = Dstar_1200;
Dv2_1200 = Dstar_1200;
DvRE1_1200 = Dstar_1200;
DvRE2_1200 = Dstar_1200;
DvMann_1200 = Dstar_1200;
DvSkin_1200 = Dstar_1200;

% Excess shear recalculated from physical stress [Pa] and Tucr [Pa].
Tv_1200 = nan(nFiles,1);
TvRE1_1200 = nan(nFiles,1);
TvRE2_1200 = nan(nFiles,1);
TvH_1200 = nan(nFiles,1);
Tv1_1200 = nan(nFiles,1);
Tv2_1200 = nan(nFiles,1);
TvMann_1200 = nan(nFiles,1);
TvSkin_1200 = nan(nFiles,1);

% Difference between the stress-based T used here and the T supplied by
% processShearStress1200.
TvDifference_1200 = nan(nFiles,1);
TvRE1Difference_1200 = nan(nFiles,1);
TvRE2Difference_1200 = nan(nFiles,1);
TvHDifference_1200 = nan(nFiles,1);
Tv1Difference_1200 = nan(nFiles,1);
Tv2Difference_1200 = nan(nFiles,1);
TvMannDifference_1200 = nan(nFiles,1);
TvSkinDifference_1200 = nan(nFiles,1);

% Critical-input consistency diagnostics.
TucrFromUcr_1200 = opt.ro .* ucr_1200.^2;
TucrFromTheta_1200 = ThetaCr_1200 .* ...
    (opt.ros - opt.ro) .* opt.g .* D50critical_m;

CriticalRelErrUcr_1200 = relativeDifference( ...
    Tucr_1200, TucrFromUcr_1200);

CriticalRelErrTheta_1200 = relativeDifference( ...
    Tucr_1200, TucrFromTheta_1200);

CriticalInputsConsistent_1200 = ...
    CriticalRelErrUcr_1200 <= opt.CriticalConsistencyTolerance & ...
    CriticalRelErrTheta_1200 <= opt.CriticalConsistencyTolerance;

% Original concentration/thickness names
Cb_1200 = nan(nFiles,1);
CbFL_1200 = nan(nFiles,1);
Dthickness_1200 = nan(nFiles,1);
Dt_FL = nan(nFiles,1);
QvR_1200 = nan(nFiles,1);
ws = nan(nFiles,1);

CbH_1200 = nan(nFiles,1);
DthicknessH_1200 = nan(nFiles,1);
QvRH_1200 = nan(nFiles,1);

Cb1_1200 = nan(nFiles,1);
Dthickness1_1200 = nan(nFiles,1);
QvR1_1200 = nan(nFiles,1);

Cb2_1200 = nan(nFiles,1);
Dthickness2_1200 = nan(nFiles,1);
QvR2_1200 = nan(nFiles,1);

DDt1200 = nan(nFiles,7);

Tv_tau = nan(nFiles,1);
nCb_new12_all = nan(nFiles,1);
nCb_new1200 = nan(nFiles,1);

% Additional concentrations for Reynolds and Manning methods
CbRE1_1200 = nan(nFiles,1);
DthicknessRE1_1200 = nan(nFiles,1);
QvRRE1_1200 = nan(nFiles,1);

CbRE2_1200 = nan(nFiles,1);
DthicknessRE2_1200 = nan(nFiles,1);
QvRRE2_1200 = nan(nFiles,1);

CbMann_1200 = nan(nFiles,1);
DthicknessMann_1200 = nan(nFiles,1);
QvRMann_1200 = nan(nFiles,1);

CbSkin_1200 = nan(nFiles,1);
DthicknessSkin_1200 = nan(nFiles,1);
QvRSkin_1200 = nan(nFiles,1);

%% Raw concentrations and validity masks for diagnostics

Cb_raw_1200 = nan(nFiles,1);
CbH_raw_1200 = nan(nFiles,1);
Cb1_raw_1200 = nan(nFiles,1);
Cb2_raw_1200 = nan(nFiles,1);
CbRE1_raw_1200 = nan(nFiles,1);
CbRE2_raw_1200 = nan(nFiles,1);
CbMann_raw_1200 = nan(nFiles,1);
CbSkin_raw_1200 = nan(nFiles,1);

CbValid_1200 = false(nFiles,1);
CbHValid_1200 = false(nFiles,1);
Cb1Valid_1200 = false(nFiles,1);
Cb2Valid_1200 = false(nFiles,1);
CbRE1Valid_1200 = false(nFiles,1);
CbRE2Valid_1200 = false(nFiles,1);
CbMannValid_1200 = false(nFiles,1);
CbSkinValid_1200 = false(nFiles,1);
CbNewValid_1200 = false(nFiles,1);

%% Kinematic transport outputs [kg/(m s)]

qs_1200 = nan(nFiles,1);
qsRE1_1200 = nan(nFiles,1);
qsRE2_1200 = nan(nFiles,1);
qsH_1200 = nan(nFiles,1);
qs1_1200 = nan(nFiles,1);
qs2_1200 = nan(nFiles,1);
qsMann_1200 = nan(nFiles,1);
qsSkin_1200 = nan(nFiles,1);
qsFL_1200 = nan(nFiles,1);
qsCbNew_1200 = nan(nFiles,1);

% Optional transport using each method-specific Van Rijn thickness
qsVR_1200 = nan(nFiles,1);
qsVRRE1_1200 = nan(nFiles,1);
qsVRRE2_1200 = nan(nFiles,1);
qsVRH_1200 = nan(nFiles,1);
qsVR1_1200 = nan(nFiles,1);
qsVR2_1200 = nan(nFiles,1);
qsVRMann_1200 = nan(nFiles,1);
qsVRSkin_1200 = nan(nFiles,1);

%% Original W* and transport-related variables moved from shear function

W_star_perDall1200 = nan(nFiles,nClasses);
tau_Mann1200vsW = nan(nFiles,nClasses);
tau_Mann1200vsWstar = nan(nFiles,nClasses);
W_i1200 = nan(nFiles,1);
W_iVRshear2_1200 = nan(nFiles,1);
W_i1200_all = cell(nFiles,1);
MTrateMean = nan(nFiles,1);
MTrateMeanALL = cell(nFiles,1);

%% Status

Message = strings(nFiles,1);

%% Main calculation loop

for i = 1:nFiles

    try

        % Measured transport summaries do not depend on mobile D50 and are
        % retained even when the calculated kinematic transport is forced
        % to zero.
        MTrateMean(i) = safeMean( ...
            [MeanTrate14(i), MeanTrate05(i)]);

        if ~isfinite(MTrateMean(i))
            MTrateMean(i) = TrateMean(i);
        end

        MTrateMeanALL_i = safeMeanMatrixColumns( ...
            QperPos_14(:,i), QperPos_05(:,i));

        MTrateMeanALL{i} = MTrateMeanALL_i;

        if NoMobileD50_1200(i)
            % No measured mobile grain size means no calculated
            % kinematic transport for this row. Keep Cb, thickness and
            % critical-shear diagnostics as NaN; set only transport to 0.
            qs_1200(i) = 0;
            qsRE1_1200(i) = 0;
            qsRE2_1200(i) = 0;
            qsH_1200(i) = 0;
            qs1_1200(i) = 0;
            qs2_1200(i) = 0;
            qsMann_1200(i) = 0;
            qsSkin_1200(i) = 0;
            qsFL_1200(i) = 0;
            qsCbNew_1200(i) = 0;

            qsVR_1200(i) = 0;
            qsVRRE1_1200(i) = 0;
            qsVRRE2_1200(i) = 0;
            qsVRH_1200(i) = 0;
            qsVR1_1200(i) = 0;
            qsVR2_1200(i) = 0;
            qsVRMann_1200(i) = 0;
            qsVRSkin_1200(i) = 0;

            Message(i) = "OK";
            continue
        end

        D50mm = D50critical_mm(i);
        D50m = D50critical_m(i);
        DvCommon = Dstar_1200(i);

        if ~all(isfinite([D50mm,D50m,DvCommon,ThetaCr_1200(i), ...
                Tucr_1200(i),ucr_1200(i)])) || ...
                D50mm <= 0 || D50m <= 0 || DvCommon <= 0 || ...
                ThetaCr_1200(i) <= 0 || Tucr_1200(i) <= 0 || ...
                ucr_1200(i) <= 0

            error(['Invalid critical-shear input from ShearStressOut. ', ...
                'Check Dstar_1200, ThetaCr_1200, Tucr_1200 and ucr_1200.']);
        end

        if ~CriticalInputsConsistent_1200(i)
            warning(['Critical inputs are not fully consistent for row %d. ', ...
                'Relative errors: Tucr-vs-ucr %.3g, Tucr-vs-theta %.3g.'], ...
                i, CriticalRelErrUcr_1200(i), ...
                CriticalRelErrTheta_1200(i));
        end

        ws(i) = settlingVelocitySoulsby( ...
            D50mm, opt.ros, opt.ro, ni(i), opt.g);

        %% Excess shear from imported physical stresses and critical stresses

        Tv_1200(i) = resolveExcessShear( ...
            tau_All1200_Pa(i), TucrLog_1200(i), ...
            Tau_All1200(i), ThetaCr_1200(i), TvInput_1200(i));

        TvRE1_1200(i) = resolveExcessShear( ...
            tau_AllRE1_1200_Pa(i), TucrRE1_1200(i), ...
            Tau_AllRE1_1200(i), ThetaCr_1200(i), TvRE1Input_1200(i));

        TvRE2_1200(i) = resolveExcessShear( ...
            tau_AllRE2_1200_Pa(i), TucrRE2_1200(i), ...
            Tau_AllRE2_1200(i), ThetaCr_1200(i), TvRE2Input_1200(i));

        TvH_1200(i) = resolveExcessShear( ...
            tau_All1200_H_Pa(i), TucrH_1200(i), ...
            Tau_All1200_H(i), ThetaCr_1200(i), TvHInput_1200(i));

        Tv1_1200(i) = resolveExcessShear( ...
            tau_All1200_1_Pa(i), Tucr1_1200(i), ...
            Tau_All1200_1(i), ThetaCr_1200(i), Tv1Input_1200(i));

        Tv2_1200(i) = resolveExcessShear( ...
            tau_All1200_2_Pa(i), Tucr2_1200(i), ...
            Tau_All1200_2(i), ThetaCr_1200(i), Tv2Input_1200(i));

        TvMann_1200(i) = resolveExcessShear( ...
            tau_Mann1200(i), TucrMann_1200(i), ...
            tau_Mann_star1200(i), ThetaCr_1200(i), TvMannInput_1200(i));

        TvSkin_1200(i) = resolveExcessShear( ...
            tau_skin_1200(i), Tucr_1200(i), ...
            theta_skin_1200(i), ThetaCr_1200(i), TvSkinInput_1200(i));

        TvDifference_1200(i) = finiteDifference( ...
            Tv_1200(i), TvInput_1200(i));

        TvRE1Difference_1200(i) = finiteDifference( ...
            TvRE1_1200(i), TvRE1Input_1200(i));

        TvRE2Difference_1200(i) = finiteDifference( ...
            TvRE2_1200(i), TvRE2Input_1200(i));

        TvHDifference_1200(i) = finiteDifference( ...
            TvH_1200(i), TvHInput_1200(i));

        Tv1Difference_1200(i) = finiteDifference( ...
            Tv1_1200(i), Tv1Input_1200(i));

        Tv2Difference_1200(i) = finiteDifference( ...
            Tv2_1200(i), Tv2Input_1200(i));

        TvMannDifference_1200(i) = finiteDifference( ...
            TvMann_1200(i), TvMannInput_1200(i));

        TvSkinDifference_1200(i) = finiteDifference( ...
            TvSkin_1200(i), TvSkinInput_1200(i));

        %% Concentration and method-specific thickness calculations

        [Cb_raw_1200(i), Dthickness_1200(i), QvR_1200(i)] = ...
            calculateVanRijnFromExcess( ...
                D50m, Tv_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [Cb_1200(i), CbValid_1200(i)] = validateCb( ...
            Cb_raw_1200(i), opt.CbMin, opt.CbMax);

        [Dthickness_1200(i), QvR_1200(i)] = ...
            invalidateAssociated( ...
                Dthickness_1200(i), ...
                QvR_1200(i), ...
                CbValid_1200(i));

        [CbRE1_raw_1200(i), DthicknessRE1_1200(i), ...
            QvRRE1_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, TvRE1_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [CbRE1_1200(i), CbRE1Valid_1200(i)] = validateCb( ...
            CbRE1_raw_1200(i), opt.CbMin, opt.CbMax);

        [DthicknessRE1_1200(i), QvRRE1_1200(i)] = ...
            invalidateAssociated( ...
                DthicknessRE1_1200(i), ...
                QvRRE1_1200(i), ...
                CbRE1Valid_1200(i));

        [CbRE2_raw_1200(i), DthicknessRE2_1200(i), ...
            QvRRE2_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, TvRE2_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [CbRE2_1200(i), CbRE2Valid_1200(i)] = validateCb( ...
            CbRE2_raw_1200(i), opt.CbMin, opt.CbMax);

        [DthicknessRE2_1200(i), QvRRE2_1200(i)] = ...
            invalidateAssociated( ...
                DthicknessRE2_1200(i), ...
                QvRRE2_1200(i), ...
                CbRE2Valid_1200(i));

        [CbH_raw_1200(i), DthicknessH_1200(i), ...
            QvRH_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, TvH_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [CbH_1200(i), CbHValid_1200(i)] = validateCb( ...
            CbH_raw_1200(i), opt.CbMin, opt.CbMax);

        [DthicknessH_1200(i), QvRH_1200(i)] = ...
            invalidateAssociated( ...
                DthicknessH_1200(i), ...
                QvRH_1200(i), ...
                CbHValid_1200(i));

        [Cb1_raw_1200(i), Dthickness1_1200(i), ...
            QvR1_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, Tv1_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [Cb1_1200(i), Cb1Valid_1200(i)] = validateCb( ...
            Cb1_raw_1200(i), opt.CbMin, opt.CbMax);

        [Dthickness1_1200(i), QvR1_1200(i)] = ...
            invalidateAssociated( ...
                Dthickness1_1200(i), ...
                QvR1_1200(i), ...
                Cb1Valid_1200(i));

        [Cb2_raw_1200(i), Dthickness2_1200(i), ...
            QvR2_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, Tv2_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [Cb2_1200(i), Cb2Valid_1200(i)] = validateCb( ...
            Cb2_raw_1200(i), opt.CbMin, opt.CbMax);

        [Dthickness2_1200(i), QvR2_1200(i)] = ...
            invalidateAssociated( ...
                Dthickness2_1200(i), ...
                QvR2_1200(i), ...
                Cb2Valid_1200(i));

        [CbMann_raw_1200(i), DthicknessMann_1200(i), ...
            QvRMann_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, TvMann_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [CbMann_1200(i), CbMannValid_1200(i)] = validateCb( ...
            CbMann_raw_1200(i), opt.CbMin, opt.CbMax);

        [DthicknessMann_1200(i), QvRMann_1200(i)] = ...
            invalidateAssociated( ...
                DthicknessMann_1200(i), ...
                QvRMann_1200(i), ...
                CbMannValid_1200(i));

        [CbSkin_raw_1200(i), DthicknessSkin_1200(i), ...
            QvRSkin_1200(i)] = calculateVanRijnFromExcess( ...
                D50m, TvSkin_1200(i), DvCommon, ...
                opt.ros, opt.ro, opt.g);

        [CbSkin_1200(i), CbSkinValid_1200(i)] = validateCb( ...
            CbSkin_raw_1200(i), opt.CbMin, opt.CbMax);

        [DthicknessSkin_1200(i), QvRSkin_1200(i)] = ...
            invalidateAssociated( ...
                DthicknessSkin_1200(i), ...
                QvRSkin_1200(i), ...
                CbSkinValid_1200(i));

        % Original FL names remain tied to the second Van Rijn shear method.
        CbFL_1200(i) = Cb2_1200(i);
        Dt_FL(i) = Dthickness2_1200(i);

        %% Original empirical active-layer-thickness array

        [Dt1, Dt2, Dt3, Dt4, Dt5, Dt6, Dt7] = ...
            Dthickness6( ...
                D50mm, ...
                ThetaCr_1200(i), ...
                Tau_All1200_1(i), ...
                opt.ros, ...
                DvCommon, ...
                Tv1_1200(i), ...
                Ushear2_1200(i), ...
                ws(i)); %#ok<ASGLU>

        % Retain all seven empirical thickness estimates independently of
        % whether a particular Cb method is valid. Final filtering and
        % fallback are applied after all rows have been processed.
        DDt1200(i,:) = [Dt1, Dt2, Dt3, Dt4, Dt5, Dt6, Dt7];

        %% Calibrated concentration from the original workflow

        if logical(opt.ApplyNewCb) && ...
                isfinite(Tv1_1200(i)) && Tv1_1200(i) > 0 && ...
                isfinite(DvCommon) && DvCommon > 0

            Tv_tau(i) = Tv1_1200(i);

            cc1 = [ ...
                7.26384427883237, ...
               -0.142516293252098, ...
                0.747997898125074, ...
               -0.164224531344291];

            nCb_new12_all(i) = ...
                cc1(1) .* ...
                (Tv_tau(i).^cc1(2) ./ ...
                DvCommon.^cc1(3)) + ...
                cc1(4);

            [nCb_new1200(i), CbNewValid_1200(i)] = ...
                validateCb( ...
                    nCb_new12_all(i), ...
                    opt.CbMin, ...
                    opt.CbMax);
        end

        %% Kinematic transport using n*D50 thickness

        qs_1200(i) = kinematicTransport( ...
            D50m, n(i), Cb_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsRE1_1200(i) = kinematicTransport( ...
            D50m, n(i), CbRE1_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsRE2_1200(i) = kinematicTransport( ...
            D50m, n(i), CbRE2_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsH_1200(i) = kinematicTransport( ...
            D50m, n(i), CbH_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qs1_1200(i) = kinematicTransport( ...
            D50m, n(i), Cb1_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qs2_1200(i) = kinematicTransport( ...
            D50m, n(i), Cb2_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsMann_1200(i) = kinematicTransport( ...
            D50m, n(i), CbMann_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsSkin_1200(i) = kinematicTransport( ...
            D50m, n(i), CbSkin_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsFL_1200(i) = kinematicTransport( ...
            D50m, n(i), CbFL_1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        qsCbNew_1200(i) = kinematicTransport( ...
            D50m, n(i), nCb_new1200(i), ...
            BTmodf_1200(i), opt.ros, Ninstrument(i));

        %% Transport using each method-specific Van Rijn thickness

        qsVR_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), Cb_1200(i), ...
            Dthickness_1200(i), opt.ros, Ninstrument(i));

        qsVRRE1_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), CbRE1_1200(i), ...
            DthicknessRE1_1200(i), opt.ros, Ninstrument(i));

        qsVRRE2_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), CbRE2_1200(i), ...
            DthicknessRE2_1200(i), opt.ros, Ninstrument(i));

        qsVRH_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), CbH_1200(i), ...
            DthicknessH_1200(i), opt.ros, Ninstrument(i));

        qsVR1_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), Cb1_1200(i), ...
            Dthickness1_1200(i), opt.ros, Ninstrument(i));

        qsVR2_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), Cb2_1200(i), ...
            Dthickness2_1200(i), opt.ros, Ninstrument(i));

        qsVRMann_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), CbMann_1200(i), ...
            DthicknessMann_1200(i), opt.ros, Ninstrument(i));

        qsVRSkin_1200(i) = transportWithThickness( ...
            BTmodf_1200(i), CbSkin_1200(i), ...
            DthicknessSkin_1200(i), opt.ros, Ninstrument(i));

        %% Original W* variables retained with transport processing

        if isfinite(uShear_Mann1200(i)) && ...
                uShear_Mann1200(i) > 0

            W_star_perDall1200(i,:) = ...
                (opt.R .* opt.g .* ...
                QperDall(i,:) ./ opt.ro) ./ ...
                (opt.ros .* ...
                PSD_percentage(i,:) .* ...
                uShear_Mann1200(i).^3);

            W_i1200(i) = ...
                (opt.R .* opt.g .* ...
                MTrateMean(i) ./ opt.ro) ./ ...
                (opt.ros .* ...
                uShear_Mann1200(i).^3);

            W_i1200_all{i} = ...
                (opt.R .* opt.g .* ...
                MTrateMeanALL_i ./ opt.ro) ./ ...
                (opt.ros .* ...
                uShear_Mann1200(i).^3);

            tau_Mann1200vsW(i,:) = ...
                tau_Mann1200(i) ./ ...
                W_star_perDall1200(i,:);

            if isfinite(D50_05(i)) && D50_05(i) > 0

                tau_Mann1200vsWstar(i,:) = ...
                    tau_Mann1200vsW(i,:) ./ ...
                    (opt.R .* opt.g .* ...
                    (D50_05(i) ./ 1000));
            end
        end

        if isfinite(Ushear2_1200(i)) && ...
                Ushear2_1200(i) > 0

            W_iVRshear2_1200(i) = ...
                (opt.R .* opt.g .* ...
                MTrateMean(i) ./ opt.ro) ./ ...
                (opt.ros .* ...
                Ushear2_1200(i).^3);
        end

        Message(i) = "OK";

    catch ME

        Message(i) = "Error: " + string(ME.message);

        warning( ...
            'Bedload concentration/transport failed for row %d, file %s:\n%s', ...
            i, ADCPfilename(i), ME.message);
    end
end

%% Final Cb from all available shear methods

CbMethodNames_1200 = [ ...
    "LogLaw", ...
    "Reynolds1", ...
    "Reynolds2", ...
    "HydraulicChezy", ...
    "VanRijn1", ...
    "VanRijn2", ...
    "Manning", ...
    "SkinFriction"];

CbCandidates_1200 = [ ...
    Cb_1200, ...
    CbRE1_1200, ...
    CbRE2_1200, ...
    CbH_1200, ...
    Cb1_1200, ...
    Cb2_1200, ...
    CbMann_1200, ...
    CbSkin_1200];

if logical(opt.IncludeCalibratedCbInFinal)
    CbMethodNames_1200(end+1) = "CalibratedCb";
    CbCandidates_1200(:,end+1) = nCb_new1200;
end

[CbFinal_1200, CbFinalMethodCount_1200, ...
    CbFinalSource_1200, CbFinalUsedFallback_1200] = ...
    aggregateBoundedCandidates( ...
        CbCandidates_1200, ...
        opt.CbMin .* ones(nFiles,1), ...
        opt.CbMax .* ones(nFiles,1), ...
        opt.CbFinalStatistic, ...
        logical(opt.UseCbFallback), ...
        opt.CbFallback .* ones(nFiles,1));

% Alias retained for workflows that use an underscore before "final".
Cb_final_1200 = CbFinal_1200;

CbFinalFromMethods_1200 = CbFinalMethodCount_1200 > 0;
CbFinalValid_1200 = isfinite(CbFinal_1200) & ...
    CbFinal_1200 >= opt.CbMin & CbFinal_1200 <= opt.CbMax;

%% Final active-layer thickness from seven empirical methods

DthicknessMethodNames_1200 = [ ...
    "Dt1", "Dt2", "Dt3", "Dt4", "Dt5", "Dt6", "Dt7"];

DthicknessCandidates_1200 = DDt1200;
DthicknessMin_1200 = opt.DthicknessMinFactor .* D50critical_m;
DthicknessMax_1200 = opt.DthicknessMaxFactor .* D50critical_m;
DthicknessFallback_1200 = ...
    opt.DthicknessFallbackFactor .* D50critical_m;

[DthicknessFinal_1200, DthicknessFinalMethodCount_1200, ...
    DthicknessFinalSource_1200, DthicknessFinalUsedFallback_1200] = ...
    aggregateBoundedCandidates( ...
        DthicknessCandidates_1200, ...
        DthicknessMin_1200, ...
        DthicknessMax_1200, ...
        opt.DthicknessFinalStatistic, ...
        logical(opt.UseDthicknessFallback), ...
        DthicknessFallback_1200);

Dthickness_final_1200 = DthicknessFinal_1200;

DthicknessFinalFromMethods_1200 = ...
    DthicknessFinalMethodCount_1200 > 0;

DthicknessFinalValid_1200 = ...
    isfinite(DthicknessFinal_1200) & ...
    DthicknessFinal_1200 >= DthicknessMin_1200 & ...
    DthicknessFinal_1200 <= DthicknessMax_1200;

%% Additional final kinematic transport using final Cb and final thickness
%
% IMPORTANT:
% The original method-specific kinematic formulations above retain
% Dthickness = n*D50mobile and are all multiplied by Ninstrument.
%
% This is one additional final value using:
%   - CbFinal_1200, already restricted to [CbMin,CbMax] or its fallback;
%   - DthicknessFinal_1200, already restricted to the requested D50 range
%     or its fallback.

qsKinematicFinal_1200 = nan(nFiles,1);

for i = 1:nFiles
    if NoMobileD50_1200(i)
        qsKinematicFinal_1200(i) = 0;
    else
        qsKinematicFinal_1200(i) = kinematicTransportFinal( ...
            CbFinal_1200(i), ...
            DthicknessFinal_1200(i), ...
            BTmodf_1200(i), ...
            opt.ros, ...
            Ninstrument(i));
    end
end

qsKinematicFinalValid_1200 = ...
    NoMobileD50_1200 | ( ...
        isfinite(qsKinematicFinal_1200) & ...
        qsKinematicFinal_1200 >= 0 & ...
        CbFinalValid_1200 & ...
        DthicknessFinalValid_1200);

% Compatibility and descriptive aliases.
qsFinalKinematic_1200 = qsKinematicFinal_1200;
qsFinal_1200 = qsKinematicFinal_1200;
qs_final_1200 = qsKinematicFinal_1200;

%% Final calculated-versus-measured transport comparison
% Convert measured rates to the same unit as qsKinematicFinal_1200.
TrateMeanComparison_1200 = ...
    TrateMean .* opt.MeasuredTrateToKgFactor;
MaxTrateComparison_1200 = ...
    MaxTrate .* opt.MeasuredTrateToKgFactor;

FinalTrateComparisonTable = makeFinalTrateComparisonTable( ...
    qsKinematicFinal_1200, ...
    TrateMeanComparison_1200, ...
    MaxTrateComparison_1200, ...
    opt.PlotFinalTrateMinN);

if logical(opt.PlotFinalTrateComparison)
    plotFinalTrateComparison( ...
        qsKinematicFinal_1200, ...
        TrateMeanComparison_1200, ...
        MaxTrateComparison_1200, ...
        logical(opt.PlotFinalTrateUseLogAxes), ...
        opt.PlotFinalTrateMinN);
end

%% Output table

BedloadTransportTable = table( ...
    InputRow, ...
    ADCPfilename, ...
    DistanceFromLeftBank, ...
    D50mobile, ...
    D50critical_mm, ...
    D50critical_m, ...
    NoMobileD50_1200, ...
    KinematicZeroDueToNoMobileD50_1200, ...
    ni, ...
    n, ...
    Ninstrument, ...
    DthicknessKin_1200, ...
    BTmodf_1200, ...
    TrateMean, ...
    MaxTrate, ...
    TrateMeanComparison_1200, ...
    MaxTrateComparison_1200, ...
    Dstar_1200, ...
    ThetaCr_1200, ...
    Tucr_1200, ...
    ucr_1200, ...
    CriticalInputsConsistent_1200, ...
    Tv_1200, ...
    TvRE1_1200, ...
    TvRE2_1200, ...
    TvH_1200, ...
    Tv1_1200, ...
    Tv2_1200, ...
    TvMann_1200, ...
    TvSkin_1200, ...
    Cb_1200, ...
    CbRE1_1200, ...
    CbRE2_1200, ...
    CbH_1200, ...
    Cb1_1200, ...
    Cb2_1200, ...
    CbMann_1200, ...
    CbSkin_1200, ...
    nCb_new1200, ...
    CbFinal_1200, ...
    CbFinalMethodCount_1200, ...
    CbFinalSource_1200, ...
    DthicknessFinal_1200, ...
    DthicknessFinalMethodCount_1200, ...
    DthicknessFinalSource_1200, ...
    qsKinematicFinal_1200, ...
    qs_1200, ...
    qsRE1_1200, ...
    qsRE2_1200, ...
    qsH_1200, ...
    qs1_1200, ...
    qs2_1200, ...
    qsMann_1200, ...
    qsSkin_1200, ...
    qsCbNew_1200, ...
    Message);

% Convert per-position series to one cell per measurement before export.
% Fixed-size grain-class matrices remain nFiles-by-nClasses summaries.
QperPos_14_cells = cell(nFiles,1);
QperPos_05_cells = cell(nFiles,1);
for iSeries = 1:nFiles
    QperPos_14_cells{iSeries} = QperPos_14(:,iSeries);
    QperPos_05_cells{iSeries} = QperPos_05(:,iSeries);
end

%% Build output structure

BedloadTransportOut = struct();

BedloadTransportOut.InputRow = InputRow;
BedloadTransportOut.ADCPfilename = ADCPfilename;
BedloadTransportOut.DistanceFromLeftBank = DistanceFromLeftBank;

BedloadTransportOut.TrateMean = TrateMean;
BedloadTransportOut.MaxTrate = MaxTrate;
BedloadTransportOut.TrateMeanComparison_1200 = ...
    TrateMeanComparison_1200;
BedloadTransportOut.MaxTrateComparison_1200 = ...
    MaxTrateComparison_1200;
BedloadTransportOut.MeasuredTrateToKgFactor = ...
    opt.MeasuredTrateToKgFactor;
BedloadTransportOut.PlotFinalTrateComparison = ...
    logical(opt.PlotFinalTrateComparison);
BedloadTransportOut.PlotFinalTrateUseLogAxes = ...
    logical(opt.PlotFinalTrateUseLogAxes);
BedloadTransportOut.PlotFinalTrateMinN = ...
    opt.PlotFinalTrateMinN;
BedloadTransportOut.MeanTrate14 = MeanTrate14;
BedloadTransportOut.MeanTrate05 = MeanTrate05;

% Required grab-sample D50 and retained aliases
BedloadTransportOut.D_50grabM = D_50grabM;
BedloadTransportOut.D50critical_mm = D50critical_mm;
BedloadTransportOut.D50critical_m = D50critical_m;
BedloadTransportOut.D50mobile = D50mobile;
BedloadTransportOut.D50mobileM = D50mobileM;
BedloadTransportOut.D50mobileInput = D50mobileInput;
BedloadTransportOut.D50mobileShear = D50mobileShear;
BedloadTransportOut.NoMobileD50_1200 = NoMobileD50_1200;
BedloadTransportOut.KinematicZeroDueToNoMobileD50_1200 = ...
    KinematicZeroDueToNoMobileD50_1200;
BedloadTransportOut.D50_14 = D50_14;
BedloadTransportOut.D50_05 = D50_05;
BedloadTransportOut.D50 = D50;

BedloadTransportOut.ni = ni;
BedloadTransportOut.n = n;
BedloadTransportOut.Ninstrument = Ninstrument;
BedloadTransportOut.CbMin = opt.CbMin;
BedloadTransportOut.CbMax = opt.CbMax;
BedloadTransportOut.CbFinalStatistic = opt.CbFinalStatistic;
BedloadTransportOut.UseCbFallback = logical(opt.UseCbFallback);
BedloadTransportOut.CbFallback = opt.CbFallback;
BedloadTransportOut.IncludeCalibratedCbInFinal = ...
    logical(opt.IncludeCalibratedCbInFinal);

BedloadTransportOut.DthicknessMinFactor = ...
    opt.DthicknessMinFactor;
BedloadTransportOut.DthicknessMaxFactor = ...
    opt.DthicknessMaxFactor;
BedloadTransportOut.DthicknessFinalStatistic = ...
    opt.DthicknessFinalStatistic;
BedloadTransportOut.UseDthicknessFallback = ...
    logical(opt.UseDthicknessFallback);
BedloadTransportOut.DthicknessFallbackFactor = ...
    opt.DthicknessFallbackFactor;

BedloadTransportOut.CriticalConsistencyTolerance = ...
    opt.CriticalConsistencyTolerance;

BedloadTransportOut.DthicknessKin_1200 = ...
    DthicknessKin_1200;

BedloadTransportOut.BTmodf_1200 = BTmodf_1200;

% Critical inputs imported from processShearStress1200
BedloadTransportOut.Dstar_1200 = Dstar_1200;
BedloadTransportOut.ThetaCr_1200 = ThetaCr_1200;
BedloadTransportOut.Tucr_1200 = Tucr_1200;
BedloadTransportOut.Tau_Crn = Tau_Crn;
BedloadTransportOut.ucr_1200 = ucr_1200;

BedloadTransportOut.TucrLog_1200 = TucrLog_1200;
BedloadTransportOut.TucrRE1_1200 = TucrRE1_1200;
BedloadTransportOut.TucrRE2_1200 = TucrRE2_1200;
BedloadTransportOut.TucrH_1200 = TucrH_1200;
BedloadTransportOut.Tucr1_1200 = Tucr1_1200;
BedloadTransportOut.Tucr2_1200 = Tucr2_1200;
BedloadTransportOut.TucrMann_1200 = TucrMann_1200;

BedloadTransportOut.TucrFromUcr_1200 = TucrFromUcr_1200;
BedloadTransportOut.TucrFromTheta_1200 = TucrFromTheta_1200;
BedloadTransportOut.CriticalRelErrUcr_1200 = ...
    CriticalRelErrUcr_1200;
BedloadTransportOut.CriticalRelErrTheta_1200 = ...
    CriticalRelErrTheta_1200;
BedloadTransportOut.CriticalInputsConsistent_1200 = ...
    CriticalInputsConsistent_1200;

% Imported shear velocities
BedloadTransportOut.UsheaAll_1200 = UsheaAll_1200;
BedloadTransportOut.U_shearRE1_1200 = U_shearRE1_1200;
BedloadTransportOut.U_shearRE2_1200 = U_shearRE2_1200;
BedloadTransportOut.u_shearH1200 = u_shearH1200;
BedloadTransportOut.Ushear1_1200 = Ushear1_1200;
BedloadTransportOut.Ushear2_1200 = Ushear2_1200;
BedloadTransportOut.uShear_Mann1200 = uShear_Mann1200;
BedloadTransportOut.uShear_skin1200 = uShear_skin1200;
BedloadTransportOut.theta_skin_1200 = theta_skin_1200;
BedloadTransportOut.tau_skin_1200 = tau_skin_1200;

% Imported normalized Shields stresses
BedloadTransportOut.Tau_All1200 = Tau_All1200;
BedloadTransportOut.Tau_AllRE1_1200 = Tau_AllRE1_1200;
BedloadTransportOut.Tau_AllRE2_1200 = Tau_AllRE2_1200;
BedloadTransportOut.Tau_All1200_H = Tau_All1200_H;
BedloadTransportOut.Tau_All1200_1 = Tau_All1200_1;
BedloadTransportOut.Tau_All1200_2 = Tau_All1200_2;
BedloadTransportOut.tau_Mann_star1200 = tau_Mann_star1200;

% Imported physical stresses [Pa]
BedloadTransportOut.tau_All1200_Pa = tau_All1200_Pa;
BedloadTransportOut.tau_AllRE1_1200_Pa = tau_AllRE1_1200_Pa;
BedloadTransportOut.tau_AllRE2_1200_Pa = tau_AllRE2_1200_Pa;
BedloadTransportOut.tau_All1200_H_Pa = tau_All1200_H_Pa;
BedloadTransportOut.tau_All1200_1_Pa = tau_All1200_1_Pa;
BedloadTransportOut.tau_All1200_2_Pa = tau_All1200_2_Pa;
BedloadTransportOut.tau_Mann1200 = tau_Mann1200;

% Excess shear used for Cb
BedloadTransportOut.Tv_1200 = Tv_1200;
BedloadTransportOut.TvRE1_1200 = TvRE1_1200;
BedloadTransportOut.TvRE2_1200 = TvRE2_1200;
BedloadTransportOut.TvH_1200 = TvH_1200;
BedloadTransportOut.Tv1_1200 = Tv1_1200;
BedloadTransportOut.Tv2_1200 = Tv2_1200;
BedloadTransportOut.TvMann_1200 = TvMann_1200;
BedloadTransportOut.TvSkin_1200 = TvSkin_1200;

% Excess shear supplied by processShearStress1200
BedloadTransportOut.TvInput_1200 = TvInput_1200;
BedloadTransportOut.TvRE1Input_1200 = TvRE1Input_1200;
BedloadTransportOut.TvRE2Input_1200 = TvRE2Input_1200;
BedloadTransportOut.TvHInput_1200 = TvHInput_1200;
BedloadTransportOut.Tv1Input_1200 = Tv1Input_1200;
BedloadTransportOut.Tv2Input_1200 = Tv2Input_1200;
BedloadTransportOut.TvMannInput_1200 = TvMannInput_1200;
BedloadTransportOut.TvSkinInput_1200 = TvSkinInput_1200;

% Differences between stress-based and supplied excess shear
BedloadTransportOut.TvDifference_1200 = TvDifference_1200;
BedloadTransportOut.TvRE1Difference_1200 = ...
    TvRE1Difference_1200;
BedloadTransportOut.TvRE2Difference_1200 = ...
    TvRE2Difference_1200;
BedloadTransportOut.TvHDifference_1200 = TvHDifference_1200;
BedloadTransportOut.Tv1Difference_1200 = Tv1Difference_1200;
BedloadTransportOut.Tv2Difference_1200 = Tv2Difference_1200;
BedloadTransportOut.TvMannDifference_1200 = ...
    TvMannDifference_1200;
BedloadTransportOut.TvSkinDifference_1200 = ...
    TvSkinDifference_1200;

% Original dimensionless grain-size aliases
BedloadTransportOut.Dv_1200 = Dv_1200;
BedloadTransportOut.DvH_1200 = DvH_1200;
BedloadTransportOut.Dv1_1200 = Dv1_1200;
BedloadTransportOut.Dv2_1200 = Dv2_1200;
BedloadTransportOut.DvRE1_1200 = DvRE1_1200;
BedloadTransportOut.DvRE2_1200 = DvRE2_1200;
BedloadTransportOut.DvMann_1200 = DvMann_1200;
BedloadTransportOut.DvSkin_1200 = DvSkin_1200;

% Concentration, thickness and Van Rijn transport
BedloadTransportOut.Cb_1200 = Cb_1200;
BedloadTransportOut.CbFL_1200 = CbFL_1200;
BedloadTransportOut.Dthickness_1200 = Dthickness_1200;
BedloadTransportOut.Dt_FL = Dt_FL;
BedloadTransportOut.QvR_1200 = QvR_1200;
BedloadTransportOut.ws = ws;

BedloadTransportOut.CbH_1200 = CbH_1200;
BedloadTransportOut.DthicknessH_1200 = DthicknessH_1200;
BedloadTransportOut.QvRH_1200 = QvRH_1200;

BedloadTransportOut.Cb1_1200 = Cb1_1200;
BedloadTransportOut.Dthickness1_1200 = Dthickness1_1200;
BedloadTransportOut.QvR1_1200 = QvR1_1200;

BedloadTransportOut.Cb2_1200 = Cb2_1200;
BedloadTransportOut.Dthickness2_1200 = Dthickness2_1200;
BedloadTransportOut.QvR2_1200 = QvR2_1200;

BedloadTransportOut.CbRE1_1200 = CbRE1_1200;
BedloadTransportOut.DthicknessRE1_1200 = ...
    DthicknessRE1_1200;
BedloadTransportOut.QvRRE1_1200 = QvRRE1_1200;

BedloadTransportOut.CbRE2_1200 = CbRE2_1200;
BedloadTransportOut.DthicknessRE2_1200 = ...
    DthicknessRE2_1200;
BedloadTransportOut.QvRRE2_1200 = QvRRE2_1200;

BedloadTransportOut.CbMann_1200 = CbMann_1200;
BedloadTransportOut.DthicknessMann_1200 = ...
    DthicknessMann_1200;
BedloadTransportOut.QvRMann_1200 = QvRMann_1200;

BedloadTransportOut.CbSkin_1200 = CbSkin_1200;
BedloadTransportOut.DthicknessSkin_1200 = ...
    DthicknessSkin_1200;
BedloadTransportOut.QvRSkin_1200 = QvRSkin_1200;

% Final concentration from all available shear methods
BedloadTransportOut.CbMethodNames_1200 = CbMethodNames_1200;
BedloadTransportOut.CbCandidates_1200 = CbCandidates_1200;
BedloadTransportOut.CbFinal_1200 = CbFinal_1200;
BedloadTransportOut.Cb_final_1200 = Cb_final_1200;
BedloadTransportOut.CbFinalMethodCount_1200 = ...
    CbFinalMethodCount_1200;
BedloadTransportOut.CbFinalSource_1200 = CbFinalSource_1200;
BedloadTransportOut.CbFinalUsedFallback_1200 = ...
    CbFinalUsedFallback_1200;
BedloadTransportOut.CbFinalFromMethods_1200 = ...
    CbFinalFromMethods_1200;
BedloadTransportOut.CbFinalValid_1200 = CbFinalValid_1200;

% Final thickness from the seven empirical methods
BedloadTransportOut.DthicknessMethodNames_1200 = ...
    DthicknessMethodNames_1200;
BedloadTransportOut.DthicknessCandidates_1200 = ...
    DthicknessCandidates_1200;
BedloadTransportOut.DthicknessMin_1200 = ...
    DthicknessMin_1200;
BedloadTransportOut.DthicknessMax_1200 = ...
    DthicknessMax_1200;
BedloadTransportOut.DthicknessFallback_1200 = ...
    DthicknessFallback_1200;
BedloadTransportOut.DthicknessFinal_1200 = ...
    DthicknessFinal_1200;
BedloadTransportOut.Dthickness_final_1200 = ...
    Dthickness_final_1200;
BedloadTransportOut.DthicknessFinalMethodCount_1200 = ...
    DthicknessFinalMethodCount_1200;
BedloadTransportOut.DthicknessFinalSource_1200 = ...
    DthicknessFinalSource_1200;
BedloadTransportOut.DthicknessFinalUsedFallback_1200 = ...
    DthicknessFinalUsedFallback_1200;
BedloadTransportOut.DthicknessFinalFromMethods_1200 = ...
    DthicknessFinalFromMethods_1200;
BedloadTransportOut.DthicknessFinalValid_1200 = ...
    DthicknessFinalValid_1200;

BedloadTransportOut.DDt1200 = DDt1200;
BedloadTransportOut.Tv_tau = Tv_tau;
BedloadTransportOut.nCb_new12_all = nCb_new12_all;
BedloadTransportOut.nCb_new1200 = nCb_new1200;

% Raw concentrations and validity masks
BedloadTransportOut.Cb_raw_1200 = Cb_raw_1200;
BedloadTransportOut.CbH_raw_1200 = CbH_raw_1200;
BedloadTransportOut.Cb1_raw_1200 = Cb1_raw_1200;
BedloadTransportOut.Cb2_raw_1200 = Cb2_raw_1200;
BedloadTransportOut.CbRE1_raw_1200 = CbRE1_raw_1200;
BedloadTransportOut.CbRE2_raw_1200 = CbRE2_raw_1200;
BedloadTransportOut.CbMann_raw_1200 = CbMann_raw_1200;
BedloadTransportOut.CbSkin_raw_1200 = CbSkin_raw_1200;

BedloadTransportOut.CbValid_1200 = CbValid_1200;
BedloadTransportOut.CbHValid_1200 = CbHValid_1200;
BedloadTransportOut.Cb1Valid_1200 = Cb1Valid_1200;
BedloadTransportOut.Cb2Valid_1200 = Cb2Valid_1200;
BedloadTransportOut.CbRE1Valid_1200 = CbRE1Valid_1200;
BedloadTransportOut.CbRE2Valid_1200 = CbRE2Valid_1200;
BedloadTransportOut.CbMannValid_1200 = CbMannValid_1200;
BedloadTransportOut.CbSkinValid_1200 = CbSkinValid_1200;
BedloadTransportOut.CbNewValid_1200 = CbNewValid_1200;

% Kinematic transport with n*D50 thickness
BedloadTransportOut.qs_1200 = qs_1200;
BedloadTransportOut.qsRE1_1200 = qsRE1_1200;
BedloadTransportOut.qsRE2_1200 = qsRE2_1200;
BedloadTransportOut.qsH_1200 = qsH_1200;
BedloadTransportOut.qs1_1200 = qs1_1200;
BedloadTransportOut.qs2_1200 = qs2_1200;
BedloadTransportOut.qsMann_1200 = qsMann_1200;
BedloadTransportOut.qsSkin_1200 = qsSkin_1200;
BedloadTransportOut.qsFL_1200 = qsFL_1200;
BedloadTransportOut.qsCbNew_1200 = qsCbNew_1200;

% Transport with method-specific Van Rijn thickness
BedloadTransportOut.qsVR_1200 = qsVR_1200;
BedloadTransportOut.qsVRRE1_1200 = qsVRRE1_1200;
BedloadTransportOut.qsVRRE2_1200 = qsVRRE2_1200;
BedloadTransportOut.qsVRH_1200 = qsVRH_1200;
BedloadTransportOut.qsVR1_1200 = qsVR1_1200;
BedloadTransportOut.qsVR2_1200 = qsVR2_1200;
BedloadTransportOut.qsVRMann_1200 = qsVRMann_1200;
BedloadTransportOut.qsVRSkin_1200 = qsVRSkin_1200;

% Additional final kinematic transport using final bounded inputs.
% The original qs_1200, qsRE1_1200, ... output names are retained and
% include the instrument-specific Ninstrument multiplier.
BedloadTransportOut.qsKinematicFinal_1200 = ...
    qsKinematicFinal_1200;
BedloadTransportOut.qsFinalKinematic_1200 = ...
    qsFinalKinematic_1200;
BedloadTransportOut.qsKinematicFinalValid_1200 = ...
    qsKinematicFinalValid_1200;

% Compatibility aliases
BedloadTransportOut.qsFinal_1200 = qsFinal_1200;
BedloadTransportOut.qs_final_1200 = qs_final_1200;

% Final calculated-versus-measured comparison diagnostics.
BedloadTransportOut.FinalTrateComparisonTable = ...
    FinalTrateComparisonTable;

% Original W* variables
BedloadTransportOut.W_star_perDall1200 = ...
    W_star_perDall1200;
BedloadTransportOut.W_i1200 = W_i1200;
BedloadTransportOut.W_iVRshear2_1200 = ...
    W_iVRshear2_1200;
BedloadTransportOut.W_i1200_all = W_i1200_all;
BedloadTransportOut.tau_Mann1200vsW = ...
    tau_Mann1200vsW;
BedloadTransportOut.tau_Mann1200vsWstar = ...
    tau_Mann1200vsWstar;
BedloadTransportOut.MTrateMean = MTrateMean;
BedloadTransportOut.MTrateMeanALL = MTrateMeanALL;
BedloadTransportOut.QperDall = QperDall;
BedloadTransportOut.PSD_percentage = PSD_percentage;
BedloadTransportOut.QperPos_14 = QperPos_14_cells;
BedloadTransportOut.QperPos_05 = QperPos_05_cells;

BedloadTransportOut.Message = Message;
BedloadTransportOut.BedloadTransportTable = ...
    BedloadTransportTable;

%% Assign to workspace

if logical(opt.AssignToCaller)
    assignStructToWorkspace(BedloadTransportOut, 'caller');
end

if logical(opt.AssignToBase)
    assignStructToWorkspace(BedloadTransportOut, 'base');
end

%% Display

if logical(opt.Verbose)

    fprintf('\n=== Bedload concentration and transport summary ===\n');
    fprintf('Input rows: %d\n', nFiles);
    fprintf('Successfully processed: %d/%d\n', ...
        sum(Message == "OK"), nFiles);
    fprintf(['Rows with missing/zero mobile D50 and forced ', ...
        'kinematic transport = 0: %d/%d.\n'], ...
        sum(KinematicZeroDueToNoMobileD50_1200), nFiles);
    fprintf('Critical shear source: ShearStressOut.\n');
    fprintf('Critical D50 source: D50mobile [mm], preferably imported from ShearStressOut.\n');
    fprintf('Cb validity range: %.4f <= Cb <= %.4f (inclusive).\n', ...
        opt.CbMin, opt.CbMax);
    fprintf('Final Cb statistic: %s. Fallback enabled: %d; value: %.4f.\n', ...
        char(opt.CbFinalStatistic), logical(opt.UseCbFallback), ...
        opt.CbFallback);
    fprintf('Final Cb from methods: %d/%d rows; fallback used: %d rows.\n', ...
        sum(CbFinalFromMethods_1200), nFiles, ...
        sum(CbFinalUsedFallback_1200));
    fprintf(['Final thickness range: %.2f*D50 to %.2f*D50; ', ...
        'statistic: %s.\n'], ...
        opt.DthicknessMinFactor, opt.DthicknessMaxFactor, ...
        char(opt.DthicknessFinalStatistic));
    fprintf('Final thickness fallback enabled: %d; factor: %.2f; used: %d rows.\n', ...
        logical(opt.UseDthicknessFallback), ...
        opt.DthicknessFallbackFactor, ...
        sum(DthicknessFinalUsedFallback_1200));
    fprintf('Consistent critical-input rows: %d/%d.\n', ...
        sum(CriticalInputsConsistent_1200), nFiles);
    fprintf(['Original method-specific kinematic outputs retain ', ...
        'Dthickness = n*D50mobile and include Ninstrument.\n']);
    fprintf('Ninstrument range: %.6g to %.6g; default/unity rows: %d/%d.\n', ...
        min(Ninstrument), max(Ninstrument), ...
        sum(abs(Ninstrument-1) <= eps(max(1,max(abs(Ninstrument))))), nFiles);
    fprintf(['Additional final kinematic output: qsKinematicFinal_1200 ', ...
        '= Ninstrument*BTmodf*CbFinal*DthicknessFinal*ros.\n']);
    fprintf('Valid final kinematic rows: %d/%d.\n', ...
        sum(qsKinematicFinalValid_1200), nFiles);
    fprintf(['Measured transport conversion to kg/(m s): factor %.6g; ', ...
        'final comparison plot=%d.\n'], ...
        opt.MeasuredTrateToKgFactor, ...
        logical(opt.PlotFinalTrateComparison));
    fprintf('Per-position series are stored only in per-measurement cells.\n');
    fprintf('No temporary value from the final processed measurement was exported.\n');
    fprintf('Kinematic transport unit: kg/(m s).\n');
    fprintf('No files were saved.\n');
    fprintf('===================================================\n');

    disp(BedloadTransportTable);
end

end

%% Helper functions

function x = getVectorFieldAny(S, fieldNames, defaultValue, n)

x = defaultValue;

for k = 1:numel(fieldNames)
    if isfield(S, fieldNames{k})
        x = S.(fieldNames{k});
        break
    end
end

x = prepareVector(x, n);

end


function x = getVectorFromStruct(S, fieldName, defaultValue, n)

if isfield(S, fieldName)
    x = S.(fieldName);
else
    x = defaultValue;
end

x = prepareVector(x, n);

end


function x = requireVectorFromStruct(S, fieldName, n)

if ~isfield(S, fieldName)
    error('ShearStressOut.%s is missing.', fieldName);
end

x = prepareVector(S.(fieldName), n);

end


function x = prepareOptionalVector(x, n)

if isempty(x)
    x = nan(n,1);
else
    x = prepareVector(x, n);
end

end


function x = prepareVector(x, n)

if isstring(x) || iscell(x)
    x = str2double(string(x));
end

x = double(x);
x = x(:);

if isempty(x)
    x = nan(n,1);
elseif numel(x) == 1 && n > 1
    x = repmat(x, n, 1);
elseif numel(x) < n
    x = [x; nan(n-numel(x),1)];
elseif numel(x) > n
    x = x(1:n);
end

end


function M = getMatrixFieldAny(S, fieldNames)

M = [];

for k = 1:numel(fieldNames)
    if isfield(S, fieldNames{k})
        M = double(S.(fieldNames{k}));
        return
    end
end

end


function M = padRows(M, nRows)

if isempty(M)
    M = nan(nRows,1);
elseif size(M,1) < nRows
    M(end+1:nRows,:) = NaN;
end

end


function BTmodf_1200 = getBTmodfVector(BTbedloadOut, nFiles)

BTmodf_1200 = getVectorFromStruct(BTbedloadOut, ...
    'BTmodf_1200', nan(nFiles,1), nFiles);

if all(isnan(BTmodf_1200)) && ...
        isfield(BTbedloadOut, 'va_final_1200')

    BTmodf_1200 = prepareVector( ...
        BTbedloadOut.va_final_1200, nFiles);
end

if all(isnan(BTmodf_1200)) && ...
        isfield(BTbedloadOut, 'BTmodf_1200ALL')

    C = BTbedloadOut.BTmodf_1200ALL;

    for i = 1:min(nFiles,numel(C))
        if iscell(C) && ~isempty(C{i})
            BTmodf_1200(i) = safeMean(C{i});
        end
    end
end

if all(isnan(BTmodf_1200)) && ...
        isfield(BTbedloadOut, 'BTmod_final_1200')

    C = BTbedloadOut.BTmod_final_1200;

    for i = 1:min(nFiles,numel(C))
        if iscell(C) && ~isempty(C{i})
            BTmodf_1200(i) = safeMean(C{i});
        end
    end
end

end


function T = resolveExcessShear( ...
    tauPa, tauCrPa, theta, thetaCr, Tinput)
% Priority:
%   1. Physical stress comparison [Pa]/[Pa].
%   2. Dimensionless Shields comparison [-]/[-].
%   3. Excess shear supplied by processShearStress1200.

T = NaN;

if all(isfinite([tauPa,tauCrPa])) && tauCrPa > 0
    T = (tauPa - tauCrPa) ./ tauCrPa;
    return
end

if all(isfinite([theta,thetaCr])) && thetaCr > 0
    T = (theta - thetaCr) ./ thetaCr;
    return
end

if isfinite(Tinput)
    T = Tinput;
end

end


function d = finiteDifference(a,b)

if isfinite(a) && isfinite(b)
    d = a-b;
else
    d = NaN;
end

end


function r = relativeDifference(a,b)

a = double(a(:));
b = double(b(:));

n = min(numel(a),numel(b));
a = a(1:n);
b = b(1:n);

r = nan(n,1);

valid = isfinite(a) & isfinite(b);

denominator = max(abs(a(valid)), eps);
r(valid) = abs(a(valid)-b(valid)) ./ denominator;

end


function [Cb, Dthickness, QvR] = ...
    calculateVanRijnFromExcess( ...
    D50m, Tv, Dv, ros, ro, g)
% D50m       [m]
% Tv         excess shear [-]
% Dv         dimensionless particle diameter [-]
% Cb         volumetric concentration [-]
% Dthickness [m]
% QvR        [kg/(m s)]

Cb = NaN;
Dthickness = NaN;
QvR = NaN;

if ~all(isfinite([D50m,Tv,Dv,ros,ro,g])) || ...
        D50m <= 0 || Dv <= 0 || ros <= ro || ...
        ro <= 0 || g <= 0 || Tv <= 0 || ~isreal(Tv)

    return
end

S = ros./ro - 1;

Dthickness = ...
    0.3 .* Dv.^0.7 .* Tv.^0.5 .* D50m;

Cb = ...
    0.18 .* Tv ./ Dv;

QvR = ...
    0.053 .* Tv.^2.1 ./ Dv.^0.3 .* ...
    sqrt(S.*g) .* D50m.^1.5 .* ros;

end


function [Cb, valid] = validateCb(CbRaw, CbMin, CbMax)

valid = ...
    isfinite(CbRaw) && ...
    isreal(CbRaw) && ...
    CbRaw >= CbMin && ...
    CbRaw <= CbMax;

if valid
    Cb = real(CbRaw);
else
    Cb = NaN;
end

end


function [Dthickness, QvR] = ...
    invalidateAssociated(Dthickness, QvR, valid)

if ~valid
    Dthickness = NaN;
    QvR = NaN;
end

end


function ws = settlingVelocitySoulsby( ...
    D50mm, ros, ro, ni, g)

s = ros./ro;
D = D50mm./1000;

if ~all(isfinite([D,s,ni,g])) || ...
        D <= 0 || s <= 1 || ni <= 0
    ws = NaN;
    return
end

Dstar = D .* ...
    ((s-1).*g./ni.^2).^(1/3);

ws = ni./D .* ...
    (sqrt(10.36.^2 + 1.049.*Dstar.^3) - 10.36);

end


function [Dt1, Dt2, Dt3, Dt4, Dt5, Dt6, Dt7] = ...
    Dthickness6( ...
    D50mm, thetaCr, Tau, rs, Dv, Tv, Ushear, ws)
% D50mm  [mm]
% thetaCr and Tau are both dimensionless Shields stresses.
% All output thicknesses are in metres.

D50 = D50mm./1000;
g = 9.81;

Dt1 = NaN;
Dt2 = NaN;
Dt3 = NaN;
Dt4 = NaN;
Dt5 = NaN;
Dt6 = NaN;
Dt7 = NaN;

if ~isfinite(D50) || D50 <= 0
    return
end

if all(isfinite([Tau,thetaCr])) && ...
        thetaCr > 0 && Tau > thetaCr

    val = ...
        14.3 .* D50 .* ...
        (Tau-thetaCr).^0.575;

    Dt1 = validThickness(val,D50);
end

if isfinite(Tau) && Tau > 0 && ...
        isfinite(rs) && rs > 0

    val = ...
        1.78 .* D50 .* ...
        (rs./1000).^0.86 .* ...
        Tau.^0.69;

    Dt2 = validThickness(val,D50);
end

if all(isfinite([Dv,Tv])) && ...
        Dv > 0 && Tv > 0

    val = ...
        0.154 .* Dv.^0.495 .* ...
        Tv.^0.73 .* D50;

    Dt3 = validThickness(val,D50);
end

if isfinite(Ushear) && Ushear >= 0

    val = ...
        0.13 .* Ushear.^2 ./ g;

    Dt4 = validThickness(val,D50);
end

if all(isfinite([Ushear,ws])) && ...
        Ushear >= 0 && ws > 0

    val = ...
        0.3 .* D50 .* ...
        (Ushear./ws).^1.5;

    Dt5 = validThickness(val,D50);

    val = ...
        D50 .* ...
        (Ushear./ws).^1.2;

    Dt6 = validThickness(val,D50);

    val = ...
        D50 .* ...
        (0.3.*(Ushear./ws).^2 + 0.1);

    Dt7 = validThickness(val,D50);
end

end


function value = validThickness(value, D50) %#ok<INUSD>
% Retain every positive finite formula result. The requested final range
% (1.5*D50 to 10*D50 by default) is applied later during aggregation.

if ~isfinite(value) || ~isreal(value) || value <= 0
    value = NaN;
else
    value = real(value);
end

end


function qs = kinematicTransport( ...
    D50, n, Cb, Btmodf, ro, Ninstrument)
%KINEMATICTRANSPORT Bedload transport [kg/(m s)].
%
% D50    [m]
% n      [-]
% Cb     [-]
% Btmodf [m/s]
% ro          sediment density [kg/m^3]
% Ninstrument instrument-specific correction coefficient [-]

Dthickness = n.*D50;

if ~all(isfinite([D50,n,Cb,Btmodf,ro,Ninstrument])) || ...
        D50 <= 0 || n <= 0 || Cb <= 0 || ...
        Btmodf < 0 || ro <= 0 || Ninstrument <= 0

    qs = NaN;
else
    qs = Ninstrument.*Btmodf.*Cb.*Dthickness.*ro;
end

end


function qs = kinematicTransportFinal( ...
    Cb, Dthickness, Btmodf, ros, Ninstrument)
%KINEMATICTRANSPORTFINAL Additional final kinematic transport [kg/(m s)].
%
% This does not replace kinematicTransport(), which retains the original
% Dthickness = n*D50mobile formulation for all method-specific outputs.
%
% Cb         final bounded/fallback volumetric concentration [-]
% Dthickness final bounded/fallback active-layer thickness [m]
% Btmodf     apparent bedload velocity [m/s]
% ros         sediment density [kg/m^3]
% Ninstrument instrument-specific correction coefficient [-]

if ~all(isfinite([Cb,Dthickness,Btmodf,ros,Ninstrument])) || ...
        Cb <= 0 || Dthickness <= 0 || ...
        Btmodf < 0 || ros <= 0 || Ninstrument <= 0

    qs = NaN;
else
    qs = Ninstrument .* Btmodf .* Cb .* Dthickness .* ros;
end

end


function qs = transportWithThickness( ...
    Btmodf, Cb, Dthickness, ros, Ninstrument)

if ~all(isfinite([Btmodf,Cb,Dthickness,ros,Ninstrument])) || ...
        Btmodf < 0 || Cb <= 0 || ...
        Dthickness <= 0 || ros <= 0 || Ninstrument <= 0

    qs = NaN;
else
    qs = Ninstrument.*Btmodf.*Cb.*Dthickness.*ros;
end

end


function y = safeMean(x)

if isempty(x)
    y = NaN;
    return
end

x = double(x(:));
x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = mean(x);
end

end


function y = safeMeanMatrixColumns(a,b)

a = double(a(:));
b = double(b(:));

nRows = max(numel(a),numel(b));

if numel(a) < nRows
    a(end+1:nRows) = NaN;
end

if numel(b) < nRows
    b(end+1:nRows) = NaN;
end

y = nanmeanLocal([a b],2);

end


function y = nanmeanLocal(A,dim)

if nargin < 2
    dim = 1;
end

A = double(A);
valid = isfinite(A);
A(~valid) = 0;

nValid = sum(valid,dim);
y = sum(A,dim)./nValid;
y(nValid == 0) = NaN;

end


function [finalValue, nUsed, source, usedFallback] = ...
    aggregateBoundedCandidates( ...
    candidates, lowerBound, upperBound, statistic, ...
    useFallback, fallbackValue)
%AGGREGATEBOUNDEDCANDIDATES Aggregate row-wise values inside inclusive bounds.

candidates = double(candidates);
nRows = size(candidates,1);

lowerBound = prepareVector(lowerBound,nRows);
upperBound = prepareVector(upperBound,nRows);
fallbackValue = prepareVector(fallbackValue,nRows);
statistic = lower(string(statistic));

finalValue = nan(nRows,1);
nUsed = zeros(nRows,1);
source = strings(nRows,1);
usedFallback = false(nRows,1);

for i = 1:nRows

    row = candidates(i,:);

    valid = ...
        isfinite(row) & ...
        imag(row) == 0 & ...
        real(row) >= lowerBound(i) & ...
        real(row) <= upperBound(i);

    values = real(row(valid));
    nUsed(i) = numel(values);

    if nUsed(i) > 0

        switch statistic
            case "mean"
                finalValue(i) = mean(values);
            case "median"
                finalValue(i) = median(values);
            otherwise
                error('Unsupported aggregation statistic: %s', statistic);
        end

        source(i) = statistic + "_of_" + string(nUsed(i)) + "_methods";

    elseif useFallback && ...
            isfinite(fallbackValue(i)) && ...
            fallbackValue(i) >= lowerBound(i) && ...
            fallbackValue(i) <= upperBound(i)

        finalValue(i) = fallbackValue(i);
        source(i) = "fallback";
        usedFallback(i) = true;

    else
        source(i) = "missing";
    end
end

end


function T = makeFinalTrateComparisonTable( ...
    calculated, measuredMean, measuredMax, minN)
%MAKEFINALTRATECOMPARISONTABLE Statistics for final calculated transport.
%
% All inputs must use the same unit, here kg/(m s). Bias is calculated as
% calculated minus measured, so positive bias indicates overprediction.

labels = ["Measured mean"; "Measured maximum"];
measured = {measuredMean, measuredMax};

N = zeros(2,1);
PearsonR = nan(2,1);
R2 = nan(2,1);
Bias = nan(2,1);
MAE = nan(2,1);
RMSE = nan(2,1);
SlopeOrigin = nan(2,1);

for k = 1:2
    [N(k),PearsonR(k),R2(k),Bias(k),MAE(k),RMSE(k), ...
        SlopeOrigin(k)] = finalComparisonStats( ...
            calculated, measured{k}, minN);
end

T = table(labels,N,PearsonR,R2,Bias,MAE,RMSE,SlopeOrigin, ...
    'VariableNames', { ...
    'MeasuredStatistic','N','PearsonR','R2','BiasCalculatedMinusMeasured', ...
    'MAE','RMSE','SlopeThroughOrigin'});

end


function plotFinalTrateComparison( ...
    calculated, measuredMean, measuredMax, useLogAxes, minN)
%PLOTFINALTRATECOMPARISON Plot final calculated rate against observations.
%
% x = measured transport rate [kg/(m s)]
% y = calculated final kinematic transport rate [kg/(m s)]

measured = {measuredMean, measuredMax};
titles = { ...
    'Final calculated vs measured mean', ...
    'Final calculated vs measured maximum'};

allValues = [calculated(:); measuredMean(:); measuredMax(:)];
allValues = allValues(isfinite(allValues));
if useLogAxes
    allValues = allValues(allValues > 0);
else
    allValues = allValues(allValues >= 0);
end

if isempty(allValues)
    warning('processBedloadTrate:NoFinalTratePlotData', ...
        'No valid data are available for the final transport comparison plot.');
    return
end

sharedLimits = expandedAgreementLimits(allValues,useLogAxes);

figure('Name','Final calculated versus measured transport','Color','w');
try
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    useTiles = true;
catch
    useTiles = false;
end

for k = 1:2
    if useTiles
        nexttile;
    else
        subplot(1,2,k);
    end
    hold on;

    obs = double(measured{k}(:));
    calc = double(calculated(:));
    nCommon = min(numel(obs),numel(calc));
    obs = obs(1:nCommon);
    calc = calc(1:nCommon);

    valid = isfinite(obs) & isfinite(calc);
    if useLogAxes
        valid = valid & obs > 0 & calc > 0;
    else
        valid = valid & obs >= 0 & calc >= 0;
    end

    plot(obs(valid),calc(valid),'o','LineStyle','none', ...
        'DisplayName','Measurements');
    plot(sharedLimits,sharedLimits,'--','DisplayName','1:1 agreement');

    xlim(sharedLimits);
    ylim(sharedLimits);
    axis square;
    grid on;
    box on;
    if useLogAxes
        set(gca,'XScale','log','YScale','log');
    end

    xlabel('Measured transport rate [kg/(m s)]');
    ylabel('Final calculated transport rate [kg/(m s)]');
    title(titles{k},'Interpreter','none');
    legend('Location','best');

    [nValid,r,r2,bias,mae,rmse,slope0] = ...
        finalComparisonStats(calc(valid),obs(valid),minN);
    addFinalComparisonText(sprintf([ ...
        'n = %d\nr = %.3f\nR^2 = %.3f\nBias = %.3g\n', ...
        'MAE = %.3g\nRMSE = %.3g\nSlope_0 = %.3f'], ...
        nValid,r,r2,bias,mae,rmse,slope0));
end

end


function [nValid,r,r2,bias,mae,rmse,slope0] = ...
    finalComparisonStats(calculated,measured,minN)
%FINALCOMPARISONSTATS Agreement statistics using paired finite values.

calculated = double(calculated(:));
measured = double(measured(:));
nCommon = min(numel(calculated),numel(measured));
calculated = calculated(1:nCommon);
measured = measured(1:nCommon);

valid = isfinite(calculated) & isfinite(measured);
calculated = calculated(valid);
measured = measured(valid);
nValid = numel(calculated);

r = NaN;
r2 = NaN;
bias = NaN;
mae = NaN;
rmse = NaN;
slope0 = NaN;

if nValid == 0
    return
end

errors = calculated-measured;
bias = mean(errors);
mae = mean(abs(errors));
rmse = sqrt(mean(errors.^2));

denSlope = sum(measured.^2);
if denSlope > 0
    slope0 = sum(measured.*calculated)./denSlope;
end

if nValid >= minN && ...
        numel(unique(measured)) >= 2 && ...
        numel(unique(calculated)) >= 2
    C = corrcoef(measured,calculated);
    r = C(1,2);
    r2 = r.^2;
end

end


function lim = expandedAgreementLimits(values,useLogAxes)
%EXPANDEDAGREEMENTLIMITS Common square limits for agreement plots.

values = double(values(:));
values = values(isfinite(values));

if useLogAxes
    values = values(values > 0);
    lo = min(values);
    hi = max(values);
    if lo == hi
        lo = lo./1.5;
        hi = hi.*1.5;
    end
    logPad = 0.05.*(log10(hi)-log10(lo));
    if logPad == 0
        logPad = 0.05;
    end
    lim = [10.^(log10(lo)-logPad),10.^(log10(hi)+logPad)];
else
    values = values(values >= 0);
    lo = min([0;values]);
    hi = max(values);
    if hi <= lo
        hi = lo + max(abs(lo).*0.1,1e-9);
    end
    pad = 0.05.*(hi-lo);
    lim = [max(0,lo-pad),hi+pad];
end

end


function addFinalComparisonText(txt)
%ADDFINALCOMPARISONTEXT Place statistics inside the current axes.

drawnow;
xl = xlim;
yl = ylim;

if strcmpi(get(gca,'XScale'),'log') && all(xl > 0)
    xpos = 10.^(log10(xl(1)) + 0.05.*diff(log10(xl)));
else
    xpos = xl(1) + 0.05.*diff(xl);
end

if strcmpi(get(gca,'YScale'),'log') && all(yl > 0)
    ypos = 10.^(log10(yl(2)) - 0.05.*diff(log10(yl)));
else
    ypos = yl(2) - 0.05.*diff(yl);
end

text(xpos,ypos,txt, ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.7 0.7 0.7], ...
    'Margin',4);

end


function assignStructToWorkspace(S,workspaceName)

f = fieldnames(S);

for i = 1:numel(f)
    assignin(workspaceName,f{i},S.(f{i}));
end

end
