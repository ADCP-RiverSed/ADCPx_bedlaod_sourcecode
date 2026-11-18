function BackscatterOut = processBackscatterM9( ...
    dataM9, ADCPin, DepthRNGout, BTvelOut, varargin)
%PROCESSBACKSCATTERM9 Process SonTek M9 five-beam bottom backscatter.
%
% BackscatterOut = processBackscatterM9( ...
%     dataM9, ADCPin, DepthRNGout, BTvelOut, 'Name',Value,...)
%
% This function converts the legacy SonTek script into the same structured,
% input-table-driven workflow used by processBackscatter. It processes:
%
%   - four 3-MHz slanted beams at 25 degrees;
%   - one 0.5-MHz vertical beam;
%   - SonTek Strength/RSSI and Contrast as separate quantities;
%   - optional water-column echo profiles for the four slanted beams;
%   - optional vertical gain-step correction from the legacy script.
%
% The Contrast-based SV values are the primary SonTek outputs (SV1_S ...
% SV4_S and SV_VB). Strength/RSSI-based alternatives are retained as
% SV1_Srssi ... SV4_Srssi and SV_VBrssi. The 3-MHz and 0.5-MHz values are
% never averaged together.
%
% INPUTS
% ------
% dataM9      : dir-style structure or folder containing processed M9 MAT
%               files. Both standardized adapter fields and legacy
%               BS_adp fields are supported.
% ADCPin      : physical/input-table structure with ADCPfilename (or an
%               M9filename/Sontekfilename alias) and optional
%               StartADCP/EndADCP.
% DepthRNGout : output from processDepthsRNG. Its depth summaries are copied
%               unchanged. Processed range series are preferred when found.
% BTvelOut    : M9 apparent-bedload-velocity output. Its exact final mask is
%               transferred to range, Strength, Contrast, pressure and echo.
%
% STANDARDIZED FIELD ALIASES
% --------------------------
% Slanted ranges:
%   BT.RNG1 ... BT.RNG4, BS_adp.RNG11 ... BS_adp.RNG44
% Strength:
%   BT.BT_Strength1 ... 4, BS_adp.BTstrength1 ... 4
% Contrast:
%   BT.BT_Contrast1 ... 4, BS_adp.BTcontrast1 ... 4
% Echo profiles:
%   ECHO.ECH1 ... 4, BS_adp.ECH1 ... 4
% Vertical beam:
%   VB.Range/RSSI/Contrast or
%   BS_adp.VB_RNGM/VB_RSSI2M/VB_CNTM
%
% A T.BS_adp or T.BT/T.ECHO/T.VB nesting level is also accepted.
%
% ACOUSTIC EQUATIONS
% ------------------
% The equations follow the attached legacy M9 script:
%
%   SL = mean(10*log10(Pressure*PressureToPa))
%
%   SV_strength = -SL + Strength + spreading - footprint
%                 + 2*TL - Offset - LegacyDbCorrection
%
%   SV_contrast = -SL + Contrast + spreading - footprint
%                 + 2*TL - Offset - LegacyDbCorrection
%
% The 3-MHz footprint uses the SonTek 1-3 and 2-4 beam-pair slopes. The
% vertical footprint follows the legacy pulse/range-squared formulation.
% Water absorption is calculated separately at 3000 and 500 kHz.
%
% IMPORTANT CALIBRATION NOTE
% --------------------------
% width3MHz, omega3MHz, width05MHz, omega05MHz, BT3MHzPL,
% BT05MHzPL, Offset3MHz and Offset05MHz were external constants in the old
% script. This function reads them from file metadata when available or from
% Name-Value options. Absolute SV is flagged unavailable when the required
% source-level, pulse-length or beam-opening inputs are missing. It does not
% fabricate RDI transmit-current/voltage variables.
%
% MAIN OPTIONS
% ------------
% Geometry/frequency:
%   theta               slanted-beam angle [deg], default 25
%   f3                   slanted-beam frequency [kHz], default 3000
%   f05                  vertical-beam frequency [kHz], default 500
%   width3MHz            slanted transducer width [m], default 0.025
%   omega3MHz            full slanted opening [rad/deg], default 0.0471733 rad
%   DistTrans3MHz        centre-to-transducer radius [m], default 0.10
%   width05MHz           vertical transducer width [m], default 0.040
%   omega05MHz           full vertical opening [rad/deg], default 0.1769 rad
%   BT3MHzPL             slanted bottom pulse length [m], default 0.21
%   BT05MHzPL            vertical bottom pulse length [m], default 0.05
%
% Acoustic:
%   alphaw3              fallback 3-MHz water absorption [dB/m]
%   alphaw05             fallback 0.5-MHz water absorption [dB/m]
%   alphaSed3/alphaSed05 sediment attenuation [dB/m]
%   PressureToPa         legacy pressure conversion, default 9806650
%   Offset3MHz           four offsets [dB], default [18.6 19.2 21.6 19.6]
%   Offset05MHz          vertical calibration offset [dB], default 11.6
%   SLaddVB              additional vertical source-level offset, default 0
%   LegacyDbCorrection   legacy common subtraction [dB], default 9
%   ExactLegacyFootprintOmega default true. Reproduces the old script's
%                        deg2rad(omegaMHz) conversion inside the footprint
%                        equation while retaining omega in radians for the
%                        SamplingVolume geometry.
%
% Scaling:
%   AmplitudeScale       numeric or 'auto'. Auto uses 0.1 for legacy
%                        BS_adp deci-dB data and otherwise inspects magnitude.
%   EchoScale            numeric or 'auto', same convention for ECH1...4.
%   RangeMaxDepth        accepted depth/range maximum [m], default 10
%
% Vertical gain correction:
%   CalibrateVB          true by default for the supplied legacy M9 data
%   VBGain               value added below the low thresholds
%   VBStrengthLowThreshold, VBStrengthHighThreshold
%   VBContrastLowThreshold, VBContrastHighThreshold
%
% Trimming/mask:
%   StartEndMode         'trimcounts' (default) or 'absolute'
%   VelocityMaskMode     'final' (default), 'boat', or 'none'
%   RequireExactVelocityMask default false. Missing masks are reported but
%                        do not invalidate otherwise usable acoustic data.
%
% Output/workspace:
%   ApplyMethod2Regression default true
%   AssignToCaller       default true
%   AssignToBase         default false
%   Verbose              default true
%
% No files are saved. All time series are stored in cells, one cell per
% measurement. Only row-wise summaries are numeric/string vectors.

%% Parse inputs

p=inputParser;
p.FunctionName=mfilename;
addRequired(p,'dataM9',@(x) isstruct(x) || ischar(x) || isstring(x));
addRequired(p,'ADCPin',@(x) isstruct(x));
addRequired(p,'DepthRNGout',@(x) isempty(x) || isstruct(x));
addRequired(p,'BTvelOut',@(x) isempty(x) || isstruct(x));

addParameter(p,'theta',25,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<90);
addParameter(p,'f3',3000,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'f05',500,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'width3MHz',0.025,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'omega3MHz',0.047173333333,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'DistTrans3MHz',0.10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'width05MHz',0.040,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'omega05MHz',0.1769,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'BT3MHzPL',0.21,@(x) isnumeric(x));
addParameter(p,'BT05MHzPL',0.05,@(x) isnumeric(x));
addParameter(p,'distanceADCPtoWaterSurface',0.05,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);

addParameter(p,'alphaw3',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x>=0));
addParameter(p,'alphaw05',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x>=0));
addParameter(p,'alphaSed3',0,@(x) isnumeric(x) && ~isempty(x));
addParameter(p,'alphaSed05',0,@(x) isnumeric(x) && ~isempty(x));
addParameter(p,'PressureToPa',9806650,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'Offset3MHz',[18.6 19.2 21.6 19.6],@(x) isnumeric(x) && ~isempty(x));
addParameter(p,'Offset05MHz',11.6,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'SLaddVB',0,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'LegacyDbCorrection',9,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'ExactLegacyFootprintOmega',true,@(x) islogical(x) || (isnumeric(x) && isscalar(x)));

addParameter(p,'AmplitudeScale','auto',@(x) isnumeric(x) || ischar(x) || isstring(x));
addParameter(p,'EchoScale','auto',@(x) isnumeric(x) || ischar(x) || isstring(x));
addParameter(p,'RangeMaxDepth',10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'Tdespike',2,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'MM',0,@(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'Win',3,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=1);
addParameter(p,'PLOT','NO',@(x) ischar(x) || isstring(x));
addParameter(p,'StrengthContrastMax',1500,@(x) isnumeric(x) && isscalar(x) && isfinite(x));

addParameter(p,'CalibrateVB',true,@(x) islogical(x) || ischar(x) || isstring(x) || isnumeric(x));
addParameter(p,'VBGain',40,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'VBStrengthLowThreshold',60,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'VBStrengthHighThreshold',100,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'VBContrastLowThreshold',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'VBContrastHighThreshold',70,@(x) isnumeric(x) && isscalar(x));

addParameter(p,'StartEndMode','trimcounts',@(x) ischar(x) || isstring(x));
addParameter(p,'UseVelocityFilterMask',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'VelocityMaskMode','final',@(x) ischar(x) || isstring(x));
addParameter(p,'RequireExactVelocityMask',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'nodataS',[],@(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(p,'ApplyMethod2Regression',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToCaller',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToBase',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'Verbose',true,@(x) islogical(x) || isnumeric(x));

parse(p,dataM9,ADCPin,DepthRNGout,BTvelOut,varargin{:});
opt=p.Results;
requiredNames={'dataM9','ADCPin','DepthRNGout','BTvelOut'};
for iReq=1:numel(requiredNames)
    if isfield(opt,requiredNames{iReq})
        opt=rmfield(opt,requiredNames{iReq});
    end
end

if isempty(opt.alphaw3)
    opt.alphaw3=0.392.*(opt.f3./1200).^2;
end
if isempty(opt.alphaw05)
    opt.alphaw05=0.392.*(opt.f05./1200).^2;
end

startEndMode=lower(strtrim(string(opt.StartEndMode)));
if ~any(startEndMode==["trimcounts","trim","absolute","indices"])
    error('StartEndMode must be trimcounts or absolute.');
end
velocityMaskMode=lower(strtrim(string(opt.VelocityMaskMode)));
if ~logical(opt.UseVelocityFilterMask), velocityMaskMode="none"; end
if ~any(velocityMaskMode==["final","boat","none"])
    error('VelocityMaskMode must be final, boat, or none.');
end
calibrateVB=isCalibrationRequested(opt.CalibrateVB);

%% Input rows and shared outputs

[ADCPfilename,ADCPfilenameField]=resolveM9FilenameInput(ADCPin);
nFiles=numel(ADCPfilename);
InputRow=(1:nFiles)';
StartADCP=getVectorFieldAny(ADCPin,{'StartADCP','Start','start'},zeros(nFiles,1),nFiles);
EndADCP=getVectorFieldAny(ADCPin,{'EndADCP','End','end'},zeros(nFiles,1),nFiles);
DistanceFromLeftBank=getVectorFieldAny(ADCPin, ...
    {'DistanceFromLeftBank','DistanceToBank'},nan(nFiles,1),nFiles);
D50_14S=getVectorFieldAny(ADCPin, ...
    {'D50_14','D50_14mm','D50mobile14'},nan(nFiles,1),nFiles);
D50_05S=getVectorFieldAny(ADCPin, ...
    {'D50_05','D50_05mm','D50mobile05'},nan(nFiles,1),nFiles);

BT3MHzPLOption=prepareVector(opt.BT3MHzPL,nFiles);
BT05MHzPLOption=prepareVector(opt.BT05MHzPL,nFiles);
alphaSed3=prepareVector(opt.alphaSed3,nFiles);
alphaSed05=prepareVector(opt.alphaSed05,nFiles);
Offset3MHz=prepareFourVector(opt.Offset3MHz);

fileList=normaliseFileList(dataM9);

BTmodf_S=getVectorFieldAny(BTvelOut, ...
    {'BTmodf_S','va_final_S','BTmodf_1200','va_final_1200'},nan(nFiles,1),nFiles);
FilteredDataN_S=getVectorFieldAny(BTvelOut, ...
    {'FilteredDataN_S','FilteredDataNS','FilteredDataN_1200'},nan(nFiles,1),nFiles);
BTstd_S=getVectorFieldAny(BTvelOut, ...
    {'BTstd_S','BTstdS','BTstd_1200'},nan(nFiles,1),nFiles);

if isempty(DepthRNGout) || ~isstruct(DepthRNGout)
    error('DepthRNGout must be the output from processDepthsRNG.');
end
Depth_S=getMatrixFieldAny(DepthRNGout, ...
    {'Depth1200','Depth_S','DepthS'},nFiles,4);
Depth_Sstd=getMatrixFieldAny(DepthRNGout, ...
    {'Depth1200std','Depth_Sstd','DepthSstd'},nFiles,4);
stdSpatial_S=getVectorFieldAny(DepthRNGout, ...
    {'stdSpatial_1200','stdSpatial_S'},nan(nFiles,1),nFiles);
stdSpatialALL_S=getCellFieldAny(DepthRNGout, ...
    {'stdSpatialALL_1200','stdSpatialALL_S'},nFiles);

% Preserve the current depth-output names for application compatibility.
Depth1200=Depth_S;
Depth1200std=Depth_Sstd;
stdSpatial_1200=stdSpatial_S;
stdSpatialALL_1200=stdSpatialALL_S;
% Legacy SonTek output alias retained for existing downstream scripts.
DepthS=Depth_S;

%% Preallocate row-wise results

FileFound=false(nFiles,1);
MatchedFile=strings(nFiles,1);
Message=strings(nFiles,1);
AcousticMessage_S=strings(nFiles,1);
DummyValue=nan(nFiles,1);
StartUsed_S=nan(nFiles,1);
EndUsed_S=nan(nFiles,1);
NensemblesRaw_S=nan(nFiles,1);
NensemblesUsed_S=nan(nFiles,1);
NfilteredBTvelocity_S=zeros(nFiles,1);
BTvelocityFilterMask_S=cell(nFiles,1);
BTvelocityFilterSource_S=strings(nFiles,1);
BTvelocityFilterExact_S=false(nFiles,1);
BTvelocityFilterAligned_S=false(nFiles,1);
BTvelocityMaskMode_S=repmat(velocityMaskMode,nFiles,1);

RangeScaleToMeters_S=nan(nFiles,1);
RangeOriginalUnit_S=strings(nFiles,1);
VBRangeScaleToMeters_S=nan(nFiles,1);
VBRangeOriginalUnit_S=strings(nFiles,1);
AmplitudeScaleUsed_S=nan(nFiles,1);
EchoScaleUsed_S=nan(nFiles,1);
RangeSource_S=strings(nFiles,4);
StrengthSource_S=strings(nFiles,4);
ContrastSource_S=strings(nFiles,4);
EchoSource_S=strings(nFiles,4);
VBRangeSource_S=strings(nFiles,1);
VBStrengthSource_S=strings(nFiles,1);
VBContrastSource_S=strings(nFiles,1);
SlantedLayout_S=strings(nFiles,1);
BackscatterLayout_S=strings(nFiles,1);
SlantedInputsAvailable_S=false(nFiles,1);
SlantedContrastAvailable_S=false(nFiles,1);
SlantedStrengthAvailable_S=false(nFiles,1);

RNG11_S=cell(nFiles,1); RNG22_S=cell(nFiles,1);
RNG33_S=cell(nFiles,1); RNG44_S=cell(nFiles,1);
V_RNG_S=cell(nFiles,1);
RNG_VB_Sm=nan(nFiles,1); RNG_VB_Sstd=nan(nFiles,1);

BT_RSSI1_S=cell(nFiles,1); BT_RSSI2_S=cell(nFiles,1);
BT_RSSI3_S=cell(nFiles,1); BT_RSSI4_S=cell(nFiles,1);
BT_CNT1_S=cell(nFiles,1); BT_CNT2_S=cell(nFiles,1);
BT_CNT3_S=cell(nFiles,1); BT_CNT4_S=cell(nFiles,1);
VB_RSSI_S=cell(nFiles,1); VB_CNT_S=cell(nFiles,1);
BT_RSSI_S=cell(nFiles,1); BT_RSSI_Sstd=cell(nFiles,1);
BT_CNT_S=cell(nFiles,1); BT_CNT_Sstd=cell(nFiles,1);
BT_RSSImean_S=nan(nFiles,4); BT_CNTmean_S=nan(nFiles,4);
VB_RSSI_Sm=nan(nFiles,1); VB_RSSI_Sstd=nan(nFiles,1);
VB_CNT_Sm=nan(nFiles,1); VB_CNT_Sstd=nan(nFiles,1);
stdSpatialRSSI_S=nan(nFiles,1); stdSpatialCNT_S=nan(nFiles,1);

WF_S=nan(nFiles,1); WSmean=nan(nFiles,1); WT_3=nan(nFiles,1);
Pressure_S=cell(nFiles,1); SL_S=nan(nFiles,1);
PressureSource_S=strings(nFiles,1);

psi1_S=nan(nFiles,1); psi2_S=nan(nFiles,1);
psi3_S=nan(nFiles,1); psi4_S=nan(nFiles,1);
alfaXmean_S=nan(nFiles,1); alfaYmean_S=nan(nFiles,1);
betaXmean_S=nan(nFiles,1); betaYmean_S=nan(nFiles,1);
ni_mean_S=nan(nFiles,1); ni_std_S=nan(nFiles,1);

R1_S=nan(nFiles,1); R2_S=nan(nFiles,1);
R3_S=nan(nFiles,1); R4_S=nan(nFiles,1);
Spread1_S=nan(nFiles,1); Spread2_S=nan(nFiles,1);
Spread3_S=nan(nFiles,1); Spread4_S=nan(nFiles,1);
SpreadMS=nan(nFiles,1); SpreadVBmean=nan(nFiles,1);

BT3MHzPL_S=nan(nFiles,1); BT05MHzPL_S=nan(nFiles,1);
Dmin_S_1=nan(nFiles,1); Dmin_S_2=nan(nFiles,1);
Dmin_S_3=nan(nFiles,1); Dmin_S_4=nan(nFiles,1);
Dmax_S_1=nan(nFiles,1); Dmax_S_2=nan(nFiles,1);
Dmax_S_3=nan(nFiles,1); Dmax_S_4=nan(nFiles,1);
Area_S_1=nan(nFiles,1); Area_S_2=nan(nFiles,1);
Area_S_3=nan(nFiles,1); Area_S_4=nan(nFiles,1);
distance_S_1=nan(nFiles,1); distance_S_2=nan(nFiles,1);
distance_S_3=nan(nFiles,1); distance_S_4=nan(nFiles,1);
Lp_S_1=nan(nFiles,1); Lp_S_2=nan(nFiles,1);
Lp_S_3=nan(nFiles,1); Lp_S_4=nan(nFiles,1);
Dmin05MHz_S=nan(nFiles,1); Dmax05MHz_S=nan(nFiles,1);
Area05MHz_S=nan(nFiles,1); distance05MHz_S=nan(nFiles,1);
Lp05MHz_S=nan(nFiles,1);
Af1meanS=nan(nFiles,1); Af2meanS=nan(nFiles,1);
Af3meanS=nan(nFiles,1); Af4meanS=nan(nFiles,1);
AfVBmean_S=nan(nFiles,1);

ECH1_S=cell(nFiles,1); ECH2_S=cell(nFiles,1);
ECH3_S=cell(nFiles,1); ECH4_S=cell(nFiles,1);
alpha1meanS=cell(nFiles,1); alpha2meanS=cell(nFiles,1);
alpha3meanS=cell(nFiles,1); alpha4meanS=cell(nFiles,1);
alpha1meanSgrad=cell(nFiles,1); alpha2meanSgrad=cell(nFiles,1);
alpha3meanSgrad=cell(nFiles,1); alpha4meanSgrad=cell(nFiles,1);
dEDR1_S=cell(nFiles,1); dEDR2_S=cell(nFiles,1);
dEDR3_S=cell(nFiles,1); dEDR4_S=cell(nFiles,1);
alphaMean3000S=nan(nFiles,4);

AW_3=nan(nFiles,1); AW_05=nan(nFiles,1);
abswat3_1=nan(nFiles,1); abswat3_2=nan(nFiles,1);
abswat3_3=nan(nFiles,1); abswat3_4=nan(nFiles,1);
abswat05=nan(nFiles,1);
TL1=nan(nFiles,1); TL2=nan(nFiles,1);
TL3=nan(nFiles,1); TL4=nan(nFiles,1);
TL1_fixS=nan(nFiles,1); TL2_fixS=nan(nFiles,1);
TL3_fixS=nan(nFiles,1); TL4_fixS=nan(nFiles,1);
TLVBMHz_fixS=nan(nFiles,1);

SV1_Srssi=cell(nFiles,1); SV2_Srssi=cell(nFiles,1);
SV3_Srssi=cell(nFiles,1); SV4_Srssi=cell(nFiles,1);
SV1_S=cell(nFiles,1); SV2_S=cell(nFiles,1);
SV3_S=cell(nFiles,1); SV4_S=cell(nFiles,1);
SV_VBrssi=cell(nFiles,1); SV_VB=cell(nFiles,1);

SVmean1rssi=nan(nFiles,1); SVmean2rssi=nan(nFiles,1);
SVmean3rssi=nan(nFiles,1); SVmean4rssi=nan(nFiles,1);
SVmean1stdRSSI=nan(nFiles,1); SVmean2stdRSSI=nan(nFiles,1);
SVmean3stdRSSI=nan(nFiles,1); SVmean4stdRSSI=nan(nFiles,1);
SVmean1=nan(nFiles,1); SVmean2=nan(nFiles,1);
SVmean3=nan(nFiles,1); SVmean4=nan(nFiles,1);
SVmean1std=nan(nFiles,1); SVmean2std=nan(nFiles,1);
SVmean3std=nan(nFiles,1); SVmean4std=nan(nFiles,1);
SV_VBmeanSrssi=nan(nFiles,1); SV_VBmeanSstdRSSI=nan(nFiles,1);
SV_VBmeanS=nan(nFiles,1); SV_VBmeanSstd=nan(nFiles,1);

stdSvSpatial_S=nan(nFiles,1);
stdSvSpatialRSSI_S=nan(nFiles,1);
stdSvSpatialALL_S=cell(nFiles,1);
stdSvSpatialRSSIALL_S=cell(nFiles,1);
S_S=nan(nFiles,1); S_Sstd=nan(nFiles,1);
S_Srssi=nan(nFiles,1); S_Sall=cell(nFiles,1);

AbsoluteBackscatterReady_S=false(nFiles,1);
StrengthBackscatterReady_S=false(nFiles,1);
VerticalBackscatterReady_S=false(nFiles,1);
VerticalStrengthBackscatterReady_S=false(nFiles,1);
VBAvailable_S=false(nFiles,1);
VBStrengthAvailable_S=false(nFiles,1);
EchoProfilesAvailable_S=false(nFiles,1);
VBcalibrationApplied_S=false(nFiles,1);

D50S=nan(nFiles,1);
qsMethod2regression_S=nan(nFiles,1);
Method2Status_S=repmat("not_requested",nFiles,1);

%% Process rows

for i=1:nFiles
    try
        if isNoDataRow(opt.nodataS,i,nFiles)
            Message(i)="Skipped: nodataS";
            continue
        end
        DummyValue(i)=1;

        matchedFile=matchADCPfile(fileList,ADCPfilename(i));
        MatchedFile(i)=matchedFile;
        FileFound(i)=matchedFile~="";
        if ~FileFound(i)
            error('No MAT file matched ADCPfilename %s.',ADCPfilename(i));
        end
        Sload=load(matchedFile);

        rangeRaw=cell(1,4); strengthRaw=cell(1,4);
        contrastRaw=cell(1,4); echoRaw=cell(1,4);
        rangeSource=strings(1,4); strengthSource=strings(1,4);
        contrastSource=strings(1,4); echoSource=strings(1,4);
        for b=1:4
            [rangeRaw{b},rangeSource(b)]=findM9BeamField(Sload,"range",b);
            [strengthRaw{b},strengthSource(b)]=findM9BeamField(Sload,"strength",b);
            [contrastRaw{b},contrastSource(b)]=findM9BeamField(Sload,"contrast",b);
            [echoRaw{b},echoSource(b)]=findM9BeamField(Sload,"echo",b);
        end
        RangeSource_S(i,:)=rangeSource;
        StrengthSource_S(i,:)=strengthSource;
        ContrastSource_S(i,:)=contrastSource;
        EchoSource_S(i,:)=echoSource;
        SlantedLayout_S(i)=summariseSourceLayout( ...
            [rangeSource,strengthSource,contrastSource]);
        [vbRangeRaw,vbStrengthRaw,vbContrastRaw, ...
            VBRangeSource_S(i),VBStrengthSource_S(i), ...
            VBContrastSource_S(i),BackscatterLayout_S(i)]= ...
            findM9VerticalBundle(Sload);
        [pressureRaw,PressureSource_S(i)]=findM9ScalarOrSeries(Sload,"pressure");

        for b=1:4
            if isempty(rangeRaw{b}) || isempty(contrastRaw{b})
                error('M9 beam %d range or primary Contrast is missing.',b);
            end
        end
        SlantedInputsAvailable_S(i)=true;
        SlantedContrastAvailable_S(i)=all(cellfun( ...
            @(x) ~isempty(x) && any(isfinite(x)),contrastRaw));
        SlantedStrengthAvailable_S(i)=all(cellfun( ...
            @(x) ~isempty(x) && any(isfinite(x)),strengthRaw));

        nCandidates=[];
        for b=1:4
            nCandidates=[nCandidates,ensembleLength(rangeRaw{b}), ... %#ok<AGROW>
                ensembleLength(contrastRaw{b})]; %#ok<AGROW>
        end
        nCandidates=nCandidates(isfinite(nCandidates) & nCandidates>0);
        if isempty(nCandidates), error('No usable M9 ensemble series.'); end
        nRaw=min(nCandidates);
        NensemblesRaw_S(i)=nRaw;

        [idxStart,idxEnd]=resolveStartEnd( ...
            StartADCP(i),EndADCP(i),nRaw,startEndMode);
        idx=idxStart:idxEnd;
        nSeg=numel(idx);
        StartUsed_S(i)=idxStart;
        EndUsed_S(i)=idxEnd;
        NensemblesUsed_S(i)=nSeg;

        [velocityRejectMask,maskSource,maskExact,maskAligned]= ...
            deriveVelocityRejectMaskM9( ...
                BTvelOut,i,nFiles,idxStart,idxEnd,nRaw,nSeg,velocityMaskMode);
        velocityRejectMask=logical(velocityRejectMask(:));
        if numel(velocityRejectMask)~=nSeg
            maskAligned=false;
            velocityRejectMask=false(nSeg,1);
            maskSource=maskSource+"_length_mismatch_disabled";
        end
        if logical(opt.RequireExactVelocityMask) && velocityMaskMode=="final" && ...
                ~(maskExact && maskAligned)
            error(['No exact aligned final M9 BT mask for row %d. ', ...
                'Rerun the M9 BT-velocity function or use ', ...
                'VelocityMaskMode=''none''. Source: %s.'],i,maskSource);
        end
        BTvelocityFilterMask_S{i}=velocityRejectMask;
        BTvelocityFilterSource_S(i)=maskSource;
        BTvelocityFilterExact_S(i)=maskExact;
        BTvelocityFilterAligned_S(i)=maskAligned;
        NfilteredBTvelocity_S(i)=sum(velocityRejectMask);

        % Prefer processed depth-series cells. Raw values are unit-detected
        % against the authoritative row summary from processDepthsRNG.
        depthCells={RNG11_S,RNG22_S,RNG33_S,RNG44_S}; %#ok<NASGU>
        rangeProcessed=cell(1,4);
        rangeM=cell(1,4);
        rawAligned=cell(1,4);
        processedRangeAvailable=false(1,4);
        for b=1:4
            rangeProcessed{b}=getAlignedSeries(DepthRNGout, ...
                rangeOutputAliases(b),i,nFiles,idxStart,idxEnd,nRaw);
            processedRangeAvailable(b)=~isempty(rangeProcessed{b}) && ...
                any(isfinite(rangeProcessed{b}));
            [rawAligned{b},rawOK]=alignSeriesFlexible( ...
                rangeRaw{b},idxStart,idxEnd,nRaw,nSeg);
            if ~processedRangeAvailable(b) && ~rawOK
                error('M9 beam %d range is not aligned.',b);
            end
        end

        referenceDepth=safeMean(Depth_S(i,:));
        combinedRaw=[];
        for b=1:4
            if ~isempty(rawAligned{b}), combinedRaw=[combinedRaw;rawAligned{b}(:)]; end %#ok<AGROW>
        end
        [~,RangeOriginalUnit_S(i),RangeScaleToMeters_S(i),~]= ...
            convertRangeToMetres(combinedRaw,referenceDepth,opt.RangeMaxDepth);

        for b=1:4
            if processedRangeAvailable(b)
                x=forceColumnLength(rangeProcessed{b},nSeg);
            else
                x=forceColumnLength(rawAligned{b},nSeg).*RangeScaleToMeters_S(i);
            end
            x(x<=0 | x>opt.RangeMaxDepth)=NaN;
            x(velocityRejectMask)=NaN;
            rangeM{b}=despikeRangeLocal(x,3);
        end
        RNG11_S{i}=rangeM{1}; RNG22_S{i}=rangeM{2};
        RNG33_S{i}=rangeM{3}; RNG44_S{i}=rangeM{4};

        vbRangeProcessed=getAlignedSeries(DepthRNGout, ...
            {'VBDepth_','V_RNG_S','VBDepth_0500'}, ...
            i,nFiles,idxStart,idxEnd,nRaw);
        [vbRangeAligned,vbRangeOK]=alignSeriesFlexible( ...
            vbRangeRaw,idxStart,idxEnd,nRaw,nSeg);
        if ~isempty(vbRangeProcessed) && any(isfinite(vbRangeProcessed))
            vbRange=forceColumnLength(vbRangeProcessed,nSeg);
            VBRangeOriginalUnit_S(i)="m";
            VBRangeScaleToMeters_S(i)=1;
        elseif vbRangeOK
            [vbRange,VBRangeOriginalUnit_S(i), ...
                VBRangeScaleToMeters_S(i),~]=convertRangeToMetres( ...
                    vbRangeAligned,referenceDepth,opt.RangeMaxDepth);
        else
            vbRange=[];
        end
        if ~isempty(vbRange)
            vbRange=forceColumnLength(vbRange,nSeg);
            vbRange(vbRange<=0 | vbRange>opt.RangeMaxDepth)=NaN;
            vbRange(velocityRejectMask)=NaN;
            V_RNG_S{i}=despikeRangeLocal(vbRange,3);
            RNG_VB_Sm(i)=safeMean(V_RNG_S{i});
            RNG_VB_Sstd(i)=safeStd(V_RNG_S{i});
        end

        % Align and scale Strength and Contrast. The same scale is used for
        % all five beams, as in the old kn multiplier.
        strengthAligned=cell(1,4); contrastAligned=cell(1,4);
        amplitudeProbe=[];
        for b=1:4
            [strengthAligned{b},okS]=alignSeriesFlexible( ...
                strengthRaw{b},idxStart,idxEnd,nRaw,nSeg);
            [contrastAligned{b},okC]=alignSeriesFlexible( ...
                contrastRaw{b},idxStart,idxEnd,nRaw,nSeg);
            if ~okC
                error('M9 beam %d primary Contrast is not aligned.',b);
            end
            if ~okS
                strengthAligned{b}=nan(nSeg,1);
                strengthSource(b)=strengthSource(b)+"_unavailable_or_unaligned";
                SlantedStrengthAvailable_S(i)=false;
            end
            amplitudeProbe=[amplitudeProbe;strengthAligned{b};contrastAligned{b}]; %#ok<AGROW>
        end
        StrengthSource_S(i,:)=strengthSource;
        sourceProbe=strjoin([strengthSource,contrastSource],";");
        ampScale=resolveDataScale(opt.AmplitudeScale,amplitudeProbe,sourceProbe);
        AmplitudeScaleUsed_S(i)=ampScale;

        beamStrength=cell(1,4); beamContrast=cell(1,4);
        colors={'r','g','b','m'};
        for b=1:4
            s=forceColumnLength(strengthAligned{b},nSeg).*ampScale;
            c=forceColumnLength(contrastAligned{b},nSeg).*ampScale;
            s(s<=0)=NaN; c(c<=0)=NaN;
            s(velocityRejectMask)=NaN; c(velocityRejectMask)=NaN;
            beamStrength{b}=despikeStrengthLocal( ...
                s,opt.Tdespike,opt.MM,opt.Win,velocityRejectMask, ...
                colors{b},opt.PLOT,opt.StrengthContrastMax);
            beamContrast{b}=despikeStrengthLocal( ...
                c,opt.Tdespike,opt.MM,opt.Win,velocityRejectMask, ...
                colors{b},opt.PLOT,opt.StrengthContrastMax);
        end
        BT_RSSI1_S{i}=beamStrength{1}; BT_RSSI2_S{i}=beamStrength{2};
        BT_RSSI3_S{i}=beamStrength{3}; BT_RSSI4_S{i}=beamStrength{4};
        BT_CNT1_S{i}=beamContrast{1}; BT_CNT2_S{i}=beamContrast{2};
        BT_CNT3_S{i}=beamContrast{3}; BT_CNT4_S{i}=beamContrast{4};

        [vbStrengthAligned,vbStrengthOK]=alignSeriesFlexible( ...
            vbStrengthRaw,idxStart,idxEnd,nRaw,nSeg);
        [vbContrastAligned,vbContrastOK]=alignSeriesFlexible( ...
            vbContrastRaw,idxStart,idxEnd,nRaw,nSeg);
        if vbContrastOK && ~isempty(V_RNG_S{i})
            if vbStrengthOK
                vbs=forceColumnLength(vbStrengthAligned,nSeg).*ampScale;
            else
                vbs=nan(nSeg,1);
            end
            vbc=forceColumnLength(vbContrastAligned,nSeg).*ampScale;
            vbs(vbs<=0)=NaN; vbc(vbc<=0)=NaN;
            vbs(velocityRejectMask)=NaN; vbc(velocityRejectMask)=NaN;
            vbs=despikeStrengthLocal(vbs,opt.Tdespike,opt.MM,opt.Win, ...
                velocityRejectMask,'k',opt.PLOT,opt.StrengthContrastMax);
            vbc=despikeStrengthLocal(vbc,opt.Tdespike,opt.MM,opt.Win, ...
                velocityRejectMask,'k',opt.PLOT,opt.StrengthContrastMax);
            if calibrateVB
                validateVBGainOptions(opt);
                low=isfinite(vbs) & vbs<opt.VBStrengthLowThreshold;
                vbs(low)=vbs(low)+opt.VBGain;
                vbs(vbs>opt.VBStrengthHighThreshold)=NaN;
                low=isfinite(vbc) & vbc<opt.VBContrastLowThreshold;
                vbc(low)=vbc(low)+opt.VBGain;
                vbc(vbc>opt.VBContrastHighThreshold)=NaN;
                VBcalibrationApplied_S(i)=true;
            end
            VB_RSSI_S{i}=vbs;
            VB_CNT_S{i}=vbc;
            VBAvailable_S(i)=any(isfinite(vbc) & isfinite(V_RNG_S{i}));
            VBStrengthAvailable_S(i)=any(isfinite(vbs) & ...
                isfinite(V_RNG_S{i}));
        end

        strengthMatrix=[beamStrength{1},beamStrength{2}, ...
            beamStrength{3},beamStrength{4}];
        contrastMatrix=[beamContrast{1},beamContrast{2}, ...
            beamContrast{3},beamContrast{4}];
        BT_RSSI_S{i}=nanmeanLocal(strengthMatrix,2);
        BT_RSSI_Sstd{i}=nanstdLocal(strengthMatrix,0,2);
        BT_CNT_S{i}=nanmeanLocal(contrastMatrix,2);
        BT_CNT_Sstd{i}=nanstdLocal(contrastMatrix,0,2);
        BT_RSSImean_S(i,:)=[safeMean(beamStrength{1}),safeMean(beamStrength{2}), ...
            safeMean(beamStrength{3}),safeMean(beamStrength{4})];
        BT_CNTmean_S(i,:)=[safeMean(beamContrast{1}),safeMean(beamContrast{2}), ...
            safeMean(beamContrast{3}),safeMean(beamContrast{4})];
        stdSpatialRSSI_S(i)=safeMean(BT_RSSI_Sstd{i});
        stdSpatialCNT_S(i)=safeMean(BT_CNT_Sstd{i});
        VB_RSSI_Sm(i)=safeMean(VB_RSSI_S{i});
        VB_RSSI_Sstd(i)=safeStd(VB_RSSI_S{i});
        VB_CNT_Sm(i)=safeMean(VB_CNT_S{i});
        VB_CNT_Sstd(i)=safeStd(VB_CNT_S{i});

        % Water-profile configuration and temperature.
        WF_S(i)=getScalarPreferred(DepthRNGout, ...
            {'WF_12','WF_1200','WF_S'},i,nFiles,NaN);
        if ~isfinite(WF_S(i))
            [wfRaw,~]=findM9ScalarOrSeries(Sload,"wf");
            WF_S(i)=convertLengthScalar(wfRaw,referenceDepth);
        end
        WSmean(i)=getScalarPreferred(DepthRNGout, ...
            {'WS_12','WS_1200','WSmean'},i,nFiles,NaN);
        if ~isfinite(WSmean(i))
            [wsRaw,~]=findM9ScalarOrSeries(Sload,"ws");
            WSmean(i)=convertCellSizeScalar(wsRaw);
        end
        WT_3(i)=getScalarPreferred(DepthRNGout, ...
            {'WT_12','Temp1','WT_3'},i,nFiles,NaN);
        if ~isfinite(WT_3(i))
            [tempRaw,~]=findM9ScalarOrSeries(Sload,"temperature");
            WT_3(i)=safeMean(tempRaw);
        end

        [pressureAligned,pressureOK]=alignSeriesFlexible( ...
            pressureRaw,idxStart,idxEnd,nRaw,nSeg);
        if pressureOK
            pressureAligned=forceColumnLength(pressureAligned,nSeg);
            pressureAligned(pressureAligned<=0)=NaN;
            pressureAligned(velocityRejectMask)=NaN;
            Pressure_S{i}=pressureAligned;
            SL_S(i)=safeMean(10.*log10(pressureAligned.*opt.PressureToPa));
        end

        % SonTek 1-3 and 2-4 beam-pair bed angles.
        r1=RNG11_S{i}; r2=RNG22_S{i}; r3=RNG33_S{i}; r4=RNG44_S{i};
        D1=r1.*tand(opt.theta); D2=r2.*tand(opt.theta);
        D3=r3.*tand(opt.theta); D4=r4.*tand(opt.theta);
        psi11=90-opt.theta-atand((r1-r3)./(D1+D3));
        psi33=90-opt.theta-atand((r3-r1)./(D1+D3));
        psi22=90-opt.theta-atand((r2-r4)./(D2+D4));
        psi44=90-opt.theta-atand((r4-r2)./(D2+D4));
        psi1_S(i)=safeMean(psi11); psi2_S(i)=safeMean(psi22);
        psi3_S(i)=safeMean(psi33); psi4_S(i)=safeMean(psi44);
        alfaX=nanmeanLocal([abs(90-opt.theta-psi11), ...
            abs(90-opt.theta-psi33)],2);
        alfaY=nanmeanLocal([abs(90-opt.theta-psi22), ...
            abs(90-opt.theta-psi44)],2);
        alfaXmean_S(i)=safeMean(alfaX);
        alfaYmean_S(i)=safeMean(alfaY);

        betaX=nan(nSeg,1); betaY=nan(nSeg,1);
        if VBAvailable_S(i)
            vb=V_RNG_S{i};
            ni1=90+acotd(D1./(-vb+r1));
            ni2=90+acotd(D2./(-vb+r2));
            ni3=90+acotd(D3./(-vb+r3));
            ni4=90+acotd(D4./(-vb+r4));
            niAbs=[abs(ni1),abs(ni2),abs(ni3),abs(ni4)];
            ni_mean_S(i)=max(nanmeanLocal(niAbs,1));
            ni_std_S(i)=max(nanstdLocal(niAbs,0,1));
            betaX=nanmeanLocal([abs(90-ni1),abs(90-ni3)],2);
            betaY=nanmeanLocal([abs(90-ni2),abs(90-ni4)],2);
            betaXmean_S(i)=safeMean(betaX);
            betaYmean_S(i)=safeMean(betaY);
        end

        R11=r1./cosd(opt.theta); R22=r2./cosd(opt.theta);
        R33=r3./cosd(opt.theta); R44=r4./cosd(opt.theta);
        R1_S(i)=safeMean(R11); R2_S(i)=safeMean(R22);
        R3_S(i)=safeMean(R33); R4_S(i)=safeMean(R44);
        Spread11=20.*log10(R11); Spread22=20.*log10(R22);
        Spread33=20.*log10(R33); Spread44=20.*log10(R44);
        Spread1_S(i)=safeMean(Spread11); Spread2_S(i)=safeMean(Spread22);
        Spread3_S(i)=safeMean(Spread33); Spread4_S(i)=safeMean(Spread44);
        SpreadMS(i)=safeMean([Spread11;Spread22;Spread33;Spread44]);
        if VBAvailable_S(i)
            SpreadVB=20.*log10(V_RNG_S{i});
            SpreadVBmean(i)=safeMean(SpreadVB);
        else
            SpreadVB=[];
        end

        BT3MHzPL_S(i)=resolveM9AcousticScalar( ...
            Sload,"pulse3",BT3MHzPLOption(i));
        BT05MHzPL_S(i)=resolveM9AcousticScalar( ...
            Sload,"pulse05",BT05MHzPLOption(i));
        width3=resolveM9AcousticScalar(Sload,"width3",opt.width3MHz);
        omega3=resolveM9AcousticScalar(Sload,"omega3",opt.omega3MHz);
        width05=resolveM9AcousticScalar(Sload,"width05",opt.width05MHz);
        omega05=resolveM9AcousticScalar(Sload,"omega05",opt.omega05MHz);
        omega3=normaliseBeamOpeningRadians(omega3);
        omega05=normaliseBeamOpeningRadians(omega05);
        if logical(opt.ExactLegacyFootprintOmega)
            omega3Footprint=deg2rad(omega3);
            omega05Footprint=deg2rad(omega05);
        else
            omega3Footprint=omega3;
            omega05Footprint=omega05;
        end

        [Dmin_S_1(i),Dmax_S_1(i),Area_S_1(i),distance_S_1(i),Lp_S_1(i)]= ...
            SamplingVolume(safeMean(r1),opt.theta,width3, ...
                rad2deg(omega3)./2,opt.DistTrans3MHz,BT3MHzPL_S(i));
        [Dmin_S_2(i),Dmax_S_2(i),Area_S_2(i),distance_S_2(i),Lp_S_2(i)]= ...
            SamplingVolume(safeMean(r2),opt.theta,width3, ...
                rad2deg(omega3)./2,opt.DistTrans3MHz,BT3MHzPL_S(i));
        [Dmin_S_3(i),Dmax_S_3(i),Area_S_3(i),distance_S_3(i),Lp_S_3(i)]= ...
            SamplingVolume(safeMean(r3),opt.theta,width3, ...
                rad2deg(omega3)./2,opt.DistTrans3MHz,BT3MHzPL_S(i));
        [Dmin_S_4(i),Dmax_S_4(i),Area_S_4(i),distance_S_4(i),Lp_S_4(i)]= ...
            SamplingVolume(safeMean(r4),opt.theta,width3, ...
                rad2deg(omega3)./2,opt.DistTrans3MHz,BT3MHzPL_S(i));
        if VBAvailable_S(i)
            [Dmin05MHz_S(i),Dmax05MHz_S(i),Area05MHz_S(i), ...
                distance05MHz_S(i),Lp05MHz_S(i)]=SamplingVolume( ...
                    RNG_VB_Sm(i),0,width05,rad2deg(omega05)./2,0, ...
                    BT05MHzPL_S(i));
        end

        Af1=slantedFootprintDb(BT3MHzPL_S(i),R11,omega3Footprint, ...
            opt.theta,alfaX,alfaY);
        Af2=slantedFootprintDb(BT3MHzPL_S(i),R22,omega3Footprint, ...
            opt.theta,alfaY,alfaX);
        Af3=slantedFootprintDb(BT3MHzPL_S(i),R33,omega3Footprint, ...
            opt.theta,alfaX,alfaY);
        Af4=slantedFootprintDb(BT3MHzPL_S(i),R44,omega3Footprint, ...
            opt.theta,alfaY,alfaX);
        Af1meanS(i)=safeMean(Af1); Af2meanS(i)=safeMean(Af2);
        Af3meanS(i)=safeMean(Af3); Af4meanS(i)=safeMean(Af4);
        if VBAvailable_S(i)
            AfVB=verticalLegacyFootprintDb( ...
                BT05MHzPL_S(i),V_RNG_S{i},omega05Footprint,betaX,betaY);
            AfVBmean_S(i)=safeMean(AfVB);
        else
            AfVB=[];
        end

        [AW_3(i),abswat3_1(i)]=callFalfaw( ...
            WT_3(i),safeMean(r1),opt.f3,opt.alphaw3);
        [~,abswat3_2(i)]=callFalfaw( ...
            WT_3(i),safeMean(r2),opt.f3,opt.alphaw3);
        [~,abswat3_3(i)]=callFalfaw( ...
            WT_3(i),safeMean(r3),opt.f3,opt.alphaw3);
        [~,abswat3_4(i)]=callFalfaw( ...
            WT_3(i),safeMean(r4),opt.f3,opt.alphaw3);
        if VBAvailable_S(i)
            [AW_05(i),abswat05(i)]=callFalfaw( ...
                WT_3(i),RNG_VB_Sm(i),opt.f05,opt.alphaw05);
        end

        % Optional four-beam echo-profile processing.
        echoCells=cell(1,4);
        alphaMeanCells=cell(1,4);
        alphaGradCells=cell(1,4);
        dEDRCells=cell(1,4);
        echoProbe=[];
        for b=1:4
            if ~isempty(echoRaw{b})
                echoCells{b}=subsetEchoMatrixM9( ...
                    echoRaw{b},idxStart,idxEnd,nRaw,nSeg);
                echoProbe=[echoProbe;echoCells{b}(:)]; %#ok<AGROW>
            end
        end
        echoScale=resolveDataScale( ...
            opt.EchoScale,echoProbe,strjoin(echoSource,";"));
        EchoScaleUsed_S(i)=echoScale;
        allEcho=true;
        for b=1:4
            if isempty(echoCells{b})
                allEcho=false;
                continue
            end
            E=double(echoCells{b}).*echoScale;
            E(:,velocityRejectMask)=NaN;
            echoCells{b}=E;
            rangeMean=safeMean(rangeM{b});
            idStop=resolveEchoStopM9(opt.distanceADCPtoWaterSurface, ...
                WF_S(i),WSmean(i),rangeMean,size(E,1));
            [alphaMeanCells{b},alphaGradCells{b},dEDRCells{b}]= ...
                processEchoAttenuationM9(E,idStop,WSmean(i));
        end
        EchoProfilesAvailable_S(i)=allEcho;
        ECH1_S{i}=echoCells{1}; ECH2_S{i}=echoCells{2};
        ECH3_S{i}=echoCells{3}; ECH4_S{i}=echoCells{4};
        alpha1meanS{i}=alphaMeanCells{1}; alpha2meanS{i}=alphaMeanCells{2};
        alpha3meanS{i}=alphaMeanCells{3}; alpha4meanS{i}=alphaMeanCells{4};
        alpha1meanSgrad{i}=alphaGradCells{1}; alpha2meanSgrad{i}=alphaGradCells{2};
        alpha3meanSgrad{i}=alphaGradCells{3}; alpha4meanSgrad{i}=alphaGradCells{4};
        dEDR1_S{i}=dEDRCells{1}; dEDR2_S{i}=dEDRCells{2};
        dEDR3_S{i}=dEDRCells{3}; dEDR4_S{i}=dEDRCells{4};
        alphaMean3000S(i,:)=[safeMean(alpha1meanS{i}),safeMean(alpha2meanS{i}), ...
            safeMean(alpha3meanS{i}),safeMean(alpha4meanS{i})];

        abswat=[abswat3_1(i),abswat3_2(i),abswat3_3(i),abswat3_4(i)];
        slantRanges={R11,R22,R33,R44};
        tlFixed=nan(1,4);
        tlProfile=nan(1,4);
        for b=1:4
            TLvec=(abswat(b)+alphaSed3(i)).*slantRanges{b};
            tlFixed(b)=safeMean(TLvec);
            if ~isempty(alphaMeanCells{b}) && isfinite(WSmean(i)) && WSmean(i)>0
                tlProfile(b)=2.*(WF_S(i)./cosd(90-opt.theta).*abswat(b) + ...
                    nansumLocal(alphaMeanCells{b}.*WSmean(i)./cosd(90-opt.theta)));
            end
        end
        TL1_fixS(i)=tlFixed(1); TL2_fixS(i)=tlFixed(2);
        TL3_fixS(i)=tlFixed(3); TL4_fixS(i)=tlFixed(4);
        TL1(i)=tlProfile(1); TL2(i)=tlProfile(2);
        TL3(i)=tlProfile(3); TL4(i)=tlProfile(4);

        validSlantedContrast=all(cellfun( ...
            @(x) any(isfinite(x)),beamContrast));
        validSlantedStrength=all(cellfun( ...
            @(x) any(isfinite(x)),beamStrength));
        SlantedContrastAvailable_S(i)=validSlantedContrast;
        SlantedStrengthAvailable_S(i)=validSlantedStrength;
        commonSlantedReady=isfinite(SL_S(i)) && isfinite(BT3MHzPL_S(i)) && ...
            BT3MHzPL_S(i)>0 && isfinite(omega3) && omega3>0 && ...
            all(isfinite([Af1meanS(i),Af2meanS(i), ...
                Af3meanS(i),Af4meanS(i)]));
        slantedReady=commonSlantedReady && validSlantedContrast;
        strengthReady=commonSlantedReady && validSlantedStrength;
        AbsoluteBackscatterReady_S(i)=slantedReady;
        StrengthBackscatterReady_S(i)=strengthReady;

        spreads={Spread11,Spread22,Spread33,Spread44};
        footprints={Af1,Af2,Af3,Af4};
        svStrength=cell(1,4); svContrast=cell(1,4);
        if slantedReady || strengthReady
            for b=1:4
                TLvec=(abswat(b)+alphaSed3(i)).*slantRanges{b};
                if strengthReady
                    svStrength{b}=-SL_S(i)+beamStrength{b}+spreads{b}- ...
                        footprints{b}+2.*TLvec-Offset3MHz(b)-opt.LegacyDbCorrection;
                end
                if slantedReady
                    svContrast{b}=-SL_S(i)+beamContrast{b}+spreads{b}- ...
                        footprints{b}+2.*TLvec-Offset3MHz(b)-opt.LegacyDbCorrection;
                end
            end
        end
        SV1_Srssi{i}=svStrength{1}; SV2_Srssi{i}=svStrength{2};
        SV3_Srssi{i}=svStrength{3}; SV4_Srssi{i}=svStrength{4};
        SV1_S{i}=svContrast{1}; SV2_S{i}=svContrast{2};
        SV3_S{i}=svContrast{3}; SV4_S{i}=svContrast{4};

        SVmean1rssi(i)=safeMean(SV1_Srssi{i});
        SVmean2rssi(i)=safeMean(SV2_Srssi{i});
        SVmean3rssi(i)=safeMean(SV3_Srssi{i});
        SVmean4rssi(i)=safeMean(SV4_Srssi{i});
        SVmean1stdRSSI(i)=safeStd(SV1_Srssi{i});
        SVmean2stdRSSI(i)=safeStd(SV2_Srssi{i});
        SVmean3stdRSSI(i)=safeStd(SV3_Srssi{i});
        SVmean4stdRSSI(i)=safeStd(SV4_Srssi{i});
        SVmean1(i)=safeMean(SV1_S{i}); SVmean2(i)=safeMean(SV2_S{i});
        SVmean3(i)=safeMean(SV3_S{i}); SVmean4(i)=safeMean(SV4_S{i});
        SVmean1std(i)=safeStd(SV1_S{i}); SVmean2std(i)=safeStd(SV2_S{i});
        SVmean3std(i)=safeStd(SV3_S{i}); SVmean4std(i)=safeStd(SV4_S{i});

        % Evaluate the optional vertical beam only after all four slanted
        % SV series and summaries have been completed. A VB-specific issue
        % therefore cannot erase otherwise valid 3-MHz results.
        verticalReady=VBAvailable_S(i) && isfinite(SL_S(i)) && ...
            isfinite(BT05MHzPL_S(i)) && BT05MHzPL_S(i)>0 && ...
            isfinite(omega05) && omega05>0 && isfinite(abswat05(i)) && ...
            isfinite(AfVBmean_S(i)) && any(isfinite(VB_CNT_S{i}));
        verticalStrengthReady=verticalReady && ...
            VBStrengthAvailable_S(i) && any(isfinite(VB_RSSI_S{i}));
        VerticalBackscatterReady_S(i)=verticalReady;
        VerticalStrengthBackscatterReady_S(i)=verticalStrengthReady;

        if verticalReady
            TL05=(abswat05(i)+alphaSed05(i)).*V_RNG_S{i};
            TLVBMHz_fixS(i)=safeMean(TL05);
            % The old script applied the mean vertical TL to every sample.
            SV_VB{i}=-opt.SLaddVB+SpreadVB-SL_S(i)+VB_CNT_S{i}- ...
                AfVB+2.*TLVBMHz_fixS(i)-opt.Offset05MHz-opt.LegacyDbCorrection;
            SV_VBmeanS(i)=safeMean(SV_VB{i});
            SV_VBmeanSstd(i)=safeStd(SV_VB{i});
            if verticalStrengthReady
                SV_VBrssi{i}=-opt.SLaddVB+SpreadVB-SL_S(i)+VB_RSSI_S{i}- ...
                    AfVB+2.*TLVBMHz_fixS(i)-opt.Offset05MHz-opt.LegacyDbCorrection;
                SV_VBmeanSrssi(i)=safeMean(SV_VBrssi{i});
                SV_VBmeanSstdRSSI(i)=safeStd(SV_VBrssi{i});
            end
        end

        if slantedReady
            svContrastMatrix=[SV1_S{i},SV2_S{i},SV3_S{i},SV4_S{i}];
            stdSvSpatialALL_S{i}=nanstdLocal(svContrastMatrix,0,2);
            stdSvSpatial_S(i)=safeMean(stdSvSpatialALL_S{i});
            S_Sall{i}=nanmeanLocal(svContrastMatrix,2);
            S_S(i)=safeMean([SVmean1(i),SVmean2(i),SVmean3(i),SVmean4(i)]);
            S_Sstd(i)=safeMean([SVmean1std(i),SVmean2std(i), ...
                SVmean3std(i),SVmean4std(i)]);
        end
        if strengthReady
            svStrengthMatrix=[SV1_Srssi{i},SV2_Srssi{i},SV3_Srssi{i},SV4_Srssi{i}];
            stdSvSpatialRSSIALL_S{i}=nanstdLocal(svStrengthMatrix,0,2);
            stdSvSpatialRSSI_S(i)=safeMean(stdSvSpatialRSSIALL_S{i});
            S_Srssi(i)=safeMean([SVmean1rssi(i),SVmean2rssi(i), ...
                SVmean3rssi(i),SVmean4rssi(i)]);
        end

        if slantedReady
            AcousticMessage_S(i)="OK: primary Contrast SV for four 3-MHz beams from "+ ...
                SlantedLayout_S(i);
            if strengthReady
                AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                    "; parallel Strength/RSSI SV ready";
            else
                AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                    "; Strength/RSSI unavailable (primary Contrast unaffected)";
            end
            if verticalReady
                AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                    " + 0.5-MHz vertical beam from "+BackscatterLayout_S(i);
            elseif VBAvailable_S(i)
                AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                    "; vertical raw data present but calibration inputs missing";
            else
                AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                    "; vertical beam not available";
            end
        else
            missing=missingAcousticInputs( ...
                SL_S(i),BT3MHzPL_S(i),omega3, ...
                [Af1meanS(i),Af2meanS(i),Af3meanS(i),Af4meanS(i)], ...
                validSlantedContrast);
            AcousticMessage_S(i)="Raw M9 backscatter processed; absolute SV unavailable: "+missing;
        end
        if velocityMaskMode~="none" && ~(maskExact && maskAligned)
            AcousticMessage_S(i)=AcousticMessage_S(i)+ ...
                "; exact BT mask unavailable, so no BT-velocity rejection was applied";
        end

        if logical(opt.ApplyMethod2Regression)
            dvals=[D50_14S(i),D50_05S(i)];
            % Do not replace a missing mobile grain-size input with the
            % other percentile. This follows the corrected physical-sample
            % handling: incomplete or zero mobile D50 information remains
            % unavailable for Method 2.
            if all(isfinite(dvals) & dvals>0)
                D50S(i)=mean(dvals)./1000;
            end
            if exist('TrateMethod2MregressionS','file')~=2
                Method2Status_S(i)="helper_missing";
            elseif ~all(isfinite([D50S(i),BTmodf_S(i),FilteredDataN_S(i), ...
                    SV_VBmeanS(i),S_S(i)]))
                Method2Status_S(i)="missing_inputs";
            else
                xSAND_S=[1.00015413749255,14.3704460790743, ...
                    9.66849951016796e-05,0.933859634629615, ...
                    -4.31752277651315,0.0452650467373799];
                xSANDgravel_S=[0.734771032680345,0.722426698958664, ...
                    0.686927942569311,0.339209927392813, ...
                    0.497412127531723,0.247711820177772];
                xGravel_S=[1.31418592253151,9.80129512608736, ...
                    0.00325463025525656,0.503791202031056, ...
                    -11.6082209460445,0.0668879652739897];
                qsMethod2regression_S(i)=TrateMethod2MregressionS( ...
                    D50S(i),xSAND_S,xSANDgravel_S,xGravel_S, ...
                    BTmodf_S(i),FilteredDataN_S(i),SV_VBmeanS(i),S_S(i));
                if isfinite(qsMethod2regression_S(i))
                    Method2Status_S(i)="OK";
                else
                    Method2Status_S(i)="helper_returned_NaN";
                end
            end
        end

        Message(i)="OK";

    catch ME
        Message(i)="Error: "+string(ME.message);
        warning('M9 backscatter failed for row %d, file %s:\n%s', ...
            i,ADCPfilename(i),ME.message);
    end
end

%% Output table and flat structure

BackscatterTable=table( ...
    InputRow,ADCPfilename,DistanceFromLeftBank,StartADCP,EndADCP, ...
    StartUsed_S,EndUsed_S,NensemblesRaw_S,NensemblesUsed_S, ...
    NfilteredBTvelocity_S,BTvelocityMaskMode_S, ...
    BTvelocityFilterSource_S,BTvelocityFilterExact_S, ...
    BTvelocityFilterAligned_S,RangeScaleToMeters_S, ...
    RangeOriginalUnit_S,VBRangeScaleToMeters_S,VBRangeOriginalUnit_S, ...
    AmplitudeScaleUsed_S,EchoScaleUsed_S,SlantedLayout_S, ...
    BackscatterLayout_S,PressureSource_S, ...
    BT_CNTmean_S,BT_RSSImean_S,SlantedContrastAvailable_S, ...
    SlantedStrengthAvailable_S, ...
    RNG_VB_Sm,RNG_VB_Sstd,VB_RSSI_Sm,VB_RSSI_Sstd, ...
    VB_CNT_Sm,VB_CNT_Sstd,SL_S,R1_S,R2_S,R3_S,R4_S, ...
    SVmean1,SVmean2,SVmean3,SVmean4,S_S,S_Srssi, ...
    SV_VBmeanS,SV_VBmeanSrssi,stdSvSpatial_S, ...
    stdSvSpatialRSSI_S,AbsoluteBackscatterReady_S, ...
    StrengthBackscatterReady_S,VerticalBackscatterReady_S, ...
    VerticalStrengthBackscatterReady_S,VBAvailable_S, ...
    VBStrengthAvailable_S,EchoProfilesAvailable_S, ...
    VBcalibrationApplied_S,AcousticMessage_S,BTmodf_S, ...
    FilteredDataN_S,BTstd_S,D50S,qsMethod2regression_S, ...
    Method2Status_S,FileFound,MatchedFile,Message);

BackscatterOut=struct();
% Keep the legacy SonTek names as aliases while using the clearer primary
% names above. The aliases point to the same per-measurement cell arrays.
SV_VBSrssi=SV_VBrssi;
SV_VBS=SV_VB;
varsToExport={ ...
    'InputRow','ADCPfilename','DistanceFromLeftBank','StartADCP','EndADCP', ...
    'StartUsed_S','EndUsed_S','NensemblesRaw_S','NensemblesUsed_S', ...
    'NfilteredBTvelocity_S','BTvelocityFilterMask_S', ...
    'BTvelocityMaskMode_S','BTvelocityFilterSource_S', ...
    'BTvelocityFilterExact_S','BTvelocityFilterAligned_S', ...
    'RangeScaleToMeters_S','RangeOriginalUnit_S', ...
    'VBRangeScaleToMeters_S','VBRangeOriginalUnit_S', ...
    'AmplitudeScaleUsed_S','EchoScaleUsed_S', ...
    'RangeSource_S','StrengthSource_S','ContrastSource_S','EchoSource_S', ...
    'VBRangeSource_S','VBStrengthSource_S','VBContrastSource_S', ...
    'SlantedLayout_S','BackscatterLayout_S','SlantedInputsAvailable_S', ...
    'SlantedContrastAvailable_S','SlantedStrengthAvailable_S', ...
    'Depth_S','DepthS','Depth_Sstd','stdSpatial_S','stdSpatialALL_S', ...
    'Depth1200','Depth1200std','stdSpatial_1200','stdSpatialALL_1200', ...
    'RNG11_S','RNG22_S','RNG33_S','RNG44_S','V_RNG_S', ...
    'RNG_VB_Sm','RNG_VB_Sstd', ...
    'BT_RSSI1_S','BT_RSSI2_S','BT_RSSI3_S','BT_RSSI4_S', ...
    'BT_CNT1_S','BT_CNT2_S','BT_CNT3_S','BT_CNT4_S', ...
    'VB_RSSI_S','VB_CNT_S','BT_RSSI_S','BT_RSSI_Sstd', ...
    'BT_CNT_S','BT_CNT_Sstd','BT_RSSImean_S','BT_CNTmean_S', ...
    'VB_RSSI_Sm','VB_RSSI_Sstd', ...
    'VB_CNT_Sm','VB_CNT_Sstd','stdSpatialRSSI_S','stdSpatialCNT_S', ...
    'WF_S','WSmean','WT_3','Pressure_S','PressureSource_S','SL_S', ...
    'psi1_S','psi2_S','psi3_S','psi4_S', ...
    'alfaXmean_S','alfaYmean_S','betaXmean_S','betaYmean_S', ...
    'ni_mean_S','ni_std_S','R1_S','R2_S','R3_S','R4_S', ...
    'Spread1_S','Spread2_S','Spread3_S','Spread4_S', ...
    'SpreadMS','SpreadVBmean','BT3MHzPL_S','BT05MHzPL_S', ...
    'Dmin_S_1','Dmin_S_2','Dmin_S_3','Dmin_S_4', ...
    'Dmax_S_1','Dmax_S_2','Dmax_S_3','Dmax_S_4', ...
    'Area_S_1','Area_S_2','Area_S_3','Area_S_4', ...
    'distance_S_1','distance_S_2','distance_S_3','distance_S_4', ...
    'Lp_S_1','Lp_S_2','Lp_S_3','Lp_S_4', ...
    'Dmin05MHz_S','Dmax05MHz_S','Area05MHz_S', ...
    'distance05MHz_S','Lp05MHz_S', ...
    'Af1meanS','Af2meanS','Af3meanS','Af4meanS','AfVBmean_S', ...
    'ECH1_S','ECH2_S','ECH3_S','ECH4_S', ...
    'alpha1meanS','alpha2meanS','alpha3meanS','alpha4meanS', ...
    'alpha1meanSgrad','alpha2meanSgrad','alpha3meanSgrad', ...
    'alpha4meanSgrad','dEDR1_S','dEDR2_S','dEDR3_S','dEDR4_S', ...
    'alphaMean3000S','AW_3','AW_05','abswat3_1','abswat3_2', ...
    'abswat3_3','abswat3_4','abswat05', ...
    'TL1','TL2','TL3','TL4','TL1_fixS','TL2_fixS', ...
    'TL3_fixS','TL4_fixS','TLVBMHz_fixS', ...
    'SV1_Srssi','SV2_Srssi','SV3_Srssi','SV4_Srssi', ...
    'SV1_S','SV2_S','SV3_S','SV4_S','SV_VBrssi','SV_VB', ...
    'SV_VBSrssi','SV_VBS', ...
    'SVmean1rssi','SVmean2rssi','SVmean3rssi','SVmean4rssi', ...
    'SVmean1stdRSSI','SVmean2stdRSSI','SVmean3stdRSSI', ...
    'SVmean4stdRSSI','SVmean1','SVmean2','SVmean3','SVmean4', ...
    'SVmean1std','SVmean2std','SVmean3std','SVmean4std', ...
    'SV_VBmeanSrssi','SV_VBmeanSstdRSSI','SV_VBmeanS', ...
    'SV_VBmeanSstd','stdSvSpatial_S','stdSvSpatialRSSI_S', ...
    'stdSvSpatialALL_S','stdSvSpatialRSSIALL_S', ...
    'S_S','S_Sstd','S_Srssi','S_Sall', ...
    'AbsoluteBackscatterReady_S','StrengthBackscatterReady_S', ...
    'VerticalBackscatterReady_S','VerticalStrengthBackscatterReady_S', ...
    'VBAvailable_S','VBStrengthAvailable_S', ...
    'EchoProfilesAvailable_S','VBcalibrationApplied_S', ...
    'AcousticMessage_S','BTmodf_S','FilteredDataN_S','BTstd_S', ...
    'D50_14S','D50_05S','D50S','qsMethod2regression_S', ...
    'Method2Status_S','DummyValue','FileFound','MatchedFile','Message', ...
    'BackscatterTable'};
for k=1:numel(varsToExport)
    name=varsToExport{k};
    BackscatterOut.(name)=eval(name);
end

BackscatterOut.Instrument="SonTek M9";
BackscatterOut.ADCPfilenameField=ADCPfilenameField;
BackscatterOut.SlantedFrequencyMHz=opt.f3./1000;
BackscatterOut.VerticalFrequencyMHz=opt.f05./1000;
BackscatterOut.PrimaryAmplitude="Contrast";
BackscatterOut.SecondaryAmplitude="Strength/RSSI";
BackscatterOut.DepthSource="DepthRNGout";
BackscatterOut.StartEndMode=startEndMode;
BackscatterOut.VelocityMaskMode=velocityMaskMode;
BackscatterOut.Options=opt;
processedRows=Message=="OK";
BackscatterOut.M9meta=struct( ...
    'SlantedFrequencyMHz',opt.f3./1000, ...
    'VerticalFrequencyMHz',opt.f05./1000, ...
    'AbsoluteBackscatterReady',any(processedRows) && ...
        all(AbsoluteBackscatterReady_S(processedRows)), ...
    'VerticalBackscatterReady',any(processedRows) && ...
        all(VerticalBackscatterReady_S(processedRows) | ...
            ~VBAvailable_S(processedRows)), ...
    'InputRows',nFiles, ...
    'ProcessedRows',sum(processedRows), ...
    'ErrorRows',sum(startsWith(Message,"Error")), ...
    'CalibrationSource',"legacy SonTek M9 configuration", ...
    'ExactLegacyFootprintOmega',logical(opt.ExactLegacyFootprintOmega), ...
    'SlantedOffsetsDB',Offset3MHz, ...
    'VerticalOffsetDB',opt.Offset05MHz, ...
    'SedimentAttenuationApplied', ...
        any(isfinite(alphaSed3) & alphaSed3~=0) || ...
        any(isfinite(alphaSed05) & alphaSed05~=0), ...
    'UsesPressureSourceLevel',true, ...
    'FabricatedRDITransmitFields',false);

if logical(opt.AssignToCaller), assignStructToWorkspace(BackscatterOut,'caller'); end
if logical(opt.AssignToBase), assignStructToWorkspace(BackscatterOut,'base'); end

if logical(opt.Verbose)
    fprintf('\n=== SonTek M9 five-beam backscatter summary ===\n');
    fprintf('Input rows: %d\n',nFiles);
    fprintf('Successfully processed: %d/%d\n',sum(Message=="OK"),nFiles);
    fprintf('Exact/aligned BT masks: %d/%d\n', ...
        sum(BTvelocityFilterExact_S & BTvelocityFilterAligned_S),nFiles);
    fprintf('3-MHz absolute SV ready: %d/%d rows.\n', ...
        sum(AbsoluteBackscatterReady_S),nFiles);
    fprintf('0.5-MHz vertical raw data available: %d/%d rows.\n', ...
        sum(VBAvailable_S),nFiles);
    fprintf('0.5-MHz vertical absolute SV ready: %d/%d rows.\n', ...
        sum(VerticalBackscatterReady_S),nFiles);
    fprintf('Four-beam echo profiles available: %d/%d rows.\n', ...
        sum(EchoProfilesAvailable_S),nFiles);
    fprintf('Vertical gain correction applied: %d/%d rows.\n', ...
        sum(VBcalibrationApplied_S),nFiles);
    fprintf('No files were saved. All time series are in cells.\n');
    fprintf('=================================================\n');
    disp(BackscatterTable);
end

end

%% Helpers

function [names,sourceField]=resolveM9FilenameInput(ADCPin)
requested={'M9filename','M9FileName','Sontekfilename','SonTekfilename', ...
    'ADCPfilenameM9','ADCPfilename_S','ADCPfilename'};
actual=fieldnames(ADCPin);
names=strings(0,1); sourceField="";
for k=1:numel(requested)
    idx=find(strcmpi(actual,requested{k}),1);
    if isempty(idx), continue; end
    candidate=string(ADCPin.(actual{idx})(:));
    if isempty(candidate) || ~any(strlength(strtrim(candidate))>0), continue; end
    names=strtrim(candidate);
    sourceField=string(actual{idx});
    return
end
error(['ADCPin must contain ADCPfilename or an M9-specific filename ', ...
    'field such as M9filename or Sontekfilename.']);
end

function fileList=normaliseFileList(inputValue)
if ischar(inputValue) || isstring(inputValue)
    inputFolder=char(string(inputValue));
    fileList=dir(fullfile(inputFolder,'*.mat'));
    % A SonTek project folder often contains the processed MAT files one or
    % more levels below the selected root. Search recursively only when the
    % selected folder itself has no MAT files.
    if isempty(fileList)
        fileList=dir(fullfile(inputFolder,'**','*.mat'));
    end
elseif isstruct(inputValue)
    fileList=inputValue;
else
    fileList=struct([]);
end
end

function matchedFile=matchADCPfile(fileList,targetName)
matchedFile="";
if isempty(fileList) || ~isfield(fileList,'name'), return; end
target=cleanFileName(targetName);
[~,targetBase,~]=fileparts(char(target));
names=string({fileList.name}');
cleanNames=strings(numel(names),1);
bases=strings(numel(names),1);
for k=1:numel(names)
    cleanNames(k)=cleanFileName(names(k));
    [~,b,~]=fileparts(char(cleanNames(k)));
    bases(k)=string(b);
end
idx=find(strcmpi(cleanNames,target),1);
if isempty(idx), idx=find(strcmpi(bases,string(targetBase)),1); end
if isempty(idx), idx=find(startsWith(lower(bases),lower(string(targetBase))),1); end
if isempty(idx), idx=find(contains(lower(bases),lower(string(targetBase))),1); end
if isempty(idx), return; end
if isfield(fileList,'folder') && ~isempty(fileList(idx).folder)
    matchedFile=string(fullfile(fileList(idx).folder,fileList(idx).name));
else
    matchedFile=string(fileList(idx).name);
end
end

function s=cleanFileName(s)
s=strtrim(string(s));
s=erase(s,"'");
s=erase(s,'"');
end

function [first,last]=resolveStartEnd(startValue,endValue,nRaw,mode)
mode=lower(strtrim(string(mode)));
if mode=="trimcounts" || mode=="trim"
    if ~isfinite(startValue) || startValue<0, nStart=0; else, nStart=round(startValue); end
    if ~isfinite(endValue) || endValue<0, nEnd=0; else, nEnd=round(endValue); end
    first=1+nStart;
    last=nRaw-nEnd;
else
    if ~isfinite(startValue) || startValue<1, first=1; else, first=round(startValue); end
    if ~isfinite(endValue) || endValue<1, last=nRaw; else, last=round(endValue); end
end
first=max(1,min(first,nRaw));
last=max(1,min(last,nRaw));
if last<first
    error('Invalid StartADCP/EndADCP interval: %d:%d.',first,last);
end
end

function tf=isNoDataRow(nodata,i,nFiles)
tf=false;
if isempty(nodata), return; end
if islogical(nodata)
    if isscalar(nodata), tf=nodata;
    elseif numel(nodata)==nFiles, tf=logical(nodata(i));
    end
else
    vals=double(nodata(:));
    if numel(vals)==nFiles && all(isfinite(vals) & ismember(vals,[0 1]))
        tf=logical(vals(i));
    else
        tf=any(round(vals(isfinite(vals)))==i);
    end
end
end

function x=getVectorFieldAny(S,names,defaultValue,n)
x=defaultValue;
if ~isempty(S) && isstruct(S)
    for k=1:numel(names)
        if isfield(S,names{k}), x=S.(names{k}); break; end
    end
end
x=prepareVector(x,n);
end

function x=prepareVector(x,n)
if isstring(x) || iscell(x), x=str2double(string(x)); end
x=double(x); x=x(:);
if isempty(x), x=nan(n,1);
elseif numel(x)==1 && n>1, x=repmat(x,n,1);
elseif numel(x)<n, x=[x;nan(n-numel(x),1)];
elseif numel(x)>n, x=x(1:n);
end
end

function x=prepareFourVector(x)
x=double(x(:).');
if isempty(x), x=zeros(1,4);
elseif isscalar(x), x=repmat(x,1,4);
elseif numel(x)<4, x=[x,repmat(x(end),1,4-numel(x))];
else, x=x(1:4);
end
end

function M=getMatrixFieldAny(S,names,nRows,nCols)
M=[];
for k=1:numel(names)
    if isfield(S,names{k}), M=double(S.(names{k})); break; end
end
if isempty(M) || ~isequal(size(M),[nRows,nCols])
    error('DepthRNGout must contain an N x 4 depth summary matrix.');
end
end

function C=getCellFieldAny(S,names,nRows)
C=[];
for k=1:numel(names)
    if isfield(S,names{k}), C=S.(names{k}); break; end
end
if ~iscell(C) || numel(C)~=nRows
    error('DepthRNGout must contain one spatial-depth cell per measurement.');
end
C=C(:);
end

function names=rangeOutputAliases(beam)
switch beam
    case 1, names={'RNG11_1200_','RNG11_S','RNG11_1200','RNG11_'};
    case 2, names={'RNG22_1200_','RNG22_S','RNG22_1200','RNG22_'};
    case 3, names={'RNG33_1200_','RNG33_S','RNG33_1200','RNG33_'};
    otherwise, names={'RNG44_1200_','RNG44_S','RNG44_1200','RNG44_'};
end
end

function val=getScalarPreferred(S,names,i,nFiles,defaultValue)
val=defaultValue;
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    if iscell(C)
        if numel(C)>=i, val=safeMean(C{i}); end
    elseif isnumeric(C) || islogical(C)
        C=double(C);
        if isvector(C) && numel(C)>=i, val=C(i);
        elseif size(C,1)==nFiles, val=safeMean(C(i,:));
        elseif size(C,2)==nFiles, val=safeMean(C(:,i));
        end
    end
    if isfinite(val), return; end
end
end

function x=getAlignedSeries(S,names,i,nFiles,idxStart,idxEnd,nRaw)
x=[];
nSeg=idxEnd-idxStart+1;
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k}); candidate=[];
    if iscell(C)
        if numel(C)>=i && ~isempty(C{i}), candidate=double(C{i}(:)); end
    elseif isnumeric(C) || islogical(C)
        C=double(C);
        if isvector(C)
            if nFiles==1 || numel(C)~=nFiles, candidate=C(:); end
        elseif size(C,1)==nFiles
            candidate=C(i,:).';
        elseif size(C,2)==nFiles
            candidate=C(:,i);
        end
    end
    if isempty(candidate), continue; end
    if numel(candidate)==nSeg, x=candidate; return;
    elseif numel(candidate)>=nRaw && numel(candidate)>=idxEnd
        x=candidate(idxStart:idxEnd); return;
    end
end
end

function [x,source]=findM9BeamField(S,kind,beam)
b=char(string(beam));
switch string(kind)
    case "range"
        leaves={['RNG' b],['RNG' b b],['Range' b]};
        prefixes={{'BS_adp'},{'T','BS_adp'},{'BT'},{'T','BT'},{} };
        requireMatrix=false;
    case "strength"
        leaves={['BT_Strength' b],['BTstrength' b],['BT_RSSI' b],['BT_AMP' b]};
        prefixes={{'BS_adp'},{'T','BS_adp'},{'BT'},{'T','BT'},{} };
        requireMatrix=false;
    case "contrast"
        leaves={['BT_Contrast' b],['BTcontrast' b],['BT_CNT' b]};
        prefixes={{'BS_adp'},{'T','BS_adp'},{'BT'},{'T','BT'},{} };
        requireMatrix=false;
    otherwise
        leaves={['ECH' b],['Echo' b]};
        prefixes={{'BS_adp'},{'T','BS_adp'},{'ECHO'},{'T','ECHO'},{} };
        requireMatrix=true;
end
[x,source]=findByPrefixesAndLeaves(S,prefixes,leaves,requireMatrix);
end

function [rangeValue,strengthValue,contrastValue, ...
    rangeSource,strengthSource,contrastSource,layout]= ...
    findM9VerticalBundle(S)
rangeValue=[]; strengthValue=[]; contrastValue=[];
rangeSource=""; strengthSource=""; contrastSource=""; layout="not_found";

% Use one internally consistent source whenever possible. A real VB
% structure is preferred. Range plus Contrast is the required primary
% bundle; Strength/RSSI is optional. Legacy files such as the supplied
% sample store the vertical quantities inside BS_adp instead.
bundlePrefixes={{'VB'},{'T','VB'},{'BS_adp'},{'T','BS_adp'},{} };
rangeLeaves={'Range','RangeRaw','Depth','StandardDepth','VBdepth', ...
    'VB_RNGM','VBDepth','VRNG1'};
strengthLeaves={'Strength','RSSI','RSSIAmp','VB_RSSI2M','VB_RSSI','VBstrength'};
contrastLeaves={'Contrast','VB_CNTM','VB_Contrast','VBcontrast'};
for p=1:numel(bundlePrefixes)
    thisPrefix={bundlePrefixes{p}};
    [r,rs]=findByPrefixesAndLeaves(S,thisPrefix,rangeLeaves,false);
    [a,as]=findByPrefixesAndLeaves(S,thisPrefix,strengthLeaves,false);
    [c,cs]=findByPrefixesAndLeaves(S,thisPrefix,contrastLeaves,false);
    if ~isempty(r) && ~isempty(c)
        rangeValue=r; strengthValue=a; contrastValue=c;
        rangeSource=rs; strengthSource=as; contrastSource=cs;
        layout=summariseSourceLayout([rs,as,cs]);
        return
    end
end

% Last-resort independent search supports partially standardized files,
% while retaining exact source names so a mixed layout is visible.
[rangeValue,rangeSource]=findM9VerticalField(S,"range");
[strengthValue,strengthSource]=findM9VerticalField(S,"strength");
[contrastValue,contrastSource]=findM9VerticalField(S,"contrast");
if ~isempty(rangeValue) || ~isempty(strengthValue) || ~isempty(contrastValue)
    layout=summariseSourceLayout( ...
        [rangeSource,strengthSource,contrastSource]);
end
end

function [x,source]=findM9VerticalField(S,kind)
switch string(kind)
    case "range"
        leaves={'Range','Depth','StandardDepth','VB_RNGM','VBdepth','VBDepth','VRNG1'};
    case "strength"
        leaves={'RSSI','Strength','RSSIAmp','VB_RSSI2M','VB_RSSI','VBstrength'};
    otherwise
        leaves={'Contrast','VB_CNTM','VB_Contrast','VBcontrast'};
end
prefixes={{'VB'},{'T','VB'},{'BS_adp'},{'T','BS_adp'},{} };
[x,source]=findByPrefixesAndLeaves(S,prefixes,leaves,false);
end

function [x,source]=findM9ScalarOrSeries(S,kind)
switch string(kind)
    case "pressure"
        leaves={'Pressure','pressure','PRESSURE'};
        prefixes={{'BS_adp'},{'T','BS_adp'},{'Data'},{'T','Data'},{} };
    case "wf"
        leaves={'WF','WaterFirstCell','BlankingDistance'};
        prefixes={{'Data'},{'T','Data'},{'BS_adp'},{'T','BS_adp'},{} };
    case "ws"
        leaves={'WS','CellSize','CellSize_m'};
        prefixes={{'Data'},{'T','Data'},{'BS_adp'},{'T','BS_adp'},{} };
    otherwise
        leaves={'Temperature','WT','Temp','temperature'};
        prefixes={{'Data'},{'T','Data'},{'BS_adp'},{'T','BS_adp'},{} };
end
[x,source]=findByPrefixesAndLeaves(S,prefixes,leaves,false);
end

function layout=summariseSourceLayout(sources)
sources=string(sources(:));
sources=sources(strlength(sources)>0);
if isempty(sources), layout="not_found"; return; end
if all(contains(sources,"BS_adp")), layout="BS_adp";
elseif all(startsWith(sources,"VB.") | startsWith(sources,"T.VB.")), layout="VB";
elseif all(startsWith(sources,"BT.") | startsWith(sources,"T.BT.")), layout="BT";
elseif all(startsWith(sources,"ECHO.") | startsWith(sources,"T.ECHO.")), layout="ECHO";
else, layout="mixed: "+strjoin(unique(sources,'stable')," + ");
end
end

function [x,source]=findByPrefixesAndLeaves(S,prefixes,leaves,requireMatrix)
x=[]; source="";
for p=1:numel(prefixes)
    for l=1:numel(leaves)
        path=[prefixes{p},{leaves{l}}];
        [found,value]=getNestedNumeric(S,path);
        if ~found, continue; end
        if requireMatrix
            usable=ismatrix(value) && size(value,1)>1 && size(value,2)>1;
        else
            usable=isvector(value) && any(isfinite(double(value(:))));
        end
        if usable
            x=double(value);
            source=string(strjoin(path,'.'));
            return
        end
    end
end
end

function [found,value]=getNestedNumeric(S,path)
found=false; value=[]; current=S;
for k=1:numel(path)
    if ~isstruct(current) || isempty(current) || ~isfield(current,path{k}), return; end
    current=current(1).(path{k});
end
if isnumeric(current) || islogical(current)
    value=current;
    found=~isempty(value);
end
end

function n=ensembleLength(x)
n=NaN;
if isempty(x), return; end
x=double(x);
if isvector(x), n=numel(x);
else, n=max(size(x));
end
end

function [x,ok]=alignSeriesExact(vals,idxStart,idxEnd,nRaw,nSeg)
x=[]; ok=false;
if isempty(vals) || ~(isnumeric(vals) || islogical(vals)), return; end
vals=double(vals(:));
if numel(vals)==nSeg, x=vals; ok=true;
elseif numel(vals)==nRaw, x=vals(idxStart:idxEnd); ok=true;
end
end

function [x,ok]=alignSeriesFlexible(vals,idxStart,idxEnd,nRaw,nSeg)
x=[]; ok=false;
if isempty(vals) || ~(isnumeric(vals) || islogical(vals)), return; end
vals=double(vals(:));
if isscalar(vals), x=repmat(vals,nSeg,1); ok=true;
elseif numel(vals)==nSeg, x=vals; ok=true;
elseif numel(vals)>=nRaw, x=vals(idxStart:idxEnd); ok=true;
end
end

function x=forceColumnLength(x,n)
x=double(x(:));
if numel(x)<n, x=[x;nan(n-numel(x),1)];
elseif numel(x)>n, x=x(1:n);
end
end

function [rangeM,unitName,scale,message]= ...
    convertRangeToMetres(rangeRaw,referenceDepth,maxDepth)
rangeRaw=double(rangeRaw(:));
rangeRaw(~isfinite(rangeRaw) | rangeRaw<=0)=NaN;
rangeM=nan(size(rangeRaw)); unitName="unknown"; scale=NaN;
message="No positive range values";
rawValid=rangeRaw(isfinite(rangeRaw));
if isempty(rawValid), return; end
scales=[1,0.1,0.01,0.001]; units=["m","dm","cm","mm"];
scores=inf(size(scales));
for k=1:numel(scales)
    c=rawValid.*scales(k); inside=c>0 & c<=maxDepth;
    if ~any(inside), continue; end
    penalty=10.*(1-mean(inside)); med=median(c(inside));
    if isfinite(referenceDepth) && referenceDepth>0
        scores(k)=abs(log(med./referenceDepth))+penalty;
    else
        scores(k)=penalty+0.01.*(k-1);
    end
end
[bestScore,best]=min(scores);
if ~isfinite(bestScore), return; end
scale=scales(best); unitName=units(best);
rangeM=rangeRaw.*scale;
rangeM(rangeM<=0 | rangeM>maxDepth)=NaN;
message="OK: interpreted as "+unitName;
end

function scale=resolveDataScale(optionValue,data,source)
if isnumeric(optionValue)
    scale=double(optionValue(1));
    if ~isfinite(scale) || scale<=0, error('Data scale must be positive.'); end
    return
end
if lower(strtrim(string(optionValue)))~="auto"
    error('AmplitudeScale/EchoScale must be numeric or auto.');
end
vals=double(data(:)); vals=vals(isfinite(vals) & vals>0);
if contains(source,"BS_adp",'IgnoreCase',true)
    scale=0.1;
elseif isempty(vals)
    scale=1;
elseif median(vals)>100
    scale=0.1;
else
    scale=1;
end
end

function tf=isCalibrationRequested(value)
if ischar(value) || isstring(value)
    textValue=lower(strtrim(string(value)));
    tf=any(contains(textValue,["correct","yes","true","on"]));
else
    tf=logical(value);
end
end

function validateVBGainOptions(opt)
required=[opt.VBGain,opt.VBStrengthLowThreshold, ...
    opt.VBStrengthHighThreshold,opt.VBContrastLowThreshold, ...
    opt.VBContrastHighThreshold];
if ~all(isfinite(required))
    error(['CalibrateVB requires finite VBGain, Strength low/high and ', ...
        'Contrast low/high thresholds.']);
end
end

function value=convertLengthScalar(raw,reference)
value=NaN;
if isempty(raw), return; end
[converted,~,~,~]=convertRangeToMetres(raw,reference,10);
value=safeMean(converted);
end

function value=convertCellSizeScalar(raw)
value=safeMean(raw);
if ~isfinite(value) || value<=0, value=NaN; return; end
if value>0.5
    candidates=value.*[0.1,0.01,0.001];
    candidates=candidates(candidates>0 & candidates<=0.5);
    if isempty(candidates), value=NaN; else, value=max(candidates); end
end
end

function value=resolveM9AcousticScalar(S,kind,fallback)
switch string(kind)
    case "pulse3", leaves={'BT3MHzPL','PulseLength3MHz','BT_PulseLength3MHz'};
    case "pulse05", leaves={'BT05MHzPL','PulseLength05MHz','VB_PulseLength'};
    case "width3", leaves={'width3MHz','TransducerWidth3MHz'};
    case "omega3", leaves={'omega3MHz','omega3MHZ','BeamOpening3MHz'};
    case "width05", leaves={'width05MHz','TransducerWidth05MHz'};
    otherwise, leaves={'omega05MHz','omega05MHZ','BeamOpening05MHz'};
end
prefixes={{'M9meta'},{'T','M9meta'},{'BS_adp'},{'T','BS_adp'},{} };
[x,~]=findByPrefixesAndLeaves(S,prefixes,leaves,false);
value=safeMean(x);
if ~isfinite(value) || value<=0, value=fallback; end
end

function omega=normaliseBeamOpeningRadians(omega)
if ~isfinite(omega) || omega<=0, omega=NaN; return; end
% M9 beam-opening values are small. Values above 0.35 are interpreted as
% degrees so that metadata such as 1.5 is not mistaken for 1.5 radians.
if omega>0.35, omega=deg2rad(omega); end
end

function y=despikeRangeLocal(x,T)
y=double(x(:)); valid=isfinite(y);
if sum(valid)<3, return; end
m=mean(y(valid)); s=std(y(valid));
if isfinite(s) && s>0, y(abs(y-m)>T.*s)=NaN; end
end

function y=despikeStrengthLocal(x,T,MM,Window,rejectMask,color,PLOT,maxValue)
y=double(x(:));
y(y>maxValue)=NaN;
if ~isempty(rejectMask) && numel(rejectMask)==numel(y)
    y(logical(rejectMask(:)))=NaN;
end
if logical(MM), y=movmean(y,round(Window),'omitnan'); end
valid=isfinite(y);
if sum(valid)>=2
    lo=mean(y(valid))-T.*std(y(valid));
    hi=mean(y(valid))+T.*std(y(valid));
    y(y<lo | y>hi)=NaN;
end
if contains(upper(string(PLOT)),"YES")
    figure('Color','w'); plot(y,color); grid on; box on;
end
end

function Af=slantedFootprintDb(pulseLength,R,omega,theta,alphaPrimary,alphaSecondary)
Af=nan(size(R));
if ~all(isfinite([pulseLength,omega])) || pulseLength<=0 || omega<=0, return; end
arg=pulseLength.*R.*omega ./ ...
    (2.*sind(theta-alphaPrimary).*cosd(alphaSecondary));
valid=isfinite(arg) & arg>0;
Af(valid)=10.*log10(arg(valid));
Af(~isreal(Af))=NaN; Af=real(Af);
end

function Af=verticalLegacyFootprintDb(pulseLength,R,omega,betaX,betaY)
Af=nan(size(R));
if ~all(isfinite([pulseLength,omega])) || pulseLength<=0 || omega<=0, return; end
arg=pulseLength.*R.^2.*omega ./ (cosd(betaX).*cosd(betaY));
valid=isfinite(arg) & arg>0;
Af(valid)=10.*log10(arg(valid));
Af(~isreal(Af))=NaN; Af=real(Af);
end

function [Dmin,Dmax,Area,distance,Lp]=SamplingVolume( ...
    distanceToBottom,slantAngle,width,beamHalfAngle,radcp,pulseLength)
Dmin=NaN; Dmax=NaN; Area=NaN; distance=NaN; Lp=NaN;
if ~all(isfinite([distanceToBottom,slantAngle,width,beamHalfAngle,radcp])) || ...
        distanceToBottom<=0 || width<=0, return; end
Dmax=((distanceToBottom+width./2.*sind(slantAngle))./ ...
    tand(90-slantAngle-beamHalfAngle) - ...
    (distanceToBottom-width./2.*sind(slantAngle))./ ...
    tand(90-slantAngle+beamHalfAngle)) + width.*cosd(slantAngle);
Dmin=2.*(distanceToBottom.*tand(beamHalfAngle)./cosd(slantAngle)+width./2);
Area=pi.*Dmax.*Dmin;
distance=distanceToBottom.*tand(slantAngle)+radcp;
Lp=safeMean(double(pulseLength(:)).*sind(90-slantAngle-beamHalfAngle));
end

function E=subsetEchoMatrixM9(Eraw,idxStart,idxEnd,nRaw,nSeg)
E=double(Eraw);
if isempty(E), return; end
if size(E,1)==nRaw || (size(E,1)>=idxEnd && size(E,2)<size(E,1))
    E=E(idxStart:idxEnd,:).';
elseif size(E,2)==nRaw || size(E,2)>=idxEnd
    E=E(:,idxStart:idxEnd);
else
    E=[];
    return
end
if size(E,2)<nSeg, E(:,end+1:nSeg)=NaN;
elseif size(E,2)>nSeg, E=E(:,1:nSeg);
end
end

function id=resolveEchoStopM9(surfaceOffset,WF,WS,rangeMean,nCells)
id=nCells;
if isfinite(WS) && WS>0
    if ~isfinite(surfaceOffset), surfaceOffset=0; end
    if ~isfinite(WF), WF=0; end
    H=surfaceOffset+WF+0.5.*WS+(0:nCells-1)'.*WS+WS./4;
    if isfinite(rangeMean)
        tmp=find(H>rangeMean,1);
        if ~isempty(tmp), id=tmp; end
    end
end
id=max(1,min(id,nCells));
end

function [alphaMean,alphaGrad,dEDR]=processEchoAttenuationM9(E,idStop,WS)
alphaMean=[]; alphaGrad=[]; dEDR=[];
if isempty(E), return; end
alpha=nan(size(E));
last=min(idStop-1,size(E,1)-1);
for j=1:last
    alpha(j+1,:)=E(j,:)-E(j+1,:);
end
alpha(alpha<0 | alpha>15)=NaN;
alphaSmooth=movmean(alpha,[10 10],2,'omitnan');
alphaMean=nanmeanLocal(alphaSmooth,2);
grad=nan(size(E));
if size(E,1)>1, grad(1:end-1,:)=diff(E,1,1); end
alphaGrad=nanmeanLocal(grad(1:max(1,idStop),:),1);
if isfinite(WS) && WS>0, dEDR=alphaMean./WS; end
end

function [AW,abswat]=callFalfaw(temp,depth,frequency,fallback)
AW=fallback; abswat=fallback;
if exist('falfaw','file')==2 && all(isfinite([temp,depth,frequency])) && depth>0
    try
        [a,b]=falfaw(temp,depth,frequency);
        if any(isfinite(a(:))), AW=safeMean(a); end
        if any(isfinite(b(:))), abswat=safeMean(b); end
    catch
    end
end
end

function [mask,source,isExact,isAligned]=deriveVelocityRejectMaskM9( ...
    BTvelOut,i,nFiles,idxStart,idxEnd,nRaw,nSeg,mode)
mask=false(nSeg,1); source="none"; isExact=mode=="none"; isAligned=mode=="none";
if mode=="none", return; end
if isempty(BTvelOut) || ~isstruct(BTvelOut)
    source="BTvelOut_missing"; return
end
if mode=="final"
    maskFields={'BTfilterMask_S','BTfilterMask_1200'};
    validFields={'BTvalidMask_S','BTvalidMask_1200'};
    for k=1:numel(maskFields)
        if ~isfield(BTvelOut,maskFields{k}), continue; end
        vals=getRowContainer(BTvelOut.(maskFields{k}),i,nFiles);
        [candidate,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg);
        if ok
            mask=candidate; source=string(maskFields{k});
            isExact=true; isAligned=true; return
        end
    end
    for k=1:numel(validFields)
        if ~isfield(BTvelOut,validFields{k}), continue; end
        vals=getRowContainer(BTvelOut.(validFields{k}),i,nFiles);
        [candidate,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg);
        if ok
            mask=~candidate; source=string(validFields{k})+"_inverted";
            isExact=true; isAligned=true; return
        end
    end
    seriesFields={'BTmod_final_S','BTmodf_SALL','BTmod_final_1200','BTmodf_1200ALL'};
    for k=1:numel(seriesFields)
        if ~isfield(BTvelOut,seriesFields{k}), continue; end
        vals=getRowContainer(BTvelOut.(seriesFields{k}),i,nFiles);
        [aligned,ok]=alignSeriesExact(vals,idxStart,idxEnd,nRaw,nSeg);
        if ok
            mask=~isfinite(aligned); source=string(seriesFields{k})+"_NaN";
            isExact=true; isAligned=true; return
        end
    end
    source="final_mask_not_found_or_not_aligned";
else
    fields={'peaksBTxy_S','peaksBTxy_1200'};
    for k=1:numel(fields)
        if ~isfield(BTvelOut,fields{k}), continue; end
        vals=getRowContainer(BTvelOut.(fields{k}),i,nFiles);
        [candidate,ok]=alignIndexOrMaskExact(vals,idxStart,idxEnd,nRaw,nSeg);
        if ok
            mask=candidate; source=string(fields{k});
            isExact=true; isAligned=true; return
        end
    end
    source="boat_mask_not_found_or_not_aligned";
end
end

function vals=getRowContainer(C,i,nFiles)
vals=[];
if iscell(C)
    if numel(C)>=i, vals=C{i}; end
elseif isnumeric(C) || islogical(C)
    if isvector(C)
        if nFiles==1 || numel(C)~=nFiles, vals=C(:); end
    elseif size(C,1)==nFiles, vals=C(i,:).';
    elseif size(C,2)==nFiles, vals=C(:,i);
    end
end
end

function [mask,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg)
mask=false(nSeg,1); ok=false;
if isempty(vals), return; end
vals=logical(vals(:));
if numel(vals)==nSeg, mask=vals; ok=true;
elseif numel(vals)==nRaw, mask=vals(idxStart:idxEnd); ok=true;
end
end

function [mask,ok]=alignIndexOrMaskExact(vals,idxStart,idxEnd,nRaw,nSeg)
mask=false(nSeg,1); ok=false;
if isempty(vals), return; end
if islogical(vals)
    [mask,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg); return
end
raw=double(vals(:)); finiteVals=raw(isfinite(raw));
if ~isempty(finiteVals) && all(ismember(finiteVals,[0 1])) && ...
        any(numel(raw)==[nSeg,nRaw])
    [mask,ok]=alignLogicalMaskExact(logical(raw),idxStart,idxEnd,nRaw,nSeg); return
end
idx=round(finiteVals);
if isempty(idx), ok=true; return; end
if all(idx>=1 & idx<=nSeg), mask(idx)=true; ok=true;
elseif all(idx>=idxStart & idx<=idxEnd), mask(idx-idxStart+1)=true; ok=true;
end
end

function text=missingAcousticInputs(SL,pulse,omega,footprintMeans,validAmplitudes)
missing=strings(0,1);
if ~isfinite(SL), missing(end+1)="pressure/source level"; end %#ok<AGROW>
if ~isfinite(pulse) || pulse<=0, missing(end+1)="BT3MHzPL"; end %#ok<AGROW>
if ~isfinite(omega) || omega<=0, missing(end+1)="omega3MHz"; end %#ok<AGROW>
if ~all(isfinite(footprintMeans)), missing(end+1)="valid ranges/footprint"; end %#ok<AGROW>
if ~validAmplitudes, missing(end+1)="valid Strength/Contrast"; end %#ok<AGROW>
if isempty(missing), text="unknown"; else, text=strjoin(missing,"+"); end
end

function y=safeMean(x)
if isempty(x), y=NaN; return; end
x=double(x(:)); x=x(isfinite(x));
if isempty(x), y=NaN; else, y=mean(x); end
end

function y=safeStd(x)
if isempty(x), y=NaN; return; end
x=double(x(:)); x=x(isfinite(x));
if numel(x)<2, y=NaN; else, y=std(x,0); end
end

function y=nanmeanLocal(A,dim)
if nargin<2, dim=1; end
A=double(A); valid=isfinite(A); A(~valid)=0;
count=sum(valid,dim); y=sum(A,dim)./count; y(count==0)=NaN;
end

function y=nanstdLocal(A,flag,dim)
if nargin<2, flag=0; end
if nargin<3, dim=1; end
A=double(A); mu=nanmeanLocal(A,dim);
rep=ones(1,ndims(A)); rep(dim)=size(A,dim);
dev=A-repmat(mu,rep); valid=isfinite(A); dev(~valid)=0;
n=sum(valid,dim);
if flag==0, den=max(n-1,0); else, den=n; end
y=sqrt(sum(dev.^2,dim)./den); y(den<=0)=NaN;
end

function y=nansumLocal(x)
x=double(x(:)); x=x(isfinite(x));
if isempty(x), y=0; else, y=sum(x); end
end

function assignStructToWorkspace(S,workspaceName)
f=fieldnames(S);
for k=1:numel(f), assignin(workspaceName,f{k},S.(f{k})); end
end
