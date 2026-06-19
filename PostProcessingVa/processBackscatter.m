function BackscatterOut = processBackscatter( ...
    dataADCP1200, ADCPin, DepthRNGout, BTvelOut, varargin)
%PROCESSBACKSCATTER Process 1200-kHz ADCP bottom backscatter.
%
% BackscatterOut = processBackscatter1200( ...
%     dataADCP1200, ADCPin, DepthRNGout, BTvelOut, 'Name',Value,...)
%
% The function converts the legacy row-by-row backscatter script into the
% same input-table-driven workflow used by the other 1200-kHz functions.
%
% REQUIRED INPUTS
% ---------------
% dataADCP1200 : dir-style structure returned by dir(), or a folder path.
% ADCPin       : input-table structure. Required field ADCPfilename.
%                StartADCP and EndADCP are used for every raw time series.
% DepthRNGout  : output from processDepthsRNG. Existing processed RNG,
%                WF, WS and temperature values are reused where available.
% BTvelOut     : output from processBTbedloadVelocity1200. BTmodf_1200,
%                FilteredDataN_1200 and BTstd_1200 are reused for Method 2.
%                The rejected/NaN velocity ensembles are also assigned NaN
%                in RNG, RSSI, power, attitude and echo time series here.
%
% IMPORTANT START/END AND FILTER ALIGNMENT
% ----------------------------------------
% By default StartADCP and EndADCP are interpreted as numbers of ensembles
% removed from the beginning and end of every raw series, respectively.
% Thus StartADCP=30 and EndADCP=30 uses ensembles 31:(N-30). This matches
% the trimming convention used by the BT-velocity and depth-processing
% functions. Set StartEndMode='absolute' only when the input values are
% literal first/last ensemble indices. In VelocityMaskMode='final', the
% function uses the explicit row-wise BTfilterMask_1200 exported by the
% updated BT-velocity function. The mask is checked against the same trimmed
% start/end interval before it is applied.
%
% OUTPUT
% ------
% BackscatterOut is the only output structure. All original variable names
% are fields directly inside BackscatterOut; there are no nested BTamp or Q
% structures. When AssignToCaller=true, the same variables are also assigned
% directly to the caller workspace using their unchanged names.
%
% MAIN OPTIONS
% ------------
% Acoustic/instrument:
%   theta1                    slant angle [deg], default 20
%   kc                        counts-to-dB coefficient, default 0.60
%   Er                        clear-water noise [counts], default 40
%   C1200                     instrument constant [dB], default -129
%   Pf, Ef                    power and efficiency factors, default 0.5
%   c                         sound speed [m/s], default 1450
%   f12                       frequency supplied to falfaw, default 1200
%   alphaw1200                fallback water absorption [dB/m], default 0.392
%   alphaSed12                sediment attenuation [dB/m], scalar/vector
%   POWER1200                 transmit-power constant [dB], default 7.8
%   percentagePulse           pulse length / bed range, default 0.20
%   width12MHz                transducer width [m], default 0.054
%   omega12MHZ                full beam opening [rad], default deg2rad(3)
%   DistTransRDI              ADCP centre-to-transducer radius [m], 0.10
%   distanceADCPtoWaterSurface transducer offset used for echo cell range,
%                              default 0.07 m
%   DtransON1m                insonified diameter at 1 m [m]. If empty,
%                              calculated by SamplingVolume.
%
% Trimming:
%   StartEndMode              'trimcounts' (default) or 'absolute'
%
% RSSI filtering:
%   Tdespike                  global-standard-deviation multiplier, default 3
%   MM                        apply moving mean before despiking, default 1
%   Win                       moving-mean window, default 3
%   PLOT                      'YES' or 'NO', default 'NO'
%   UseVelocityFilterMask     backward-compatible switch, default true
%   VelocityMaskMode          'final' (default), 'boat', or 'none'
%       final : use the explicit per-file BTfilterMask_1200 exported by
%               processBTbedloadVelocity1200. This is the recommended mode
%               when BT velocity and backscatter must use identical ensembles.
%       boat  : use only peaksBTxy_1200 from the boat-motion threshold.
%       none  : do not transfer any BT-velocity mask to backscatter.
%   RequireExactVelocityMask  default true. In final mode, stop with an
%                              informative error if the per-file mask or
%                              trim alignment is unavailable/mismatched.
%   RSSIMax                   values above this count become NaN, default 1500
%
% Method 2:
%   ApplyMethod2Regression    default true. If the external helper
%                             TrateMethod2Mregression1200 is unavailable,
%                             qsMethod2regression_1200 remains NaN.
%
% Workspace:
%   AssignToCaller            default true
%   AssignToBase              default false
%   Verbose                   default true
%
% EXAMPLE
% -------
% BackscatterOut = processBackscatter1200( ...
%     dataADCP1200, PhSamplesin, DepthRNGout, BTvelOut, ...
%     'kc',0.60, ...
%     'Er',40, ...
%     'alphaSed12',0, ...
%     'UseVelocityFilterMask',true, ...
%     'VelocityMaskMode','final', ...
%     'RequireExactVelocityMask',true, ...
%     'AssignToCaller',true, ...
%     'Verbose',true);
%
% No files are saved.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'dataADCP1200',@(x) isstruct(x) || ischar(x) || isstring(x));
addRequired(p,'ADCPin',@(x) isstruct(x));
addRequired(p,'DepthRNGout',@(x) isempty(x) || isstruct(x));
addRequired(p,'BTvelOut',@(x) isempty(x) || isstruct(x));

addParameter(p,'theta1',20,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<90);
addParameter(p,'kc',0.60,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'Er',40,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'C1200',-129,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'Pf',0.5,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'Ef',0.5,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'c',1450,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'f12',1200,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'alphaw1200',0.392,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'alphaSed12',0,@(x) isnumeric(x) && ~isempty(x));
addParameter(p,'POWER1200',7.8,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'percentagePulse',0.20,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'width12MHz',0.054,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'omega12MHZ',deg2rad(3),@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'DistTransRDI',0.10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'distanceADCPtoWaterSurface',0.07,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'DtransON1m',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x>0));
addParameter(p,'StartEndMode','trimcounts',@(x) ischar(x) || isstring(x));

addParameter(p,'Tdespike',3,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'MM',1,@(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'Win',3,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=1);
addParameter(p,'PLOT','NO',@(x) ischar(x) || isstring(x));
addParameter(p,'UseVelocityFilterMask',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'VelocityMaskMode','final',@(x) ischar(x) || isstring(x));
addParameter(p,'RequireExactVelocityMask',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'RSSIMax',1500,@(x) isnumeric(x) && isscalar(x) && isfinite(x));

addParameter(p,'nodata1200',[],@(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(p,'ApplyMethod2Regression',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToCaller',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToBase',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'Verbose',true,@(x) islogical(x) || isnumeric(x));

parse(p,dataADCP1200,ADCPin,DepthRNGout,BTvelOut,varargin{:});
opt = p.Results;
% Keep only configurable options in the exported Options structure.
requiredNames={'dataADCP1200','ADCPin','DepthRNGout','BTvelOut'};
for iReq=1:numel(requiredNames)
    if isfield(opt,requiredNames{iReq}), opt=rmfield(opt,requiredNames{iReq}); end
end
startEndMode = lower(strtrim(string(opt.StartEndMode)));
if ~any(startEndMode == ["trimcounts","trim","absolute","indices"])
    error('StartEndMode must be trimcounts or absolute.');
end
velocityMaskMode = lower(strtrim(string(opt.VelocityMaskMode)));
if ~logical(opt.UseVelocityFilterMask)
    velocityMaskMode = "none";
end
if ~any(velocityMaskMode == ["final","boat","none"])
    error('VelocityMaskMode must be final, boat, or none.');
end

%% Input rows and file list
if ~isfield(ADCPin,'ADCPfilename')
    error('ADCPin.ADCPfilename is missing.');
end

ADCPfilename = string(ADCPin.ADCPfilename(:));
nFiles = numel(ADCPfilename);
InputRow = (1:nFiles)';
StartADCP = getVectorFieldAny(ADCPin,{'StartADCP','Start','start'},ones(nFiles,1),nFiles);
EndADCP = getVectorFieldAny(ADCPin,{'EndADCP','End','end'},nan(nFiles,1),nFiles);
DistanceFromLeftBank = getVectorFieldAny(ADCPin, ...
    {'DistanceFromLeftBank','DistanceToBank'},nan(nFiles,1),nFiles);
D50_14 = getVectorFieldAny(ADCPin,{'D50_14','D50_14mm','D50mobile14'},nan(nFiles,1),nFiles);
D50_05 = getVectorFieldAny(ADCPin,{'D50_05','D50_05mm','D50mobile05'},nan(nFiles,1),nFiles);
D50 = getVectorFieldAny(ADCPin,{'D50','D_50'},nan(nFiles,1),nFiles);
D_50grabM = getVectorFieldAny(ADCPin,{'D_50grabM','D50grabM','D50_grabM','D50grab'},nan(nFiles,1),nFiles);

fileList = normaliseFileList(dataADCP1200);

alphaSed12 = prepareVector(opt.alphaSed12,nFiles);

if isempty(opt.DtransON1m)
    [Dmin1m,~,~,~,~] = SamplingVolume(1,opt.theta1,opt.width12MHz, ...
        rad2deg(opt.omega12MHZ)./2,opt.DistTransRDI,opt.percentagePulse);
    DtransON1m = Dmin1m;
else
    DtransON1m = opt.DtransON1m;
end

%% Reuse BT-velocity summary outputs
BTmodf_1200 = getVectorFieldAny(BTvelOut, ...
    {'BTmodf_1200','va_final_1200'},nan(nFiles,1),nFiles);
FilteredDataN_1200 = getVectorFieldAny(BTvelOut, ...
    {'FilteredDataN_1200','FilteredDataN1200'},nan(nFiles,1),nFiles);
BTstd_1200 = getVectorFieldAny(BTvelOut, ...
    {'BTstd_1200','BTstd1200'},nan(nFiles,1),nFiles);

%% Preallocate row-wise outputs
FileFound = false(nFiles,1);
MatchedFile = strings(nFiles,1);
Message = strings(nFiles,1);
DummyValue = nan(nFiles,1);
StartUsed_1200 = nan(nFiles,1);
EndUsed_1200 = nan(nFiles,1);
NensemblesRaw_1200 = nan(nFiles,1);
NensemblesUsed_1200 = nan(nFiles,1);
NfilteredBTvelocity_1200 = zeros(nFiles,1);
BTvelocityFilterMask_1200 = cell(nFiles,1);
BTvelocityFilterSource_1200 = strings(nFiles,1);
BTvelocityFilterExact_1200 = false(nFiles,1);
BTvelocityFilterAligned_1200 = false(nFiles,1);
BTvelocityMaskMode_1200 = repmat(velocityMaskMode,nFiles,1);

WF_12 = nan(nFiles,1);
WS_12 = nan(nFiles,1);
WT_12 = nan(nFiles,1);
Temp1 = nan(nFiles,1);

Depth1200_cell = cell(nFiles,1);
Depth1200std_cell = cell(nFiles,1);
stdSpatial_1200 = nan(nFiles,1);
stdSpatialALL_1200 = cell(nFiles,1);

RNG11_1200 = cell(nFiles,1);
RNG22_1200 = cell(nFiles,1);
RNG33_1200 = cell(nFiles,1);
RNG44_1200 = cell(nFiles,1);

BT_RSSI1_1200 = cell(nFiles,1);
BT_RSSI2_1200 = cell(nFiles,1);
BT_RSSI3_1200 = cell(nFiles,1);
BT_RSSI4_1200 = cell(nFiles,1);
BT_RSSI1200_cell = cell(nFiles,1);
stdSpatialRSSI_1200 = nan(nFiles,1);

psi1_1200 = nan(nFiles,1);
psi2_1200 = nan(nFiles,1);
psi3_1200 = nan(nFiles,1);
psi4_1200 = nan(nFiles,1);
alfaXmean_1200 = nan(nFiles,1);
alfaYmean_1200 = nan(nFiles,1);

PtransW1200 = cell(nFiles,1);
PtransWmean1200 = nan(nFiles,1);
Pdwb1200 = cell(nFiles,1);
STDPowerdb1200 = nan(nFiles,1);
STDPW1200 = nan(nFiles,1);
XmitCC1200 = cell(nFiles,1);
XmitVV1200 = cell(nFiles,1);
PitchH1200 = cell(nFiles,1);
RollL1200 = cell(nFiles,1);

R1_1200 = nan(nFiles,1);
R2_1200 = nan(nFiles,1);
R3_1200 = nan(nFiles,1);
R4_1200 = nan(nFiles,1);
Spread1_1200 = nan(nFiles,1);
Spread2_1200 = nan(nFiles,1);
Spread3_1200 = nan(nFiles,1);
Spread4_1200 = nan(nFiles,1);
SpreadM1200 = nan(nFiles,1);
LdbMean1200 = nan(nFiles,1);

Dmin_1200_1 = nan(nFiles,1); Dmin_1200_2 = nan(nFiles,1);
Dmin_1200_3 = nan(nFiles,1); Dmin_1200_4 = nan(nFiles,1);
Dmax_1200_1 = nan(nFiles,1); Dmax_1200_2 = nan(nFiles,1);
Dmax_1200_3 = nan(nFiles,1); Dmax_1200_4 = nan(nFiles,1);
Area_1200_1 = nan(nFiles,1); Area_1200_2 = nan(nFiles,1);
Area_1200_3 = nan(nFiles,1); Area_1200_4 = nan(nFiles,1);
distance_1200_1 = nan(nFiles,1); distance_1200_2 = nan(nFiles,1);
distance_1200_3 = nan(nFiles,1); distance_1200_4 = nan(nFiles,1);
Lp_1200_1 = nan(nFiles,1); Lp_1200_2 = nan(nFiles,1);
Lp_1200_3 = nan(nFiles,1); Lp_1200_4 = nan(nFiles,1);

Af1mean1200 = nan(nFiles,1);
Af2mean1200 = nan(nFiles,1);
Af3mean1200 = nan(nFiles,1);
Af4mean1200 = nan(nFiles,1);

SL1200 = cell(nFiles,1);

alpha1mean = cell(nFiles,1);
alpha2mean = cell(nFiles,1);
alpha3mean = cell(nFiles,1);
alpha4mean = cell(nFiles,1);
alpha1mean1200grad = cell(nFiles,1);
alpha2mean1200grad = cell(nFiles,1);
alpha3mean1200grad = cell(nFiles,1);
alpha4mean1200grad = cell(nFiles,1);
dEDR1_1200 = cell(nFiles,1);
dEDR2_1200 = cell(nFiles,1);
dEDR3_1200 = cell(nFiles,1);
dEDR4_1200 = cell(nFiles,1);
alphaMean1200 = nan(nFiles,4);

AW_12 = nan(nFiles,1);
abswat12_1 = nan(nFiles,1);
abswat12_2 = nan(nFiles,1);
abswat12_3 = nan(nFiles,1);
abswat12_4 = nan(nFiles,1);

TL1_1200 = nan(nFiles,1); TL2_1200 = nan(nFiles,1);
TL3_1200 = nan(nFiles,1); TL4_1200 = nan(nFiles,1);
TL1_fix1200 = nan(nFiles,1); TL2_fix1200 = nan(nFiles,1);
TL3_fix1200 = nan(nFiles,1); TL4_fix1200 = nan(nFiles,1);

SV1_1200 = cell(nFiles,1); SV2_1200 = cell(nFiles,1);
SV3_1200 = cell(nFiles,1); SV4_1200 = cell(nFiles,1);
SV11_1200 = cell(nFiles,1); SV22_1200 = cell(nFiles,1);
SV33_1200 = cell(nFiles,1); SV44_1200 = cell(nFiles,1);

SVBmean1_1200 = nan(nFiles,1); SVBmean2_1200 = nan(nFiles,1);
SVBmean3_1200 = nan(nFiles,1); SVBmean4_1200 = nan(nFiles,1);
SVB1mean1_1200 = nan(nFiles,1); SVB1mean2_1200 = nan(nFiles,1);
SVB1mean3_1200 = nan(nFiles,1); SVB1mean4_1200 = nan(nFiles,1);
SVBstd1_1200 = nan(nFiles,1); SVBstd2_1200 = nan(nFiles,1);
SVBstd3_1200 = nan(nFiles,1); SVBstd4_1200 = nan(nFiles,1);
SVB1std1_1200 = nan(nFiles,1); SVB1std2_1200 = nan(nFiles,1);
SVB1std3_1200 = nan(nFiles,1); SVB1std4_1200 = nan(nFiles,1);
stdSvSpatial_1200 = nan(nFiles,1);
stdSvSpatialALL_1200 = cell(nFiles,1);
S_1200 = nan(nFiles,1);
D50_12 = nan(nFiles,1);
qsMethod2regression_1200 = nan(nFiles,1);
Method2Status_1200 = strings(nFiles,1);
Method2Status_1200(:) = "not_requested";

%% Main row loop
for i = 1:nFiles
    try
        if isNoDataRow(opt.nodata1200,i,nFiles)
            DummyValue(i) = NaN;
            Message(i) = "Skipped: nodata1200";
            continue
        end
        DummyValue(i) = 1;

        matchedFile = matchADCPfile(fileList,ADCPfilename(i));
        MatchedFile(i) = matchedFile;
        FileFound(i) = matchedFile ~= "";
        if ~FileFound(i)
            error('No MAT file matched ADCPfilename %s.',ADCPfilename(i));
        end

        Sload = load(matchedFile);
        requiredTop = {'BT','Data','ECHO'};
        for k = 1:numel(requiredTop)
            if ~isfield(Sload,requiredTop{k})
                error('Required structure %s is missing in %s.',requiredTop{k},matchedFile);
            end
        end
        BT = Sload.BT;
        Data = Sload.Data;
        ECHO = Sload.ECHO;

        requiredBT = {'RNG1','RNG2','RNG3','RNG4', ...
            'BT_RSSI1','BT_RSSI2','BT_RSSI3','BT_RSSI4'};
        for k = 1:numel(requiredBT)
            if ~isfield(BT,requiredBT{k})
                error('BT.%s is missing in %s.',requiredBT{k},matchedFile);
            end
        end
        requiredData = {'XmitC','XmitV','Pitch','Roll'};
        for k = 1:numel(requiredData)
            if ~isfield(Data,requiredData{k})
                error('Data.%s is missing in %s.',requiredData{k},matchedFile);
            end
        end
        requiredEcho = {'ECH1','ECH2','ECH3','ECH4'};
        for k = 1:numel(requiredEcho)
            if ~isfield(ECHO,requiredEcho{k})
                error('ECHO.%s is missing in %s.',requiredEcho{k},matchedFile);
            end
        end

        rawLengths = [numel(BT.RNG1),numel(BT.RNG2),numel(BT.RNG3),numel(BT.RNG4), ...
            numel(BT.BT_RSSI1),numel(BT.BT_RSSI2),numel(BT.BT_RSSI3),numel(BT.BT_RSSI4), ...
            numel(Data.XmitC),numel(Data.XmitV),numel(Data.Pitch),numel(Data.Roll)];
        nRaw = min(rawLengths(rawLengths>0));
        if isempty(nRaw) || nRaw < 1
            error('No usable ensemble data in %s.',matchedFile);
        end
        NensemblesRaw_1200(i) = nRaw;

        [idxStart,idxEnd] = resolveStartEnd(StartADCP(i),EndADCP(i),nRaw,startEndMode);
        StartUsed_1200(i) = idxStart;
        EndUsed_1200(i) = idxEnd;
        idx = idxStart:idxEnd;
        nSeg = numel(idx);
        NensemblesUsed_1200(i) = nSeg;

        [velocityRejectMask,maskSource,maskExact,maskAligned] = ...
            deriveVelocityRejectMask( ...
                BTvelOut,i,nFiles,idxStart,idxEnd,nRaw,nSeg,velocityMaskMode);

        velocityRejectMask = logical(velocityRejectMask(:));
        if numel(velocityRejectMask) ~= nSeg
            maskAligned = false;
            if logical(opt.RequireExactVelocityMask) && velocityMaskMode == "final"
                error(['BT velocity mask length mismatch for row %d. Expected %d ', ...
                    'trimmed ensembles but received %d. Rerun the updated ', ...
                    'processBTbedloadVelocity1200 first.'], ...
                    i,nSeg,numel(velocityRejectMask));
            end
            velocityRejectMask = false(nSeg,1);
            maskSource = maskSource + "_length_mismatch_disabled";
        end

        if logical(opt.RequireExactVelocityMask) && velocityMaskMode == "final" && ...
                ~(maskExact && maskAligned)
            error(['No exact aligned final BT mask was available for row %d. ', ...
                'Rerun the updated processBTbedloadVelocity1200 and pass its ', ...
                'output as BTvelOut, or use VelocityMaskMode=''none''. ', ...
                'Detected source: %s.'],i,maskSource);
        end

        BTvelocityFilterMask_1200{i} = velocityRejectMask;
        BTvelocityFilterSource_1200(i) = maskSource;
        BTvelocityFilterExact_1200(i) = maskExact;
        BTvelocityFilterAligned_1200(i) = maskAligned;
        NfilteredBTvelocity_1200(i) = sum(velocityRejectMask);

        % Prefer already processed depth/range series where they align;
        % otherwise use the raw fields from the matched MAT file.
        RNG11 = getAlignedSeries(DepthRNGout, ...
            {'RNG11_','RNG11','RNG11_1200','RNG11_1200_all'},i,nFiles,idxStart,idxEnd,nRaw);
        RNG22 = getAlignedSeries(DepthRNGout, ...
            {'RNG22_','RNG22','RNG22_1200','RNG22_1200_all'},i,nFiles,idxStart,idxEnd,nRaw);
        RNG33 = getAlignedSeries(DepthRNGout, ...
            {'RNG33_','RNG33','RNG33_1200','RNG33_1200_all'},i,nFiles,idxStart,idxEnd,nRaw);
        RNG44 = getAlignedSeries(DepthRNGout, ...
            {'RNG44_','RNG44','RNG44_1200','RNG44_1200_all'},i,nFiles,idxStart,idxEnd,nRaw);
        if isempty(RNG11), RNG11 = double(BT.RNG1(idx)); end
        if isempty(RNG22), RNG22 = double(BT.RNG2(idx)); end
        if isempty(RNG33), RNG33 = double(BT.RNG3(idx)); end
        if isempty(RNG44), RNG44 = double(BT.RNG4(idx)); end
        RNG11 = forceColumnLength(RNG11,nSeg);
        RNG22 = forceColumnLength(RNG22,nSeg);
        RNG33 = forceColumnLength(RNG33,nSeg);
        RNG44 = forceColumnLength(RNG44,nSeg);
        RNG11(RNG11==0) = NaN; RNG22(RNG22==0) = NaN;
        RNG33(RNG33==0) = NaN; RNG44(RNG44==0) = NaN;
        RNG11(velocityRejectMask) = NaN; RNG22(velocityRejectMask) = NaN;
        RNG33(velocityRejectMask) = NaN; RNG44(velocityRejectMask) = NaN;

        RNG11_1200{i} = despikeRangeLocal(RNG11,2);
        RNG22_1200{i} = despikeRangeLocal(RNG22,2);
        RNG33_1200{i} = despikeRangeLocal(RNG33,2);
        RNG44_1200{i} = despikeRangeLocal(RNG44,2);

        rangeMatrix = [RNG11_1200{i},RNG22_1200{i},RNG33_1200{i},RNG44_1200{i}];
        Depth1200_cell{i} = nanmeanLocal(rangeMatrix,2);
        Depth1200std_cell{i} = nanstdLocal(rangeMatrix,0,2);
        stdSpatialALL_1200{i} = nanstdLocal(rangeMatrix,0,2);
        stdSpatial_1200(i) = safeMean(stdSpatialALL_1200{i});

        WF_12(i) = getScalarPreferred(DepthRNGout,{'WF_12','WF_1200'},i,nFiles,NaN);
        if ~isfinite(WF_12(i)) && isfield(Data,'WF')
            WF_12(i) = safeMean(Data.WF)./100;
        end
        WS_12(i) = getScalarPreferred(DepthRNGout,{'WS_12','WS_1200'},i,nFiles,NaN);
        if ~isfinite(WS_12(i)) && isfield(Data,'WS')
            WS_12(i) = safeMean(Data.WS)./100;
        end
        WT_12(i) = getScalarPreferred(DepthRNGout,{'WT_12','WT12','Temp1'},i,nFiles,NaN);
        if ~isfinite(WT_12(i)) && isfield(Data,'WT')
            WT_12(i) = safeMean(Data.WT);
        end
        Temp1(i) = WT_12(i);

        % Apply the exact BT-velocity rejection mask before RSSI despiking.
        rssi1 = forceColumnLength(double(BT.BT_RSSI1(idx)),nSeg);
        rssi2 = forceColumnLength(double(BT.BT_RSSI2(idx)),nSeg);
        rssi3 = forceColumnLength(double(BT.BT_RSSI3(idx)),nSeg);
        rssi4 = forceColumnLength(double(BT.BT_RSSI4(idx)),nSeg);
        rssi1(velocityRejectMask)=NaN; rssi2(velocityRejectMask)=NaN;
        rssi3(velocityRejectMask)=NaN; rssi4(velocityRejectMask)=NaN;

        BT_RSSI1_1200{i} = despikeBTstrength(rssi1,opt.Tdespike,opt.MM, ...
            opt.Win,velocityRejectMask,'r',opt.PLOT,opt.RSSIMax);
        BT_RSSI2_1200{i} = despikeBTstrength(rssi2,opt.Tdespike,opt.MM, ...
            opt.Win,velocityRejectMask,'g',opt.PLOT,opt.RSSIMax);
        BT_RSSI3_1200{i} = despikeBTstrength(rssi3,opt.Tdespike,opt.MM, ...
            opt.Win,velocityRejectMask,'b',opt.PLOT,opt.RSSIMax);
        BT_RSSI4_1200{i} = despikeBTstrength(rssi4,opt.Tdespike,opt.MM, ...
            opt.Win,velocityRejectMask,'m',opt.PLOT,opt.RSSIMax);
        rssiMatrix = [BT_RSSI1_1200{i},BT_RSSI2_1200{i}, ...
            BT_RSSI3_1200{i},BT_RSSI4_1200{i}];
        BT_RSSI1200_cell{i} = nanmeanLocal(rssiMatrix,2);
        stdSpatialRSSI_1200(i) = safeMean(nanstdLocal(rssiMatrix,0,2));

        Pitch = forceColumnLength(double(Data.Pitch(idx)),nSeg);
        Roll = forceColumnLength(double(Data.Roll(idx)),nSeg);
        Pitch(velocityRejectMask)=NaN;
        Roll(velocityRejectMask)=NaN;
        PitchH1200{i}=Pitch;
        RollL1200{i}=Roll;

        % Local bed-angle correction.
        D1=RNG11_1200{i}.*tand(opt.theta1);
        D2=RNG22_1200{i}.*tand(opt.theta1);
        D3=RNG33_1200{i}.*tand(opt.theta1);
        D4=RNG44_1200{i}.*tand(opt.theta1);
        psi11=90-opt.theta1-atand((RNG11_1200{i}-RNG22_1200{i})./(D1+D2));
        psi22=90-opt.theta1-atand((RNG22_1200{i}-RNG11_1200{i})./(D2+D1));
        psi33=90-opt.theta1-atand((RNG33_1200{i}-RNG44_1200{i})./(D4+D3));
        psi44=90-opt.theta1-atand((RNG44_1200{i}-RNG33_1200{i})./(D3+D4));
        psi1_1200(i)=safeMean(psi11+Pitch);
        psi2_1200(i)=safeMean(psi22-Pitch);
        psi3_1200(i)=safeMean(psi33+Roll);
        psi4_1200(i)=safeMean(psi44-Roll);

        XmitCfull=double(Data.XmitC).*0.011451;
        XmitVfull=double(Data.XmitV).*0.253765;
        XmitC=forceColumnLength(XmitCfull(idx),nSeg);
        XmitV=forceColumnLength(XmitVfull(idx),nSeg);
        XmitC(velocityRejectMask)=NaN;
        XmitV(velocityRejectMask)=NaN;
        XmitCC1200{i}=XmitC;
        XmitVV1200{i}=XmitV;

        PtransW1200{i}=10.*log10(XmitC.*XmitV);
        PtransWmean1200(i)=safeMean(PtransW1200{i});
        Pdwb1200{i}=opt.POWER1200 + 20.*log10(XmitV./12);
        STDPowerdb1200(i)=safeStd(Pdwb1200{i});
        STDPW1200(i)=safeStd(PtransW1200{i});

        R11=RNG11_1200{i}./cosd(opt.theta1);
        R22=RNG22_1200{i}./cosd(opt.theta1);
        R33=RNG33_1200{i}./cosd(opt.theta1);
        R44=RNG44_1200{i}./cosd(opt.theta1);
        R1_1200(i)=safeMean(R11); R2_1200(i)=safeMean(R22);
        R3_1200(i)=safeMean(R33); R4_1200(i)=safeMean(R44);

        Spread11=20.*log10(R11); Spread22=20.*log10(R22);
        Spread33=20.*log10(R33); Spread44=20.*log10(R44);
        Spread1_1200(i)=safeMean(Spread11); Spread2_1200(i)=safeMean(Spread22);
        Spread3_1200(i)=safeMean(Spread33); Spread4_1200(i)=safeMean(Spread44);
        SpreadM1200(i)=safeMean([Spread11;Spread22;Spread33;Spread44]);

        Ldbm1=opt.percentagePulse.*RNG11_1200{i};
        Ldbm2=opt.percentagePulse.*RNG22_1200{i};
        Ldbm3=opt.percentagePulse.*RNG33_1200{i};
        Ldbm4=opt.percentagePulse.*RNG44_1200{i};
        LdbMean1200(i)=safeMean([Ldbm1;Ldbm2;Ldbm3;Ldbm4]);

        [Dmin_1200_1(i),Dmax_1200_1(i),Area_1200_1(i), ...
            distance_1200_1(i),Lp_1200_1(i)] = SamplingVolume( ...
            safeMean(RNG11_1200{i}),opt.theta1,opt.width12MHz, ...
            rad2deg(opt.omega12MHZ)./2,opt.DistTransRDI,Ldbm1);
        [Dmin_1200_2(i),Dmax_1200_2(i),Area_1200_2(i), ...
            distance_1200_2(i),Lp_1200_2(i)] = SamplingVolume( ...
            safeMean(RNG22_1200{i}),opt.theta1,opt.width12MHz, ...
            rad2deg(opt.omega12MHZ)./2,opt.DistTransRDI,Ldbm2);
        [Dmin_1200_3(i),Dmax_1200_3(i),Area_1200_3(i), ...
            distance_1200_3(i),Lp_1200_3(i)] = SamplingVolume( ...
            safeMean(RNG33_1200{i}),opt.theta1,opt.width12MHz, ...
            rad2deg(opt.omega12MHZ)./2,opt.DistTransRDI,Ldbm3);
        [Dmin_1200_4(i),Dmax_1200_4(i),Area_1200_4(i), ...
            distance_1200_4(i),Lp_1200_4(i)] = SamplingVolume( ...
            safeMean(RNG44_1200{i}),opt.theta1,opt.width12MHz, ...
            rad2deg(opt.omega12MHZ)./2,opt.DistTransRDI,Ldbm4);

        alfaX=nanmeanLocal([abs(90-opt.theta1-psi11), ...
            abs(90-opt.theta1-psi22)],2);
        alfaY=nanmeanLocal([abs(90-opt.theta1-psi33), ...
            abs(90-opt.theta1-psi44)],2);
        alfaXmean_1200(i)=safeMean(alfaX);
        alfaYmean_1200(i)=safeMean(alfaY);

        Af1=surfaceScatteringArea(Ldbm1,R11,opt.omega12MHZ, ...
            opt.theta1,alfaX,alfaY);
        Af2=surfaceScatteringArea(Ldbm2,R22,opt.omega12MHZ, ...
            opt.theta1,alfaY,alfaX);
        Af3=surfaceScatteringArea(Ldbm3,R33,opt.omega12MHZ, ...
            opt.theta1,alfaX,alfaY);
        Af4=surfaceScatteringArea(Ldbm4,R44,opt.omega12MHZ, ...
            opt.theta1,alfaY,alfaX);
        Af1mean1200(i)=safeMean(Af1); Af2mean1200(i)=safeMean(Af2);
        Af3mean1200(i)=safeMean(Af3); Af4mean1200(i)=safeMean(Af4);

        II=4.*XmitC.*XmitV.*opt.Pf.*opt.Ef./(pi.*DtransON1m.^2);
        PP1200=sqrt(1000.*opt.c.*II);
        % Preserve the original natural-log source-level formulation.
        SL1200{i}=10.*log(PP1200);

        % Subset echo matrices to the same StartADCP:EndADCP interval and
        % assign NaN to the same velocity-rejected ensemble columns.
        E1=subsetEchoMatrix(ECHO.ECH1,idxStart,idxEnd,nRaw,nSeg);
        E2=subsetEchoMatrix(ECHO.ECH2,idxStart,idxEnd,nRaw,nSeg);
        E3=subsetEchoMatrix(ECHO.ECH3,idxStart,idxEnd,nRaw,nSeg);
        E4=subsetEchoMatrix(ECHO.ECH4,idxStart,idxEnd,nRaw,nSeg);
        E1(:,velocityRejectMask)=NaN; E2(:,velocityRejectMask)=NaN;
        E3(:,velocityRejectMask)=NaN; E4(:,velocityRejectMask)=NaN;

        nCells=max([size(E1,1),size(E2,1),size(E3,1),size(E4,1)]);
        H=opt.distanceADCPtoWaterSurface + WF_12(i) + ...
            0.5.*(WS_12(i)+(0:nCells-1)'.*WS_12(i)) + WS_12(i)./2;
        IDstop1=resolveEchoStop(H,safeMean(RNG11_1200{i}),size(E1,1));
        IDstop2=resolveEchoStop(H,safeMean(RNG22_1200{i}),size(E2,1));
        IDstop3=resolveEchoStop(H,safeMean(RNG33_1200{i}),size(E3,1));
        IDstop4=resolveEchoStop(H,safeMean(RNG44_1200{i}),size(E4,1));

        algrad1=gradient(E1); algrad2=gradient(E2);
        algrad3=gradient(E3); algrad4=gradient(E4);
        alpha1=echoAttenuation(E1,IDstop1,opt.kc);
        alpha2=echoAttenuation(E2,IDstop2,opt.kc);
        alpha3=echoAttenuation(E3,IDstop3,opt.kc);
        alpha4=echoAttenuation(E4,IDstop4,opt.kc);
        alpha11=movmean(alpha1,[10 10],2,'omitnan');
        alpha22=movmean(alpha2,[10 10],2,'omitnan');
        alpha33=movmean(alpha3,[10 10],2,'omitnan');
        alpha44=movmean(alpha4,[10 10],2,'omitnan');

        if isfinite(WS_12(i)) && WS_12(i)>0
            dEDR1_1200{i}=nanmeanLocal(alpha11./WS_12(i),2);
            dEDR2_1200{i}=nanmeanLocal(alpha22./WS_12(i),2);
            dEDR3_1200{i}=nanmeanLocal(alpha33./WS_12(i),2);
            dEDR4_1200{i}=nanmeanLocal(alpha44./WS_12(i),2);
        end
        alpha1mean{i}=nanmeanLocal(alpha11,2);
        alpha2mean{i}=nanmeanLocal(alpha22,2);
        alpha3mean{i}=nanmeanLocal(alpha33,2);
        alpha4mean{i}=nanmeanLocal(alpha44,2);
        alpha1mean1200grad{i}=nanmeanLocal(algrad1(1:IDstop1,:),1);
        alpha2mean1200grad{i}=nanmeanLocal(algrad2(1:IDstop2,:),1);
        alpha3mean1200grad{i}=nanmeanLocal(algrad3(1:IDstop3,:),1);
        alpha4mean1200grad{i}=nanmeanLocal(algrad4(1:IDstop4,:),1);
        alphaMean1200(i,:)=[safeMean(alpha1mean{i}),safeMean(alpha2mean{i}), ...
            safeMean(alpha3mean{i}),safeMean(alpha4mean{i})];

        [AW_12(i),abswat12_1(i)]=callFalfaw(Temp1(i),safeMean(RNG11_1200{i}), ...
            opt.f12,opt.alphaw1200);
        [~,abswat12_2(i)]=callFalfaw(Temp1(i),safeMean(RNG22_1200{i}), ...
            opt.f12,opt.alphaw1200);
        [~,abswat12_3(i)]=callFalfaw(Temp1(i),safeMean(RNG33_1200{i}), ...
            opt.f12,opt.alphaw1200);
        [~,abswat12_4(i)]=callFalfaw(Temp1(i),safeMean(RNG44_1200{i}), ...
            opt.f12,opt.alphaw1200);

        alphaSed12_=alphaSed12(i);
        TL1_=(abswat12_1(i)+alphaSed12_).*R11;
        TL2_=(abswat12_2(i)+alphaSed12_).*R22;
        TL3_=(abswat12_3(i)+alphaSed12_).*R33;
        TL4_=(abswat12_4(i)+alphaSed12_).*R44;
        TL1_fix1200(i)=safeMean(TL1_); TL2_fix1200(i)=safeMean(TL2_);
        TL3_fix1200(i)=safeMean(TL3_); TL4_fix1200(i)=safeMean(TL4_);

        if isfinite(WS_12(i)) && WS_12(i)>0
            TL1_1200(i)=2.*(WF_12(i)./cosd(70).*abswat12_1(i) + ...
                nansumLocal(alpha1mean{i}.*WS_12(i)./cosd(70)));
            TL2_1200(i)=2.*(WF_12(i)./cosd(70).*abswat12_2(i) + ...
                nansumLocal(alpha2mean{i}.*WS_12(i)./cosd(70)));
            TL3_1200(i)=2.*(WF_12(i)./cosd(70).*abswat12_3(i) + ...
                nansumLocal(alpha3mean{i}.*WS_12(i)./cosd(70)));
            TL4_1200(i)=2.*(WF_12(i)./cosd(70).*abswat12_4(i) + ...
                nansumLocal(alpha4mean{i}.*WS_12(i)./cosd(70)));
        end

        SV1_1200{i}=opt.C1200 + Spread11 - Pdwb1200{i} - Af1 + ...
            2.*TL1_ + opt.kc.*(BT_RSSI1_1200{i}-opt.Er);
        SV2_1200{i}=opt.C1200 + Spread22 - Pdwb1200{i} - Af2 + ...
            2.*TL2_ + opt.kc.*(BT_RSSI2_1200{i}-opt.Er);
        SV3_1200{i}=opt.C1200 + Spread33 - Pdwb1200{i} - Af3 + ...
            2.*TL3_ + opt.kc.*(BT_RSSI3_1200{i}-opt.Er);
        SV4_1200{i}=opt.C1200 + Spread44 - Pdwb1200{i} - Af4 + ...
            2.*TL4_ + opt.kc.*(BT_RSSI4_1200{i}-opt.Er);

        SV11_1200{i}=-SL1200{i} + opt.kc.*BT_RSSI1_1200{i} + Spread11 + 2.*TL1_;
        SV22_1200{i}=-SL1200{i} + opt.kc.*BT_RSSI2_1200{i} + Spread22 + 2.*TL2_;
        SV33_1200{i}=-SL1200{i} + opt.kc.*BT_RSSI3_1200{i} + Spread33 + 2.*TL3_;
        SV44_1200{i}=-SL1200{i} + opt.kc.*BT_RSSI4_1200{i} + Spread44 + 2.*TL4_;

        SVBmean1_1200(i)=safeMean(SV1_1200{i});
        SVBmean2_1200(i)=safeMean(SV2_1200{i});
        SVBmean3_1200(i)=safeMean(SV3_1200{i});
        SVBmean4_1200(i)=safeMean(SV4_1200{i});
        SVB1mean1_1200(i)=safeMean(SV11_1200{i});
        SVB1mean2_1200(i)=safeMean(SV22_1200{i});
        SVB1mean3_1200(i)=safeMean(SV33_1200{i});
        SVB1mean4_1200(i)=safeMean(SV44_1200{i});
        SVBstd1_1200(i)=safeStd(SV1_1200{i});
        SVBstd2_1200(i)=safeStd(SV2_1200{i});
        SVBstd3_1200(i)=safeStd(SV3_1200{i});
        SVBstd4_1200(i)=safeStd(SV4_1200{i});
        SVB1std1_1200(i)=safeStd(SV11_1200{i});
        SVB1std2_1200(i)=safeStd(SV22_1200{i});
        SVB1std3_1200(i)=safeStd(SV33_1200{i});
        SVB1std4_1200(i)=safeStd(SV44_1200{i});

        svMatrix=[SV1_1200{i},SV2_1200{i},SV3_1200{i},SV4_1200{i}];
        stdSvSpatialALL_1200{i}=nanstdLocal(svMatrix,0,2);
        stdSvSpatial_1200(i)=safeMean(stdSvSpatialALL_1200{i});
        S_1200(i)=safeMean([SVBmean1_1200(i),SVBmean2_1200(i), ...
            SVBmean3_1200(i),SVBmean4_1200(i)]);

        if logical(opt.ApplyMethod2Regression)
            D50vals=[D50_14(i),D50_05(i)];
            D50vals=D50vals(isfinite(D50vals) & D50vals>0);
            D50_12(i)=safeMean(D50vals);
            if ~isfinite(D50_12(i)) || D50_12(i)<=0
                if isfinite(D50(i)) && D50(i)>0
                    D50_12(i)=D50(i);
                elseif isfinite(D_50grabM(i)) && D_50grabM(i)>0
                    D50_12(i)=D_50grabM(i);
                end
            end
            Var1=BTmodf_1200(i);
            Var2=FilteredDataN_1200(i);
            Var4=S_1200(i);
            Var5=BTstd_1200(i);
            if exist('TrateMethod2Mregression1200','file')~=2
                Method2Status_1200(i)="helper_missing";
            elseif ~all(isfinite([D50_12(i),Var1,Var2,Var4,Var5]))
                missingNames=strings(0,1);
                if ~isfinite(D50_12(i)), missingNames(end+1)="D50_12"; end %#ok<AGROW>
                if ~isfinite(Var1), missingNames(end+1)="BTmodf_1200"; end %#ok<AGROW>
                if ~isfinite(Var2), missingNames(end+1)="FilteredDataN_1200"; end %#ok<AGROW>
                if ~isfinite(Var4), missingNames(end+1)="S_1200"; end %#ok<AGROW>
                if ~isfinite(Var5), missingNames(end+1)="BTstd_1200"; end %#ok<AGROW>
                Method2Status_1200(i)="missing_"+strjoin(missingNames,"+");
            else
                xSAND_12=[-0.757403896527827,0.116047708004703, ...
                    0.553042745645718,0.101474754190850, ...
                   -0.198010090921993,1.85394968562851];
                xSANDgravel_12=[15.7618140988780,1.56257214251156, ...
                   -0.00162826490796587,3.77328682594328, ...
                   -0.00909604444541630,-1.16152304889946];
                xGravel_12=[-1.89356670251487,0.890393231016732, ...
                    0.264482313118918,0.394699875620202,60.6118408818348, ...
                    22.3354445706887,-20.9300150845519,2.03824955829931, ...
                   -0.178535192064460];
                qsMethod2regression_1200(i)=TrateMethod2Mregression1200( ...
                    D50_12(i),xSAND_12,xSANDgravel_12,xGravel_12, ...
                    Var1,Var2,Var4,Var5);
                if isfinite(qsMethod2regression_1200(i))
                    Method2Status_1200(i)="OK";
                else
                    Method2Status_1200(i)="helper_returned_NaN";
                end
            end
        end

        Message(i)="OK";

    catch ME
        Message(i)="Error: "+string(ME.message);
        warning('Backscatter processing failed for row %d, file %s:\n%s', ...
            i,ADCPfilename(i),ME.message);
    end
end

%% Convert variable-length cell series into padded legacy matrices
Depth1200 = cellsToPaddedMatrix(Depth1200_cell,nFiles);
Depth1200std = cellsToPaddedMatrix(Depth1200std_cell,nFiles);
BT_RSSI1200 = cellsToPaddedMatrix(BT_RSSI1200_cell,nFiles);
BT_RSSI_1200 = BT_RSSI1200;

%% Single output structure and summary table
BackscatterTable = table(InputRow,ADCPfilename,DistanceFromLeftBank, ...
    StartADCP,EndADCP,StartUsed_1200,EndUsed_1200,NensemblesRaw_1200, ...
    NensemblesUsed_1200,NfilteredBTvelocity_1200, ...
    BTvelocityMaskMode_1200,BTvelocityFilterSource_1200, ...
    BTvelocityFilterExact_1200,BTvelocityFilterAligned_1200, ...
    WF_12,WS_12,Temp1,stdSpatial_1200, ...
    stdSpatialRSSI_1200,R1_1200,R2_1200,R3_1200,R4_1200, ...
    SVBmean1_1200,SVBmean2_1200,SVBmean3_1200,SVBmean4_1200, ...
    S_1200,stdSvSpatial_1200,BTmodf_1200,FilteredDataN_1200, ...
    BTstd_1200,D50_12,qsMethod2regression_1200,Method2Status_1200, ...
    FileFound,MatchedFile,Message);

BackscatterOut = struct();
varsToExport = { ...
    'InputRow','ADCPfilename','DistanceFromLeftBank','StartADCP','EndADCP', ...
    'StartUsed_1200','EndUsed_1200','NensemblesRaw_1200', ...
    'NensemblesUsed_1200','NfilteredBTvelocity_1200', ...
    'BTvelocityFilterMask_1200','BTvelocityMaskMode_1200', ...
    'BTvelocityFilterSource_1200','BTvelocityFilterExact_1200', ...
    'BTvelocityFilterAligned_1200', ...
    'WF_12','WS_12','WT_12','Temp1', ...
    'Depth1200','Depth1200std','stdSpatial_1200','stdSpatialALL_1200', ...
    'RNG11_1200','RNG22_1200','RNG33_1200','RNG44_1200', ...
    'BT_RSSI1200','BT_RSSI_1200','BT_RSSI1_1200','BT_RSSI2_1200', ...
    'BT_RSSI3_1200','BT_RSSI4_1200','stdSpatialRSSI_1200', ...
    'psi1_1200','psi2_1200','psi3_1200','psi4_1200', ...
    'alfaXmean_1200','alfaYmean_1200', ...
    'PtransW1200','PtransWmean1200','Pdwb1200','STDPowerdb1200', ...
    'STDPW1200','XmitCC1200','XmitVV1200','PitchH1200','RollL1200', ...
    'R1_1200','R2_1200','R3_1200','R4_1200', ...
    'Spread1_1200','Spread2_1200','Spread3_1200','Spread4_1200', ...
    'SpreadM1200','LdbMean1200', ...
    'Dmin_1200_1','Dmin_1200_2','Dmin_1200_3','Dmin_1200_4', ...
    'Dmax_1200_1','Dmax_1200_2','Dmax_1200_3','Dmax_1200_4', ...
    'Area_1200_1','Area_1200_2','Area_1200_3','Area_1200_4', ...
    'distance_1200_1','distance_1200_2','distance_1200_3','distance_1200_4', ...
    'Lp_1200_1','Lp_1200_2','Lp_1200_3','Lp_1200_4', ...
    'Af1mean1200','Af2mean1200','Af3mean1200','Af4mean1200','SL1200', ...
    'alpha1mean','alpha2mean','alpha3mean','alpha4mean', ...
    'alpha1mean1200grad','alpha2mean1200grad', ...
    'alpha3mean1200grad','alpha4mean1200grad','alphaMean1200', ...
    'dEDR1_1200','dEDR2_1200','dEDR3_1200','dEDR4_1200', ...
    'AW_12','abswat12_1','abswat12_2','abswat12_3','abswat12_4', ...
    'TL1_1200','TL2_1200','TL3_1200','TL4_1200', ...
    'TL1_fix1200','TL2_fix1200','TL3_fix1200','TL4_fix1200', ...
    'SV1_1200','SV2_1200','SV3_1200','SV4_1200', ...
    'SV11_1200','SV22_1200','SV33_1200','SV44_1200', ...
    'SVBmean1_1200','SVBmean2_1200','SVBmean3_1200','SVBmean4_1200', ...
    'SVB1mean1_1200','SVB1mean2_1200','SVB1mean3_1200','SVB1mean4_1200', ...
    'SVBstd1_1200','SVBstd2_1200','SVBstd3_1200','SVBstd4_1200', ...
    'SVB1std1_1200','SVB1std2_1200','SVB1std3_1200','SVB1std4_1200', ...
    'stdSvSpatial_1200','stdSvSpatialALL_1200','S_1200', ...
    'BTmodf_1200','FilteredDataN_1200','BTstd_1200', ...
    'D50_14','D50_05','D50','D_50grabM','D50_12', ...
    'qsMethod2regression_1200','Method2Status_1200', ...
    'DummyValue','FileFound','MatchedFile','Message','BackscatterTable'};
for k = 1:numel(varsToExport)
    name = varsToExport{k};
    BackscatterOut.(name) = eval(name);
end
BackscatterOut.DtransON1m = DtransON1m;
BackscatterOut.StartEndMode = startEndMode;
BackscatterOut.VelocityMaskMode = velocityMaskMode;
BackscatterOut.Options = opt;

%% Assign to workspace
if logical(opt.AssignToCaller)
    assignStructToWorkspace(BackscatterOut,'caller');
end
if logical(opt.AssignToBase)
    assignStructToWorkspace(BackscatterOut,'base');
end

%% Display
if logical(opt.Verbose)
    fprintf('\n=== 1200-kHz backscatter processing summary ===\n');
    fprintf('Input rows: %d\n',nFiles);
    fprintf('Successfully processed: %d/%d\n',sum(Message=="OK"),nFiles);
    fprintf('Rows skipped by nodata1200: %d\n',sum(startsWith(Message,"Skipped")));
    fprintf('Velocity mask mode: %s.\n',velocityMaskMode);
    fprintf('Exact/aligned masks: %d/%d rows.\n', ...
        sum(BTvelocityFilterExact_1200 & BTvelocityFilterAligned_1200),nFiles);
    fprintf('Velocity-filtered ensembles assigned NaN: %d\n',sum(NfilteredBTvelocity_1200));
    fprintf('Start/end mode: %s.\n',startEndMode);
    fprintf('Method 2 regression helper available: %d\n', ...
        exist('TrateMethod2Mregression1200','file')==2);
    fprintf('Finite Method 2 outputs: %d/%d\n', ...
        sum(isfinite(qsMethod2regression_1200)),nFiles);
    fprintf('No files were saved.\n');
    fprintf('================================================\n');
    disp(BackscatterTable);
end

end

%% ------------------------------------------------------------------------
function fileList = normaliseFileList(inputValue)
if ischar(inputValue) || isstring(inputValue)
    folder = char(string(inputValue));
    fileList = dir(fullfile(folder,'*.mat'));
elseif isstruct(inputValue)
    fileList = inputValue;
else
    fileList = struct([]);
end
end

function matchedFile = matchADCPfile(fileList,targetName)
matchedFile = "";
if isempty(fileList) || ~isfield(fileList,'name'), return; end

target = cleanFileName(targetName);
[~,targetBase,~] = fileparts(char(target));
names = string({fileList.name}');
cleanNames = strings(numel(names),1);
bases = strings(numel(names),1);
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
% By default StartADCP and EndADCP are trim counts, not absolute indices.
mode=lower(strtrim(string(mode)));
if mode=="trimcounts" || mode=="trim"
    if ~isfinite(startValue) || startValue<0, nTrimStart=0; else, nTrimStart=round(startValue); end
    if ~isfinite(endValue) || endValue<0, nTrimEnd=0; else, nTrimEnd=round(endValue); end
    first=1+nTrimStart;
    last=nRaw-nTrimEnd;
else
    if ~isfinite(startValue) || startValue<1, first=1; else, first=round(startValue); end
    if ~isfinite(endValue) || endValue<1, last=nRaw; else, last=round(endValue); end
end
first=max(1,min(first,nRaw));
last=max(1,min(last,nRaw));
if last<first
    error(['The selected StartADCP/EndADCP values leave no ensembles. ', ...
        'Resolved interval is %d:%d for a raw series of length %d.'], ...
        first,last,nRaw);
end
end

function tf=isNoDataRow(spec,i,nFiles)
tf=false;
if isempty(spec), return; end
if islogical(spec)
    spec=spec(:);
    if numel(spec)==nFiles, tf=spec(i); end
else
    tf=any(round(spec(:))==i);
end
end

function x=getVectorFieldAny(S,names,defaultValue,n)
x=defaultValue;
if isempty(S) || ~isstruct(S)
    x=prepareVector(x,n); return
end
for k=1:numel(names)
    if isfield(S,names{k})
        x=S.(names{k}); break
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
elseif numel(x)>n, x=x(1:n); end
end

function val=getScalarPreferred(S,names,i,nFiles,defaultValue)
val=defaultValue;
if isempty(S) || ~isstruct(S), return; end
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    if iscell(C)
        if numel(C)>=i, val=safeMean(C{i}); end
    elseif isnumeric(C)
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
if isempty(S) || ~isstruct(S), return; end
nSeg=idxEnd-idxStart+1;
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    candidate=[];
    if iscell(C)
        if numel(C)>=i && ~isempty(C{i}), candidate=double(C{i}(:)); end
    elseif isnumeric(C)
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
    if numel(candidate)==nSeg
        x=candidate; return
    elseif numel(candidate)>=idxEnd && numel(candidate)>=nRaw
        x=candidate(idxStart:idxEnd); return
    elseif numel(candidate)>=idxEnd
        x=candidate(idxStart:idxEnd); return
    end
end
end

function x=forceColumnLength(x,n)
x=double(x(:));
if numel(x)<n, x=[x;nan(n-numel(x),1)];
elseif numel(x)>n, x=x(1:n); end
end

function [mask,source,isExact,isAligned]=deriveVelocityRejectMask( ...
    BTvelOut,i,nFiles,idxStart,idxEnd,nRaw,nSeg,mode)
%DERIVEVELOCITYREJECTMASK Return one mask aligned with the trimmed row.
%
% mode='final' uses only the explicit final mask or final per-file BT time
% series. Intermediate and last-successful-row aliases are deliberately not
% combined because doing so can apply a mask from the wrong file.
% mode='boat' uses only peaksBTxy_1200.
% mode='none' returns an all-false mask.

mask=false(nSeg,1);
source="none";
isExact=false;
isAligned=false;
mode=lower(strtrim(string(mode)));

if mode=="none"
    source="disabled";
    isExact=true;
    isAligned=true;
    return
end

if isempty(BTvelOut) || ~isstruct(BTvelOut)
    source="BTvelOut_missing";
    return
end

% Verify that the BT function used the same raw ensemble interval.
btStart=getScalarFieldAny(BTvelOut, ...
    {'BTfilterMaskStart_1200','startAll'},i,nFiles);
btEnd=getScalarFieldAny(BTvelOut, ...
    {'BTfilterMaskEnd_1200','finishAll'},i,nFiles);
if isfinite(btStart) && round(btStart)~=idxStart
    source="trim_start_mismatch";
    return
end
if isfinite(btEnd) && round(btEnd)~=idxEnd
    source="trim_end_mismatch";
    return
end

switch mode
    case "final"
        % Preferred explicit final rejection mask from the updated BT
        % processing function. true means rejected/invalid.
        if isfield(BTvelOut,'BTfilterMask_1200')
            vals=getRowContainer(BTvelOut.BTfilterMask_1200,i,nFiles);
            [candidate,ok]=alignLogicalMaskExact( ...
                vals,idxStart,idxEnd,nRaw,nSeg);
            if ok
                mask=candidate;
                source="BTfilterMask_1200";
                isExact=true;
                isAligned=true;
                return
            end
        end

        % Equivalent retained-sample mask, if supplied.
        if isfield(BTvelOut,'BTvalidMask_1200')
            vals=getRowContainer(BTvelOut.BTvalidMask_1200,i,nFiles);
            [candidate,ok]=alignLogicalMaskExact( ...
                vals,idxStart,idxEnd,nRaw,nSeg);
            if ok
                mask=~candidate;
                source="BTvalidMask_1200_inverted";
                isExact=true;
                isAligned=true;
                return
            end
        end

        % Backward-compatible per-file final series. These are safe because
        % they are row-wise cells; last-successful-row aliases such as
        % BTE_1200_ and BTN_1200_ are intentionally excluded.
        finalSeriesFields={'BTmod_final_1200','BTmodf_1200ALL'};
        for k=1:numel(finalSeriesFields)
            if ~isfield(BTvelOut,finalSeriesFields{k}), continue; end
            vals=getRowContainer(BTvelOut.(finalSeriesFields{k}),i,nFiles);
            [aligned,ok]=alignSeriesExact(vals,idxStart,idxEnd,nRaw,nSeg);
            if ok
                mask=~isfinite(aligned);
                source=string(finalSeriesFields{k})+"_NaN";
                isExact=true;
                isAligned=true;
                return
            end
        end

        % Final x/y components are another exact representation.
        if isfield(BTvelOut,'BTx_final_1200') && ...
                isfield(BTvelOut,'BTy_final_1200')
            xvals=getRowContainer(BTvelOut.BTx_final_1200,i,nFiles);
            yvals=getRowContainer(BTvelOut.BTy_final_1200,i,nFiles);
            [x,okx]=alignSeriesExact(xvals,idxStart,idxEnd,nRaw,nSeg);
            [y,oky]=alignSeriesExact(yvals,idxStart,idxEnd,nRaw,nSeg);
            if okx && oky
                mask=~isfinite(x) | ~isfinite(y);
                source="BTx_final_1200+BTy_final_1200_NaN";
                isExact=true;
                isAligned=true;
                return
            end
        end

        source="final_mask_not_found_or_not_aligned";

    case "boat"
        if isfield(BTvelOut,'peaksBTxy_1200')
            vals=getRowContainer(BTvelOut.peaksBTxy_1200,i,nFiles);
            [candidate,ok]=alignIndexOrMaskExact( ...
                vals,idxStart,idxEnd,nRaw,nSeg);
            if ok
                mask=candidate;
                source="peaksBTxy_1200";
                isExact=true;
                isAligned=true;
                return
            end
        end
        source="boat_mask_not_found_or_not_aligned";
end
end

function value=getScalarFieldAny(S,names,i,nFiles)
value=NaN;
for k=1:numel(names)
    if ~isfield(S,names{k}), continue; end
    C=S.(names{k});
    if iscell(C)
        if numel(C)>=i && ~isempty(C{i}), value=safeMean(C{i}); end
    elseif isnumeric(C) || islogical(C)
        C=double(C);
        if isvector(C) && numel(C)>=i
            value=C(i);
        elseif size(C,1)==nFiles && size(C,2)>=1
            value=C(i,1);
        elseif size(C,2)==nFiles && size(C,1)>=1
            value=C(1,i);
        end
    end
    if isfinite(value), return; end
end
end

function [mask,ok]=alignIndexOrMaskExact(vals,idxStart,idxEnd,nRaw,nSeg)
mask=false(nSeg,1);
ok=false;
if isempty(vals), return; end
if islogical(vals)
    [mask,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg);
    return
end
valsRaw=double(vals(:));
finiteVals=valsRaw(isfinite(valsRaw));
if ~isempty(finiteVals) && all(ismember(finiteVals,[0 1])) && ...
        any(numel(valsRaw)==[nSeg,nRaw])
    [mask,ok]=alignLogicalMaskExact(logical(valsRaw), ...
        idxStart,idxEnd,nRaw,nSeg);
    return
end
vals=round(finiteVals);
if isempty(vals)
    ok=true;
    return
end
if all(vals>=1 & vals<=nSeg)
    mask(vals)=true;
    ok=true;
elseif all(vals>=idxStart & vals<=idxEnd)
    mask(vals-idxStart+1)=true;
    ok=true;
end
end

function [mask,ok]=alignLogicalMaskExact(vals,idxStart,idxEnd,nRaw,nSeg)
mask=false(nSeg,1);
ok=false;
if isempty(vals), return; end
vals=logical(vals(:));
if numel(vals)==nSeg
    mask=vals;
    ok=true;
elseif numel(vals)==nRaw
    mask=vals(idxStart:idxEnd);
    ok=true;
end
end

function [x,ok]=alignSeriesExact(vals,idxStart,idxEnd,nRaw,nSeg)
x=[];
ok=false;
if isempty(vals) || ~(isnumeric(vals) || islogical(vals)), return; end
vals=double(vals(:));
if numel(vals)==nSeg
    x=vals;
    ok=true;
elseif numel(vals)==nRaw
    x=vals(idxStart:idxEnd);
    ok=true;
end
end

function vals=getRowContainer(C,i,nFiles)
vals=[];
if iscell(C)
    if numel(C)>=i, vals=C{i}; end
elseif isnumeric(C) || islogical(C)
    if isvector(C)
        if nFiles==1 || numel(C)~=nFiles, vals=C(:); end
    elseif size(C,1)==nFiles
        vals=C(i,:).';
    elseif size(C,2)==nFiles
        vals=C(:,i);
    end
end
end

function y=despikeRangeLocal(x,T)
y=double(x(:));
valid=isfinite(y);
if sum(valid)<3, return; end
m=mean(y(valid)); s=std(y(valid));
if isfinite(s) && s>0
    y(abs(y-m)>T.*s)=NaN;
end
end

function BT=despikeBTstrength(BTv,T,MM,Window,ValuesOUt,Color,PLOT,RSSIMax)
%DESPIKEBTSTRENGTH Despike bottom-track RSSI while retaining legacy name.
A=double(BTv(:));
A(A>RSSIMax)=NaN;
if ~isempty(ValuesOUt)
    if islogical(ValuesOUt)
        m=logical(ValuesOUt(:));
        if numel(m)==numel(A), A(m)=NaN; end
    else
        idx=round(double(ValuesOUt(:)));
        idx=idx(isfinite(idx) & idx>=1 & idx<=numel(A));
        A(idx)=NaN;
    end
end
if logical(MM)
    A=movmean(A,round(Window),'omitnan');
end
valid=isfinite(A);
if sum(valid)>=2
    T1=mean(A(valid))-T.*std(A(valid));
    To=mean(A(valid))+T.*std(A(valid));
    A(A>To | A<T1)=NaN;
end
if contains(upper(string(PLOT)),"YES")
    figure('Color','w'); hold on;
    bad=isnan(A);
    plot(find(bad),BTv(bad),'o','MarkerSize',8);
    plot(A,Color);
    grid on; box on;
end
BT=A;
end

function [Dmin,Dmax,Area,distance,Lp]=SamplingVolume( ...
    distanceToBottom,SlantingAngle,transducerWidth,beamopen,radcp,PLength)
%SAMPLINGVOLUME Legacy sampling-geometry formulation.
Dmin=NaN; Dmax=NaN; Area=NaN; distance=NaN; Lp=NaN;
if ~all(isfinite([distanceToBottom,SlantingAngle,transducerWidth,beamopen,radcp])) || ...
        distanceToBottom<=0 || transducerWidth<=0
    return
end
Dmax=((distanceToBottom+transducerWidth/2*sind(SlantingAngle))/ ...
    tand(90-SlantingAngle-beamopen) - ...
    (distanceToBottom-transducerWidth/2*sind(SlantingAngle))/ ...
    tand(90-SlantingAngle+beamopen)) + ...
    2*transducerWidth/2*cosd(SlantingAngle);
Dmin=2*(distanceToBottom*tand(beamopen)/cosd(SlantingAngle)+transducerWidth/2);
Area=pi*Dmax*Dmin;
distance=distanceToBottom*tand(SlantingAngle)+radcp;
Lp=safeMean(double(PLength(:)).*sind(90-SlantingAngle-beamopen));
end

function Af=surfaceScatteringArea(Ldbm,R,omega,theta,alphaPrimary,alphaSecondary)
Af=nan(size(R));
if isempty(R), return; end
num=Ldbm.*R.*omega;
den=2.*sind(theta-alphaPrimary).*cosd(alphaSecondary);
arg=num./den;
valid=isfinite(arg) & arg>0;
Af(valid)=10.*log10(arg(valid));
Af(~isreal(Af))=NaN;
Af=real(Af);
end

function E=subsetEchoMatrix(Eraw,idxStart,idxEnd,nRaw,nSeg)
E=double(Eraw);
if isempty(E), E=nan(1,nSeg); return; end
if size(E,2)>=idxEnd && (size(E,2)>=nRaw || size(E,2)>=idxEnd)
    E=E(:,idxStart:idxEnd);
elseif size(E,1)>=idxEnd && size(E,2)<idxEnd
    E=E(idxStart:idxEnd,:).';
else
    if size(E,2)<nSeg
        E(:,end+1:nSeg)=NaN;
    elseif size(E,2)>nSeg
        E=E(:,1:nSeg);
    end
end
end

function id=resolveEchoStop(H,rangeMean,nRows)
id=nRows;
if isfinite(rangeMean)
    tmp=find(H(1:min(numel(H),nRows))>rangeMean,1);
    if ~isempty(tmp), id=tmp; end
end
id=max(1,min(id,nRows));
end

function alpha=echoAttenuation(E,idStop,kc)
alpha=nan(size(E));
last=min(idStop-1,size(E,1)-1);
for jj=2:last
    alpha(jj,:)=E(jj,:).*kc-E(jj+1,:).*kc;
end
alpha(alpha<0 | alpha>20)=NaN;
end

function [AW,abswat]=callFalfaw(Temp,Depth,f12,fallback)
AW=fallback; abswat=fallback;
if exist('falfaw','file')==2 && all(isfinite([Temp,Depth,f12])) && Depth>0
    try
        [AWtmp,abstmp]=falfaw(Temp,Depth,f12);
        if any(isfinite(AWtmp(:))), AW=safeMean(AWtmp); end
        if any(isfinite(abstmp(:))), abswat=safeMean(abstmp); end
    catch
    end
end
end

function M=cellsToPaddedMatrix(C,nRows)
maxLen=0;
for i=1:min(nRows,numel(C)), maxLen=max(maxLen,numel(C{i})); end
M=nan(nRows,maxLen);
for i=1:min(nRows,numel(C))
    x=double(C{i}(:));
    M(i,1:numel(x))=x.';
end
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
A=double(A);
mu=nanmeanLocal(A,dim);
rep=ones(1,ndims(A)); rep(dim)=size(A,dim);
dev=A-repmat(mu,rep);
valid=isfinite(A); dev(~valid)=0;
n=sum(valid,dim);
if flag==0, den=max(n-1,0); else, den=n; end
y=sqrt(sum(dev.^2,dim)./den);
y(den<=0)=NaN;
end

function y=nansumLocal(x)
x=double(x(:)); x=x(isfinite(x));
if isempty(x), y=0; else, y=sum(x); end
end

function assignStructToWorkspace(S,workspaceName)
f=fieldnames(S);
for k=1:numel(f)
    assignin(workspaceName,f{k},S.(f{k}));
end
end
