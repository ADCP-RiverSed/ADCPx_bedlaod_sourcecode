function ValidationReport = validateFeatureExtractionStructuresM9( ...
    tratesFile, shearFile, backscatterFile, varargin)
%VALIDATEFEATUREEXTRACTIONSTRUCTURESM9 Validate saved MAT files against
% the exact structures and fields used by the supplied feature-extraction
% script.
%
% ValidationReport = validateFeatureExtractionStructuresM9( ...
%     tratesFile,shearFile,backscatterFile,'ThrowOnMissing',true)
%
% Required top-level variables:
%   tratesFile:      Data_field, BTvel, Q
%   shearFile:       ShearData
%   backscatterFile: BTamp, DataBTamp
%
% The three file names are checked for the same leading YYYYMMDD
% prefix and the same optional location token, for example:
%   20170322_Tomorrica_TratesData.mat
%
% The function checks field presence, placeholder values, row counts and
% the equality of BTamp and DataBTamp. It also confirms that no saved field
% name still contains 1200 or an isolated abbreviated 12 instrument token.
% A 12 embedded in a larger number, such as QPer0125D, is permitted.
% Optional M9 vertical-beam data are checked when VBAvailable_S is true.
% The function does not modify the files.

p=inputParser;
addRequired(p,'tratesFile',@(x)ischar(x)||isstring(x));
addRequired(p,'shearFile',@(x)ischar(x)||isstring(x));
addRequired(p,'backscatterFile',@(x)ischar(x)||isstring(x));
addParameter(p,'ThrowOnMissing',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOnPlaceholder',false,@(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOn1200Fields',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOnInvalidVerticalBackscatter',false, ...
    @(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOnInvalidAcousticBackscatter',false, ...
    @(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOnDateMismatch',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'ExpectedLocationName',"",@(x)ischar(x)||isstring(x));
addParameter(p,'ThrowOnLocationMismatch',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'Verbose',true,@(x)islogical(x)||isnumeric(x));
parse(p,tratesFile,shearFile,backscatterFile,varargin{:});
opt=p.Results;

tratesFile=char(string(tratesFile));
shearFile=char(string(shearFile));
backscatterFile=char(string(backscatterFile));

assert(exist(tratesFile,'file')==2,'Trates MAT file not found: %s',tratesFile);
assert(exist(shearFile,'file')==2,'Shear MAT file not found: %s',shearFile);
assert(exist(backscatterFile,'file')==2,'Backscatter MAT file not found: %s',backscatterFile);

D=load(tratesFile);
S=load(shearFile);
B=load(backscatterFile);

topRequired = { ...
    'tratesFile','Data_field'; ...
    'tratesFile','BTvel'; ...
    'tratesFile','Q'; ...
    'shearFile','ShearData'; ...
    'backscatterFile','BTamp'; ...
    'backscatterFile','DataBTamp'};

TopFileRole=strings(size(topRequired,1),1);
TopVariable=strings(size(topRequired,1),1);
TopExists=false(size(topRequired,1),1);
for i=1:size(topRequired,1)
    role=topRequired{i,1}; varName=topRequired{i,2};
    TopFileRole(i)=string(role); TopVariable(i)=string(varName);
    switch role
        case 'tratesFile', C=D;
        case 'shearFile', C=S;
        otherwise, C=B;
    end
    TopExists(i)=isfield(C,varName);
end
TopLevel=table(TopFileRole,TopVariable,TopExists);

if isfield(D,'Data_field'), Data_field=D.Data_field; else, Data_field=struct(); end
if isfield(S,'ShearData'), ShearData=S.ShearData; else, ShearData=struct(); end
if isfield(D,'BTvel'), BTvel=D.BTvel; else, BTvel=struct(); end
if isfield(D,'Q'), Q=D.Q; else, Q=struct(); end %#ok<NASGU>
if isfield(B,'BTamp'), BTamp=B.BTamp; else, BTamp=struct(); end
if isfield(B,'DataBTamp'), DataBTamp=B.DataBTamp; else, DataBTamp=struct(); end

featureDataFields = { ...
    'MeanTrate','MaxTrate','MinTrate', ...
    'QPer6D','QPer4D','QPer2D','QPer0125D','QPer063D', ...
    'Q84_14','Q84_05','Q16_14','Q16_05','Dmass05','Dmass14', ...
    'D_50grabM','D_50grabMAX','D_90grabM','D_90grabMAX', ...
    'D_16grabM','D_16grabMAX','D_84grabM','D_84grabMAX', ...
    'D84_14','D16_14','D84_05','D16_05','LeftBankDist14','D90_05'};
featureShearFields = { ...
    'Depth_S','VmodAll_S', ...
    'BT1_D_Sstd','BT2_D_Sstd','BT3_D_Sstd','BT4_D_Sstd', ...
    'Tau_All_S_1','Tau_All_S_2','Ushear1_S','Ushear2_S', ...
    'UsheaAll_S','Cb2_S','nCb_new_S','ucr_S','Tucr_S', ...
    'TKE_prof_S','tau_rm_S','tau_Mann_star_S','W_i_S', ...
    'W_iVRshear2_S','Tau_shear2VR_S','Dv1_S','Tv1_S', ...
    'HEightDUNE_S','Vmod1low_SstdMOD','Vmod_SstdMOD', ...
    'WT_S','u_shearH_S'};
featureBTvelFields = { ...
    'BTmodf_S','BTstd_Sa','BTerror_Sm', ...
    'FilteredDataN_S','FilteredDataE_S', ...
    'BT_stdf_S','BTstd_S','DeltaT_S'};
featureBTampFields = { ...
    'Ks2_S','D50_05','Depth_S','Depth_Sstd', ...
    'Af1mean_S','Af2mean_S','Af3mean_S','Af4mean_S', ...
    'Spread1_S','Spread2_S','Spread3_S','Spread4_S', ...
    'SVBmean1_S','SVBmean2_S','SVBmean3_S','SVBmean4_S', ...
    'SVB1mean1_S','SVB1mean2_S','SVB1mean3_S','SVB1mean4_S', ...
    'SVBstd1_S','SVBstd2_S','SVBstd3_S','SVBstd4_S', ...
    'SV1_S','SV2_S','SV3_S','SV4_S', ...
    'SVB1std1_S','SVB1std2_S','SVB1std3_S','SVB1std4_S', ...
    'BT_RSSI_S','alphaMean_S','stdSpatialRSSI_S', ...
    'stdSvSpatial_S','alpha1mean_S','alpha2mean_S', ...
    'alpha3mean_S','alpha4mean_S'};

nRows=inferRows(Data_field,ShearData,BTvel,BTamp);
Fields=[checkFields(Data_field,'Data_field',featureDataFields,nRows); ...
        checkFields(ShearData,'ShearData',featureShearFields,nRows); ...
        checkFields(BTvel,'BTvel',featureBTvelFields,nRows); ...
        checkFields(BTamp,'BTamp',featureBTampFields,nRows)];

BTampAliasEqual=isstruct(BTamp)&&isstruct(DataBTamp)&&isequaln(BTamp,DataBTamp);
RowSummary=table(nRows,fieldLength(Data_field,'MeanTrate'), ...
    fieldLength(ShearData,'Depth_S'),fieldLength(BTvel,'BTmodf_S'), ...
    fieldLength(BTamp,'SVBmean1_S'), ...
    'VariableNames',{'ExpectedRows','DataRows','ShearRows','BTvelRows','BackscatterRows'});

RemainingLegacyInstrumentFields=unique([ ...
    collectLegacyInstrumentNames(Data_field,"Data_field"); ...
    collectLegacyInstrumentNames(ShearData,"ShearData"); ...
    collectLegacyInstrumentNames(BTvel,"BTvel"); ...
    collectLegacyInstrumentNames(Q,"Q"); ...
    collectLegacyInstrumentNames(BTamp,"BTamp"); ...
    collectLegacyInstrumentNames(DataBTamp,"DataBTamp")],'stable');

Measurement=(1:nRows).';
VBAvailable=logicalRows(BTamp,'VBAvailable_S',nRows,false);
VerticalReady=logicalRows(BTamp,'VerticalBackscatterReady_S',nRows,false);
VBContrastFinite=finiteRows(BTamp,'SV_VBmeanS',nRows);
VBStrengthFinite=finiteRows(BTamp,'SV_VBmeanSrssi',nRows);
VBValid=~VBAvailable | (VerticalReady & ...
    (VBContrastFinite | VBStrengthFinite));
M9VerticalBackscatter=table(Measurement,VBAvailable,VerticalReady, ...
    VBContrastFinite,VBStrengthFinite,VBValid);
InvalidVerticalRows=Measurement(~VBValid);

[AcousticBackscatterSummary,AcousticBackscatterFields] = ...
    checkM9AcousticsAndBackscatter(BTamp,any(VBAvailable));
InvalidAcousticBackscatter = ...
    any(AcousticBackscatterSummary.Required & ...
    ~AcousticBackscatterSummary.OverallReady);

[TratesDatePrefix,TratesBaseName] = leadingDatePrefix(tratesFile);
[ShearDatePrefix,ShearBaseName] = leadingDatePrefix(shearFile);
[BackscatterDatePrefix,BackscatterBaseName] = leadingDatePrefix(backscatterFile);

TratesLocationToken = extractLocationToken(TratesBaseName,TratesDatePrefix,"TratesData");
ShearLocationToken = extractLocationToken(ShearBaseName,ShearDatePrefix,"ShearData");
BackscatterLocationToken = extractLocationToken( ...
    BackscatterBaseName,BackscatterDatePrefix,"BackscatterData");

FileIdentity = table( ...
    ["Trates";"Shear";"Backscatter"], ...
    [string(tratesFile);string(shearFile);string(backscatterFile)], ...
    [TratesBaseName;ShearBaseName;BackscatterBaseName], ...
    [TratesDatePrefix;ShearDatePrefix;BackscatterDatePrefix], ...
    [TratesLocationToken;ShearLocationToken;BackscatterLocationToken], ...
    'VariableNames',{'FileRole','FilePath','BaseName','DatePrefix','LocationToken'});
DatePrefixes = FileIdentity(:,{'FileRole','FilePath','BaseName','DatePrefix'});
DatePrefixesPresent = all(FileIdentity.DatePrefix ~= "");
DatePrefixesMatch = DatePrefixesPresent && ...
    numel(unique(FileIdentity.DatePrefix))==1;
LocationTokensPresent = all(FileIdentity.LocationToken ~= "");
LocationTokensMatch = numel(unique(FileIdentity.LocationToken))==1;
ExpectedLocationToken = sanitizeLocationToken(opt.ExpectedLocationName);
if ExpectedLocationToken == ""
    ExpectedLocationMatches = true;
else
    ExpectedLocationMatches = all(FileIdentity.LocationToken == ExpectedLocationToken);
end

ValidationReport=struct();
ValidationReport.TopLevel=TopLevel;
ValidationReport.FileIdentity=FileIdentity;
ValidationReport.DatePrefixes=DatePrefixes;
ValidationReport.DatePrefixesPresent=DatePrefixesPresent;
ValidationReport.DatePrefixesMatch=DatePrefixesMatch;
ValidationReport.LocationTokensPresent=LocationTokensPresent;
ValidationReport.LocationTokensMatch=LocationTokensMatch;
ValidationReport.ExpectedLocationToken=ExpectedLocationToken;
ValidationReport.ExpectedLocationMatches=ExpectedLocationMatches;
ValidationReport.Fields=Fields;
ValidationReport.RowSummary=RowSummary;
ValidationReport.BTampAliasEqual=BTampAliasEqual;
ValidationReport.RemainingLegacyInstrumentFields= ...
    RemainingLegacyInstrumentFields;
% Retained as a compatibility alias; it now reports both legacy tokens.
ValidationReport.Remaining1200Fields=RemainingLegacyInstrumentFields;
ValidationReport.AllInstrumentFieldsUseS= ...
    isempty(RemainingLegacyInstrumentFields);
ValidationReport.M9VerticalBackscatter=M9VerticalBackscatter;
ValidationReport.InvalidVerticalRows=InvalidVerticalRows;
ValidationReport.AcousticBackscatterSummary=AcousticBackscatterSummary;
ValidationReport.AcousticBackscatterFields=AcousticBackscatterFields;
ValidationReport.InvalidAcousticBackscatter=InvalidAcousticBackscatter;
ValidationReport.MissingTopLevel=TopVariable(~TopExists);
ValidationReport.MissingFields=Fields.StructureField(~Fields.Exists);
ValidationReport.PlaceholderFields=Fields.StructureField(Fields.IsPlaceholder);
ValidationReport.RowIncompatibleFields=Fields.StructureField(~Fields.RowCompatible);
ValidationReport.CriticalUnresolved=intersect(ValidationReport.PlaceholderFields,[ ...
    "BTvel.DeltaT_S";"BTamp.Ks2_S"; ...
    "ShearData.Tau_shear2VR_S"; ...
    "ShearData.Vmod_SstdMOD"; ...
    "ShearData.Vmod1low_SstdMOD"]);

if logical(opt.ThrowOnMissing)
    if any(~TopExists)
        error('Missing required top-level MAT variables: %s', ...
            strjoin(TopVariable(~TopExists),', '));
    end
    if any(~Fields.Exists)
        error('Missing required structure fields: %s', ...
            strjoin(Fields.StructureField(~Fields.Exists),', '));
    end
end
if logical(opt.ThrowOn1200Fields) && ...
        ~isempty(RemainingLegacyInstrumentFields)
    error(['M9 structure field names must use _S. Fields still containing ', ...
        '1200 or an isolated 12 token: %s'], ...
        strjoin(RemainingLegacyInstrumentFields,', '));
end
if logical(opt.ThrowOnInvalidVerticalBackscatter) && ...
        ~isempty(InvalidVerticalRows)
    error(['Vertical beam is marked available but has no valid processed ', ...
        'backscatter at measurement row(s): %s'], ...
        strjoin(string(InvalidVerticalRows),', '));
end
if logical(opt.ThrowOnInvalidAcousticBackscatter) && ...
        InvalidAcousticBackscatter
    bad=AcousticBackscatterSummary.BeamSystem( ...
        AcousticBackscatterSummary.Required & ...
        ~AcousticBackscatterSummary.OverallReady);
    error(['Required M9 acoustic/backscatter data are incomplete for: %s. ', ...
        'Inspect ValidationReport.AcousticBackscatterFields.'], ...
        strjoin(bad,', '));
end
if logical(opt.ThrowOnDateMismatch) && ~DatePrefixesMatch
    error(['Saved MAT filenames must start with the same YYYYMMDD prefix. ', ...
        'Observed prefixes: %s'],strjoin(DatePrefixes.DatePrefix,', '));
end

if logical(opt.ThrowOnLocationMismatch) && ...
        (~LocationTokensMatch || ~ExpectedLocationMatches)
    error(['Saved MAT filenames must use the same location token after ', ...
        'the date. Observed tokens: %s | Expected: %s'], ...
        strjoin(FileIdentity.LocationToken,', '),ExpectedLocationToken);
end

if logical(opt.ThrowOnPlaceholder) && any(Fields.IsPlaceholder)
    error('Required structure fields still contain placeholders: %s', ...
        strjoin(Fields.StructureField(Fields.IsPlaceholder),', '));
end

if logical(opt.Verbose)
    fprintf('\n=== M9 feature-extraction structure validation ===\n');
    fprintf('Top-level variables present: %d/%d\n',sum(TopExists),numel(TopExists));
    fprintf('Date prefixes present and matching: %d\n',DatePrefixesMatch);
    fprintf('Location tokens present: %d | matching: %d\n', ...
        LocationTokensPresent,LocationTokensMatch);
    if ExpectedLocationToken ~= ""
        fprintf('Expected location token: %s | matched: %d\n', ...
            ExpectedLocationToken,ExpectedLocationMatches);
    end
    disp(FileIdentity);
    fprintf('Required fields present: %d/%d\n',sum(Fields.Exists),height(Fields));
    fprintf('Required fields containing placeholders: %d\n',sum(Fields.IsPlaceholder));
    fprintf('Row-compatible required fields: %d/%d\n',sum(Fields.RowCompatible),height(Fields));
    fprintf('BTamp and DataBTamp are identical: %d\n',BTampAliasEqual);
    fprintf('All instrument-specific field names use _S: %d\n', ...
        isempty(RemainingLegacyInstrumentFields));
    fprintf('Vertical beam available: %d row(s) | invalid: %d row(s)\n', ...
        sum(VBAvailable),numel(InvalidVerticalRows));
    disp(AcousticBackscatterSummary);
    badAcousticFields=AcousticBackscatterFields.Required & ...
        ~AcousticBackscatterFields.HasUsableValue;
    if any(badAcousticFields)
        fprintf('Missing/invalid required acoustic or backscatter groups:\n');
        disp(AcousticBackscatterFields(badAcousticFields,:));
    end
    disp(RowSummary);
    if ~isempty(ValidationReport.CriticalUnresolved)
        fprintf('Critical unresolved quantities used by the reader:\n');
        disp(ValidationReport.CriticalUnresolved);
    end
    fprintf('================================================\n');
end
end

function [prefix,baseName] = leadingDatePrefix(filePath)
[~,base,~]=fileparts(char(string(filePath)));
baseName=string(base);
m=regexp(base,'^((?:19|20)\d{6})','tokens','once');
if isempty(m)
    prefix="";
else
    prefix=string(m{1});
end
end

function token = extractLocationToken(baseName,datePrefix,roleSuffix)
baseName=string(baseName);
datePrefix=string(datePrefix);
roleSuffix=string(roleSuffix);
token="";
if datePrefix==""
    return
end
prefix=datePrefix+"_";
suffix="_"+roleSuffix;
if baseName == datePrefix+"_"+roleSuffix
    token="";
    return
end
if ~startsWith(baseName,prefix) || ~endsWith(baseName,suffix)
    return
end
middle=extractBetween(baseName,strlength(prefix)+1, ...
    strlength(baseName)-strlength(suffix));
if isempty(middle)
    token="";
else
    token=string(middle);
end
end

function token = sanitizeLocationToken(value)
value=strtrim(string(value));
value=erase(value,"'");
value=erase(value,'"');
if ismissing(value) || value=="" || lower(value)=="nan"
    token="";
    return
end
token=string(regexprep(char(value),'[^A-Za-z0-9]+','_'));
token=string(regexprep(char(token),'^_+|_+$',''));
token=string(regexprep(char(token),'_+','_'));
end

function T=checkFields(S,structureName,fields,nRows)
n=numel(fields);
Structure=repmat(string(structureName),n,1);
Field=string(fields(:));
StructureField=Structure+"."+Field;
Exists=false(n,1); IsPlaceholder=false(n,1);
RowCompatible=false(n,1); SizeText=strings(n,1);
for i=1:n
    f=fields{i}; Exists(i)=isstruct(S)&&isfield(S,f);
    if Exists(i)
        x=S.(f); IsPlaceholder(i)=placeholder(x);
        SizeText(i)=strjoin(string(size(x)),'x');
        RowCompatible(i)=rowCompatible(x,nRows,f);
    else
        IsPlaceholder(i)=true; SizeText(i)="missing";
    end
end
T=table(Structure,Field,StructureField,Exists,IsPlaceholder,RowCompatible,SizeText);
end

function n=inferRows(varargin)
n=1;
keys={'MeanTrate','Depth_S','BTmodf_S','SVBmean1_S'};
for i=1:nargin
    S=varargin{i}; if ~isstruct(S), continue; end
    for k=1:numel(keys)
        if isfield(S,keys{k})
            m=fieldLength(S,keys{k}); if m>n, n=m; end
        end
    end
end
end

function n=fieldLength(S,f)
n=0; if ~isstruct(S)||~isfield(S,f), return; end
x=S.(f);
if iscell(x)||isstring(x), n=numel(x);
elseif isnumeric(x)||islogical(x)
    if isvector(x), n=numel(x); else, n=size(x,1); end
end
end

function tf=rowCompatible(x,nRows,f)
if isempty(x), tf=false; return; end
if iscell(x)||isstring(x), tf=numel(x)==nRows; return; end
if isnumeric(x)||islogical(x)
    if any(strcmp(f,{'Depth_S','Depth_Sstd','BT_RSSI_S','alphaMean_S'}))
        tf=size(x,1)==nRows;
    elseif isvector(x), tf=numel(x)==nRows;
    else, tf=size(x,1)==nRows;
    end
else, tf=true;
end
end

function tf=placeholder(x)
if isempty(x), tf=true;
elseif isnumeric(x), tf=all(~isfinite(x(:)));
elseif isstring(x), tf=all(x=="missing"|ismissing(x));
elseif iscell(x), tf=all(cellfun(@isempty,x(:)));
elseif istable(x), tf=height(x)==0;
else, tf=false;
end
end

function [Summary,Fields] = checkM9AcousticsAndBackscatter(BTamp,vbRequired)
%CHECKM9ACOUSTICSANDBACKSCATTER Check saved 3 MHz and optional 0.5 MHz data.

groups3 = { ...
    "Acoustic","frequency parameter",{'f3'}; ...
    "Acoustic","wavelength",{'lambda3'}; ...
    "Acoustic","wave number",{'k3'}; ...
    "Acoustic","form function",{'F3'}; ...
    "Acoustic","Kappa",{'Kappa3'}; ...
    "Acoustic","mean form function",{'FMean3'}; ...
    "Acoustic","scattering coefficient",{'Ks_3','Ks2_3'}; ...
    "Backscatter","beam 1 range",{'R1_S'}; ...
    "Backscatter","beam 2 range",{'R2_S'}; ...
    "Backscatter","beam 3 range",{'R3_S'}; ...
    "Backscatter","beam 4 range",{'R4_S'}; ...
    "Backscatter","beam 1 amplitude",{'BT_CNT1_S','BT_RSSI1_S'}; ...
    "Backscatter","beam 2 amplitude",{'BT_CNT2_S','BT_RSSI2_S'}; ...
    "Backscatter","beam 3 amplitude",{'BT_CNT3_S','BT_RSSI3_S'}; ...
    "Backscatter","beam 4 amplitude",{'BT_CNT4_S','BT_RSSI4_S'}; ...
    "Backscatter","beam 1 footprint",{'Af1meanS','Af1mean_S'}; ...
    "Backscatter","beam 2 footprint",{'Af2meanS','Af2mean_S'}; ...
    "Backscatter","beam 3 footprint",{'Af3meanS','Af3mean_S'}; ...
    "Backscatter","beam 4 footprint",{'Af4meanS','Af4mean_S'}; ...
    "Backscatter","3 MHz attenuation",{'alphaMean3000S','alphaMean_S'}; ...
    "Backscatter","source level",{'SL_S'}; ...
    "Backscatter","beam 1 Sv profile",{'SV1_S'}; ...
    "Backscatter","beam 2 Sv profile",{'SV2_S'}; ...
    "Backscatter","beam 3 Sv profile",{'SV3_S'}; ...
    "Backscatter","beam 4 Sv profile",{'SV4_S'}; ...
    "Backscatter","beam 1 mean Sv",{'SVmean1','SVBmean1_S'}; ...
    "Backscatter","beam 2 mean Sv",{'SVmean2','SVBmean2_S'}; ...
    "Backscatter","beam 3 mean Sv",{'SVmean3','SVBmean3_S'}; ...
    "Backscatter","beam 4 mean Sv",{'SVmean4','SVBmean4_S'}; ...
    "Backscatter","spatial Sv variability",{'stdSvSpatial_S'}};

groups05 = { ...
    "Acoustic","frequency parameter",{'f05'}; ...
    "Acoustic","wavelength",{'lambda05'}; ...
    "Acoustic","wave number",{'k05'}; ...
    "Acoustic","form function",{'F05'}; ...
    "Acoustic","Kappa",{'Kappa05'}; ...
    "Acoustic","mean form function",{'FMean05'}; ...
    "Acoustic","scattering coefficient",{'Ks_05','Ks2_05'}; ...
    "Backscatter","vertical range/depth",{'RNG_VB_Sm'}; ...
    "Backscatter","vertical amplitude",{'VB_CNT_Sm','VB_RSSI_Sm'}; ...
    "Backscatter","vertical footprint",{'AfVBmean_S'}; ...
    "Backscatter","0.5 MHz water absorption",{'abswat05','AW_05'}; ...
    "Backscatter","vertical spreading",{'SpreadVBmean'}; ...
    "Backscatter","vertical transmission loss",{'TLVBMHz_fixS'}; ...
    "Backscatter","vertical Sv profile", ...
        {'SV_VB','SV_VBrssi','SV_VBS','SV_VBSrssi'}; ...
    "Backscatter","vertical mean Sv",{'SV_VBmeanS','SV_VBmeanSrssi'}; ...
    "Backscatter","vertical readiness",{'VerticalBackscatterReady_S'}};

Fields=[evaluateFieldGroups(BTamp,"3 MHz slanted",groups3,true); ...
        evaluateFieldGroups(BTamp,"0.5 MHz vertical",groups05,vbRequired)];

[f3,src3]=resolveFrequencyMHz(BTamp,'SlantedFrequencyMHz','f3');
[f05,src05]=resolveFrequencyMHz(BTamp,'VerticalFrequencyMHz','f05');
frequency=[f3;f05]; expected=[3;0.5];
frequencyMatch=isfinite(frequency) & ...
    abs(frequency-expected)<=max(0.01,0.02.*expected);
required=[true;logical(vbRequired)];

acousticReady=false(2,1); backscatterReady=false(2,1);
systems=["3 MHz slanted";"0.5 MHz vertical"];
for k=1:2
    rows=Fields.BeamSystem==systems(k);
    acousticRows=rows & Fields.Category=="Acoustic";
    backscatterRows=rows & Fields.Category=="Backscatter";
    acousticReady(k)=all(Fields.HasUsableValue(acousticRows));
    backscatterReady(k)=all(Fields.HasUsableValue(backscatterRows));
end
overall=~required | (frequencyMatch & acousticReady & backscatterReady);
status=strings(2,1);
for k=1:2
    if ~required(k)
        status(k)="not_required_vertical_beam_unavailable";
    elseif overall(k)
        status(k)="ready";
    elseif ~frequencyMatch(k)
        status(k)="frequency_missing_or_incorrect";
    elseif ~acousticReady(k) && ~backscatterReady(k)
        status(k)="acoustics_and_backscatter_incomplete";
    elseif ~acousticReady(k)
        status(k)="acoustics_incomplete";
    else
        status(k)="backscatter_incomplete";
    end
end

BeamSystem=systems;
ExpectedFrequencyMHz=expected;
ObservedFrequencyMHz=frequency;
FrequencySource=[src3;src05];
FrequencyMatches=frequencyMatch;
Required=required;
AcousticParametersSaved=acousticReady;
BackscatterSaved=backscatterReady;
OverallReady=overall;
Status=status;
Summary=table(BeamSystem,ExpectedFrequencyMHz,ObservedFrequencyMHz, ...
    FrequencySource,FrequencyMatches,Required,AcousticParametersSaved, ...
    BackscatterSaved,OverallReady,Status);
end

function T=evaluateFieldGroups(S,beamSystem,groups,isRequired)
n=size(groups,1);
BeamSystem=repmat(string(beamSystem),n,1);
Category=strings(n,1); Check=strings(n,1); Required=false(n,1);
Candidates=strings(n,1); FoundField=strings(n,1);
Exists=false(n,1); HasUsableValue=false(n,1); SizeText=strings(n,1);
for k=1:n
    Category(k)=string(groups{k,1});
    Check(k)=string(groups{k,2});
    candidateNames=groups{k,3};
    Candidates(k)=strjoin(string(candidateNames),' | ');
    Required(k)=logical(isRequired);
    for j=1:numel(candidateNames)
        name=candidateNames{j};
        if isstruct(S)&&isfield(S,name)
            if FoundField(k)==""
                FoundField(k)=string(name);
                Exists(k)=true;
                SizeText(k)=strjoin(string(size(S.(name))),'x');
            end
            if hasUsableValue(S.(name))
                FoundField(k)=string(name);
                Exists(k)=true;
                HasUsableValue(k)=true;
                SizeText(k)=strjoin(string(size(S.(name))),'x');
                break
            end
        end
    end
    if FoundField(k)=="", SizeText(k)="missing"; end
end
T=table(BeamSystem,Category,Check,Required,Candidates,FoundField, ...
    Exists,HasUsableValue,SizeText);
end

function tf=hasUsableValue(x)
if isempty(x)
    tf=false;
elseif isnumeric(x)
    tf=any(isfinite(double(x(:))));
elseif islogical(x)
    tf=any(x(:));
elseif isstring(x)
    t=strtrim(x(:));
    tf=any(~ismissing(t) & t~="" & lower(t)~="missing" & lower(t)~="nan");
elseif ischar(x)
    tf=~isempty(strtrim(x));
elseif iscell(x)
    tf=false;
    for k=1:numel(x)
        if hasUsableValue(x{k}), tf=true; return; end
    end
elseif istable(x)
    tf=height(x)>0;
elseif isstruct(x)
    tf=~isempty(fieldnames(x));
else
    tf=true;
end
end

function [frequencyMHz,source]=resolveFrequencyMHz(S,metadataField,rawField)
frequencyMHz=NaN; source="missing";
candidates={ ...
    {metadataField}, ...
    {'M9meta',metadataField}, ...
    {rawField}, ...
    {'Options',rawField}};
for k=1:numel(candidates)
    [found,value]=nestedFieldValue(S,candidates{k});
    if ~found || ~isnumeric(value), continue; end
    value=double(value(:)); value=value(isfinite(value));
    if isempty(value), continue; end
    raw=median(value);
    if raw>=1e5
        frequencyMHz=raw/1e6;
    elseif raw>=100
        frequencyMHz=raw/1000;
    else
        frequencyMHz=raw;
    end
    source=strjoin(string(candidates{k}),'.');
    return
end
end

function [found,value]=nestedFieldValue(S,path)
found=false; value=[]; current=S;
for k=1:numel(path)
    name=path{k};
    if ~isstruct(current)||~isscalar(current)||~isfield(current,name)
        return
    end
    current=current.(name);
end
found=true; value=current;
end

function y=logicalRows(S,fieldName,nRows,defaultValue)
y=repmat(logical(defaultValue),nRows,1);
if ~isstruct(S)||~isfield(S,fieldName), return; end
x=S.(fieldName);
if islogical(x)||isnumeric(x)
    x=logical(x(:));
    y(1:min(nRows,numel(x)))=x(1:min(nRows,numel(x)));
end
end

function y=finiteRows(S,fieldName,nRows)
y=false(nRows,1);
if ~isstruct(S)||~isfield(S,fieldName), return; end
x=S.(fieldName);
if isnumeric(x)
    if isvector(x)
        v=isfinite(x(:));
    elseif size(x,1)==nRows
        v=any(isfinite(x),2);
    elseif size(x,2)==nRows
        v=any(isfinite(x),1).';
    else
        v=false(0,1);
    end
    y(1:min(nRows,numel(v)))=v(1:min(nRows,numel(v)));
elseif iscell(x)
    for k=1:min(nRows,numel(x))
        item=x{k};
        y(k)=isnumeric(item)&&any(isfinite(double(item(:))));
    end
end
end

function namesFound=collectLegacyInstrumentNames(value,path)
namesFound=strings(0,1);
if isstruct(value)
    fields=fieldnames(value);
    for k=1:numel(fields)
        thisName=fields{k};
        if isLegacyInstrumentName(thisName)
            namesFound(end+1,1)=path+"."+string(thisName); %#ok<AGROW>
        end
        for n=1:numel(value)
            if numel(value)==1
                childPath=path+"."+string(thisName);
            else
                childPath=path+"("+string(n)+")."+string(thisName);
            end
            namesFound=[namesFound;collectLegacyInstrumentNames( ...
                value(n).(thisName),childPath)]; %#ok<AGROW>
        end
    end
elseif istable(value)||isa(value,'timetable')
    variables=value.Properties.VariableNames;
    for k=1:numel(variables)
        thisName=variables{k};
        if isLegacyInstrumentName(thisName)
            namesFound(end+1,1)=path+"."+string(thisName); %#ok<AGROW>
        end
        namesFound=[namesFound;collectLegacyInstrumentNames( ...
            value.(thisName), ...
            path+"."+string(thisName))]; %#ok<AGROW>
    end
elseif iscell(value)
    for n=1:numel(value)
        namesFound=[namesFound;collectLegacyInstrumentNames(value{n}, ...
            path+"{"+string(n)+"}")]; %#ok<AGROW>
    end
end
end

function tf=isLegacyInstrumentName(name)
% Do not match 12 when it is part of a larger number such as 0125.
tf=contains(name,'1200') || ...
    ~isempty(regexp(name,'(?<!\d)12(?!\d)','once'));
end
