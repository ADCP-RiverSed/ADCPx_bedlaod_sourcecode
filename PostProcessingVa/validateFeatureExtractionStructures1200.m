function ValidationReport = validateFeatureExtractionStructures1200( ...
    tratesFile, shearFile, backscatterFile, varargin)
%VALIDATEFEATUREEXTRACTIONSTRUCTURES1200 Validate saved MAT files against
% the exact structures and fields used by the supplied feature-extraction
% script.
%
% ValidationReport = validateFeatureExtractionStructures1200( ...
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
% the equality of BTamp and DataBTamp. It does not modify the files.

p=inputParser;
addRequired(p,'tratesFile',@(x)ischar(x)||isstring(x));
addRequired(p,'shearFile',@(x)ischar(x)||isstring(x));
addRequired(p,'backscatterFile',@(x)ischar(x)||isstring(x));
addParameter(p,'ThrowOnMissing',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'ThrowOnPlaceholder',false,@(x)islogical(x)||isnumeric(x));
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
    'Depth_1200','VmodAll_1200', ...
    'BT1_D_1200std','BT2_D_1200std','BT3_D_1200std','BT4_D_1200std', ...
    'Tau_All1200_1','Tau_All1200_2','Ushear1_1200','Ushear2_1200', ...
    'UsheaAll_1200','Cb2_1200','nCb_new1200','ucr_1200','Tucr_1200', ...
    'TKE_prof1200','tau_rm1200','tau_Mann_star1200','W_i1200', ...
    'W_iVRshear2_1200','Tau_shear2VR_1200','Dv1_1200','Tv1_1200', ...
    'HEightDUNE1200','Vmod1low_1200stdMOD','Vmod_1200stdMOD', ...
    'WT1200','u_shearH1200'};
featureBTvelFields = { ...
    'BTmodf_1200','BTstd1200a','BTerror1200m', ...
    'FilteredDataN_1200','FilteredDataE_1200', ...
    'BT_stdf_1200','BTstd_1200','DeltaT1200'};
featureBTampFields = { ...
    'Ks2_12','D50_05','Depth1200','Depth1200std', ...
    'Af1mean1200','Af2mean1200','Af3mean1200','Af4mean1200', ...
    'Spread1_1200','Spread2_1200','Spread3_1200','Spread4_1200', ...
    'SVBmean1_1200','SVBmean2_1200','SVBmean3_1200','SVBmean4_1200', ...
    'SVB1mean1_1200','SVB1mean2_1200','SVB1mean3_1200','SVB1mean4_1200', ...
    'SVBstd1_1200','SVBstd2_1200','SVBstd3_1200','SVBstd4_1200', ...
    'SV1_1200','SV2_1200','SV3_1200','SV4_1200', ...
    'SVB1std1_1200','SVB1std2_1200','SVB1std3_1200','SVB1std4_1200', ...
    'BT_RSSI_1200','alphaMean1200','stdSpatialRSSI_1200', ...
    'stdSvSpatial_1200','alpha1mean1200','alpha2mean1200', ...
    'alpha3mean1200','alpha4mean1200'};

nRows=inferRows(Data_field,ShearData,BTvel,BTamp);
Fields=[checkFields(Data_field,'Data_field',featureDataFields,nRows); ...
        checkFields(ShearData,'ShearData',featureShearFields,nRows); ...
        checkFields(BTvel,'BTvel',featureBTvelFields,nRows); ...
        checkFields(BTamp,'BTamp',featureBTampFields,nRows)];

BTampAliasEqual=isstruct(BTamp)&&isstruct(DataBTamp)&&isequaln(BTamp,DataBTamp);
RowSummary=table(nRows,fieldLength(Data_field,'MeanTrate'), ...
    fieldLength(ShearData,'Depth_1200'),fieldLength(BTvel,'BTmodf_1200'), ...
    fieldLength(BTamp,'SVBmean1_1200'), ...
    'VariableNames',{'ExpectedRows','DataRows','ShearRows','BTvelRows','BackscatterRows'});

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
ValidationReport.MissingTopLevel=TopVariable(~TopExists);
ValidationReport.MissingFields=Fields.StructureField(~Fields.Exists);
ValidationReport.PlaceholderFields=Fields.StructureField(Fields.IsPlaceholder);
ValidationReport.RowIncompatibleFields=Fields.StructureField(~Fields.RowCompatible);
ValidationReport.CriticalUnresolved=intersect(ValidationReport.PlaceholderFields,[ ...
    "BTvel.DeltaT1200";"BTamp.Ks2_12"; ...
    "ShearData.Tau_shear2VR_1200"; ...
    "ShearData.Vmod_1200stdMOD"; ...
    "ShearData.Vmod1low_1200stdMOD"]);

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
    fprintf('\n=== Feature-extraction structure validation ===\n');
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
keys={'MeanTrate','Depth_1200','BTmodf_1200','SVBmean1_1200'};
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
    if any(strcmp(f,{'Depth1200','Depth1200std','BT_RSSI_1200','alphaMean1200'}))
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
