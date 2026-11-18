function [Data_field, ShearData, BTvel, Q, BTamp, StructureReport] = ...
    buildProcessingStructuresM9(ADCPin, WaterProfileOut, DepthRNGout, ...
    ShearStressOut, BTvelOut, BedloadTransportOut, BackscatterOut, varargin)
%BUILDPROCESSINGSTRUCTURESM9 Build legacy structures using M9 _S fields.
%
% [Data_field, ShearData, BTvel, Q, BTamp, StructureReport] = ...
%     buildProcessingStructuresM9(ADCPin, WaterProfileOut, DepthRNGout, ...
%     ShearStressOut, BTvelOut, BedloadTransportOut, BackscatterOut, ...
%     'Name',Value,...)
%
% PURPOSE
% -------
% This is the final packaging step after the modular SonTek M9 processing
% functions. It creates the structures used by the historical multi-river
% analysis scripts while retaining their M9-compatible legacy field names:
%
%   Data_field   field samples, grain sizes and measured transport
%   ShearData    water profile, shear, concentration and thickness results
%   BTvel        apparent bedload-velocity results and exact filter masks
%   Q            calculated transport-rate results
%   BTamp        bottom-backscatter and acoustic results
%
% The function does not recalculate the physical results. It copies values
% from the modern output structures, creates clear legacy aliases where the
% equivalence is known, and creates a placeholder for every requested legacy
% field that is unavailable. Missing fields are listed in StructureReport.
%
% INPUTS
% ------
% ADCPin              input-table structure, normally PhSamplesin
% WaterProfileOut     output from processWaterProfile, or []
% DepthRNGout         output from processDepthsRNG, or []
% ShearStressOut      output from processShearStress1200, or []
% BTvelOut            output from processBTbedloadVelocity1200, or []
% BedloadTransportOut output from processBedloadTrate, or []
% BackscatterOut      output from processBackscatterM9, or []
% 'Acoustic'           output from computeADCPacousticParameters. Every
%                      field is copied directly into BTamp/DataBTamp.
%
% IMPORTANT BEHAVIOUR
% -------------------
% 1. All current fields from each modern output are retained in the most
%    relevant legacy structure. New final Cb, thickness, transport, exact
%    BT masks and backscatter diagnostics are therefore not lost.
% 2. Every historical field explicitly listed in the old saving script is
%    present even when it cannot be resolved. Numeric row variables become
%    NaN vectors, time-series containers become cell arrays, status fields
%    become "missing", and missing tables become empty tables.
% 3. Folder-specific Colin/Massimo/Rheine branches are replaced by a union
%    of all fields. This prevents a field from disappearing merely because
%    the folder name changed.
% 4. DataBTamp is created as a compatibility alias of BTamp when structures
%    are assigned or saved, because older analysis scripts use both names.
%
% OPTIONS
% -------
% 'FolderSontek'        optional folder text retained in StructureReport
% 'FolderM9'            source folder used to extract the YYYYMMDD prefix
% 'FolderTrates'        destination folder for Data_field, BTvel and Q
% 'FolderShear'         destination folder for ShearData
% 'FolderBS'            destination folder for BTamp and DataBTamp
% 'FileDatePrefix'      optional YYYYMMDD override; normally extracted from
%                       FolderM9, FolderSontek or source filenames
% 'LocationName'        optional location name inserted after the date and
%                       underscore, e.g. 20170322_Tomorrica_TratesData.mat
% 'dataADCPS'          raw dir-style file list to store in Data_field and
%                       ShearData
% 'dataADCPSboat'      optional boat-file list
% 'Acoustic'            output structure from computeADCPacousticParameters.
%                       Acoustic fields are flattened directly into BTamp;
%                       no nested Acoustic structure is saved.
% 'ExtraVariables'      backward-compatible supplemental-variable structure.
%                       Its fields are routed into Data_field, ShearData,
%                       BTvel, Q or BTamp; it is never saved as a separate
%                       output structure. Direct name-value variables are
%                       preferred, for example 'DeltaT_S',DeltaT_S.
% 'FlipPSDFields'       apply the legacy flipud to PSD_percentage and
%                       QperDall; default true
% 'UseCallerWorkspace'  search the caller workspace for unresolved legacy
%                       variables; default true
% 'AssignToCaller'      assign Data_field, ShearData, BTvel, Q, BTamp,
%                       DataBTamp and StructureReport; default true
% 'AssignToBase'        assign the same structures to base; default false
% 'SaveDataFile'        optional explicit MAT file containing Data_field,
%                       BTvel and Q; date/location prefix is added automatically
% 'SaveShearFile'       optional explicit MAT file containing ShearData only;
%                       date/location prefix is added automatically
% 'SaveBackscatterFile' optional explicit MAT file containing BTamp and
%                       DataBTamp; date/location prefix is added automatically
% 'SaveCombinedFile'    MAT file containing all structures and report
% 'SaveVersion'         '-v7.3' (default) or another MATLAB save version
% 'EnforceFeatureReaderFields' ensure every field referenced by the supplied
%                       feature-extraction script exists; default true
% 'FailOnMissingRequiredValues' error if any required reader field remains
%                       a placeholder; default false
% 'FailOn1200FieldNames' backward-compatible option that errors before
%                       saving if any output field name still contains
%                       1200 or an isolated 12 token; default true
% 'WarnOnPhysicalSampleNaN' issue warnings when measured transport or
%                       grain-size fields contain NaN/missing values;
%                       default true
% 'PhysicalSampleFields' optional list of physical-sample fields to check.
%                       The default includes measured transport statistics,
%                       D10/D16/D50/D84/D90 variants, mobile/grab grain sizes,
%                       DmeaNUMpsd and SSC.
% Any unrecognized name-value pair is treated as a supplemental variable
% and routed to its existing destination structure by its field name. No
% supplemental structure is created or saved.
% 'Verbose'             print a summary; default true
%
% RECOMMENDED CALL
% ----------------
% [Data_field,ShearData,BTvel,Q,BTamp,StructureReport] = ...
%     buildProcessingStructuresM9( ...
%         PhSamplesin,WaterProfileOut,DepthRNGout,ShearStressOut, ...
%         BTvelOut,BedloadTransportOut,BackscatterOut, ...
%         'FolderSontek',FolderSontek, ...
%         'FolderM9',FolderM9, ...
%         'FolderTrates',FolderTrates, ...
%         'FolderShear',FolderShear, ...
%         'FolderBS',FolderBS, ...
%         'LocationName','Tomorrica', ...
%         'dataADCPS',dataADCPS, ...
%         'Acoustic',Acoustic, ...
%         'UseCallerWorkspace',true, ...
%         'AssignToCaller',true);
%
% With the three destination folders supplied, the function saves:
%   FolderTrates/YYYYMMDD_LOCATION_TratesData.mat  -> Data_field, BTvel, Q
%   FolderShear/YYYYMMDD_LOCATION_ShearData.mat    -> ShearData
%   FolderBS/YYYYMMDD_LOCATION_BackscatterData.mat -> BTamp, DataBTamp
% If LocationName is empty, the legacy YYYYMMDD_* names are used.
%
% Field names containing 1200 or an isolated abbreviated 12 instrument
% token are normalized to _S. A 12 embedded in a larger number, such as
% QPer0125D, is not changed. All stored values are retained unchanged.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;
p.KeepUnmatched = true;  % direct supplemental variables, routed later

addRequired(p,'ADCPin',@(x) isempty(x) || isstruct(x));
addRequired(p,'WaterProfileOut',@(x) isempty(x) || isstruct(x));
addRequired(p,'DepthRNGout',@(x) isempty(x) || isstruct(x));
addRequired(p,'ShearStressOut',@(x) isempty(x) || isstruct(x));
addRequired(p,'BTvelOut',@(x) isempty(x) || isstruct(x));
addRequired(p,'BedloadTransportOut',@(x) isempty(x) || isstruct(x));
addRequired(p,'BackscatterOut',@(x) isempty(x) || isstruct(x));

addParameter(p,'FolderSontek',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderM9',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderTrates',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderShear',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderBS',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FileDatePrefix',"",@(x) ischar(x) || isstring(x));
addParameter(p,'LocationName',"",@(x) ischar(x) || isstring(x));
addParameter(p,'dataADCPS',[],@(x) true);
addParameter(p,'dataADCPSboat',[],@(x) true);
addParameter(p,'Acoustic',struct(),@(x) isempty(x) || isstruct(x));
addParameter(p,'ExtraVariables',struct(),@(x) isempty(x) || isstruct(x));
addParameter(p,'FlipPSDFields',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'UseCallerWorkspace',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToCaller',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'AssignToBase',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'SaveDataFile',"",@(x) ischar(x) || isstring(x));
addParameter(p,'SaveShearFile',"",@(x) ischar(x) || isstring(x));
addParameter(p,'SaveBackscatterFile',"",@(x) ischar(x) || isstring(x));
addParameter(p,'SaveCombinedFile',"",@(x) ischar(x) || isstring(x));
addParameter(p,'SaveVersion','-v7.3',@(x) ischar(x) || isstring(x));
addParameter(p,'EnforceFeatureReaderFields',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'FailOnMissingRequiredValues',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'FailOn1200FieldNames',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'WarnOnPhysicalSampleNaN',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'PhysicalSampleFields',defaultPhysicalSampleFieldNames(), ...
    @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(p,'Verbose',true,@(x) islogical(x) || isnumeric(x));

parse(p,ADCPin,WaterProfileOut,DepthRNGout,ShearStressOut, ...
    BTvelOut,BedloadTransportOut,BackscatterOut,varargin{:});
opt = p.Results;
DirectVariables = scalarStructOrEmpty(p.Unmatched);

ADCPin = scalarStructOrEmpty(ADCPin);
WaterProfileOut = scalarStructOrEmpty(WaterProfileOut);
DepthRNGout = scalarStructOrEmpty(DepthRNGout);
ShearStressOut = scalarStructOrEmpty(ShearStressOut);
BTvelOut = scalarStructOrEmpty(BTvelOut);
BedloadTransportOut = scalarStructOrEmpty(BedloadTransportOut);
BackscatterOut = scalarStructOrEmpty(BackscatterOut);
Acoustic = scalarStructOrEmpty(opt.Acoustic);
ExtraVariables = scalarStructOrEmpty(opt.ExtraVariables);

% Normalize any RDI-style names that remain in supplied structures. Values
% are not changed. If an _S target already exists, a finite/non-placeholder
% value is preferred and the duplicate 1200/12 field is removed.
InputFieldRenaming = emptyFieldRenameReport();
[ADCPin,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    ADCPin,"ADCPin","ADCPin",InputFieldRenaming);
[WaterProfileOut,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    WaterProfileOut,"WaterProfileOut","WaterProfileOut",InputFieldRenaming);
[DepthRNGout,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    DepthRNGout,"DepthRNGout","DepthRNGout",InputFieldRenaming);
[ShearStressOut,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    ShearStressOut,"ShearStressOut","ShearStressOut",InputFieldRenaming);
[BTvelOut,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    BTvelOut,"BTvelOut","BTvelOut",InputFieldRenaming);
[BedloadTransportOut,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    BedloadTransportOut,"BedloadTransportOut", ...
    "BedloadTransportOut",InputFieldRenaming);
[BackscatterOut,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    BackscatterOut,"BackscatterOut","BackscatterOut",InputFieldRenaming);
[Acoustic,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    Acoustic,"Acoustic","Acoustic",InputFieldRenaming);
[ExtraVariables,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    ExtraVariables,"ExtraVariables","ExtraVariables",InputFieldRenaming);
[DirectVariables,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    DirectVariables,"DirectVariables","DirectVariables",InputFieldRenaming);
SupplementalVariables = mergeStructs(ExtraVariables,DirectVariables);

% The acoustic calculation normally returns [Acoustic,ADCPin]. Prefer the
% explicit Acoustic input. For backwards compatibility, recover a caller
% variable named Acoustic or copy recognized acoustic fields from the
% updated ADCPin. Fields are later flattened directly into BTamp.
acousticFieldCandidates = defaultAcousticFieldNames();
if isempty(fieldnames(Acoustic)) && logical(opt.UseCallerWorkspace)
    try
        if evalin('caller','exist(''Acoustic'',''var'')')
            Acoustic = scalarStructOrEmpty(evalin('caller','Acoustic'));
        end
    catch
        % Explicit input and ADCPin fallback remain available.
    end
end
% Always merge recognized acoustic fields carried by the updated ADCPin.
% The explicit Acoustic output takes precedence when the same field exists.
AcousticFromADCPin = selectFieldsPresent(ADCPin,acousticFieldCandidates);
Acoustic = mergeStructs(AcousticFromADCPin,Acoustic);
[Acoustic,InputFieldRenaming] = normalizeLegacyInstrumentFieldsToS( ...
    Acoustic,"Acoustic","Acoustic",InputFieldRenaming);

nRows = inferRowCount({ADCPin,WaterProfileOut,DepthRNGout,ShearStressOut, ...
    BTvelOut,BedloadTransportOut,BackscatterOut,Acoustic,SupplementalVariables});

%% Historical field inventories
% These are the union of the fields in the supplied legacy saving script.
dataLegacyFields = { ...
    'D65_4','D65_14','D50_05','D50_14','D90_05','D90_14', ...
    'D84_05','D84_14','D16_05','D16_14', ...
    'Q50_05','Q50_14','Q90_05','Q90_14','Q84_05','Q84_14', ...
    'Q16_05','Q16_14','QperPos','Q50s','Q50g','Q_std_14','Q_std_05', ...
    'FsS','FsG','D50s_05','D50g_05','D50s_14','D50g_14', ...
    'QPer8D','QPer6D','QPer4D','QPer2D','QPer1D','QPer063D', ...
    'QPer05D','QPer025D','QPer0125D','QPer00635D', ...
    'PSD_percentage','QperDall', ...
    'Fs8D','Fs6D','Fs4D','Fs2D','Fs1D','Fs063D','Fs05D', ...
    'Fs025D','Fs0125D','Fs00635D', ...
    'Dmass14','Dmass05','D_g_mean_05','D_g_std_05', ...
    'D_g_mean_14','D_g_std_14','DnumMean05','DnumMean14', ...
    'Rep','vs','D_','dataADCPS','dataADCPSboat', ...
    'D_50grabM','D_50grabSTD','D_50grabMAX','D_50grabMIN', ...
    'D_90grabM','D_90grabSTD','D_90grabMAX','D_90grabMIN', ...
    'D_10grabM','D_10grabSTD','D_10grabMAX','D_10grabMIN', ...
    'D_16grabM','D_16grabSTD','D_16grabMAX','D_16grabMIN', ...
    'D_84grabM','D_84grabSTD','D_84grabMAX','D_84grabMIN', ...
    'D_mgrabM','D_mgrabSTD','D_mgrabMAX','D_mgrabMIN', ...
    'Max05','Min05','MeanTrate05','Max14','Min14','MeanTrate14', ...
    'MMaxTrate','MMinTrate','MMeanTrate', ...
    'MMax05','MMin05','MMeanTrate05', ...
    'MeanMassTrateGravel_05','MeanMassTrateSand_05', ...
    'MMax14','MMin14','MMeanTrate14', ...
    'MeanMassTrateGravel_14','MeanMassTrateSand_14', ...
    'LeftBankDist14','ks','Position2ADCPs','Position1ADCPs', ...
    'MaxTrate','MinTrate','MeanTrate','Temp1','Dstd_14','Dstd_05', ...
    'D84','D16','D50','D90'};

shearLegacyFields = { ...
    'WT_S','dataADCPS','VmodAll_S','Vmod_SstdMOD', ...
    'Vmod1low_SstdMOD','u_shearH_S','Tau_All_S_H', ...
    'Tau_All_S','CbH_S','TvH_S','DvH_S', ...
    'DthicknessH_S','QvRH_S','ucr_S','Fr_S', ...
    'Tausk_S','Tucr_S','Ushear1_S','Tau_All_S_1', ...
    'Cb1_S','Tv1_S','Dv1_S','Tv2_S','Dv2_S', ...
    'DDt_S','Ushear2_S','Tau_All_S_2','Cb2_S', ...
    'nCb_new_S','UsheaAll_S','Cb_S','HEightDUNE_S', ...
    'Depth_S','BT1_D_Sstd','BT2_D_Sstd', ...
    'BT3_D_Sstd','BT4_D_Sstd', ...
    'U_shearRE1_S','U_shearRE2_S', ...
    'Tau_AllRE1_S','Tau_AllRE2_S','WT_S', ...
    'tau_Mann_S','Tau_shear2VR_S','tau_Mann_star_S', ...
    'tau_rm_S','uShear_Mann_S','tau_rm_starGaumann_S', ...
    'tau_rm_star_S','W_star_perDall_S','tau_Mann_SvsW', ...
    'tau_Mann_SvsWstar','W_i_S','W_iVRshear2_S', ...
    'tauRef_S','tauRefVRs_S','tauRef_SVRstar', ...
    'tauRef_Sstar','tauRef_SPerfrac', ...
    'tau_turb_prof_S','TKE_prof_S','uw_mean_S'};

btLegacyFields = { ...
    'BTE_S','BTN_S','BT_x_S','BTx_Sf','BTx_S', ...
    'BT_y_S','BTy_Sf','BTy_S','BTx_Sa','BTy_Sa', ...
    'BTerror_Sm','BTerror_Sall','BTstd_Sa','BTmod_Sa', ...
    'BTmod_SaW','BTstd_SaW','BTmodM1_S','BTmod_S', ...
    'BTmodf_S','FilteredDataN_S','FilteredDataE_S', ...
    'NANsbeforeX_Sperc','NANsbeforeY_Sperc','BTstd_S', ...
    'BT_stdf_S','ErrorVel_Svalid','ErrorV_S', ...
    'alphaV_S','alphaBTf_S','VmodALL_S','BTmodf_SALL', ...
    'Diff_S','maxDiff_S','meanDiff_S','MeanDepth_S', ...
    'beamDiff_S','meanDiffbeams_S','NNinhomogenious_S', ...
    'MesString_S','CV1_S','CV2_S', ...
    'VelHomString1_S','VelHomString2_S','DeltaT_S'};

qLegacyFields = { ...
    'qsNew_S','qsNew_SstdBT','qs_S','qs_SstdBT', ...
    'qs1_S','qs1_SstdBT','qs_SLambda', ...
    'qs_SstdBTLambda','qs2_SstdBT','qs2_S', ...
    'qsH_SstdBT','qsH_S','qs1_SDD', ...
    'qsMethod2regression_S'};

btampLegacyFields = { ...
    'Depth_S','stdSpatial_S','Depth_Sstd','stdSpatialALL_S', ...
    'BT_RSSI_S','alphaMean_S', ...
    'alpha1mean_Sgrad','alpha2mean_Sgrad', ...
    'alpha3mean_Sgrad','alpha4mean_Sgrad', ...
    'alpha1mean_S','alpha2mean_S','alpha3mean_S','alpha4mean_S', ...
    'dEDR1_S','dEDR2_S','dEDR3_S','dEDR4_S', ...
    'R1_S','R2_S','R3_S','R4_S', ...
    'Spread1_S','Spread2_S','Spread3_S','Spread4_S', ...
    'SpreadM_S','SVBmean1_S','SVBmean2_S', ...
    'SVBmean3_S','SVBmean4_S', ...
    'SVB1mean1_S','SVB1mean2_S','SVB1mean3_S','SVB1mean4_S', ...
    'SVBstd1_S','SVBstd2_S','SVBstd3_S','SVBstd4_S', ...
    'SVB1std1_S','SVB1std2_S','SVB1std3_S','SVB1std4_S', ...
    'SV11_S','SV22_S','SV33_S','SV44_S', ...
    'SV1_S','SV2_S','SV3_S','SV4_S', ...
    'stdSvSpatial_S','stdSvSpatialALL_S', ...
    'BT_RSSI1_S','BT_RSSI2_S','BT_RSSI3_S','BT_RSSI4_S', ...
    'stdSpatialRSSI_S','TL1_S','TL2_S','TL3_S','TL4_S', ...
    'TL1_fix_S','TL2_fix_S','TL3_fix_S','TL4_fix_S', ...
    'Dmax_S_1','Dmax_S_2','Dmax_S_3','Dmax_S_4', ...
    'Area_S_1','Area_S_2','Area_S_3','Area_S_4', ...
    'Af1mean_S','Af2mean_S','Af3mean_S','Af4mean_S', ...
    'distance_S_1','distance_S_2','distance_S_3','distance_S_4', ...
    'Lp_S_1','Lp_S_2','Lp_S_3','Lp_S_4', ...
    'Pdwb_S','SL_S','STDPowerdb_S','STDPW_S', ...
    'XmitVV_S','XmitCC_S','PitchH_S','RollL_S', ...
    'Ks2_S','D50_05','D50_14'};


%% Exact contract used by the supplied feature-extraction/loading script
% These inventories are intentionally explicit. The three saved files must
% contain:
%   Trates file:      Data_field, BTvel, Q
%   Shear file:       ShearData
%   Backscatter file: BTamp, DataBTamp
%
% Every field below is accessed directly by the reader script. Missing
% quantities are still created as placeholders so load/index operations do
% not fail. StructureReport.ReaderCompatibility identifies placeholders.
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
    'Ks2_S','D50_05', ...
    'Depth_S','Depth_Sstd', ...
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

% Explicitly retained modern variables that are essential for interpreting
% the final methods even when an older version of an output structure is used.
newShearFields = { ...
    'D50critical_mm','D50critical_m','D50mobile','Rep_S','Dstar_S', ...
    'ThetaCr_S','ThetaCrVanRijn_S','ThetaCrBrownlie_S', ...
    'ThetaCrSoulsby_S','CriticalInputsConsistent_S', ...
    'CbCandidates_S','CbMethodNames_S','CbFinal_S', ...
    'CbFinalMethodCount_S','CbFinalSource_S','CbFinalUsedFallback_S', ...
    'DthicknessCandidates_S','DthicknessMethodNames_S', ...
    'DthicknessFinal_S','DthicknessFinalMethodCount_S', ...
    'DthicknessFinalSource_S','DthicknessFinalUsedFallback_S', ...
    'DthicknessKin_S','Ninstrument','uShear_skin_S', ...
    'theta_skin_S','tau_skin_S'};

newBTFields = { ...
    'BTfilterMask_S','BTvalidMask_S','BTfilterMaskLength_S', ...
    'NfilteredFinal_S','NvalidFinal_S','BTfilterMaskConvention', ...
    'BTfilterMaskStart_S','BTfilterMaskEnd_S', ...
    'BTx_final_S','BTy_final_S','BTmod_final_S','va_final_S'};

newQFields = { ...
    'qsKinematicFinal_S','qsFinalKinematic_S','qsFinal_S', ...
    'qs_final_S','qsKinematicFinalValid_S', ...
    'qsRE1_S','qsRE2_S','qsMann_S','qsSkin_S', ...
    'qsVR_S','qsVRRE1_S','qsVRRE2_S','qsVRH_S', ...
    'qsVR1_S','qsVR2_S','qsVRMann_S','qsVRSkin_S', ...
    'QvR_S','QvRRE1_S','QvRRE2_S','QvRH_S', ...
    'QvR1_S','QvR2_S','QvRMann_S','QvRSkin_S', ...
    'CbFinal_S','CbFinalSource_S','DthicknessFinal_S', ...
    'DthicknessFinalSource_S','BTmodf_S','Ninstrument', ...
    'FinalTrateComparisonTable'};

newBTampFields = { ...
    'S_S','qsMethod2regression_S','Method2Status_S', ...
    'BTvelocityFilterMask_S','BTvelocityMaskMode_S', ...
    'BTvelocityFilterSource_S','BTvelocityFilterExact_S', ...
    'BTvelocityFilterAligned_S','NfilteredBTvelocity_S', ...
    'StartUsed_S','EndUsed_S','NensemblesUsed_S', ...
    'BackscatterTable'};

%% Alias definitions
aliases = struct();

% Data_field aliases and compatibility fallbacks.
aliases.D65_4 = {'D65_4','D65_14','D65'};
aliases.D65_14 = {'D65_14','D65'};
aliases.D50_05 = {'D50_05','D50'};
aliases.D50_14 = {'D50_14','D50'};
aliases.D90_05 = {'D90_05','D90'};
aliases.D90_14 = {'D90_14','D90'};
aliases.D84_05 = {'D84_05','D84'};
aliases.D84_14 = {'D84_14','D84'};
aliases.D16_05 = {'D16_05','D16'};
aliases.D16_14 = {'D16_14','D16'};
aliases.QperPos = {'QperPos','QperPos_14','QperPos_05'};
aliases.Rep = {'Rep','Rep_S'};
aliases.vs = {'vs','ws'};
aliases.LeftBankDist14 = {'LeftBankDist14','Position2ADCPs','DistanceFromLeftBank'};
aliases.Position2ADCPs = {'Position2ADCPs','DistanceFromLeftBank'};
aliases.Position1ADCPs = {'Position1ADCPs'};
aliases.MeanTrate = {'MeanTrate','TrateMean'};
aliases.Dstd_14 = {'Dstd_14','D_g_std_14'};
aliases.Dstd_05 = {'Dstd_05','D_g_std_05'};
aliases.D_50grabMAX = {'D_50grabMAX','D_50grabM'};
aliases.D_90grabMAX = {'D_90grabMAX','D_90grabM'};
aliases.D_84grabM = {'D_84grabM','D84'};
aliases.D_84grabMAX = {'D_84grabMAX','D_84grabM','D84'};
aliases.D_16grabM = {'D_16grabM','D16'};
aliases.D_16grabMAX = {'D_16grabMAX','D_16grabM','D16'};
aliases.D_mgrabM = {'D_mgrabM','D_16grabM','D16'};
aliases.Temp1 = {'Temp1','WT_S','WT_3','WT_Sm9','WT_12'};
aliases.alpha1mean_S = {'alpha1mean_S','alpha1meanS','alpha1mean'};
aliases.alpha2mean_S = {'alpha2mean_S','alpha2meanS','alpha2mean'};
aliases.alpha3mean_S = {'alpha3mean_S','alpha3meanS','alpha3mean'};
aliases.alpha4mean_S = {'alpha4mean_S','alpha4meanS','alpha4mean'};


aliases.MaxTrate = {'MaxTrate','TrateMax','MTrateMax','MaxTrateAll'};
aliases.MinTrate = {'MinTrate','TrateMin','MTrateMin','MinTrateAll'};
aliases.Depth_S = {'Depth_S','DepthS'};
aliases.Depth_Sstd = {'Depth_Sstd','DepthSstd'};
aliases.DeltaT_S = {'DeltaT_S','DeltaT','dt_S'};
aliases.Ks2_S = {'Ks2_S','Ks2_3','Ks_3','ks2_S','Ks_S', ...
    'Ks2_12','ks2_12','Ks_12'};
aliases.Vmod_SstdMOD = {'Vmod_SstdMOD','VmodAll_Sstd'};
aliases.Vmod1low_SstdMOD = {'Vmod1low_SstdMOD'};
% Final concentration aliases. CbFinal_S is the retained canonical name
% in ShearData and Q; older capitalization variants are accepted as input.
aliases.CbFinal_S = {'CbFinal_S','Cbfinal_S','CbFinal','Cbfinal'};
aliases.CbFinalSource_S = {'CbFinalSource_S','CbfinalSource_S', ...
    'CbFinalSource','CbfinalSource'};

% ShearData aliases.
aliases.WT_S = {'WT_S','WT_3','WT_Sm9','Temp1','WT_12'};
aliases.BT1_D_Sstd = {'BT1_D_Sstd','RNG11std_S','RNG11_Sstd'};
aliases.BT2_D_Sstd = {'BT2_D_Sstd','RNG22std_S','RNG22_Sstd'};
aliases.BT3_D_Sstd = {'BT3_D_Sstd','RNG33std_S','RNG33_Sstd'};
aliases.BT4_D_Sstd = {'BT4_D_Sstd','RNG44std_S','RNG44_Sstd'};

% BTvel aliases.
aliases.BTmod_S = {'BTmod_S','BTmod_Sorg'};
aliases.MeanDepth_S = {'MeanDepth_S','Depth_S','DepthS'};

% Q aliases.
aliases.qsNew_S = {'qsNew_S','qsCbNew_S'};
aliases.qs1_SDD = {'qs1_SDD','qsVR1_S'};

% SonTek M9 backscatter aliases. The modern M9 processor retains its
% physically explicit native names; these aliases also populate the legacy
% feature-reader names carrying the _S instrument suffix.
aliases.alphaMean_S = {'alphaMean_S','alphaMean3000S'};
aliases.alpha1mean_Sgrad = {'alpha1mean_Sgrad','alpha1meanSgrad'};
aliases.alpha2mean_Sgrad = {'alpha2mean_Sgrad','alpha2meanSgrad'};
aliases.alpha3mean_Sgrad = {'alpha3mean_Sgrad','alpha3meanSgrad'};
aliases.alpha4mean_Sgrad = {'alpha4mean_Sgrad','alpha4meanSgrad'};
aliases.Af1mean_S = {'Af1mean_S','Af1meanS'};
aliases.Af2mean_S = {'Af2mean_S','Af2meanS'};
aliases.Af3mean_S = {'Af3mean_S','Af3meanS'};
aliases.Af4mean_S = {'Af4mean_S','Af4meanS'};
aliases.SpreadM_S = {'SpreadM_S','SpreadMS'};
aliases.SVBmean1_S = {'SVBmean1_S','SVmean1'};
aliases.SVBmean2_S = {'SVBmean2_S','SVmean2'};
aliases.SVBmean3_S = {'SVBmean3_S','SVmean3'};
aliases.SVBmean4_S = {'SVBmean4_S','SVmean4'};
aliases.SVBstd1_S = {'SVBstd1_S','SVmean1std'};
aliases.SVBstd2_S = {'SVBstd2_S','SVmean2std'};
aliases.SVBstd3_S = {'SVBstd3_S','SVmean3std'};
aliases.SVBstd4_S = {'SVBstd4_S','SVmean4std'};
aliases.SVB1mean1_S = {'SVB1mean1_S','SVmean1rssi'};
aliases.SVB1mean2_S = {'SVB1mean2_S','SVmean2rssi'};
aliases.SVB1mean3_S = {'SVB1mean3_S','SVmean3rssi'};
aliases.SVB1mean4_S = {'SVB1mean4_S','SVmean4rssi'};
aliases.SVB1std1_S = {'SVB1std1_S','SVmean1stdRSSI'};
aliases.SVB1std2_S = {'SVB1std2_S','SVmean2stdRSSI'};
aliases.SVB1std3_S = {'SVB1std3_S','SVmean3stdRSSI'};
aliases.SVB1std4_S = {'SVB1std4_S','SVmean4stdRSSI'};
aliases.SV11_S = {'SV11_S','SV1_Srssi'};
aliases.SV22_S = {'SV22_S','SV2_Srssi'};
aliases.SV33_S = {'SV33_S','SV3_Srssi'};
aliases.SV44_S = {'SV44_S','SV4_Srssi'};
aliases.TL1_S = {'TL1_S','TL1'};
aliases.TL2_S = {'TL2_S','TL2'};
aliases.TL3_S = {'TL3_S','TL3'};
aliases.TL4_S = {'TL4_S','TL4'};
aliases.TL1_fix_S = {'TL1_fix_S','TL1_fixS'};
aliases.TL2_fix_S = {'TL2_fix_S','TL2_fixS'};
aliases.TL3_fix_S = {'TL3_fix_S','TL3_fixS'};
aliases.TL4_fix_S = {'TL4_fix_S','TL4_fixS'};

%% Special sources and caller-workspace capture
special = struct();
if ~isempty(opt.dataADCPS)
    special.dataADCPS = opt.dataADCPS;
end
if ~isempty(opt.dataADCPSboat)
    special.dataADCPSboat = opt.dataADCPSboat;
end
special.FolderSontek = string(opt.FolderSontek);
special.FolderM9 = string(opt.FolderM9);

allRequested = unique([dataLegacyFields,shearLegacyFields,btLegacyFields, ...
    qLegacyFields,btampLegacyFields,newShearFields,newBTFields, ...
    newQFields,newBTampFields,featureDataFields,featureShearFields, ...
    featureBTvelFields,featureBTampFields,acousticFieldCandidates, ...
    allAliasNames(aliases)],'stable');
callerVariables = struct();
if logical(opt.UseCallerWorkspace)
    % This lookup is intentionally performed in the main function. Here,
    % evalin('caller',...) refers to the workspace that called
    % buildProcessingStructuresM9 rather than to a local helper.
    for iCaller=1:numel(allRequested)
        callerName = allRequested{iCaller};
        if ~isvarname(callerName), continue; end
        try
            callerExists = evalin('caller', ...
                sprintf('exist(''%s'',''var'')',callerName));
            if callerExists
                callerVariables.(callerName) = evalin('caller',callerName);
            end
        catch
            % Caller lookup is optional; unresolved fields become placeholders.
        end
    end
end

%% Route explicitly supplied supplemental variables
% Supplemental values are never stored in a sixth output structure. Each
% field is inserted directly into its established legacy destination.
[extraDataFields,extraShearFields,extraBTvelFields,extraQFields, ...
    extraBTampFields,unroutedSupplementalFields] = ...
    routeSupplementalVariables(SupplementalVariables, ...
        unique([dataLegacyFields,featureDataFields],'stable'), ...
        unique([shearLegacyFields,newShearFields,featureShearFields],'stable'), ...
        unique([btLegacyFields,newBTFields,featureBTvelFields],'stable'), ...
        unique([qLegacyFields,newQFields],'stable'), ...
        unique([btampLegacyFields,newBTampFields,featureBTampFields, ...
        acousticFieldCandidates],'stable'));

%% Build Data_field
% Preserve every input-table field first. Explicit data supplements are
% merged last so a deliberately supplied value replaces a placeholder.
Data_field = mergeStructs(ADCPin,extraDataFields);
dataSources = {special,extraDataFields,ADCPin,BedloadTransportOut, ...
    ShearStressOut,BackscatterOut,DepthRNGout,WaterProfileOut,BTvelOut, ...
    callerVariables};
dataSourceNames = {'special options','supplemental Data_field variables','ADCPin', ...
    'BedloadTransportOut','ShearStressOut','BackscatterOut', ...
    'DepthRNGout','WaterProfileOut','BTvelOut','caller workspace'};
[Data_field,dataReport] = ensureFields(Data_field,dataLegacyFields, ...
    aliases,dataSources,dataSourceNames,nRows,'Data_field');
if logical(opt.EnforceFeatureReaderFields)
    [Data_field,dataReaderReport] = ensureFields(Data_field,featureDataFields, ...
        aliases,dataSources,dataSourceNames,nRows,'Data_field');
    dataReport = appendReports(dataReport,dataReaderReport);
end

% Preserve the original orientation convention used by the legacy saving
% script for the full PSD matrices.
if logical(opt.FlipPSDFields)
    Data_field.PSD_percentage = flipud(Data_field.PSD_percentage);
    Data_field.QperDall = flipud(Data_field.QperDall);
end

%% Build ShearData
% Merge current outputs so all new shear/concentration/thickness fields are
% retained. Transport output is merged last because it contains the final
% consistency-checked Cb and thickness variables.
ShearData = mergeStructs(WaterProfileOut,DepthRNGout,ShearStressOut, ...
    BedloadTransportOut,extraShearFields);
shearSources = {special,extraShearFields,BedloadTransportOut,ShearStressOut, ...
    WaterProfileOut,DepthRNGout,ADCPin,BackscatterOut,BTvelOut, ...
    callerVariables};
shearSourceNames = {'special options','supplemental ShearData variables', ...
    'BedloadTransportOut','ShearStressOut','WaterProfileOut', ...
    'DepthRNGout','ADCPin','BackscatterOut','BTvelOut', ...
    'caller workspace'};
[ShearData,shearReportLegacy] = ensureFields(ShearData,shearLegacyFields, ...
    aliases,shearSources,shearSourceNames,nRows,'ShearData');
[ShearData,shearReportNew] = ensureFields(ShearData,newShearFields, ...
    aliases,shearSources,shearSourceNames,nRows,'ShearData');
shearReport = appendReports(shearReportLegacy,shearReportNew);
if logical(opt.EnforceFeatureReaderFields)
    [ShearData,shearReaderReport] = ensureFields(ShearData,featureShearFields, ...
        aliases,shearSources,shearSourceNames,nRows,'ShearData');
    shearReport = appendReports(shearReport,shearReaderReport);
end

% The old reader expects one temporal standard deviation per bottom-track
% beam. The current shear function exports these diagnostics together as
% D0BeamStd_S. Use those columns when the individual legacy names were
% not available.
[ShearData,shearReport] = resolveLegacyBeamStd( ...
    ShearData,ShearStressOut,DepthRNGout,shearReport,nRows);

% Explicitly guarantee that the final concentration from processBedloadTrate
% is retained in ShearData. The canonical saved field name remains
% CbFinal_S. This block also accepts older capitalization variants.
[ShearData,shearReport,CbFinalSaveReport] = ensureCbFinalInShearData( ...
    ShearData,BedloadTransportOut,ShearStressOut,callerVariables, ...
    shearReport,nRows);
if ~CbFinalSaveReport.PresentAndFinite
    warning('buildProcessingStructuresM9:CbFinalMissing', ...
        ['ShearData.CbFinal_S was created, but it contains no finite ', ...
         'values. Check BedloadTransportOut.CbFinal_S and the Cb ', ...
         'fallback settings before using the saved shear data.']);
end

%% Build BTvel
BTvel = mergeStructs(BTvelOut,extraBTvelFields);
btSources = {extraBTvelFields,BTvelOut,DepthRNGout,WaterProfileOut, ...
    BackscatterOut,ADCPin,callerVariables};
btSourceNames = {'supplemental BTvel variables','BTvelOut','DepthRNGout', ...
    'WaterProfileOut','BackscatterOut','ADCPin','caller workspace'};
[BTvel,btReportLegacy] = ensureFields(BTvel,btLegacyFields,aliases, ...
    btSources,btSourceNames,nRows,'BTvel');
[BTvel,btReportNew] = ensureFields(BTvel,newBTFields,aliases, ...
    btSources,btSourceNames,nRows,'BTvel');
btReport = appendReports(btReportLegacy,btReportNew);
if logical(opt.EnforceFeatureReaderFields)
    [BTvel,btReaderReport] = ensureFields(BTvel,featureBTvelFields,aliases, ...
        btSources,btSourceNames,nRows,'BTvel');
    btReport = appendReports(btReport,btReaderReport);
end

% Derive MeanDepth_S only when it was not supplied meaningfully.
if isPlaceholderValue(BTvel.MeanDepth_S) && isfield(DepthRNGout,'Depth_S')
    BTvel.MeanDepth_S = rowMeanFinite(DepthRNGout.Depth_S,nRows);
    btReport = updateReportResolution(btReport,'MeanDepth_S', ...
        'DepthRNGout','Depth_S (row mean)');
end

%% Build Q
Q = mergeStructs(selectTransportFields(BedloadTransportOut),extraQFields);
qSources = {extraQFields,BedloadTransportOut,BackscatterOut,ShearStressOut, ...
    BTvelOut,ADCPin,callerVariables};
qSourceNames = {'supplemental Q variables','BedloadTransportOut','BackscatterOut', ...
    'ShearStressOut','BTvelOut','ADCPin','caller workspace'};
[Q,qReportLegacy] = ensureFields(Q,qLegacyFields,aliases,qSources, ...
    qSourceNames,nRows,'Q');
[Q,qReportNew] = ensureFields(Q,newQFields,aliases,qSources, ...
    qSourceNames,nRows,'Q');
qReport = appendReports(qReportLegacy,qReportNew);

%% Build BTamp
BTamp = mergeStructs(BackscatterOut,Acoustic,extraBTampFields);
btampSources = {special,extraBTampFields,Acoustic,BackscatterOut, ...
    DepthRNGout,ADCPin,BTvelOut,BedloadTransportOut,callerVariables};
btampSourceNames = {'special options','supplemental BTamp variables', ...
    'Acoustic','BackscatterOut','DepthRNGout','ADCPin','BTvelOut', ...
    'BedloadTransportOut','caller workspace'};
[BTamp,btampReportLegacy] = ensureFields(BTamp,btampLegacyFields,aliases, ...
    btampSources,btampSourceNames,nRows,'BTamp');
[BTamp,btampReportNew] = ensureFields(BTamp,newBTampFields,aliases, ...
    btampSources,btampSourceNames,nRows,'BTamp');
btampReport = appendReports(btampReportLegacy,btampReportNew);
if logical(opt.EnforceFeatureReaderFields)
    [BTamp,btampReaderReport] = ensureFields(BTamp,featureBTampFields, ...
        aliases,btampSources,btampSourceNames,nRows,'BTamp');
    btampReport = appendReports(btampReport,btampReaderReport);
end

% The reader obtains D50 from DataBTamp rather than Data_field. Make the
% compatibility copy explicit and deterministic.
if isPlaceholderValue(BTamp.D50_05) && ~isPlaceholderValue(Data_field.D50_05)
    BTamp.D50_05 = Data_field.D50_05;
    btampReport = updateReportResolution(btampReport,'D50_05', ...
        'Data_field','D50_05');
end
if isfield(Data_field,'D50_14') && ...
        (~isfield(BTamp,'D50_14') || isPlaceholderValue(BTamp.D50_14))
    BTamp.D50_14 = Data_field.D50_14;
    btampReport = updateReportResolution(btampReport,'D50_14', ...
        'Data_field','D50_14');
end

% Normalize only the fields consumed directly by the supplied reader.
Data_field = normalizeReaderStructure(Data_field,featureDataFields,nRows,'Data_field');
ShearData = normalizeReaderStructure(ShearData,featureShearFields,nRows,'ShearData');
BTvel = normalizeReaderStructure(BTvel,featureBTvelFields,nRows,'BTvel');
BTamp = normalizeReaderStructure(BTamp,featureBTampFields,nRows,'BTamp');

% Historical downstream scripts sometimes load DataBTamp instead of BTamp.
DataBTamp = BTamp; %#ok<NASGU>

% Check the measured physical-sample inputs before saving. The check uses
% ADCPin first so an all-NaN source column is not mistaken for an ordinary
% legacy placeholder created later in Data_field.
physicalSampleFields = normalizeFieldNameList(opt.PhysicalSampleFields);
PhysicalSampleNaNReport = makePhysicalSampleNaNReport( ...
    ADCPin,Data_field,physicalSampleFields,nRows);
if logical(opt.WarnOnPhysicalSampleNaN)
    emitPhysicalSampleNaNWarnings(PhysicalSampleNaNReport,ADCPin);
end

% Record the exact acoustic fields flattened into BTamp. The Acoustic
% structure itself is intentionally not nested or saved separately.
acousticNames = string(fieldnames(Acoustic));
if isempty(acousticNames)
    acousticToBTampReport = table(strings(0,1),false(0,1), ...
        'VariableNames',{'Field','PresentInBTamp'});
else
    acousticPresent = false(numel(acousticNames),1);
    for iAcoustic = 1:numel(acousticNames)
        acousticPresent(iAcoustic) = isfield(BTamp,char(acousticNames(iAcoustic)));
    end
    acousticToBTampReport = table(acousticNames,acousticPresent, ...
        'VariableNames',{'Field','PresentInBTamp'});
end

%% Build report
StructureReport = struct();
StructureReport.CreatedBy = mfilename;
StructureReport.NumberOfRows = nRows;
StructureReport.FolderSontek = string(opt.FolderSontek);
StructureReport.FolderM9 = string(opt.FolderM9);
StructureReport.FolderTrates = string(opt.FolderTrates);
StructureReport.FolderShear = string(opt.FolderShear);
StructureReport.FolderBS = string(opt.FolderBS);
StructureReport.Data_field = dataReport;
StructureReport.ShearData = shearReport;
StructureReport.BTvel = btReport;
StructureReport.Q = qReport;
StructureReport.BTamp = btampReport;
StructureReport.CbFinalInShearData = CbFinalSaveReport;
StructureReport.PhysicalSampleNaNReport = PhysicalSampleNaNReport;
StructureReport.AcousticToBTamp = acousticToBTampReport;
StructureReport.InputFieldRenaming = InputFieldRenaming;
StructureReport.SupplementalRouting = table( ...
    [repmat("Data_field",numel(fieldnames(extraDataFields)),1); ...
     repmat("ShearData",numel(fieldnames(extraShearFields)),1); ...
     repmat("BTvel",numel(fieldnames(extraBTvelFields)),1); ...
     repmat("Q",numel(fieldnames(extraQFields)),1); ...
     repmat("BTamp",numel(fieldnames(extraBTampFields)),1)], ...
    string([fieldnames(extraDataFields);fieldnames(extraShearFields); ...
     fieldnames(extraBTvelFields);fieldnames(extraQFields); ...
     fieldnames(extraBTampFields)]), ...
    'VariableNames',{'DestinationStructure','Field'});
StructureReport.UnroutedSupplementalFields = string(unroutedSupplementalFields(:));
StructureReport.ReaderCompatibility = makeReaderCompatibilityReport( ...
    Data_field,ShearData,BTvel,Q,BTamp, ...
    featureDataFields,featureShearFields,featureBTvelFields,featureBTampFields,nRows);
StructureReport.ReaderTopLevelVariables = table( ...
    ["Trates file";"Trates file";"Trates file";"Shear file"; ...
     "Backscatter file";"Backscatter file"], ...
    ["Data_field";"BTvel";"Q";"ShearData";"BTamp";"DataBTamp"], ...
    'VariableNames',{'FileRole','RequiredVariable'});

RemainingLegacyInstrumentFieldNames = unique([ ...
    collectLegacyInstrumentFieldNames(Data_field,"Data_field"); ...
    collectLegacyInstrumentFieldNames(ShearData,"ShearData"); ...
    collectLegacyInstrumentFieldNames(BTvel,"BTvel"); ...
    collectLegacyInstrumentFieldNames(Q,"Q"); ...
    collectLegacyInstrumentFieldNames(BTamp,"BTamp")],'stable');
StructureReport.RemainingLegacyInstrumentFieldNames = ...
    RemainingLegacyInstrumentFieldNames;
% Retained as a compatibility alias; it now reports both legacy tokens.
StructureReport.Remaining1200FieldNames = ...
    RemainingLegacyInstrumentFieldNames;
StructureReport.AllInstrumentFieldsUseS = ...
    isempty(RemainingLegacyInstrumentFieldNames);
if logical(opt.FailOn1200FieldNames) && ...
        ~isempty(RemainingLegacyInstrumentFieldNames)
    error('buildProcessingStructuresM9:NamesStillContainLegacyToken', ...
        ['M9 packaging stopped because output field names still contain ', ...
         '1200 or an isolated 12 token: %s'], ...
        strjoin(RemainingLegacyInstrumentFieldNames,', '));
end
StructureReport.MissingFields = unique([ ...
    string(dataReport.Field(dataReport.WasPlaceholder)); ...
    string(shearReport.Field(shearReport.WasPlaceholder)); ...
    string(btReport.Field(btReport.WasPlaceholder)); ...
    string(qReport.Field(qReport.WasPlaceholder)); ...
    string(btampReport.Field(btampReport.WasPlaceholder))]);
StructureReport.Notes = [ ...
    "All legacy field names from the supplied saving script are created."; ...
    "Every instrument-specific legacy field uses the normalized _S token; input fields containing 1200 or an isolated 12 token are renamed before packaging."; ...
    "All current fields from the modern output structures are retained."; ...
    "PSD_percentage and QperDall use the legacy flipud orientation when FlipPSDFields=true."; ...
    "Fr_S is taken from the actual Fr_S output; the old accidental assignment Fr_S=Tausk_S is not repeated."; ...
    "Tau_shear2VR_S is not inferred from another stress variable because its legacy dimensional definition is ambiguous."; ...
    "Ks2_S is stored directly in BTamp/DataBTamp and may be resolved from the M9 3 MHz fields Ks2_3 or Ks_3."; ...
    "Legacy transport uncertainty fields ending in stdBT are stored directly in Q when supplied; otherwise they remain placeholders."; ...
    "CbFinal_S is explicitly retained in ShearData from BedloadTransportOut, with older Cbfinal capitalization variants accepted as input."; ...
    "Measured transport and grain-size input fields are checked for NaN/missing values before saving; details are stored in PhysicalSampleNaNReport."; ...
    "DataBTamp is created as an exact compatibility alias of BTamp when assigning or saving."; ...
    "Every field from computeADCPacousticParameters output Acoustic is flattened directly into BTamp and DataBTamp."; ...
    "When Acoustic is omitted, recognized acoustic fields are recovered from the updated ADCPin when available."; ...
    "Legacy Ks2_12 and Ks_12 inputs are accepted, normalized, and saved only as Ks2_S or Ks_S."; ...
    "The Trates MAT file contains Data_field, BTvel and Q; the Shear MAT file contains only ShearData."; ...
    "All automatically generated filenames start with the same YYYYMMDD source-date prefix."; ...
    "When LocationName is supplied, it is sanitized and inserted after the date in all saved filenames."; ...
    "ReaderCompatibility lists every field referenced directly by the supplied feature-extraction script."; ...
    "Direct supplemental variables are routed into the existing five structures; no ExtraVariables structure is assigned or saved."; ...
    "DeltaT_S routes to BTvel; Ks2_S routes to BTamp/DataBTamp; Tau_shear2VR_S, Vmod_SstdMOD and Vmod1low_SstdMOD route to ShearData."; ...
    "Unavailable values are not fabricated; they remain explicit placeholders and are flagged." ...
    ];

if logical(opt.FailOnMissingRequiredValues)
    bad = StructureReport.ReaderCompatibility.IsPlaceholder;
    if any(bad)
        missingText = strjoin( ...
            StructureReport.ReaderCompatibility.Structure(bad) + "." + ...
            StructureReport.ReaderCompatibility.Field(bad),", ");
        error('Required feature-reader values remain unresolved: %s',missingText);
    end
end

%% Resolve source date and destination file paths
needsDatePrefix = strlength(string(opt.FolderTrates))>0 || ...
    strlength(string(opt.FolderShear))>0 || ...
    strlength(string(opt.FolderBS))>0 || ...
    strlength(string(opt.SaveDataFile))>0 || ...
    strlength(string(opt.SaveShearFile))>0 || ...
    strlength(string(opt.SaveBackscatterFile))>0 || ...
    strlength(string(opt.SaveCombinedFile))>0;

[FileDatePrefix,FileDateSource] = resolveSourceDatePrefix( ...
    opt.FileDatePrefix,opt.FolderM9,opt.FolderSontek,opt.dataADCPS,ADCPin);

[LocationName,LocationToken] = normalizeLocationName(opt.LocationName);
FileNamePrefix = composeFileNamePrefix(FileDatePrefix,LocationToken);

if needsDatePrefix && FileDatePrefix == ""
    error(['Could not extract an 8-digit source date from FolderM9, ', ...
        'FolderSontek, dataADCPS or ADCPin.ADCPfilename. ', ...
        'Provide ''FileDatePrefix'',''YYYYMMDD''.']);
end

saveDataPath = resolveOutputFilePath( ...
    opt.SaveDataFile,opt.FolderTrates,FileDatePrefix,LocationToken,"TratesData.mat");
saveShearPath = resolveOutputFilePath( ...
    opt.SaveShearFile,opt.FolderShear,FileDatePrefix,LocationToken,"ShearData.mat");
saveBackscatterPath = resolveOutputFilePath( ...
    opt.SaveBackscatterFile,opt.FolderBS,FileDatePrefix,LocationToken,"BackscatterData.mat");
saveCombinedPath = resolveExplicitSavePath( ...
    opt.SaveCombinedFile,FileDatePrefix,LocationToken);

StructureReport.FileDatePrefix = FileDatePrefix;
StructureReport.FileDateSource = FileDateSource;
StructureReport.LocationName = LocationName;
StructureReport.LocationToken = LocationToken;
StructureReport.FileNamePrefix = FileNamePrefix;
StructureReport.SavePlan = table( ...
    ["Trates";"Shear";"Backscatter";"Combined"], ...
    [saveDataPath;saveShearPath;saveBackscatterPath;saveCombinedPath], ...
    ["Data_field, BTvel, Q";"ShearData"; ...
     "BTamp, DataBTamp";"all structures"], ...
    'VariableNames',{'FileRole','FilePath','TopLevelVariables'});

%% Assign structures
if logical(opt.AssignToCaller)
    assignin('caller','Data_field',Data_field);
    assignin('caller','ShearData',ShearData);
    assignin('caller','BTvel',BTvel);
    assignin('caller','Q',Q);
    assignin('caller','BTamp',BTamp);
    assignin('caller','DataBTamp',DataBTamp);
    assignin('caller','StructureReport',StructureReport);
end
if logical(opt.AssignToBase)
    assignin('base','Data_field',Data_field);
    assignin('base','ShearData',ShearData);
    assignin('base','BTvel',BTvel);
    assignin('base','Q',Q);
    assignin('base','BTamp',BTamp);
    assignin('base','DataBTamp',DataBTamp);
    assignin('base','StructureReport',StructureReport);
end

%% Optional MAT-file saving
saveVersion = char(string(opt.SaveVersion));
if saveDataPath ~= ""
    fileName = prepareSavePath(saveDataPath);
    saveThreeWithVersion(fileName,saveVersion, ...
        'Data_field',Data_field,'BTvel',BTvel,'Q',Q);
end
if saveShearPath ~= ""
    fileName = prepareSavePath(saveShearPath);
    saveWithVersion(fileName,saveVersion,'ShearData',ShearData);
end
if saveBackscatterPath ~= ""
    fileName = prepareSavePath(saveBackscatterPath);
    saveTwoWithVersion(fileName,saveVersion,'BTamp',BTamp,'DataBTamp',DataBTamp);
end
if saveCombinedPath ~= ""
    fileName = prepareSavePath(saveCombinedPath);
    saveCombined(fileName,saveVersion,Data_field,ShearData,BTvel,Q,BTamp, ...
        DataBTamp,StructureReport);
end

%% Summary
if logical(opt.Verbose)
    fprintf('\n=== Legacy/project structure packaging summary ===\n');
    fprintf('Rows inferred: %d\n',nRows);
    fprintf('Data_field fields: %d | placeholders: %d\n', ...
        numel(fieldnames(Data_field)),sum(dataReport.WasPlaceholder));
    fprintf('ShearData fields:  %d | placeholders: %d\n', ...
        numel(fieldnames(ShearData)),sum(shearReport.WasPlaceholder));
    fprintf('BTvel fields:      %d | placeholders: %d\n', ...
        numel(fieldnames(BTvel)),sum(btReport.WasPlaceholder));
    fprintf('Q fields:          %d | placeholders: %d\n', ...
        numel(fieldnames(Q)),sum(qReport.WasPlaceholder));
    fprintf('BTamp fields:      %d | placeholders: %d\n', ...
        numel(fieldnames(BTamp)),sum(btampReport.WasPlaceholder));
    fprintf('DataBTamp compatibility alias created.\n');
    fprintf('Input fields renamed from 1200 to _S: %d\n',height(InputFieldRenaming));
    fprintf('Acoustic fields copied directly into BTamp/DataBTamp: %d\n', ...
        numel(fieldnames(Acoustic)));
    fprintf('CbFinal_S retained in ShearData: %s\n', ...
        char(string(CbFinalSaveReport.PresentAndFinite)));
    if isempty(PhysicalSampleNaNReport)
        fprintf('Physical-sample NaN check: no checked fields contained NaN/missing values.\n');
    else
        fprintf('Physical-sample NaN check: %d field(s) contain NaN/missing values.\n', ...
            height(PhysicalSampleNaNReport));
    end
    if FileDatePrefix ~= ""
        fprintf('Source date prefix: %s (from %s)\n',FileDatePrefix,FileDateSource);
        if LocationToken ~= ""
            fprintf('Location token: %s (from LocationName=%s)\n',LocationToken,LocationName);
        end
        fprintf('Common filename prefix: %s\n',FileNamePrefix);
    end
    if saveDataPath ~= "", fprintf('Trates file: %s\n',saveDataPath); end
    if saveShearPath ~= "", fprintf('Shear file: %s\n',saveShearPath); end
    if saveBackscatterPath ~= "", fprintf('Backscatter file: %s\n',saveBackscatterPath); end
    fprintf('Supplemental fields routed directly: Data=%d, Shear=%d, BTvel=%d, Q=%d, BTamp=%d.\n', ...
        numel(fieldnames(extraDataFields)),numel(fieldnames(extraShearFields)), ...
        numel(fieldnames(extraBTvelFields)),numel(fieldnames(extraQFields)), ...
        numel(fieldnames(extraBTampFields)));
    if ~isempty(unroutedSupplementalFields)
        fprintf('Unrouted supplemental fields: %s\n', ...
            strjoin(string(unroutedSupplementalFields),', '));
    end
    fprintf('Feature-reader fields checked: %d | placeholders: %d\n', ...
        height(StructureReport.ReaderCompatibility), ...
        sum(StructureReport.ReaderCompatibility.IsPlaceholder));
    if ~isempty(StructureReport.MissingFields)
        fprintf('Unresolved legacy fields are preserved as placeholders.\n');
        fprintf('Inspect StructureReport and especially StructureReport.MissingFields.\n');
    end
    fprintf('==================================================\n');
end
end

%% Local helpers
function names = defaultPhysicalSampleFieldNames()
%DEFAULTPHYSICALSAMPLEFIELDNAMES Physical measurements checked for NaN.
% Only fields that actually exist in ADCPin, or contain a non-placeholder
% value in Data_field, are checked and reported.
names = { ...
    'TrateMean','MeanTrate','Tmean','Tmea', ...
    'MaxTrate','TrateMax','Tmax', ...
    'MinTrate','TrateMin','Tmin', ...
    'MeanTrate05','Max05','Min05','MeanTrate14','Max14','Min14', ...
    'MMeanTrate','MMaxTrate','MMinTrate', ...
    'D10','D16','D50','D65','D84','D90', ...
    'D10_05','D16_05','D50_05','D65_05','D84_05','D90_05', ...
    'D10_14','D16_14','D50_14','D65_14','D84_14','D90_14', ...
    'D50mobile','D50grab','D90grab', ...
    'D_10grabM','D_10grabMAX','D_16grabM','D_16grabMAX', ...
    'D_50grabM','D_50grabMAX','D_84grabM','D_84grabMAX', ...
    'D_90grabM','D_90grabMAX','DmeaNUMpsd','SSC'};
end

function names = normalizeFieldNameList(value)
%NORMALIZEFIELDNAMELIST Convert a char/string/cell list to valid unique names.
if isempty(value)
    names = {};
elseif ischar(value)
    names = {value};
elseif isstring(value)
    names = cellstr(value(:));
elseif iscell(value)
    names = cell(size(value(:)));
    for i=1:numel(value)
        names{i}=char(string(value{i}));
    end
else
    names = {};
end
names = names(:).';
keep = false(size(names));
for i=1:numel(names)
    names{i}=strtrim(names{i});
    keep(i)=~isempty(names{i}) && isvarname(names{i});
end
names=unique(names(keep),'stable');
end

function [ShearData,report,status] = ensureCbFinalInShearData( ...
    ShearData,BedloadTransportOut,ShearStressOut,callerVariables,report,nRows)
%ENSURECBFINALINSHEARDATA Retain final Cb in the saved ShearData structure.
canonical='CbFinal_S';
candidates={'CbFinal_S','Cbfinal_S','CbFinal','Cbfinal'};
sources={BedloadTransportOut,ShearStressOut,callerVariables};
sourceNames={'BedloadTransportOut','ShearStressOut','caller workspace'};
sourceContainer="existing ShearData";
sourceField=canonical;
if isfield(ShearData,canonical) && ~isPlaceholderValue(ShearData.(canonical))
    value=ShearData.(canonical); %#ok<NASGU>
else
    [found,value,sourceContainer,sourceField]=resolveField( ...
        candidates,sources,sourceNames);
    if found
        ShearData.(canonical)=value;
        report=updateReportResolution(report,canonical,sourceContainer,sourceField);
    else
        if ~isfield(ShearData,canonical)
            ShearData.(canonical)=nan(nRows,1);
        end
        sourceContainer="placeholder";
        sourceField="missing";
    end
end
if isfield(ShearData,canonical)
    x=ShearData.(canonical);
    if isnumeric(x) && isvector(x)
        ShearData.(canonical)=fitLength(double(x(:)),nRows);
    end
end
present=isfield(ShearData,canonical);
finiteCount=0;
if present && isnumeric(ShearData.(canonical))
    cbValue=ShearData.(canonical);
    finiteCount=sum(isfinite(cbValue(:)));
end
status=table(string(canonical),logical(present),finiteCount>0,finiteCount, ...
    string(sourceContainer),string(sourceField), ...
    'VariableNames',{'Field','Present','PresentAndFinite','FiniteValueCount', ...
    'SourceContainer','SourceField'});
end

function T = makePhysicalSampleNaNReport(ADCPin,Data_field,fields,nRows)
%MAKEPHYSICALSAMPLENANREPORT Report NaN/missing values in measured inputs.
% Prefer exact fields in ADCPin. This avoids duplicate warnings from legacy
% aliases subsequently created in Data_field (for example TrateMean and
% MeanTrate referring to the same input column). Data_field is used only as
% a fallback when none of the requested fields exists in ADCPin.
Field=strings(0,1); Source=strings(0,1); NaNCount=zeros(0,1);
AffectedRowCount=zeros(0,1); AffectedRows=strings(0,1);
CheckedSize=strings(0,1);
hasAnyADCPinPhysicalField=false;
for j=1:numel(fields)
    if isfield(ADCPin,fields{j})
        hasAnyADCPinPhysicalField=true;
        break
    end
end
for i=1:numel(fields)
    f=fields{i};
    source=""; x=[];
    if isfield(ADCPin,f)
        x=ADCPin.(f); source="ADCPin";
    elseif ~hasAnyADCPinPhysicalField && isfield(Data_field,f) && ...
            ~isPlaceholderValue(Data_field.(f))
        x=Data_field.(f); source="Data_field fallback";
    else
        continue
    end
    [rowMask,nMissing,sizeText,wasCheckable]=missingMaskByRow(x,nRows);
    if ~wasCheckable || nMissing==0, continue; end
    rows=find(rowMask);
    Field(end+1,1)=string(f); %#ok<AGROW>
    Source(end+1,1)=source; %#ok<AGROW>
    NaNCount(end+1,1)=nMissing; %#ok<AGROW>
    AffectedRowCount(end+1,1)=numel(rows); %#ok<AGROW>
    AffectedRows(end+1,1)=join(string(rows(:).'),", "); %#ok<AGROW>
    CheckedSize(end+1,1)=sizeText; %#ok<AGROW>
end
T=table(Field,Source,NaNCount,AffectedRowCount,AffectedRows,CheckedSize);
end

function [rowMask,nMissing,sizeText,checkable] = missingMaskByRow(x,nRows)
%MISSINGMASKBYROW Detect numeric NaN/nonfinite or text/cell missing values.
rowMask=false(nRows,1); nMissing=0; checkable=true;
try
    sizeText=strjoin(string(size(x)),'x');
catch
    sizeText="unknown";
end
if isnumeric(x)
    A=double(x);
    if isempty(A), checkable=false; return; end
    missing=~isfinite(A);
    nMissing=sum(missing(:));
    if isvector(A)
        rowMask=fitLogicalLength(missing(:),nRows);
    elseif size(A,1)==nRows
        rowMask=any(missing,2);
    elseif size(A,2)==nRows
        rowMask=any(missing,1).';
    else
        rowMask(:)=any(missing(:));
    end
elseif isstring(x)
    v=x(:);
    missing=ismissing(v) | strlength(strtrim(v))==0 | lower(strtrim(v))=="nan";
    nMissing=sum(missing);
    rowMask=fitLogicalLength(missing,nRows);
elseif iscell(x)
    v=x(:); missing=false(numel(v),1);
    for k=1:numel(v)
        item=v{k};
        if isempty(item)
            missing(k)=true;
        elseif isnumeric(item)
            missing(k)=all(~isfinite(double(item(:))));
        elseif isstring(item) || ischar(item)
            t=strtrim(string(item));
            missing(k)=all(ismissing(t) | t=="" | lower(t)=="nan");
        end
    end
    nMissing=sum(missing);
    rowMask=fitLogicalLength(missing,nRows);
else
    checkable=false;
end
end

function y = fitLogicalLength(x,n)
x=logical(x(:));
if numel(x)==1 && n>1
    y=repmat(x,n,1);
elseif numel(x)<n
    y=[x;false(n-numel(x),1)];
else
    y=x(1:n);
end
end

function emitPhysicalSampleNaNWarnings(T,ADCPin)
%EMITPHYSICALSAMPLENANWARNINGS Issue one concise warning per affected field.
if isempty(T), return; end
fileNames=strings(0,1);
if isstruct(ADCPin) && isfield(ADCPin,'ADCPfilename')
    fileNames=string(ADCPin.ADCPfilename(:));
end
for i=1:height(T)
    rows=str2double(strtrim(split(T.AffectedRows(i),",")));
    rows=rows(isfinite(rows));
    fileText="";
    if ~isempty(fileNames) && ~isempty(rows)
        rows=rows(rows>=1 & rows<=numel(fileNames));
        if ~isempty(rows)
            listed=fileNames(rows);
            maxShown=min(5,numel(listed));
            fileText=" Files: "+join(listed(1:maxShown),", ");
            if numel(listed)>maxShown, fileText=fileText+", ..."; end
        end
    end
    warning('buildProcessingStructuresM9:PhysicalSampleNaN', ...
        ['Physical sample field %s contains %d NaN/missing value(s) ', ...
         'affecting row(s) %s.%s'], ...
        char(T.Field(i)),T.NaNCount(i),char(T.AffectedRows(i)),char(fileText));
end
end

function names = defaultAcousticFieldNames()
% Fields produced by computeADCPacousticParameters and earlier input-table
% acoustic calculations. Legacy 1200/12 names use the canonical _S token.
names = { ...
    'c','ro','rs','k','f3','f1','f05','f_S','f06', ...
    'lambda3','k3','x3','F3','Kappa3','FMean3','Ks_3','Ks2_3', ...
    'lambda1','k1','x1','F1','Kappa1','FMean1','Ks_1','Ks2_1', ...
    'lambda05','k05','x05','F05','Kappa05','FMean05', ...
    'Ks_05','Ks2_05', ...
    'lambda_S','k_S','x_S','F_S','Kappa_S','FMean_S', ...
    'Ks_S','Ks2_S', ...
    'lambda06','k06','x06','F06','Kappa06','FMean06', ...
    'Ks_06','Ks2_06', ...
    'psi3','psi1','psi05','psi06','psi_S','DmeaNUMpsd'};
end

function out = selectFieldsPresent(S,names)
out = struct();
if isempty(S) || ~isstruct(S) || ~isscalar(S), return; end
for i = 1:numel(names)
    f = names{i};
    if isfield(S,f)
        out.(f) = S.(f);
    end
end
end

function S = scalarStructOrEmpty(S)
if isempty(S)
    S = struct();
elseif ~isstruct(S)
    error('Expected a structure or empty input.');
elseif ~isscalar(S)
    error('Modern output inputs must be scalar structures.');
end
end

function T = emptyFieldRenameReport()
T = table('Size',[0 6], ...
    'VariableTypes',repmat("string",1,6), ...
    'VariableNames',{'RootStructure','ContainerPath','OldName', ...
    'NewName','ObjectType','Action'});
end

function [value,report] = normalizeLegacyInstrumentFieldsToS( ...
    value,rootName,path,report)
% Recursively normalize 1200 and isolated 12 tokens without data edits.

if isstruct(value)
    names = fieldnames(value);
    if isempty(names), return; end

    if isempty(value)
        converted = cellfun(@m9FieldName,names,'UniformOutput',false);
        for k=1:numel(names)
            if ~strcmp(names{k},converted{k})
                report = appendFieldRename(report,rootName,path,names{k}, ...
                    converted{k},"struct_field","renamed_empty_structure_field");
            end
        end
        [converted,keep] = unique(converted,'stable');
        if numel(keep) < numel(names)
            % Duplicate targets in an empty structure contain no values to
            % reconcile; retaining one target is deterministic.
        end
        oldSize = size(value);
        template = cell2struct(cell(numel(converted),1),converted(:),1);
        value = repmat(template,oldSize);
        return
    end

    % Process nested contents before renaming their parent fields.
    for k=1:numel(names)
        oldName = names{k};
        for n=1:numel(value)
            if numel(value)==1
                childPath = path+"."+string(oldName);
            else
                childPath = path+"("+string(n)+")."+string(oldName);
            end
            [value(n).(oldName),report] = normalizeLegacyInstrumentFieldsToS( ...
                value(n).(oldName),rootName,childPath,report);
        end
    end

    names = fieldnames(value);
    for k=1:numel(names)
        oldName = names{k};
        if ~isLegacyInstrumentFieldName(oldName) || ...
                ~isfield(value,oldName), continue; end
        newName = m9FieldName(oldName);
        if isfield(value,newName)
            useOld = structFieldIsPlaceholder(value,newName) && ...
                ~structFieldIsPlaceholder(value,oldName);
            if useOld
                for n=1:numel(value)
                    value(n).(newName)=value(n).(oldName);
                end
                action="collision_replaced_placeholder_s";
            elseif structFieldsEqualLocal(value,oldName,newName)
                action="identical_collision_removed_legacy";
            else
                action="collision_kept_existing_s";
            end
            value=rmfield(value,oldName);
        else
            for n=1:numel(value)
                value(n).(newName)=value(n).(oldName);
            end
            value=rmfield(value,oldName);
            action="renamed";
        end
        report = appendFieldRename(report,rootName,path,oldName,newName, ...
            "struct_field",action);
    end

elseif istable(value) || isa(value,'timetable')
    names=value.Properties.VariableNames;
    for k=1:numel(names)
        oldName=names{k};
        [value.(oldName),report]=normalizeLegacyInstrumentFieldsToS( ...
            value.(oldName), ...
            rootName,path+"."+string(oldName),report);
    end
    names=value.Properties.VariableNames;
    for k=1:numel(names)
        oldName=names{k};
        if ~isLegacyInstrumentFieldName(oldName) || ...
                ~ismember(oldName,value.Properties.VariableNames)
            continue
        end
        newName=m9FieldName(oldName);
        if ismember(newName,value.Properties.VariableNames)
            if isPlaceholderValue(value.(newName)) && ...
                    ~isPlaceholderValue(value.(oldName))
                value.(newName)=value.(oldName);
                action="collision_replaced_placeholder_s";
            elseif isequaln(value.(oldName),value.(newName))
                action="identical_collision_removed_legacy";
            else
                action="collision_kept_existing_s";
            end
            value(:,{oldName})=[];
        else
            idx=strcmp(value.Properties.VariableNames,oldName);
            value.Properties.VariableNames{idx}=newName;
            action="renamed";
        end
        report=appendFieldRename(report,rootName,path,oldName,newName, ...
            "table_variable",action);
    end

elseif iscell(value)
    for n=1:numel(value)
        [value{n},report]=normalizeLegacyInstrumentFieldsToS( ...
            value{n},rootName, ...
            path+"{"+string(n)+"}",report);
    end
end
end

function newName = m9FieldName(oldName)
% Normalize full 1200 and abbreviated isolated 12 instrument tokens.
newName=regexprep(oldName,'_?1200','_S');
newName=regexprep(newName,'_?(?<!\d)12(?!\d)','_S');
end

function tf = isLegacyInstrumentFieldName(name)
% Do not match 12 when it is part of a larger number such as 0125.
tf = contains(name,'1200') || ...
    ~isempty(regexp(name,'(?<!\d)12(?!\d)','once'));
end

function tf = structFieldIsPlaceholder(S,fieldName)
tf=true;
for n=1:numel(S)
    if ~isPlaceholderValue(S(n).(fieldName))
        tf=false;
        return
    end
end
end

function tf = structFieldsEqualLocal(S,fieldA,fieldB)
tf=true;
for n=1:numel(S)
    if ~isequaln(S(n).(fieldA),S(n).(fieldB))
        tf=false;
        return
    end
end
end

function report = appendFieldRename(report,rootName,path,oldName,newName, ...
    objectType,action)
row=table(string(rootName),string(path),string(oldName),string(newName), ...
    string(objectType),string(action), ...
    'VariableNames',report.Properties.VariableNames);
report=[report;row]; %#ok<AGROW>
end

function namesFound = collectLegacyInstrumentFieldNames(value,path)
namesFound=strings(0,1);
if isstruct(value)
    names=fieldnames(value);
    for k=1:numel(names)
        thisName=names{k};
        if isLegacyInstrumentFieldName(thisName)
            namesFound(end+1,1)=path+"."+string(thisName); %#ok<AGROW>
        end
        for n=1:numel(value)
            if numel(value)==1
                childPath=path+"."+string(thisName);
            else
                childPath=path+"("+string(n)+")."+string(thisName);
            end
            namesFound=[namesFound;collectLegacyInstrumentFieldNames( ...
                value(n).(thisName),childPath)]; %#ok<AGROW>
        end
    end
elseif istable(value) || isa(value,'timetable')
    names=value.Properties.VariableNames;
    for k=1:numel(names)
        thisName=names{k};
        if isLegacyInstrumentFieldName(thisName)
            namesFound(end+1,1)=path+"."+string(thisName); %#ok<AGROW>
        end
        namesFound=[namesFound;collectLegacyInstrumentFieldNames( ...
            value.(thisName), ...
            path+"."+string(thisName))]; %#ok<AGROW>
    end
elseif iscell(value)
    for n=1:numel(value)
        namesFound=[namesFound;collectLegacyInstrumentFieldNames(value{n}, ...
            path+"{"+string(n)+"}")]; %#ok<AGROW>
    end
end
end

function nRows = inferRowCount(sources)
nRows = 0;
preferred = {'ADCPfilename','InputRow','DistanceFromLeftBank','D50_14', ...
    'BTmodf_S','Depth_S','S_S'};
for i = 1:numel(sources)
    S = sources{i};
    if ~isstruct(S) || isempty(fieldnames(S)), continue; end
    for j = 1:numel(preferred)
        f = preferred{j};
        if isfield(S,f)
            x = S.(f);
            nCandidate = inferLength(x);
            if nCandidate>nRows, nRows=nCandidate; end
        end
    end
end
if nRows<1, nRows=1; end
end

function n = inferLength(x)
n = 0;
if isempty(x), return; end
if istable(x), n=height(x); return; end
if iscell(x) || isstring(x) || ischar(x)
    if ischar(x), n=1; else, n=numel(x); end
    return
end
if isnumeric(x) || islogical(x)
    if isvector(x), n=numel(x); else, n=size(x,1); end
    return
end
if isstruct(x), n=numel(x); end
end

function out = mergeStructs(varargin)
out = struct();
for i = 1:nargin
    S = varargin{i};
    if isempty(S) || ~isstruct(S), continue; end
    if ~isscalar(S), continue; end
    f = fieldnames(S);
    for j = 1:numel(f)
        out.(f{j}) = S.(f{j});
    end
end
end

function names = allAliasNames(aliases)
names = {};
f = fieldnames(aliases);
for i=1:numel(f)
    a = aliases.(f{i});
    names = [names,a]; %#ok<AGROW>
end
end

function [target,report] = ensureFields(target,fields,aliases,sources, ...
    sourceNames,nRows,structureName)
n = numel(fields);
Structure = repmat(string(structureName),n,1);
Field = string(fields(:));
SourceContainer = strings(n,1);
SourceField = strings(n,1);
WasPlaceholder = false(n,1);

for i=1:n
    fieldName = fields{i};
    if isfield(target,fieldName) && ~isPlaceholderValue(target.(fieldName))
        SourceContainer(i) = "existing merged output";
        SourceField(i) = string(fieldName);
        continue
    end

    aliasList = {fieldName};
    if isfield(aliases,fieldName)
        aliasList = aliases.(fieldName);
    end

    [found,value,sourceContainer,sourceField] = resolveField( ...
        aliasList,sources,sourceNames);
    if found
        target.(fieldName) = value;
        SourceContainer(i) = string(sourceContainer);
        SourceField(i) = string(sourceField);
    else
        target.(fieldName) = makePlaceholder(fieldName,nRows);
        SourceContainer(i) = "placeholder";
        SourceField(i) = "missing";
        WasPlaceholder(i) = true;
    end
end

report = table(Structure,Field,SourceContainer,SourceField,WasPlaceholder);
end

function [found,value,sourceContainer,sourceField] = resolveField( ...
    aliases,sources,sourceNames)
found = false;
value = [];
sourceContainer = "";
sourceField = "";
for i=1:numel(sources)
    S = sources{i};
    if isempty(S) || ~isstruct(S) || ~isscalar(S), continue; end
    for j=1:numel(aliases)
        f = aliases{j};
        if isfield(S,f)
            candidate = S.(f);
            if isPlaceholderValue(candidate)
                continue
            end
            found = true;
            value = candidate;
            sourceContainer = sourceNames{i};
            sourceField = f;
            return
        end
    end
end
end

function value = makePlaceholder(fieldName,nRows)
cellFields = { ...
    'dataADCPS','dataADCPSboat','TKE_prof_S', ...
    'tau_turb_prof_S','uw_mean_S','BTerror_Sall', ...
    'ErrorVel_Svalid','alphaV_S','alphaBTf_S', ...
    'VmodALL_S','BTmodf_SALL','BTE_S','BTN_S', ...
    'BTx_S','BTy_S','BTx_Sa','BTy_Sa', ...
    'SV1_S','SV2_S','SV3_S','SV4_S', ...
    'SV11_S','SV22_S','SV33_S','SV44_S', ...
    'stdSpatialALL_S','stdSvSpatialALL_S', ...
    'BT_RSSI1_S','BT_RSSI2_S','BT_RSSI3_S','BT_RSSI4_S', ...
    'alpha1mean_S','alpha2mean_S','alpha3mean_S','alpha4mean_S', ...
    'dEDR1_S','dEDR2_S','dEDR3_S','dEDR4_S', ...
    'Pdwb_S','SL_S','XmitVV_S','XmitCC_S','PitchH_S','RollL_S', ...
    'BTfilterMask_S','BTvalidMask_S','BTx_final_S', ...
    'BTy_final_S','BTmod_final_S','BTvelocityFilterMask_S'};
stringFields = { ...
    'MesString_S','VelHomString1_S','VelHomString2_S', ...
    'CbFinalSource_S','DthicknessFinalSource_S', ...
    'Method2Status_S','BTvelocityMaskMode_S', ...
    'BTvelocityFilterSource_S','BTfilterMaskConvention'};
logicalFields = { ...
    'CriticalInputsConsistent_S','CbFinalUsedFallback_S', ...
    'DthicknessFinalUsedFallback_S','qsKinematicFinalValid_S', ...
    'BTvelocityFilterExact_S','BTvelocityFilterAligned_S'};

if any(strcmp(fieldName,cellFields))
    value = cell(nRows,1);
elseif any(strcmp(fieldName,stringFields)) || endsWith(fieldName,'Source_S') || ...
        endsWith(fieldName,'Status_S') || strcmp(fieldName,'Message')
    value = repmat("missing",nRows,1);
elseif any(strcmp(fieldName,logicalFields))
    value = false(nRows,1);
elseif endsWith(fieldName,'Table')
    value = table();
elseif strcmp(fieldName,'DDt_S') || strcmp(fieldName,'DthicknessCandidates_S')
    value = nan(nRows,7);
elseif strcmp(fieldName,'CbCandidates_S')
    value = nan(nRows,8);
elseif strcmp(fieldName,'D0BeamStd_S')
    value = nan(nRows,4);
elseif any(strcmp(fieldName,{'Depth_S','Depth_Sstd','BT_RSSI_S', ...
        'BT_RSSI_S','alphaMean_S'}))
    value = nan(nRows,4);
elseif any(strcmp(fieldName,{'PSD_percentage','QperDall'}))
    value = nan(nRows,0);
elseif startsWith(fieldName,'FileFound')
    value = false(nRows,1);
else
    value = nan(nRows,1);
end
end

function report = appendReports(A,B)
if isempty(A), report=B; return; end
if isempty(B), report=A; return; end
report = [A;B];
[~,idx] = unique(report.Field,'stable');
report = report(idx,:);
end

function report = updateReportResolution(report,fieldName,containerName,sourceField)
idx = find(report.Field==string(fieldName),1,'first');
if isempty(idx), return; end
report.SourceContainer(idx) = string(containerName);
report.SourceField(idx) = string(sourceField);
report.WasPlaceholder(idx) = false;
end

function tf = isPlaceholderValue(x)
tf = isempty(x);
if isnumeric(x)
    tf = isempty(x) || all(~isfinite(x(:)));
elseif isstring(x)
    tf = isempty(x) || all(x=="missing" | ismissing(x));
elseif iscell(x)
    if isempty(x)
        tf = true;
    else
        tf = all(cellfun(@isempty,x(:)));
    end
elseif istable(x)
    tf = isempty(x) || height(x)==0;
end
end

function y = rowMeanFinite(A,nRows)
A = double(A);
if isempty(A)
    y = nan(nRows,1);
elseif isvector(A)
    y = A(:);
    y = fitLength(y,nRows);
elseif size(A,1)==nRows
    y = finiteMean(A,2);
elseif size(A,2)==nRows
    y = finiteMean(A,1).';
else
    y = nan(nRows,1);
end
end

function y = finiteMean(A,dim)
valid = isfinite(A);
B = A;
B(~valid)=0;
n = sum(valid,dim);
y = sum(B,dim)./n;
y(n==0)=NaN;
end

function x = fitLength(x,n)
x=x(:);
if numel(x)==1 && n>1
    x=repmat(x,n,1);
elseif numel(x)<n
    x=[x;nan(n-numel(x),1)];
elseif numel(x)>n
    x=x(1:n);
end
end


function [dataExtra,shearExtra,btExtra,qExtra,ampExtra,unrouted] = ...
    routeSupplementalVariables(S,dataFields,shearFields,btFields,qFields,ampFields)
%ROUTESUPPLEMENTALVARIABLES Put additional values into existing structures.
% A field may intentionally be copied to more than one destination when the
% legacy reader expects the same value in two structures (for example
% D50_05 in Data_field and DataBTamp/BTamp).
dataExtra=struct(); shearExtra=struct(); btExtra=struct();
qExtra=struct(); ampExtra=struct(); unrouted={};
if isempty(S) || ~isstruct(S) || ~isscalar(S), return; end
f=fieldnames(S);
for i=1:numel(f)
    name=f{i}; value=S.(name); matched=false;
    if any(strcmp(name,dataFields))
        dataExtra.(name)=value; matched=true;
    end
    if any(strcmp(name,shearFields))
        shearExtra.(name)=value; matched=true;
    end
    if any(strcmp(name,btFields))
        btExtra.(name)=value; matched=true;
    end
    if any(strcmp(name,qFields))
        qExtra.(name)=value; matched=true;
    end
    if any(strcmp(name,ampFields))
        ampExtra.(name)=value; matched=true;
    end
    if matched, continue; end

    % Conservative prefix routing for genuinely new variables. More
    % specific acoustic and transport prefixes are checked first.
    if startsWith(name,'qs') || startsWith(name,'QvR')
        qExtra.(name)=value; matched=true;
    elseif startsWith(name,'BT_RSSI') || startsWith(name,'SV') || ...
            startsWith(name,'Spread') || startsWith(name,'TL') || ...
            startsWith(name,'Af') || startsWith(name,'dEDR') || ...
            startsWith(name,'Pdwb') || startsWith(name,'SL') || ...
            startsWith(name,'Xmit') || startsWith(name,'PitchH') || ...
            startsWith(name,'RollL') || startsWith(name,'stdSv') || ...
            startsWith(name,'stdSpatialRSSI')
        ampExtra.(name)=value; matched=true;
    elseif startsWith(name,'BTx') || startsWith(name,'BTy') || ...
            startsWith(name,'BTE') || startsWith(name,'BTN') || ...
            startsWith(name,'BTmod') || startsWith(name,'BTstd') || ...
            startsWith(name,'FilteredData') || startsWith(name,'NANs') || ...
            startsWith(name,'alphaBT') || startsWith(name,'alphaV') || ...
            startsWith(name,'WaterBias') || startsWith(name,'DeltaT') || ...
            startsWith(name,'BTfilter') || startsWith(name,'BTvalid')
        btExtra.(name)=value; matched=true;
    elseif startsWith(name,'Tau') || startsWith(name,'tau') || ...
            startsWith(name,'Ushear') || startsWith(name,'U_shear') || ...
            startsWith(name,'u_shear') || startsWith(name,'uShear') || ...
            startsWith(name,'Cb') || startsWith(name,'Tv') || ...
            startsWith(name,'Dv') || startsWith(name,'W_') || ...
            startsWith(name,'ucr') || startsWith(name,'Tucr') || ...
            startsWith(name,'Theta') || startsWith(name,'TKE') || ...
            startsWith(name,'Vmod') || startsWith(name,'HEightDUNE') || ...
            startsWith(name,'Fr_') || startsWith(name,'Tausk')
        shearExtra.(name)=value; matched=true;
    end
    if ~matched
        unrouted{end+1,1}=name; %#ok<AGROW>
    end
end
end

function Q = selectTransportFields(S)
Q = struct();
if isempty(S) || ~isstruct(S), return; end
f = fieldnames(S);
keyExact = {'Ninstrument','BTmodf_S','CbFinal_S','CbFinalSource_S', ...
    'CbFinalMethodCount_S','DthicknessFinal_S', ...
    'DthicknessFinalSource_S','DthicknessFinalMethodCount_S', ...
    'FinalTrateComparisonTable','TrateMeanComparison_S', ...
    'MaxTrateComparison_S','MeasuredTrateToKgFactor'};
for i=1:numel(f)
    name=f{i};
    if startsWith(name,'qs') || startsWith(name,'QvR') || any(strcmp(name,keyExact))
        Q.(name)=S.(name);
    end
end
end


function [ShearData,report] = resolveLegacyBeamStd( ...
    ShearData,ShearStressOut,DepthRNGout,report,nRows)
names = {'BT1_D_Sstd','BT2_D_Sstd', ...
    'BT3_D_Sstd','BT4_D_Sstd'};

% Preferred source: the four-column diagnostic already generated by the
% current shear function.
if isfield(ShearStressOut,'D0BeamStd_S')
    M = double(ShearStressOut.D0BeamStd_S);
    if size(M,1)==4 && size(M,2)==nRows, M=M.'; end
    if size(M,1)==nRows && size(M,2)>=4
        for b=1:4
            if isPlaceholderValue(ShearData.(names{b}))
                ShearData.(names{b}) = M(:,b);
                report = updateReportResolution(report,names{b}, ...
                    'ShearStressOut','D0BeamStd_S column');
            end
        end
    end
end

% Secondary source: calculate a temporal standard deviation from each
% processed beam-range cell. This does not change the range processing.
rawNames = { ...
    {'RNG11_S','RNG11_','RNG11'}, ...
    {'RNG22_S','RNG22_','RNG22'}, ...
    {'RNG33_S','RNG33_','RNG33'}, ...
    {'RNG44_S','RNG44_','RNG44'}};
for b=1:4
    if ~isPlaceholderValue(ShearData.(names{b})), continue; end
    [ok,v,src] = rowStdFromPossibleCells(DepthRNGout,rawNames{b},nRows);
    if ok
        ShearData.(names{b}) = v;
        report = updateReportResolution(report,names{b}, ...
            'DepthRNGout',src);
    end
end
end

function [ok,v,source] = rowStdFromPossibleCells(S,names,nRows)
ok=false; v=nan(nRows,1); source="";
if isempty(S) || ~isstruct(S), return; end
for k=1:numel(names)
    f=names{k};
    if ~isfield(S,f), continue; end
    C=S.(f);
    if iscell(C)
        vv=nan(nRows,1);
        for i=1:min(nRows,numel(C))
            x=double(C{i}(:)); x=x(isfinite(x));
            if numel(x)>=2, vv(i)=std(x,0); end
        end
        if any(isfinite(vv))
            ok=true; v=vv; source=string(f)+" (row temporal std)";
            return
        end
    end
end
end

function S = normalizeReaderStructure(S,fields,nRows,structureName)
if isempty(S) || ~isstruct(S), return; end

matrix4Fields = {'Depth_S','Depth_Sstd','BT_RSSI_S','alphaMean_S'};
cellRowFields = { ...
    'TKE_prof_S','SV1_S','SV2_S','SV3_S','SV4_S', ...
    'alpha1mean_S','alpha2mean_S','alpha3mean_S','alpha4mean_S'};

for k=1:numel(fields)
    f=fields{k};
    if ~isfield(S,f), continue; end
    x=S.(f);

    if any(strcmp(f,matrix4Fields))
        if isnumeric(x)
            if size(x,1)==4 && size(x,2)==nRows
                x=x.';
            elseif isvector(x) && numel(x)==4 && nRows==1
                x=reshape(x,1,4);
            end
        end
    elseif any(strcmp(f,cellRowFields))
        if ~iscell(x)
            if isnumeric(x) && nRows==1
                x={x};
            else
                x=cell(nRows,1);
            end
        else
            x=x(:);
            if numel(x)<nRows, x=[x;cell(nRows-numel(x),1)]; end
            if numel(x)>nRows, x=x(1:nRows); end
        end
    elseif isnumeric(x) || islogical(x) || isstring(x)
        if isvector(x)
            x=x(:);
            if numel(x)==1 && nRows>1
                x=repmat(x,nRows,1);
            elseif numel(x)<nRows
                if isnumeric(x), x=[x;nan(nRows-numel(x),1)];
                elseif islogical(x), x=[x;false(nRows-numel(x),1)];
                else, x=[x;repmat("missing",nRows-numel(x),1)];
                end
            elseif numel(x)>nRows
                x=x(1:nRows);
            end
        end
    end
    S.(f)=x;
end

% The reader concatenates these fields vertically, so enforce columns.
columnFields = {'LeftBankDist14','Ks2_S','D50_05','D50_14'};
if any(strcmp(structureName,{'Data_field','BTamp'}))
    for k=1:numel(columnFields)
        f=columnFields{k};
        if isfield(S,f) && isvector(S.(f))
            tmpValue=S.(f);
            S.(f)=tmpValue(:);
        end
    end
end
end

function T = makeReaderCompatibilityReport( ...
    Data_field,ShearData,BTvel,Q,BTamp, ...
    dataFields,shearFields,btFields,btampFields,nRows)
structures = {Data_field,ShearData,BTvel,BTamp};
structureNames = {'Data_field','ShearData','BTvel','BTamp'};
fieldGroups = {dataFields,shearFields,btFields,btampFields};

Structure=strings(0,1); Field=strings(0,1);
Exists=false(0,1); IsPlaceholder=false(0,1);
SizeText=strings(0,1); RowCompatible=false(0,1);

for g=1:numel(structures)
    S=structures{g}; fields=fieldGroups{g};
    for k=1:numel(fields)
        f=fields{k};
        Structure(end+1,1)=string(structureNames{g}); %#ok<AGROW>
        Field(end+1,1)=string(f); %#ok<AGROW>
        ex=isfield(S,f);
        Exists(end+1,1)=ex; %#ok<AGROW>
        if ex
            x=S.(f);
            IsPlaceholder(end+1,1)=isPlaceholderValue(x); %#ok<AGROW>
            sz=size(x);
            SizeText(end+1,1)=strjoin(string(sz),'x'); %#ok<AGROW>
            RowCompatible(end+1,1)=readerRowCompatible(x,nRows,f); %#ok<AGROW>
        else
            IsPlaceholder(end+1,1)=true; %#ok<AGROW>
            SizeText(end+1,1)="missing"; %#ok<AGROW>
            RowCompatible(end+1,1)=false; %#ok<AGROW>
        end
    end
end

% Q is a required top-level variable, but this reader does not access a Q
% field directly. Record that the structure itself exists.
Structure(end+1,1)="Q";
Field(end+1,1)="<top-level structure>";
Exists(end+1,1)=isstruct(Q);
IsPlaceholder(end+1,1)=~isstruct(Q);
SizeText(end+1,1)="struct";
RowCompatible(end+1,1)=isstruct(Q);

T=table(Structure,Field,Exists,IsPlaceholder,SizeText,RowCompatible);
end

function tf = readerRowCompatible(x,nRows,fieldName)
if isempty(x), tf=false; return; end
if istable(x), tf=true; return; end
if iscell(x) || isstring(x)
    tf=numel(x)==nRows;
    return
end
if isnumeric(x) || islogical(x)
    if any(strcmp(fieldName,{'Depth_S','Depth_Sstd','BT_RSSI_S','alphaMean_S'}))
        tf=size(x,1)==nRows;
    elseif isvector(x)
        tf=numel(x)==nRows;
    else
        tf=size(x,1)==nRows;
    end
    return
end
tf=true;
end

function [datePrefix,dateSource] = resolveSourceDatePrefix( ...
    override,FolderM9,FolderSontek,dataADCPS,ADCPin)
datePrefix=""; dateSource="";

if strlength(string(override))>0
    token = extractDateToken(string(override));
    if token==""
        error('FileDatePrefix must contain a valid YYYYMMDD date.');
    end
    datePrefix=token;
    dateSource="FileDatePrefix";
    return
end

pathCandidates = [string(FolderM9); string(FolderSontek)];
pathLabels = ["FolderM9";"FolderSontek"];
for i=1:numel(pathCandidates)
    txt=strtrim(pathCandidates(i));
    if txt=="", continue; end
    [token,part]=extractDateFromPathEnd(txt);
    if token~=""
        datePrefix=token;
        dateSource=pathLabels(i)+":"+part;
        return
    end
end

if isstruct(dataADCPS) && isfield(dataADCPS,'name')
    names=string({dataADCPS.name});
    for i=1:numel(names)
        token=extractDateToken(names(i));
        if token~=""
            datePrefix=token;
            dateSource="dataADCPS.name:"+names(i);
            return
        end
    end
end

if isstruct(ADCPin) && isfield(ADCPin,'ADCPfilename')
    names=string(ADCPin.ADCPfilename(:));
    for i=1:numel(names)
        token=extractDateToken(names(i));
        if token~=""
            datePrefix=token;
            dateSource="ADCPin.ADCPfilename:"+names(i);
            return
        end
    end
end
end

function [token,matchedPart] = extractDateFromPathEnd(pathText)
token=""; matchedPart="";
parts=regexp(char(pathText),'[\\/]+','split');
for i=numel(parts):-1:1
    if isempty(parts{i}), continue; end
    candidate=extractDateToken(string(parts{i}));
    if candidate~=""
        token=candidate;
        matchedPart=string(parts{i});
        return
    end
end
end

function token = extractDateToken(textValue)
token="";
textValue=char(string(textValue));
m=regexp(textValue,'(?<!\d)((?:19|20)\d{6})(?!\d)','tokens','once');
if isempty(m), return; end
candidate=string(m{1});
try
    datetime(candidate,'InputFormat','yyyyMMdd');
    token=candidate;
catch
    token="";
end
end

function filePath = resolveOutputFilePath(explicitFile,destinationFolder,datePrefix,locationToken,defaultSuffix)
if strlength(string(explicitFile))>0
    filePath=resolveExplicitSavePath(explicitFile,datePrefix,locationToken);
elseif strlength(string(destinationFolder))>0
    prefix=composeFileNamePrefix(datePrefix,locationToken);
    filePath=string(fullfile(char(string(destinationFolder)), ...
        char(prefix+"_"+string(defaultSuffix))));
else
    filePath="";
end
end

function filePath = resolveExplicitSavePath(explicitFile,datePrefix,locationToken)
if strlength(string(explicitFile))==0
    filePath="";
    return
end
filePath=string(explicitFile);
[parent,base,ext]=fileparts(char(filePath));
if isempty(ext), ext='.mat'; end
baseString=applyFileNamePrefix(string(base),datePrefix,locationToken);
filePath=string(fullfile(parent,char(baseString+string(ext))));
end

function prefix = composeFileNamePrefix(datePrefix,locationToken)
datePrefix=string(datePrefix);
locationToken=string(locationToken);
if datePrefix==""
    prefix="";
elseif locationToken==""
    prefix=datePrefix;
else
    prefix=datePrefix+"_"+locationToken;
end
end

function baseString = applyFileNamePrefix(baseString,datePrefix,locationToken)
baseString=string(baseString);
datePrefix=string(datePrefix);
locationToken=string(locationToken);
prefix=composeFileNamePrefix(datePrefix,locationToken);
if prefix==""
    return
end
if startsWith(baseString,prefix+"_") || baseString==prefix
    return
end
if locationToken~="" && startsWith(baseString,datePrefix+"_")
    tail=extractAfter(baseString,strlength(datePrefix)+1);
    baseString=prefix+"_"+tail;
elseif baseString==datePrefix
    baseString=prefix;
else
    baseString=prefix+"_"+baseString;
end
end

function [locationName,locationToken] = normalizeLocationName(locationValue)
locationName=strtrim(string(locationValue));
locationName=erase(locationName,"'");
locationName=erase(locationName,'"');
if ismissing(locationName) || lower(locationName)=="nan"
    locationName="";
end
if locationName==""
    locationToken="";
    return
end
locationToken=string(regexprep(char(locationName),'[^A-Za-z0-9]+','_'));
locationToken=string(regexprep(char(locationToken),'^_+|_+$',''));
locationToken=string(regexprep(char(locationToken),'_+','_'));
if locationToken==""
    error('LocationName must contain at least one letter or number.');
end
end

function fileName = prepareSavePath(fileName)
fileName = char(string(fileName));
[parent,~,~] = fileparts(fileName);
if ~isempty(parent) && exist(parent,'dir')~=7
    mkdir(parent);
end
end

function saveWithVersion(fileName,saveVersion,name1,value1)
S=struct(); S.(name1)=value1;
save(fileName,'-struct','S',saveVersion);
end

function saveTwoWithVersion(fileName,saveVersion,name1,value1,name2,value2)
S=struct(); S.(name1)=value1; S.(name2)=value2;
save(fileName,'-struct','S',saveVersion);
end

function saveThreeWithVersion(fileName,saveVersion,name1,value1,name2,value2,name3,value3)
S=struct(); S.(name1)=value1; S.(name2)=value2; S.(name3)=value3;
save(fileName,'-struct','S',saveVersion);
end

function saveCombined(fileName,saveVersion,Data_field,ShearData,BTvel,Q,BTamp, ...
    DataBTamp,StructureReport)
S=struct();
S.Data_field=Data_field;
S.ShearData=ShearData;
S.BTvel=BTvel;
S.Q=Q;
S.BTamp=BTamp;
S.DataBTamp=DataBTamp;
S.StructureReport=StructureReport;
save(fileName,'-struct','S',saveVersion);
end
