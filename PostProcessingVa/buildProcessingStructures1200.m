function [Data_field, ShearData, BTvel, Q, BTamp, StructureReport] = ...
    buildProcessingStructures1200(ADCPin, WaterProfileOut, DepthRNGout, ...
    ShearStressOut, BTvelOut, BedloadTransportOut, BackscatterOut, varargin)
%BUILDPROCESSINGSTRUCTURES1200 Build the legacy project output structures.
%
% [Data_field, ShearData, BTvel, Q, BTamp, StructureReport] = ...
%     buildProcessingStructures1200(ADCPin, WaterProfileOut, DepthRNGout, ...
%     ShearStressOut, BTvelOut, BedloadTransportOut, BackscatterOut, ...
%     'Name',Value,...)
%
% PURPOSE
% -------
% This is the final packaging step after the modular 1200-kHz processing
% functions. It creates the structures used by the historical multi-river
% analysis scripts while retaining the original field names:
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
% WaterProfileOut     output from processWaterProfile1200, or []
% DepthRNGout         output from processDepthsRNG, or []
% ShearStressOut      output from processShearStress1200, or []
% BTvelOut            output from processBTbedloadVelocity1200, or []
% BedloadTransportOut output from processBedloadTrate, or []
% BackscatterOut      output from processBackscatter1200, or []
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
% 'Folder1200'          source folder used to extract the YYYYMMDD prefix
% 'FolderTrates'        destination folder for Data_field, BTvel and Q
% 'FolderShear'         destination folder for ShearData
% 'FolderBS'            destination folder for BTamp and DataBTamp
% 'FileDatePrefix'      optional YYYYMMDD override; normally extracted from
%                       Folder1200, FolderSontek or source filenames
% 'LocationName'        optional location name inserted after the date and
%                       underscore, e.g. 20170322_Tomorrica_TratesData.mat
% 'dataADCP1200'        raw dir-style file list to store in Data_field and
%                       ShearData
% 'dataADCP1200boat'    optional boat-file list
% 'Acoustic'            output structure from computeADCPacousticParameters.
%                       Acoustic fields are flattened directly into BTamp;
%                       no nested Acoustic structure is saved.
% 'ExtraVariables'      backward-compatible supplemental-variable structure.
%                       Its fields are routed into Data_field, ShearData,
%                       BTvel, Q or BTamp; it is never saved as a separate
%                       output structure. Direct name-value variables are
%                       preferred, for example 'DeltaT1200',DeltaT1200.
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
% Any unrecognized name-value pair is treated as a supplemental variable
% and routed to its existing destination structure by its field name. No
% supplemental structure is created or saved.
% 'Verbose'             print a summary; default true
%
% RECOMMENDED CALL
% ----------------
% [Data_field,ShearData,BTvel,Q,BTamp,StructureReport] = ...
%     buildProcessingStructures1200( ...
%         PhSamplesin,WaterProfileOut,DepthRNGout,ShearStressOut, ...
%         BTvelOut,BedloadTransportOut,BackscatterOut, ...
%         'FolderSontek',FolderSontek, ...
%         'Folder1200',Folder1200, ...
%         'FolderTrates',FolderTrates, ...
%         'FolderShear',FolderShear, ...
%         'FolderBS',FolderBS, ...
%         'LocationName','Tomorrica', ...
%         'dataADCP1200',dataADCP1200, ...
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
% No existing field name is renamed or removed.

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
addParameter(p,'Folder1200',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderTrates',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderShear',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FolderBS',"",@(x) ischar(x) || isstring(x));
addParameter(p,'FileDatePrefix',"",@(x) ischar(x) || isstring(x));
addParameter(p,'LocationName',"",@(x) ischar(x) || isstring(x));
addParameter(p,'dataADCP1200',[],@(x) true);
addParameter(p,'dataADCP1200boat',[],@(x) true);
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
    'Rep','vs','D_','dataADCP1200','dataADCP1200boat', ...
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
    'WT1200','dataADCP1200','VmodAll_1200','Vmod_1200stdMOD', ...
    'Vmod1low_1200stdMOD','u_shearH1200','Tau_All1200_H', ...
    'Tau_All1200','CbH_1200','TvH_1200','DvH_1200', ...
    'DthicknessH_1200','QvRH_1200','ucr_1200','Fr_1200', ...
    'Tausk_1200','Tucr_1200','Ushear1_1200','Tau_All1200_1', ...
    'Cb1_1200','Tv1_1200','Dv1_1200','Tv2_1200','Dv2_1200', ...
    'DDt1200','Ushear2_1200','Tau_All1200_2','Cb2_1200', ...
    'nCb_new1200','UsheaAll_1200','Cb_1200','HEightDUNE1200', ...
    'Depth_1200','BT1_D_1200std','BT2_D_1200std', ...
    'BT3_D_1200std','BT4_D_1200std', ...
    'U_shearRE1_1200','U_shearRE2_1200', ...
    'Tau_AllRE1_1200','Tau_AllRE2_1200','WT_1200', ...
    'tau_Mann1200','Tau_shear2VR_1200','tau_Mann_star1200', ...
    'tau_rm1200','uShear_Mann1200','tau_rm_starGaumann1200', ...
    'tau_rm_star1200','W_star_perDall1200','tau_Mann1200vsW', ...
    'tau_Mann1200vsWstar','W_i1200','W_iVRshear2_1200', ...
    'tauRef1200','tauRefVRs1200','tauRef1200VRstar', ...
    'tauRef1200star','tauRef1200Perfrac', ...
    'tau_turb_prof1200','TKE_prof1200','uw_mean1200'};

btLegacyFields = { ...
    'BTE_1200','BTN_1200','BT_x_1200','BTx_1200f','BTx_1200', ...
    'BT_y_1200','BTy_1200f','BTy_1200','BTx_1200a','BTy_1200a', ...
    'BTerror1200m','BTerror1200all','BTstd1200a','BTmod1200a', ...
    'BTmod1200aW','BTstd1200aW','BTmodM1_12','BTmod1200', ...
    'BTmodf_1200','FilteredDataN_1200','FilteredDataE_1200', ...
    'NANsbeforeX1200perc','NANsbeforeY1200perc','BTstd_1200', ...
    'BT_stdf_1200','ErrorVel1200valid','ErrorV1200', ...
    'alphaV_1200','alphaBTf_1200','VmodALL_1200','BTmodf_1200ALL', ...
    'Diff_1200','maxDiff_1200','meanDiff_1200','MeanDepth1200', ...
    'beamDiff1200','meanDiffbeams_1200','NNinhomogenious1200', ...
    'MesString1200','CV1_1200','CV2_1200', ...
    'VelHomString1_1200','VelHomString2_1200','DeltaT1200'};

qLegacyFields = { ...
    'qsNew_1200','qsNew_1200stdBT','qs_1200','qs_1200stdBT', ...
    'qs1_1200','qs1_1200stdBT','qs_1200Lambda', ...
    'qs_1200stdBTLambda','qs2_1200stdBT','qs2_1200', ...
    'qsH_1200stdBT','qsH_1200','qs1_1200DD', ...
    'qsMethod2regression_1200'};

btampLegacyFields = { ...
    'Depth1200','stdSpatial_1200','Depth1200std','stdSpatialALL_1200', ...
    'BT_RSSI_1200','alphaMean1200', ...
    'alpha1mean1200grad','alpha2mean1200grad', ...
    'alpha3mean1200grad','alpha4mean1200grad', ...
    'alpha1mean1200','alpha2mean1200','alpha3mean1200','alpha4mean1200', ...
    'dEDR1_1200','dEDR2_1200','dEDR3_1200','dEDR4_1200', ...
    'R1_1200','R2_1200','R3_1200','R4_1200', ...
    'Spread1_1200','Spread2_1200','Spread3_1200','Spread4_1200', ...
    'SpreadM1200','SVBmean1_1200','SVBmean2_1200', ...
    'SVBmean3_1200','SVBmean4_1200', ...
    'SVB1mean1_1200','SVB1mean2_1200','SVB1mean3_1200','SVB1mean4_1200', ...
    'SVBstd1_1200','SVBstd2_1200','SVBstd3_1200','SVBstd4_1200', ...
    'SVB1std1_1200','SVB1std2_1200','SVB1std3_1200','SVB1std4_1200', ...
    'SV11_1200','SV22_1200','SV33_1200','SV44_1200', ...
    'SV1_1200','SV2_1200','SV3_1200','SV4_1200', ...
    'stdSvSpatial_1200','stdSvSpatialALL_1200', ...
    'BT_RSSI1_1200','BT_RSSI2_1200','BT_RSSI3_1200','BT_RSSI4_1200', ...
    'stdSpatialRSSI_1200','TL1_1200','TL2_1200','TL3_1200','TL4_1200', ...
    'TL1_fix1200','TL2_fix1200','TL3_fix1200','TL4_fix1200', ...
    'Dmax_1200_1','Dmax_1200_2','Dmax_1200_3','Dmax_1200_4', ...
    'Area_1200_1','Area_1200_2','Area_1200_3','Area_1200_4', ...
    'Af1mean1200','Af2mean1200','Af3mean1200','Af4mean1200', ...
    'distance_1200_1','distance_1200_2','distance_1200_3','distance_1200_4', ...
    'Lp_1200_1','Lp_1200_2','Lp_1200_3','Lp_1200_4', ...
    'Pdwb1200','SL1200','STDPowerdb1200','STDPW1200', ...
    'XmitVV1200','XmitCC1200','PitchH1200','RollL1200', ...
    'Ks2_12','D50_05','D50_14'};


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
    'Ks2_12','D50_05', ...
    'Depth1200','Depth1200std', ...
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

% Explicitly retained modern variables that are essential for interpreting
% the final methods even when an older version of an output structure is used.
newShearFields = { ...
    'D50critical_mm','D50critical_m','D50mobile','Rep_1200','Dstar_1200', ...
    'ThetaCr_1200','ThetaCrVanRijn_1200','ThetaCrBrownlie_1200', ...
    'ThetaCrSoulsby_1200','CriticalInputsConsistent_1200', ...
    'CbCandidates_1200','CbMethodNames_1200','CbFinal_1200', ...
    'CbFinalMethodCount_1200','CbFinalSource_1200','CbFinalUsedFallback_1200', ...
    'DthicknessCandidates_1200','DthicknessMethodNames_1200', ...
    'DthicknessFinal_1200','DthicknessFinalMethodCount_1200', ...
    'DthicknessFinalSource_1200','DthicknessFinalUsedFallback_1200', ...
    'DthicknessKin_1200','Ninstrument','uShear_skin1200', ...
    'theta_skin_1200','tau_skin_1200'};

newBTFields = { ...
    'BTfilterMask_1200','BTvalidMask_1200','BTfilterMaskLength_1200', ...
    'NfilteredFinal_1200','NvalidFinal_1200','BTfilterMaskConvention', ...
    'BTfilterMaskStart_1200','BTfilterMaskEnd_1200', ...
    'BTx_final_1200','BTy_final_1200','BTmod_final_1200','va_final_1200'};

newQFields = { ...
    'qsKinematicFinal_1200','qsFinalKinematic_1200','qsFinal_1200', ...
    'qs_final_1200','qsKinematicFinalValid_1200', ...
    'qsRE1_1200','qsRE2_1200','qsMann_1200','qsSkin_1200', ...
    'qsVR_1200','qsVRRE1_1200','qsVRRE2_1200','qsVRH_1200', ...
    'qsVR1_1200','qsVR2_1200','qsVRMann_1200','qsVRSkin_1200', ...
    'QvR_1200','QvRRE1_1200','QvRRE2_1200','QvRH_1200', ...
    'QvR1_1200','QvR2_1200','QvRMann_1200','QvRSkin_1200', ...
    'CbFinal_1200','CbFinalSource_1200','DthicknessFinal_1200', ...
    'DthicknessFinalSource_1200','BTmodf_1200','Ninstrument', ...
    'FinalTrateComparisonTable'};

newBTampFields = { ...
    'S_1200','qsMethod2regression_1200','Method2Status_1200', ...
    'BTvelocityFilterMask_1200','BTvelocityMaskMode_1200', ...
    'BTvelocityFilterSource_1200','BTvelocityFilterExact_1200', ...
    'BTvelocityFilterAligned_1200','NfilteredBTvelocity_1200', ...
    'StartUsed_1200','EndUsed_1200','NensemblesUsed_1200', ...
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
aliases.Rep = {'Rep','Rep_1200'};
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
aliases.Temp1 = {'Temp1','WT_12','WT1200','WT_1200'};
aliases.alpha1mean1200 = {'alpha1mean1200','alpha1mean'};
aliases.alpha2mean1200 = {'alpha2mean1200','alpha2mean'};
aliases.alpha3mean1200 = {'alpha3mean1200','alpha3mean'};
aliases.alpha4mean1200 = {'alpha4mean1200','alpha4mean'};


aliases.MaxTrate = {'MaxTrate','TrateMax','MTrateMax','MaxTrateAll'};
aliases.MinTrate = {'MinTrate','TrateMin','MTrateMin','MinTrateAll'};
aliases.Depth_1200 = {'Depth_1200','Depth1200'};
aliases.DeltaT1200 = {'DeltaT1200','DeltaT_1200','DeltaT','dt1200','dt_1200'};
aliases.Ks2_12 = {'Ks2_12','Ks2_1200','ks2_12','Ks_12'};
aliases.Vmod_1200stdMOD = {'Vmod_1200stdMOD','VmodAll_1200std','VmodAll1200std'};
aliases.Vmod1low_1200stdMOD = {'Vmod1low_1200stdMOD','Vmod1low1200stdMOD'};

% ShearData aliases.
aliases.WT1200 = {'WT1200','WT_1200','WT_12','Temp1'};
aliases.WT_1200 = {'WT_1200','WT1200','WT_12','Temp1'};
aliases.BT1_D_1200std = {'BT1_D_1200std','RNG11std_1200','RNG11_1200std'};
aliases.BT2_D_1200std = {'BT2_D_1200std','RNG22std_1200','RNG22_1200std'};
aliases.BT3_D_1200std = {'BT3_D_1200std','RNG33std_1200','RNG33_1200std'};
aliases.BT4_D_1200std = {'BT4_D_1200std','RNG44std_1200','RNG44_1200std'};

% BTvel aliases.
aliases.BTmod1200 = {'BTmod1200','BTmod1200org'};
aliases.MeanDepth1200 = {'MeanDepth1200','MeanDepth_1200'};

% Q aliases.
aliases.qsNew_1200 = {'qsNew_1200','qsCbNew_1200'};
aliases.qs1_1200DD = {'qs1_1200DD','qsVR1_1200'};

%% Special sources and caller-workspace capture
special = struct();
if ~isempty(opt.dataADCP1200)
    special.dataADCP1200 = opt.dataADCP1200;
end
if ~isempty(opt.dataADCP1200boat)
    special.dataADCP1200boat = opt.dataADCP1200boat;
end
special.FolderSontek = string(opt.FolderSontek);
special.Folder1200 = string(opt.Folder1200);

allRequested = unique([dataLegacyFields,shearLegacyFields,btLegacyFields, ...
    qLegacyFields,btampLegacyFields,newShearFields,newBTFields, ...
    newQFields,newBTampFields,featureDataFields,featureShearFields, ...
    featureBTvelFields,featureBTampFields,acousticFieldCandidates, ...
    allAliasNames(aliases)],'stable');
callerVariables = struct();
if logical(opt.UseCallerWorkspace)
    % This lookup is intentionally performed in the main function. Here,
    % evalin('caller',...) refers to the workspace that called
    % buildProcessingStructures1200 rather than to a local helper.
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
% D0BeamStd_1200. Use those columns when the individual legacy names were
% not available.
[ShearData,shearReport] = resolveLegacyBeamStd( ...
    ShearData,ShearStressOut,DepthRNGout,shearReport,nRows);

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

% Derive MeanDepth1200 only when it was not supplied meaningfully.
if isPlaceholderValue(BTvel.MeanDepth1200) && isfield(DepthRNGout,'Depth1200')
    BTvel.MeanDepth1200 = rowMeanFinite(DepthRNGout.Depth1200,nRows);
    btReport = updateReportResolution(btReport,'MeanDepth1200', ...
        'DepthRNGout','Depth1200 (row mean)');
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
StructureReport.Folder1200 = string(opt.Folder1200);
StructureReport.FolderTrates = string(opt.FolderTrates);
StructureReport.FolderShear = string(opt.FolderShear);
StructureReport.FolderBS = string(opt.FolderBS);
StructureReport.Data_field = dataReport;
StructureReport.ShearData = shearReport;
StructureReport.BTvel = btReport;
StructureReport.Q = qReport;
StructureReport.BTamp = btampReport;
StructureReport.AcousticToBTamp = acousticToBTampReport;
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
StructureReport.MissingFields = unique([ ...
    string(dataReport.Field(dataReport.WasPlaceholder)); ...
    string(shearReport.Field(shearReport.WasPlaceholder)); ...
    string(btReport.Field(btReport.WasPlaceholder)); ...
    string(qReport.Field(qReport.WasPlaceholder)); ...
    string(btampReport.Field(btampReport.WasPlaceholder))]);
StructureReport.Notes = [ ...
    "All legacy field names from the supplied saving script are created."; ...
    "All current fields from the modern output structures are retained."; ...
    "PSD_percentage and QperDall use the legacy flipud orientation when FlipPSDFields=true."; ...
    "Fr_1200 is taken from the actual Fr_1200 output; the old accidental assignment Fr_1200=Tausk_1200 is not repeated."; ...
    "Tau_shear2VR_1200 is not inferred from another stress variable because its legacy dimensional definition is ambiguous."; ...
    "Ks2_12 is stored directly in BTamp/DataBTamp when supplied directly or found in the caller workspace."; ...
    "Legacy transport uncertainty fields ending in stdBT are stored directly in Q when supplied; otherwise they remain placeholders."; ...
    "DataBTamp is created as an exact compatibility alias of BTamp when assigning or saving."; ...
    "Every field from computeADCPacousticParameters output Acoustic is flattened directly into BTamp and DataBTamp."; ...
    "When Acoustic is omitted, recognized acoustic fields are recovered from the updated ADCPin when available."; ...
    "Legacy Ks2_12 is resolved from Ks_12 when that is the name returned by the acoustic calculation."; ...
    "The Trates MAT file contains Data_field, BTvel and Q; the Shear MAT file contains only ShearData."; ...
    "All automatically generated filenames start with the same YYYYMMDD source-date prefix."; ...
    "When LocationName is supplied, it is sanitized and inserted after the date in all saved filenames."; ...
    "ReaderCompatibility lists every field referenced directly by the supplied feature-extraction script."; ...
    "Direct supplemental variables are routed into the existing five structures; no ExtraVariables structure is assigned or saved."; ...
    "DeltaT1200 routes to BTvel; Ks2_12 routes to BTamp/DataBTamp; Tau_shear2VR_1200, Vmod_1200stdMOD and Vmod1low_1200stdMOD route to ShearData."; ...
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
    opt.FileDatePrefix,opt.Folder1200,opt.FolderSontek,opt.dataADCP1200,ADCPin);

[LocationName,LocationToken] = normalizeLocationName(opt.LocationName);
FileNamePrefix = composeFileNamePrefix(FileDatePrefix,LocationToken);

if needsDatePrefix && FileDatePrefix == ""
    error(['Could not extract an 8-digit source date from Folder1200, ', ...
        'FolderSontek, dataADCP1200 or ADCPin.ADCPfilename. ', ...
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
    fprintf('Acoustic fields copied directly into BTamp/DataBTamp: %d\n', ...
        numel(fieldnames(Acoustic)));
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
function names = defaultAcousticFieldNames()
% Fields produced by computeADCPacousticParameters and earlier input-table
% acoustic calculations. All names are retained unchanged in BTamp.
names = { ...
    'c','ro','rs','k','f3','f1','f05','f12','f06', ...
    'lambda3','k3','x3','F3','Kappa3','FMean3','Ks_3', ...
    'lambda1','k1','x1','F1','Kappa1','FMean1','Ks_1', ...
    'lambda05','k05','x05','F05','Kappa05','FMean05','Ks_05', ...
    'lambda12','k12','x12','F12','Kappa12','FMean12','Ks_12', ...
    'lambda06','k06','x06','F06','Kappa06','FMean06','Ks_06', ...
    'psi06','psi12','DmeaNUMpsd'};
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

function nRows = inferRowCount(sources)
nRows = 0;
preferred = {'ADCPfilename','InputRow','DistanceFromLeftBank','D50_14', ...
    'BTmodf_1200','Depth_1200','S_1200'};
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
    'dataADCP1200','dataADCP1200boat','TKE_prof1200', ...
    'tau_turb_prof1200','uw_mean1200','BTerror1200all', ...
    'ErrorVel1200valid','alphaV_1200','alphaBTf_1200', ...
    'VmodALL_1200','BTmodf_1200ALL','BTE_1200','BTN_1200', ...
    'BTx_1200','BTy_1200','BTx_1200a','BTy_1200a', ...
    'SV1_1200','SV2_1200','SV3_1200','SV4_1200', ...
    'SV11_1200','SV22_1200','SV33_1200','SV44_1200', ...
    'stdSpatialALL_1200','stdSvSpatialALL_1200', ...
    'BT_RSSI1_1200','BT_RSSI2_1200','BT_RSSI3_1200','BT_RSSI4_1200', ...
    'alpha1mean1200','alpha2mean1200','alpha3mean1200','alpha4mean1200', ...
    'dEDR1_1200','dEDR2_1200','dEDR3_1200','dEDR4_1200', ...
    'Pdwb1200','SL1200','XmitVV1200','XmitCC1200','PitchH1200','RollL1200', ...
    'BTfilterMask_1200','BTvalidMask_1200','BTx_final_1200', ...
    'BTy_final_1200','BTmod_final_1200','BTvelocityFilterMask_1200'};
stringFields = { ...
    'MesString1200','VelHomString1_1200','VelHomString2_1200', ...
    'CbFinalSource_1200','DthicknessFinalSource_1200', ...
    'Method2Status_1200','BTvelocityMaskMode_1200', ...
    'BTvelocityFilterSource_1200','BTfilterMaskConvention'};
logicalFields = { ...
    'CriticalInputsConsistent_1200','CbFinalUsedFallback_1200', ...
    'DthicknessFinalUsedFallback_1200','qsKinematicFinalValid_1200', ...
    'BTvelocityFilterExact_1200','BTvelocityFilterAligned_1200'};

if any(strcmp(fieldName,cellFields))
    value = cell(nRows,1);
elseif any(strcmp(fieldName,stringFields)) || endsWith(fieldName,'Source_1200') || ...
        endsWith(fieldName,'Status_1200') || strcmp(fieldName,'Message')
    value = repmat("missing",nRows,1);
elseif any(strcmp(fieldName,logicalFields))
    value = false(nRows,1);
elseif endsWith(fieldName,'Table')
    value = table();
elseif strcmp(fieldName,'DDt1200') || strcmp(fieldName,'DthicknessCandidates_1200')
    value = nan(nRows,7);
elseif strcmp(fieldName,'CbCandidates_1200')
    value = nan(nRows,8);
elseif strcmp(fieldName,'D0BeamStd_1200')
    value = nan(nRows,4);
elseif any(strcmp(fieldName,{'Depth1200','Depth1200std','BT_RSSI_1200', ...
        'BT_RSSI1200','alphaMean1200'}))
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
keyExact = {'Ninstrument','BTmodf_1200','CbFinal_1200','CbFinalSource_1200', ...
    'CbFinalMethodCount_1200','DthicknessFinal_1200', ...
    'DthicknessFinalSource_1200','DthicknessFinalMethodCount_1200', ...
    'FinalTrateComparisonTable','TrateMeanComparison_1200', ...
    'MaxTrateComparison_1200','MeasuredTrateToKgFactor'};
for i=1:numel(f)
    name=f{i};
    if startsWith(name,'qs') || startsWith(name,'QvR') || any(strcmp(name,keyExact))
        Q.(name)=S.(name);
    end
end
end


function [ShearData,report] = resolveLegacyBeamStd( ...
    ShearData,ShearStressOut,DepthRNGout,report,nRows)
names = {'BT1_D_1200std','BT2_D_1200std', ...
    'BT3_D_1200std','BT4_D_1200std'};

% Preferred source: the four-column diagnostic already generated by the
% current shear function.
if isfield(ShearStressOut,'D0BeamStd_1200')
    M = double(ShearStressOut.D0BeamStd_1200);
    if size(M,1)==4 && size(M,2)==nRows, M=M.'; end
    if size(M,1)==nRows && size(M,2)>=4
        for b=1:4
            if isPlaceholderValue(ShearData.(names{b}))
                ShearData.(names{b}) = M(:,b);
                report = updateReportResolution(report,names{b}, ...
                    'ShearStressOut','D0BeamStd_1200 column');
            end
        end
    end
end

% Secondary source: calculate a temporal standard deviation from each
% processed beam-range cell. This does not change the range processing.
rawNames = { ...
    {'RNG11_1200','RNG11_','RNG11'}, ...
    {'RNG22_1200','RNG22_','RNG22'}, ...
    {'RNG33_1200','RNG33_','RNG33'}, ...
    {'RNG44_1200','RNG44_','RNG44'}};
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

matrix4Fields = {'Depth1200','Depth1200std','BT_RSSI_1200','alphaMean1200'};
cellRowFields = { ...
    'TKE_prof1200','SV1_1200','SV2_1200','SV3_1200','SV4_1200', ...
    'alpha1mean1200','alpha2mean1200','alpha3mean1200','alpha4mean1200'};

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
columnFields = {'LeftBankDist14','Ks2_12','D50_05','D50_14'};
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
    if any(strcmp(fieldName,{'Depth1200','Depth1200std','BT_RSSI_1200','alphaMean1200'}))
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
    override,Folder1200,FolderSontek,dataADCP1200,ADCPin)
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

pathCandidates = [string(Folder1200); string(FolderSontek)];
pathLabels = ["Folder1200";"FolderSontek"];
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

if isstruct(dataADCP1200) && isfield(dataADCP1200,'name')
    names=string({dataADCP1200.name});
    for i=1:numel(names)
        token=extractDateToken(names(i));
        if token~=""
            datePrefix=token;
            dateSource="dataADCP1200.name:"+names(i);
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
