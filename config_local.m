function cfg = config_local()
%CONFIG_LOCAL Local configuration for the RDI 1200 kHz ADCP workflow.
%
% This file contains computer-specific paths and project settings.
% Keep it in the repository root when running the workflow locally.
% DO NOT commit this file to GitHub. It should be listed in .gitignore.
%
% The public repository should contain config_example.m instead.

%% Project identity
cfg.mainDir      = 'C:\Users\slavenc';
cfg.locationName = 'Tomorrica';
cfg.datasetDate  = '20161202';   % Documentation only; filenames use the date extracted from folder1200.
cfg.ADCP_used    = 'RDI1200';

%% Input folders
cfg.phSamplesDir = fullfile(cfg.mainDir, ...
    'OneDrive - NTNU', 'Slaven_PhD', 'RiverPro_BedLoadAlbania', ...
    'DataForprocessing2026', '20161202_Tomorrice_Albania', ...
    'inputtable');

cfg.folder1200 = fullfile(cfg.mainDir, ...
    'OneDrive - NTNU', 'Slaven_PhD', 'RiverPro_BedLoadAlbania', ...
    'DataForprocessing2026', '20161202_Tomorrice_Albania', ...
    '20161202Dec2', 'processedADCPdata');

%% Output folders
% buildProcessingStructures1200 creates these folders when necessary.
cfg.folderShear = fullfile(cfg.mainDir, ...
    'OneDrive - NTNU', 'Slaven_PhD', 'FieldDataProcessing', ...
    'DATA', 'shearDATAnewFilters21boat_T');

cfg.folderTrates = fullfile(cfg.mainDir, ...
    'OneDrive - NTNU', 'Slaven_PhD', 'FieldDataProcessing', ...
    'DATA', 'TransportDATnewFilters21boat_T');

cfg.folderBS = fullfile(cfg.mainDir, ...
    'OneDrive - NTNU', 'Slaven_PhD', 'FieldDataProcessing', ...
    'DATA', 'BSstrengtgh_VBfilt21boat_T');

%% Instrument geometry and acoustic frequency
cfg.submergeDistance = 0.05;    % m, transducer below water surface
cfg.f12               = 1.0e6;  % Hz, nominal RDI frequency
cfg.beamAngleDeg      = 20;     % degrees
cfg.lowerFraction     = 0.50;   % fraction of water profile used for lower-profile processing

%% Fluid and sediment properties
cfg.ro         = 1000;  % kg/m^3, water density
cfg.rs         = 2650;  % kg/m^3, sediment density
cfg.vonKarman  = 0.40;  % von Karman constant
cfg.defaultNi  = 1.2e-6; % m^2/s, fallback kinematic viscosity

%% Hydraulic parameters
% IMPORTANT: processShearStress1200 converts SO only when SOIsPermille=true.
% SO=0.01 with SOIsPermille=false represents a 1% slope as 0.01 m/m.
cfg.SO             = 0.01;
cfg.SOIsPermille   = false;
cfg.CZ             = [];
cfg.aVanRijn       = 1;
cfg.ksD90Factor    = 2;      % ks = ksD90Factor * D90
cfg.criticalShieldsMethod = 'both';
cfg.applySkinFrictionPartition = true;

%% Bottom-track velocity filtering
cfg.coordinateSystem = 'instrument';
cfg.filtOUTwin       = 3;
cfg.filtOUTstd       = 3;
cfg.movVelTresh      = 0.20;  % m/s
cfg.RDICutOFFtrim    = 0.009; % m/s
cfg.Ws_min           = 0.005; % m
cfg.extraFilter      = 'no';
cfg.BTerrorMultiplier = 4;

cfg.applySTDDespike                = true;
cfg.applyBoatMotionFilter          = true;
cfg.applyMovingMedianFilter        = true;
cfg.applyExtraPercentileFilter     = false;
cfg.applyErrorVelocityFilter       = true;
cfg.applyGPSFilter                 = false;
cfg.applyAttitudeFilter            = true;
cfg.applyWaterBiasFilter           = true;
cfg.applyDirectionFilter           = true;
cfg.applyInstrumentDirectionFilter = true;

cfg.waterBiasReference         = 'last2cells';
cfg.waterBiasComparison        = 'component';
cfg.waterBiasUseBothComponents = true;
cfg.waterBiasMultiplier        = 1.0;

%% Kinematic bedload-transport model
cfg.nActiveLayer = 5;

% Current project value used in the Tomorrica workflow.
% Set to 1 for an uncalibrated/general application. Values other than 1
% should be supported by instrument- or dataset-specific calibration.
cfg.Ninstrument = 3;

cfg.CbMin          = 0.015;
cfg.CbMax          = 0.55;
cfg.CbFallback     = 0.10;
cfg.CbFinalStatistic = 'mean';

cfg.DthicknessMinFactor      = 1;
cfg.DthicknessMaxFactor      = 12;
cfg.DthicknessFallbackFactor = 2;
cfg.DthicknessFinalStatistic = 'mean';

% Measured TrateMean/MaxTrate/MinTrate are stored in g/(m s), whereas the
% calculated transport outputs are in kg/(m s).
cfg.MeasuredTrateToKgFactor = 0.001;

%% Acoustic and backscatter parameters
cfg.kc         = 0.60;  % counts-to-dB conversion coefficient
cfg.Er         = 40;    % reference/noise level in counts
cfg.C1200      = -129;  % instrument constant, dB
cfg.alphaSed12 = 0;     % sediment attenuation contribution
cfg.Tdespike   = 2;
cfg.MM         = 1;
cfg.Win        = 3;

% Apply the exact final mask exported by processBTbedloadVelocity1200.
cfg.velocityMaskMode         = 'final';
cfg.useVelocityFilterMask    = true;
cfg.requireExactVelocityMask = true;
cfg.applyMethod2Regression   = true;

%% Plot controls
cfg.makePlots       = true;
cfg.plotRowsDepth   = [1 3 5];
cfg.plotRowsVelocity = [1 3 7];
cfg.plotVaTrate     = true;
cfg.plotFinalTrateComparison = true;
cfg.plotFinalTrateUseLogAxes = false;

%% Workspace and console behaviour
% AssignToCaller is retained for compatibility with the established
% feature-extraction and structure-building workflow.
cfg.assignToCaller = true;
cfg.assignToBase   = false;
cfg.verbose        = true;

%% Validation behaviour
cfg.throwOnMissing          = true;
cfg.throwOnPlaceholder      = false;
cfg.throwOnDateMismatch     = true;
cfg.throwOnLocationMismatch = true;

end
