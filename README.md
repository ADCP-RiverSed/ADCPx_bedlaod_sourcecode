# ADCPx_bedlaod_sourcecode
This code process ADCP - stationary bedload data and compares (the standalone app does only processing, and it is all automatized) to measured data which is added in the input table as given in the intructions 
# ADCP Bedload Processing Toolbox

MATLAB workflow for processing stationary RDI (RiverPro and RioGrande) kHz ADCP measurements and estimating riverbed mobility and bedload transport from bottom-track velocity, hydraulic variables, and acoustic backscatter.

The workflow was developed for multi-river research datasets and currently includes a configured example for the Tomorrica River dataset in Albania. Raw data and private project paths are intentionally kept outside the repository.

## Main capabilities

The toolbox processes:

- ADCP bottom ranges and four-beam depth statistics;
- water-column velocity profiles;
- apparent bottom-track bedload velocity;
- velocity filtering and per-ensemble validity masks;
- shear velocity, bed shear stress, and critical mobility parameters;
- kinematic and empirical bedload-transport estimates;
- acoustic parameters and bottom-backscatter strength;
- standardized MAT-file structures for feature extraction and machine learning;
- automatic validation of the saved structures and filenames.

## Research-code status

This is research software under active development. Results must be checked against independent measurements and instrument-specific calibration before scientific, engineering, or operational use.

The code contains dataset-specific coefficients and empirical assumptions. In particular, `Ninstrument`, the acoustic constants, roughness assumptions, and regression models must not be transferred blindly between instruments or rivers.

## Repository layout

```text
adcp-bedload-processing/
├─ README.md
├─ GITHUB_UPLOAD.md
├─ .gitignore
├─ config_example.m
├─ run_RDI1200_processing.m
├─ src/
│  ├─ readInputTable_PhSamples.m
│  ├─ processDepthsRNG.m
│  ├─ processWaterProfile.m
│  ├─ processBTbedloadVelocity1200.m
│  ├─ processShearStress1200.m
│  ├─ processBedloadTrate.m
│  ├─ computeADCPacousticParameters.m
│  ├─ processBackscatter1200.m
│  ├─ buildProcessingStructures1200.m
│  ├─ validateFeatureExtractionStructures1200.m
│  └─ supporting helper functions
├─ examples/
├─ data/
└─ results/
```

`config_local.m` is used only on the local computer and must remain excluded by `.gitignore`.

## Requirements

- MATLAB.
- All project functions and helper functions included in `src/` or otherwise available on the MATLAB path.
- Statistics and Machine Learning Toolbox is recommended because the workflow uses functions such as `filloutliers`, `isoutlier`, `prctile`, and correlation statistics.
- Additional toolboxes may be required by individual fitting or optimization helpers.

Check the installed MATLAB products with:

```matlab
ver
license('test', 'Statistics_Toolbox')
```

The MATLAB release used for each published result should be recorded in the release notes or associated manuscript.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/adcp-bedload-processing.git
```

Open MATLAB and change to the repository folder.

### 2. Create the local configuration

Copy:

```text
config_example.m
```

to:

```text
config_local.m
```

Edit the local paths and dataset-specific parameters. The supplied private `config_local.m` for the Tomorrica workflow contains the current local paths, but it must not be committed.

### 3. Confirm that the functions are available

```matlab
which readInputTable_PhSamples -all
which processDepthsRNG -all
which processWaterProfile -all
which processBTbedloadVelocity1200 -all
which processShearStress1200 -all
which processBedloadTrate -all
which computeADCPacousticParameters -all
which processBackscatter1200 -all
which buildProcessingStructures1200 -all
which validateFeatureExtractionStructures1200 -all
```

Each command should point to the intended repository version. Remove or rename stale duplicate function files elsewhere on the MATLAB path.

### 4. Run the workflow

```matlab
run_RDI1200_processing
```

## Processing sequence

The main script performs the following steps:

1. Loads `config_local.m`.
2. Reads the field-data input table.
3. Matches input rows to preprocessed ADCP MAT files.
4. Processes four-beam bottom ranges and depths.
5. Processes the water-column velocity profile.
6. Filters apparent bottom-track bedload velocity.
7. Exports an exact final validity mask for every ADCP file.
8. Calculates shear stress and critical-mobility parameters.
9. Calculates bedload transport using kinematic and empirical methods.
10. Calculates acoustic parameters.
11. Processes bottom backscatter using the same trimmed interval and, optionally, the same final BT mask.
12. Builds the legacy-compatible structures used by downstream feature-extraction scripts.
13. Saves date- and location-prefixed MAT files.
14. Validates the saved filenames, top-level variables, required fields, and placeholders.

## Input data

### Preprocessed ADCP files

The workflow expects one MAT file per stationary measurement position. The filenames must match `ADCPfilename` in the input table.

Typical raw structures used from each MAT file include:

- `BT` for bottom-track data;
- `Data` for configuration, attitude, temperature, and transmit variables;
- `ECHO` for echo-intensity profiles;
- water-profile velocity structures used by `processWaterProfile`.

### Input table

The project-specific input-table reader returns `PhSamplesin` or the updated `ADCPin` structure.

Important fields include:

- `ADCPfilename`;
- `StartADCP`;
- `EndADCP`;
- `DistanceFromLeftBank`;
- `D50`, `D90`, and related grain-size statistics;
- `TrateMean`, `MaxTrate`, and `MinTrate`;
- fractional transport and PSD variables used by the downstream feature-extraction workflow.

`StartADCP` and `EndADCP` are interpreted as **numbers of ensembles removed from the beginning and end**, not as absolute ensemble indices.

For a file containing `N` ensembles:

```matlab
firstUsed = 1 + StartADCP;
lastUsed  = N - EndADCP;
```

## Important units

| Quantity | Expected unit |
|---|---:|
| Bottom-track and water velocity | m/s |
| Depth and range | m |
| Grain sizes in the input table | normally mm |
| Grain sizes used internally in hydraulic equations | m |
| Roughness height `ks` | m |
| Water density | kg/m³ |
| Sediment density | kg/m³ |
| Measured transport rate | g/(m·s) |
| Calculated transport rate | kg/(m·s) |
| Frequency | Hz |
| Bed shear stress | Pa |
| Kinematic viscosity | m²/s |

The current measured transport conversion is:

```matlab
MeasuredTrateToKgFactor = 0.001;
```

Therefore:

```matlab
Trate_kg_ms = Trate_g_ms .* 0.001;
```

## Hydraulic slope convention

`processShearStress1200` uses:

```matlab
if SOIsPermille
    S0 = SO ./ 1000;
else
    S0 = SO;
end
```

Consequently, a 1% slope expressed as a dimensionless ratio must be configured as:

```matlab
cfg.SO = 0.01;
cfg.SOIsPermille = false;
```

Do not use `SO = 1` with `SOIsPermille = false` unless the intended slope is actually 1 m/m.

## Bottom-track velocity filtering

`processBTbedloadVelocity1200` applies selectable filtering steps, including:

- standard-deviation despiking;
- boat-motion threshold filtering;
- moving-median filtering;
- optional percentile filtering;
- error-velocity filtering;
- optional GPS filtering;
- heading, pitch, and roll filtering;
- water-bias filtering using the lower water-profile cells;
- direction-consistency filtering.

The function retains the original variable names and also exports one exact final mask per file:

```matlab
BTvelOut.BTfilterMask_1200{i}
BTvelOut.BTvalidMask_1200{i}
```

Mask convention:

```matlab
BTfilterMask_1200{i} == true   % rejected or invalid ensemble
BTvalidMask_1200{i}  == true   % accepted ensemble
```

The mask length should equal:

```matlab
finishAll(i) - startAll(i) + 1
```

## Backscatter-mask alignment

For direct comparability between bottom-track velocity and backscatter, the recommended settings are:

```matlab
'VelocityMaskMode', 'final', ...
'UseVelocityFilterMask', true, ...
'RequireExactVelocityMask', true
```

This applies the exact final BT mask to the corresponding range, RSSI, attitude, transmit-power, and echo time series.

To use the same trimmed time interval without removing the BT-rejected ensembles:

```matlab
'VelocityMaskMode', 'none'
```

The output diagnostics include the mask source, length, alignment, and number of filtered ensembles.

## Bedload-transport model

The final kinematic estimate uses:

```text
qs = Ninstrument × BTmodf_1200 × CbFinal_1200
     × DthicknessFinal_1200 × sediment density
```

where:

- `BTmodf_1200` is the filtered apparent bottom-track velocity;
- `CbFinal_1200` is selected from the available concentration methods and bounded by configured limits;
- `DthicknessFinal_1200` is selected from available active-layer-thickness methods;
- fallback values are used when no candidate satisfies the configured limits;
- `Ninstrument` is a calibration coefficient.

The Tomorrica local configuration currently uses:

```matlab
cfg.Ninstrument = 3;
```

For an uncalibrated or transferable application, use:

```matlab
cfg.Ninstrument = 1;
```

and document any subsequent calibration separately.

## Acoustic processing

`computeADCPacousticParameters` calculates instrument and sediment acoustic variables from the input-table grain sizes, sound speed, sediment density, and frequency.

All fields in `Acoustic` are copied directly into:

```matlab
BTamp
DataBTamp
```

The original names are retained. A legacy alias is also generated when necessary:

```matlab
BTamp.Ks2_12 = BTamp.Ks_12;
```

The backscatter calculation includes configurable values such as:

- `kc`, counts-to-dB conversion coefficient;
- `Er`, reference/noise level;
- `C1200`, instrument constant;
- water and sediment attenuation;
- pulse length and sampling geometry;
- surface-scattering corrections;
- transmit-power corrections.

These values are instrument-specific and should be documented for each deployment.

## Saved files and structures

The date is extracted from the final or preceding components of `Folder1200`. The configured location token is appended after the date.

Example:

```text
20161202_Tomorrica_TratesData.mat
20161202_Tomorrica_ShearData.mat
20161202_Tomorrica_BackscatterData.mat
```

### Transport folder

```text
<date>_<location>_TratesData.mat
```

contains:

```matlab
Data_field
BTvel
Q
```

### Shear folder

```text
<date>_<location>_ShearData.mat
```

contains:

```matlab
ShearData
```

### Backscatter folder

```text
<date>_<location>_BackscatterData.mat
```

contains:

```matlab
BTamp
DataBTamp
```

`DataBTamp` is retained because the established feature-extraction scripts use both names.

## Structure responsibilities

### `Data_field`

Measured and sampled field information, including:

- grain-size distributions;
- fractional transport rates;
- grab-sample statistics;
- mean, minimum, and maximum measured transport rates;
- cross-section positions.

### `ShearData`

Hydraulic and shear-related variables, including:

- depth and water velocity;
- shear velocities;
- total and skin-friction stresses;
- critical stress and critical velocity;
- excess-shear parameters;
- concentration and active-layer diagnostics;
- turbulence quantities when available.

### `BTvel`

Bottom-track velocity and filtering variables, including:

- raw and filtered components;
- filtered velocity magnitude;
- error velocity;
- filtered fractions;
- all intermediate filtering stages;
- final validity masks and alignment diagnostics.

### `Q`

Calculated bedload-transport variables, including:

- original method-specific transport estimates;
- empirical capacity estimates;
- final kinematic estimate;
- regression-based backscatter estimate when available;
- comparison statistics against measured mean and maximum rates.

### `BTamp` and `DataBTamp`

Acoustic and backscatter variables, including:

- processed range and depth statistics;
- RSSI for all four beams;
- spreading and transmission losses;
- volume-backscattering strength;
- attenuation profiles;
- surface-scattering corrections;
- transmit-power and sampling-geometry variables;
- all fields from `Acoustic`.

## Validation

The main script validates the exact files reported by:

```matlab
StructureReport.FileNamePrefix
```

Example:

```matlab
filePrefix = char(StructureReport.FileNamePrefix);

TratesFile = fullfile(cfg.folderTrates, ...
    [filePrefix '_TratesData.mat']);
ShearFile = fullfile(cfg.folderShear, ...
    [filePrefix '_ShearData.mat']);
BackscatterFile = fullfile(cfg.folderBS, ...
    [filePrefix '_BackscatterData.mat']);

ValidationReport = validateFeatureExtractionStructures1200( ...
    TratesFile, ...
    ShearFile, ...
    BackscatterFile, ...
    'ExpectedLocationName', cfg.locationName, ...
    'ThrowOnMissing', cfg.throwOnMissing, ...
    'ThrowOnPlaceholder', cfg.throwOnPlaceholder, ...
    'ThrowOnDateMismatch', cfg.throwOnDateMismatch, ...
    'ThrowOnLocationMismatch', cfg.throwOnLocationMismatch, ...
    'Verbose', cfg.verbose);
```

Inspect unresolved fields with:

```matlab
disp(ValidationReport.MissingFields)
disp(ValidationReport.PlaceholderFields)
disp(ValidationReport.CriticalUnresolved)
disp(ValidationReport.RowSummary)
```

During development, `ThrowOnPlaceholder=false` is recommended so optional legacy fields can be inspected without terminating the run.

## Downstream feature-extraction compatibility

The structure builder retains the names and orientations expected by the existing feature-extraction scripts. Important examples include:

```matlab
Data_field.MeanTrate
Data_field.MaxTrate
Data_field.QPer6D
Data_field.D_50grabM
Data_field.LeftBankDist14

ShearData.Depth_1200
ShearData.VmodAll_1200
ShearData.Tau_All1200_1
ShearData.Ushear1_1200
ShearData.TKE_prof1200

BTvel.BTmodf_1200
BTvel.FilteredDataN_1200
BTvel.BT_stdf_1200
BTvel.DeltaT1200

BTamp.Depth1200
BTamp.SVBmean1_1200
BTamp.SV1_1200
BTamp.BT_RSSI_1200
BTamp.alpha1mean1200
DataBTamp.Ks2_12
```

Unavailable legacy variables are saved as documented placeholders and listed in the validation report rather than silently omitted.

## Common checks before trusting a run

### Confirm the correct function versions

```matlab
which processBTbedloadVelocity1200 -all
which processBackscatter1200 -all
which buildProcessingStructures1200 -all
```

### Confirm trim lengths

```matlab
table( ...
    BTvelOut.InputRow, ...
    BTvelOut.startAll, ...
    BTvelOut.finishAll, ...
    BTvelOut.BTfilterMaskLength_1200, ...
    BTvelOut.NvalidFinal_1200)
```

### Confirm exact backscatter-mask alignment

```matlab
table( ...
    BackscatterOut.InputRow, ...
    BackscatterOut.NensemblesUsed_1200, ...
    BackscatterOut.NfilteredBTvelocity_1200, ...
    BackscatterOut.BTvelocityFilterExact_1200, ...
    BackscatterOut.BTvelocityFilterAligned_1200)
```

### Confirm acoustic fields

```matlab
requiredAcoustic = { ...
    'c','rs','f12','lambda12','k12', ...
    'x12','F12','Kappa12','FMean12','Ks_12','Ks2_12'};

assert(all(cellfun(@(f) isfield(BTamp,f), requiredAcoustic)), ...
    'One or more required acoustic fields are missing from BTamp.');
```

### Clear stale functions after replacing files

```matlab
clear functions
rehash
```

## Data and privacy policy

Do not commit:

- raw ADCP MAT files;
- measured transport datasets;
- input spreadsheets containing restricted field data;
- generated MAT, FIG, or large binary result files;
- private OneDrive or institutional-storage paths;
- `config_local.m`;
- confidential client or hydropower-project information.

The repository `.gitignore` should include at least:

```gitignore
config_local.m
*.mat
*.fig
*.xlsx
*.xls
*.csv
data/
results/
output/
```

A small anonymized test dataset may be included only when redistribution is permitted.

## Contributing

Use feature branches for changes:

```bash
git checkout -b feature/descriptive-name
```

Before merging:

1. run the complete workflow on a known dataset;
2. run structure validation;
3. compare important outputs against the previous version;
4. document parameter or unit changes;
5. avoid changing established output-variable names without a compatibility alias.

## Citation

An associated manuscript is currently under review:

> Conevski, S., Guerrero, M., Tabesh, M., Rennie, C. D., Winterscheid, A., Schippa, L., and Rüther, N. *Towards Reducing Uncertainty in ADCP Bedload Estimates Using Machine Learning: A Multi-River Study*. Manuscript under review at the *Journal of Hydraulic Research*.

Update this section with the final DOI and bibliographic information after publication.

## Contact

**Slaven Conevski**  
Department of Civil and Environmental Engineering, Norwegian University of Science and Technology  
and Water Resources and Environment, Multiconsult, Norway

For public releases, use the repository issue tracker for bug reports and reproducible examples.

## License

No license should be added until ownership and redistribution rights have been confirmed for every included function and dependency. Without an explicit license, others do not automatically receive permission to copy, modify, or redistribute the code.
