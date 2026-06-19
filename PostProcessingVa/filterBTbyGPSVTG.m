function [BTx_f,BTy_f,BTmag_f,maskFiltered,filteredIdx,pctFiltered,stats] = ...
    filterBTbyGPSVTG(GPSlong,GPSlat,deltaT,BTx,BTy,VTG,varargin)
% FILTERBTBYGPSVTG
%
% Filter ADCP bottom-track (BT) velocities based on platform motion.
%
% -------------------------------------------------------------------------
% MAIN IDEA
% -------------------------------------------------------------------------
% We remove BT values when the boat/platform motion is too large.
%
% The filtering is based on:
%   1) GPS-derived velocity from consecutive GPS positions
%   2) optional VTG velocity
%   3) optional EH / Pitch / Roll attitude filtering
%
% BT itself is NOT compared to GPS speed anymore.
%
% Instead, BT values are deleted when:
%
%   GPSspeed > GPSThreshold
%   and/or
%   abs(VTG) > VTGThreshold
%
% This means:
%   - if the platform moves too fast, BT is considered unreliable
%   - all corresponding BT values at those indices are replaced by NaN
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% GPSlong, GPSlat
%   GPS coordinate vectors of length N
%
%   If CoordType = 'utm':
%       GPSlong = Easting  [m]
%       GPSlat  = Northing [m]
%
%   If CoordType = 'deg':
%       GPSlong = longitude [degrees]
%       GPSlat  = latitude  [degrees]
%
% deltaT
%   Time step in seconds
%
%   Supported forms:
%   - scalar
%   - vector length N
%   - vector length N-1
%
% BTx, BTy
%   Bottom-track velocity components in x-y or east-north directions
%
% VTG
%   Optional VTG velocity vector, same length as BT
%   If VTG is not available, pass:
%       VTG = []
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE INPUTS
% -------------------------------------------------------------------------
% 'CoordType'          : 'utm' (default) or 'deg'
%
% 'GPSThreshold'       : GPS speed threshold
%                        default = inf
%                        BT is removed where:
%                            GPSspeed > GPSThreshold
%
% 'VTGThreshold'       : VTG speed threshold
%                        default = []
%                        If empty, VTGThreshold = GPSThreshold
%                        BT is removed where:
%                            abs(VTG) > VTGThreshold
%
% 'UseAttitudeFilter'  : true/false, default = false
%
% 'EH'                 : optional EH vector
% 'Pitch'              : optional Pitch vector
% 'Roll'               : optional Roll vector
%
% 'AttitudeMethod'     : 'quartiles', 'threshold', or 'both'
%                        default = 'quartiles'
%
% 'EHThreshold'        : threshold for abs(EH), optional
% 'PitchThreshold'     : threshold for abs(Pitch), optional
% 'RollThreshold'      : threshold for abs(Roll), optional
%
% 'TreatZeroGPSasInvalid' : true/false, default = true
%                        If true, GPS coordinates (0,0) are treated as
%                        invalid / missing GPS and ignored in GPS-speed
%                        calculation.
%
% -------------------------------------------------------------------------
% HOW THE FILTER WORKS
% -------------------------------------------------------------------------
% Step 1:
%   Compute distance between consecutive GPS points
%
% Step 2:
%   Compute GPS speed:
%       GPSspeed = distance / deltaT
%
% Step 3:
%   Compute BT magnitude:
%       BTmag = sqrt(BTx.^2 + BTy.^2)
%
% Step 4:
%   Build GPS mask:
%       maskGPS = GPSspeed > GPSThreshold
%
% Step 5:
%   Build VTG mask:
%       maskVTG = abs(VTG) > VTGThreshold
%
% Step 6 (optional):
%   Build EH / Pitch / Roll masks
%
% Step 7:
%   Combine all masks:
%       maskFiltered = maskGPS | maskVTG | maskAttitude
%
% Step 8:
%   Replace BT values at these indices with NaN
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% BTx_f, BTy_f
%   Filtered BT components
%
% BTmag_f
%   Filtered BT magnitude
%
% maskFiltered
%   Final logical mask used to filter BT
%
% filteredIdx
%   Indices of filtered values
%
% pctFiltered
%   Percentage of filtered BT samples relative to originally valid BT
%
% stats
%   Diagnostic structure with fields:
%       .gpsStepDistance
%       .gpsSpeed
%       .BTmag_raw
%       .BTmag_filtered
%       .maskGPS
%       .maskVTG
%       .maskEH_quart
%       .maskPitch_quart
%       .maskRoll_quart
%       .maskEH_thr
%       .maskPitch_thr
%       .maskRoll_thr
%       .maskEH
%       .maskPitch
%       .maskRoll
%       .maskAttitude
%       .maskFiltered
%       .filteredIdx
%       .nFiltered
%       .pctFiltered
%       .GPSThreshold
%       .VTGThreshold
%       .hasVTG
%
% -------------------------------------------------------------------------
% EXAMPLE 1: UTM coordinates, GPS threshold only
% -------------------------------------------------------------------------
% [BTx_f,BTy_f,BTmag_f,maskFiltered,filteredIdx,pctFiltered,stats] = ...
%     filterBTbyGPSVTG(GPSx_utm,GPSy_utm,deltaT,BTx,BTy,[], ...
%     'CoordType','utm','GPSThreshold',0.5);
%
% This removes BT values where GPS speed > 0.5 m/s
%
% -------------------------------------------------------------------------
% EXAMPLE 2: Degrees coordinates, GPS + VTG threshold
% -------------------------------------------------------------------------
% [BTx_f,BTy_f,BTmag_f,maskFiltered,filteredIdx,pctFiltered,stats] = ...
%     filterBTbyGPSVTG(GPSlon,GPSlat,deltaT,BTx,BTy,VTG, ...
%     'CoordType','deg','GPSThreshold',0.5,'VTGThreshold',0.5);
%
% This removes BT values where:
%   GPSspeed > 0.5
%   or abs(VTG) > 0.5
%
% -------------------------------------------------------------------------
% EXAMPLE 3: Add attitude filtering
% -------------------------------------------------------------------------
% [BTx_f,BTy_f,BTmag_f,maskFiltered,filteredIdx,pctFiltered,stats] = ...
%     filterBTbyGPSVTG(GPSlon,GPSlat,deltaT,BTx,BTy,VTG, ...
%     'CoordType','deg','GPSThreshold',0.5,'VTGThreshold',0.5, ...
%     'EH',EH,'Pitch',Pitch,'Roll',Roll, ...
%     'UseAttitudeFilter',true,'AttitudeMethod','both', ...
%     'PitchThreshold',5,'RollThreshold',5,'EHThreshold',10);
%

    % ---------------------------------------------------------------------
    % 0) Basic checks
    % ---------------------------------------------------------------------
    if ~isvector(GPSlong) || ~isvector(GPSlat) || ~isvector(BTx) || ~isvector(BTy)
        error('GPSlong, GPSlat, BTx, and BTy must all be vectors.');
    end

    if ~isempty(VTG) && ~isvector(VTG)
        error('VTG must be a vector or empty [].');
    end

    GPSlong = GPSlong(:);
    GPSlat  = GPSlat(:);
    BTx     = BTx(:);
    BTy     = BTy(:);

    if ~isempty(VTG)
        VTG = VTG(:);
    end

    N = numel(BTx);

    if numel(GPSlong) ~= N || numel(GPSlat) ~= N || numel(BTy) ~= N
        error('GPSlong, GPSlat, BTx, and BTy must have the same length.');
    end

    if ~isempty(VTG) && numel(VTG) ~= N
        error('VTG must have the same length as BTx, or be empty [].');
    end

    % ---------------------------------------------------------------------
    % 1) Parse options
    % ---------------------------------------------------------------------
    p = inputParser;
    addParameter(p,'CoordType','utm',@(x)ischar(x) || isstring(x));

    addParameter(p,'GPSThreshold',inf,@(x)isnumeric(x) && isscalar(x) && x>=0);
    addParameter(p,'VTGThreshold',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));

    addParameter(p,'UseAttitudeFilter',false,@(x)islogical(x) || isnumeric(x));
    addParameter(p,'EH',[],@(x)isnumeric(x) || islogical(x));
    addParameter(p,'Pitch',[],@(x)isnumeric(x) || islogical(x));
    addParameter(p,'Roll',[],@(x)isnumeric(x) || islogical(x));

    addParameter(p,'AttitudeMethod','quartiles',@(x)ischar(x) || isstring(x));
    addParameter(p,'EHThreshold',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
    addParameter(p,'PitchThreshold',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
    addParameter(p,'RollThreshold',[],@(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));

    addParameter(p,'TreatZeroGPSasInvalid',true,@(x)islogical(x) || isnumeric(x));
    parse(p,varargin{:});

    coordType             = lower(char(p.Results.CoordType));
    GPSThreshold          = p.Results.GPSThreshold;
    VTGThreshold          = p.Results.VTGThreshold;
    useAttitudeFilter     = logical(p.Results.UseAttitudeFilter);

    EH                    = p.Results.EH;
    Pitch                 = p.Results.Pitch;
    Roll                  = p.Results.Roll;
    attitudeMethod        = lower(char(p.Results.AttitudeMethod));

    EHThreshold           = p.Results.EHThreshold;
    PitchThreshold        = p.Results.PitchThreshold;
    RollThreshold         = p.Results.RollThreshold;

    treatZeroGPSasInvalid = logical(p.Results.TreatZeroGPSasInvalid);

    if isempty(VTGThreshold)
        VTGThreshold = GPSThreshold;
    end

    if ~ismember(attitudeMethod,{'quartiles','threshold','both'})
        error('AttitudeMethod must be ''quartiles'', ''threshold'', or ''both''.');
    end

    % ---------------------------------------------------------------------
    % 2) Detect invalid GPS samples
    % ---------------------------------------------------------------------
    invalidGPS = ~isfinite(GPSlong) | ~isfinite(GPSlat) | isnan(GPSlong) | isnan(GPSlat);

    if treatZeroGPSasInvalid
        invalidGPS = invalidGPS | (GPSlong == 0 & GPSlat == 0);
    end

    % ---------------------------------------------------------------------
    % 3) Distance between consecutive GPS points
    % ---------------------------------------------------------------------
    gpsStepDistance = NaN(N,1);

    switch coordType
        case 'utm'
            for i = 2:N
                if ~invalidGPS(i) && ~invalidGPS(i-1)
                    dX = GPSlong(i) - GPSlong(i-1);
                    dY = GPSlat(i)  - GPSlat(i-1);
                    gpsStepDistance(i) = hypot(dX,dY);
                end
            end

        case 'deg'
            for i = 2:N
                if ~invalidGPS(i) && ~invalidGPS(i-1)
                    gpsStepDistance(i) = haversineDistance( ...
                        GPSlat(i-1), GPSlong(i-1), ...
                        GPSlat(i),   GPSlong(i));
                end
            end

        otherwise
            error('CoordType must be ''utm'' or ''deg''.');
    end

    % ---------------------------------------------------------------------
    % 4) Time interval per GPS step
    % ---------------------------------------------------------------------
    dtStep = makeStepDeltaT(deltaT,N);

    % ---------------------------------------------------------------------
    % 5) GPS speed
    % ---------------------------------------------------------------------
    gpsSpeed = NaN(N,1);
    gpsSpeed(2:end) = gpsStepDistance(2:end) ./ dtStep;

    % ---------------------------------------------------------------------
    % 6) BT magnitude
    % ---------------------------------------------------------------------
    BTmag_raw = hypot(BTx,BTy);

    % ---------------------------------------------------------------------
    % 7) GPS threshold mask
    % ---------------------------------------------------------------------
    % BT is removed wherever GPS speed is larger than GPSThreshold
    maskGPS = false(N,1);
    validGPSspeed = ~isnan(gpsSpeed) & isfinite(gpsSpeed);
    maskGPS(validGPSspeed) = gpsSpeed(validGPSspeed) > GPSThreshold;

    % ---------------------------------------------------------------------
    % 8) VTG threshold mask (optional)
    % ---------------------------------------------------------------------
    % BT is removed wherever abs(VTG) is larger than VTGThreshold
    maskVTG = false(N,1);

    if ~isempty(VTG)
        validVTG = ~isnan(VTG) & isfinite(VTG);
        maskVTG(validVTG) = abs(VTG(validVTG)) > VTGThreshold;
    end

    % ---------------------------------------------------------------------
    % 9) Optional EH / Pitch / Roll masks
    % ---------------------------------------------------------------------
    maskEH_quart    = false(N,1);
    maskPitch_quart = false(N,1);
    maskRoll_quart  = false(N,1);

    maskEH_thr      = false(N,1);
    maskPitch_thr   = false(N,1);
    maskRoll_thr    = false(N,1);

    maskEH          = false(N,1);
    maskPitch       = false(N,1);
    maskRoll        = false(N,1);

    if useAttitudeFilter

        if ~isempty(EH)
            EH = EH(:);
            if numel(EH) ~= N
                error('EH must have the same length as BTx.');
            end

            if strcmp(attitudeMethod,'quartiles') || strcmp(attitudeMethod,'both')
                maskEH_quart = quartileOutlierMask(EH);
            end

            if strcmp(attitudeMethod,'threshold') || strcmp(attitudeMethod,'both')
                if ~isempty(EHThreshold)
                    maskEH_thr = abs(EH) > EHThreshold;
                    maskEH_thr(~isfinite(EH)) = false;
                end
            end

            maskEH = maskEH_quart | maskEH_thr;
        end

        if ~isempty(Pitch)
            Pitch = Pitch(:);
            if numel(Pitch) ~= N
                error('Pitch must have the same length as BTx.');
            end

            if strcmp(attitudeMethod,'quartiles') || strcmp(attitudeMethod,'both')
                maskPitch_quart = quartileOutlierMask(Pitch);
            end

            if strcmp(attitudeMethod,'threshold') || strcmp(attitudeMethod,'both')
                if ~isempty(PitchThreshold)
                    maskPitch_thr = abs(Pitch) > PitchThreshold;
                    maskPitch_thr(~isfinite(Pitch)) = false;
                end
            end

            maskPitch = maskPitch_quart | maskPitch_thr;
        end

        if ~isempty(Roll)
            Roll = Roll(:);
            if numel(Roll) ~= N
                error('Roll must have the same length as BTx.');
            end

            if strcmp(attitudeMethod,'quartiles') || strcmp(attitudeMethod,'both')
                maskRoll_quart = quartileOutlierMask(Roll);
            end

            if strcmp(attitudeMethod,'threshold') || strcmp(attitudeMethod,'both')
                if ~isempty(RollThreshold)
                    maskRoll_thr = abs(Roll) > RollThreshold;
                    maskRoll_thr(~isfinite(Roll)) = false;
                end
            end

            maskRoll = maskRoll_quart | maskRoll_thr;
        end
    end

    maskAttitude = maskEH | maskPitch | maskRoll;

    % ---------------------------------------------------------------------
    % 10) Final combined mask
    % ---------------------------------------------------------------------
    maskFiltered = maskGPS | maskVTG | maskAttitude;

    % Only count/filter originally valid BT values
    validBT0 = ~isnan(BTx) & ~isnan(BTy) & isfinite(BTx) & isfinite(BTy);
    maskFiltered(~validBT0) = false;

    % ---------------------------------------------------------------------
    % 11) Apply mask
    % ---------------------------------------------------------------------
    BTx_f   = BTx;
    BTy_f   = BTy;
    BTmag_f = BTmag_raw;

    BTx_f(maskFiltered)   = NaN;
    BTy_f(maskFiltered)   = NaN;
    BTmag_f(maskFiltered) = NaN;

    % ---------------------------------------------------------------------
    % 12) Diagnostics
    % ---------------------------------------------------------------------
    filteredIdx = find(maskFiltered);
    nFiltered   = nnz(maskFiltered);
    nValidBT    = nnz(validBT0);

    if nValidBT > 0
        pctFiltered = 100 * nFiltered / nValidBT;
    else
        pctFiltered = NaN;
    end

    stats = struct();
    stats.coordType            = coordType;
    stats.attitudeMethod       = attitudeMethod;
    stats.gpsStepDistance      = gpsStepDistance;
    stats.gpsSpeed             = gpsSpeed;
    stats.BTmag_raw            = BTmag_raw;
    stats.BTmag_filtered       = BTmag_f;

    stats.maskGPS              = maskGPS;
    stats.maskVTG              = maskVTG;

    stats.maskEH_quart         = maskEH_quart;
    stats.maskPitch_quart      = maskPitch_quart;
    stats.maskRoll_quart       = maskRoll_quart;

    stats.maskEH_thr           = maskEH_thr;
    stats.maskPitch_thr        = maskPitch_thr;
    stats.maskRoll_thr         = maskRoll_thr;

    stats.maskEH               = maskEH;
    stats.maskPitch            = maskPitch;
    stats.maskRoll             = maskRoll;
    stats.maskAttitude         = maskAttitude;
    stats.maskFiltered         = maskFiltered;

    stats.filteredIdx          = filteredIdx;
    stats.nFiltered            = nFiltered;
    stats.nValidBT             = nValidBT;
    stats.pctFiltered          = pctFiltered;

    stats.GPSThreshold         = GPSThreshold;
    stats.VTGThreshold         = VTGThreshold;
    stats.hasVTG               = ~isempty(VTG);
    stats.invalidGPS           = invalidGPS;
    stats.treatZeroGPSasInvalid = treatZeroGPSasInvalid;
end

% =========================================================================
% Helper: create N-1 time-step vector
% =========================================================================
function dtStep = makeStepDeltaT(deltaT,N)
    if isscalar(deltaT)
        dtStep = repmat(deltaT,N-1,1);
    else
        deltaT = deltaT(:);

        if numel(deltaT) == N
            dtStep = deltaT(2:end);
        elseif numel(deltaT) == N-1
            dtStep = deltaT;
        else
            error('deltaT must be scalar, length N, or length N-1.');
        end
    end

    if any(~isfinite(dtStep) | dtStep <= 0)
        error('deltaT values must be finite and > 0.');
    end
end

% =========================================================================
% Helper: quartile outlier mask
% =========================================================================
function TF = quartileOutlierMask(x)
    TF = false(size(x));
    valid = isfinite(x) & ~isnan(x);

    if nnz(valid) >= 3
        TF(valid) = isoutlier(x(valid),'quartiles');
    end
end

% =========================================================================
% Helper: haversine distance [m] for lat/lon in degrees
% =========================================================================
function d = haversineDistance(lat1,lon1,lat2,lon2)
    R = 6371000; % Earth radius [m]

    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    dLat = lat2 - lat1;
    dLon = lon2 - lon1;

    a = sin(dLat/2).^2 + cos(lat1).*cos(lat2).*sin(dLon/2).^2;
    c = 2 * atan2(sqrt(a),sqrt(1-a));
    d = R * c;
end