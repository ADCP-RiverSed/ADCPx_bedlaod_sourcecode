function c = soundSpeedSeawater(T, S, z)
% soundSpeedSeawater  Speed of sound in seawater (Mackenzie, 1981)
%
%   c = soundSpeedSeawater(T, S, z)
%
%   Inputs:
%       T : temperature in °C      (can be scalar, vector, or matrix)
%       S : salinity in PSU       (same size as T or scalar)
%       z : depth in m            (same size as T or scalar)
%
%   Output:
%       c : sound speed in m/s

T = double(T);
S = double(S);
z = double(z);

c = 1448.96 ...
    + 4.591.*T ...
    - 5.304e-2.*T.^2 ...
    + 2.374e-4.*T.^3 ...
    + 1.340.*(S - 35) ...
    + 1.630e-2.*z ...
    + 1.675e-7.*z.^2 ...
    - 1.025e-2.*T.*(S - 35) ...
    - 7.139e-13.*T.*z.^3;
end
