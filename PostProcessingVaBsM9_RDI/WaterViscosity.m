
function[ni]= WaterViscosity(WT)
ni=interp1([0:5:30],[1.9 1.5 1.3 1.2 1.0 0.9 0.8],WT)*10^-6;        % water viscosity
end
