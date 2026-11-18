% calculate the vicous cattering coeficient
function psiVisc06 = viscScat(f06,k06, D_SS1, ro_S, ni)
%f06- frequency o fthe instrument ,
% k06 - shape number, D_SS1 - particle diameter , ro_S- density of
% sediments, ni -viscosity 
gamma06=sqrt((3.14.*f06)/ni); SS06=(9./(2*gamma06.*D_SS1)).*(1+(2./(gamma06.*D_SS1)));
Tss06= 0.5+9./(2*gamma06.*D_SS1); %HANS and Thorne 2012 /Urick 1948
psiVisc06=k06./(2.*ro_S).*(ro_S./1000-1).^2.*(SS06./(SS06.^2+(ro_S./1000+Tss06).^2));%viscous attenuation coefficient [m^2/kg]
end