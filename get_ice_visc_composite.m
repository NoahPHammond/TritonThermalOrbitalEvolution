function [ viscosity , mechanism ] = get_ice_visc_composite(temp,d,stress)
%get_ice_visc Returns the viscosity of water ice for volume diffusion
%rheology

%Goldsby and Kolhstedt composite flow law constants
A=[1.2e-10,2.2e-7,6.2e-14,4e-19];
n=[1,2.4,1.8,4];
p=[2,0,1.4,0];
Q=[5.94e4,6e4,4.9e4,6e4];
R=8.31;

%Calculate strain rate for each mechanism
strain_rate=A.*(stress.^n./d.^p).*exp(-Q./(R.*temp));

strain_rate_diff=strain_rate(1);
strain_rate_basal=strain_rate(2);
strain_rate_gbs=strain_rate(3);
strain_rate_disl=strain_rate(4);

strain_rate_total=strain_rate_diff+strain_rate_disl+((1/strain_rate_gbs)+(1/strain_rate_basal))^-1;

viscosity=(stress/(2*strain_rate_total));

[v2,mechanism]=max(strain_rate); %Returns 1 for diffusion, 2 for gbs, 3 for basal and 4 for dislocation creep

end

