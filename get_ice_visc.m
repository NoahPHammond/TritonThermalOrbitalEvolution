function [ viscosity ] = get_ice_visc(temp,grainsize)
%get_ice_visc Returns the viscosity of water ice for volume diffusion
%rheology

gas_const = 8.314;  % J/mol K
activation_energy = 5.94e4;  % J/mol
diffusion_coeff = 9.1e-4;  % m2/s
molar_volume = 1.95e-5;  % m3/mol
viscosity = exp(activation_energy/(gas_const*temp))*gas_const*temp*(grainsize^2)/(14*molar_volume*diffusion_coeff);

end

