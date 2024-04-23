%%%%%%% Initialize basic parameters %%%%%%%   
%Triton Parameters
grav_const = 6.6738e-11;  % m3/kg s2
triton_radius = 1.3534e6;  % m; Stern et al. 2015
triton_mass = 2.139e22;  % kg; Brozovic et al. 2015
triton_gravity = grav_const*triton_mass/(triton_radius^2);  % m/s2
triton_volume = (4.0*pi/3.0)*triton_radius^3;  % m3
triton_surface_area = 4.0*pi*triton_radius^2;  % m2
triton_density = triton_mass/triton_volume;  % kg/m3
water_density = 1000.0;  % kg/m3
ice_density = 917;  %kg/m3
rock_density = 3000.0;  % kg/m3
rock_radius = (((3.0*triton_mass/(4.0*pi)) - (water_density*(triton_radius^3)))/(rock_density - water_density))^(1/3);
rock_volume = (4.0*pi/3)*rock_radius^3;  % m3
rock_mass = rock_density*rock_volume;  % kg
water_depth = triton_radius - rock_radius;
h2_real=0.15; %Initial h_2
obliquity=0;

%Neptune parameters
neptune_mass = 1.024e26;
neptune_radius = 2.4633e7;
neptune_Q = 1e4;
neptune_k2 =0.19;


%Physical Constants
ln2 = log(2);
gas_const = 8.314;  % J/mol K
water_heat_capacity = 4185.5;  % J/(kg*K)
core_heat_capacity = 1000;  % J/kg K
convecting_ice_temp=260;
ice_heat_capacity = 74.11 + (7.56*convecting_ice_temp(1));  % J/kg K; Choukroun and Grasset J. Chem. Phys. 2010 eq. 4
ice_heat_of_fusion = 3.2e5;  % J/kg
ice_therm_conductivity = 3.5;  % W/m K  approximate value for cold ice, NOTE: should replace with temperature-dependent equation? but then would need to break down temps within conductive lid
core_therm_conductivity = 3.5;  % W/m K %%%%%%%%%%%%%%% Approximate value chosen %%%%%%%%%%%%%%%
ice_thermal_expansion = 1.6e-4;  % 1/K
ice_youngs_modulus = 9e9;
%ice_thermal_expansion = 5e-5;  % 1/K
ice_kappa = ice_therm_conductivity/(water_density*ice_heat_capacity);  % m2/s %%%%%%%%%%%%%%% Water_Density? %%%%%%%%%%%%%%%
core_kappa = core_therm_conductivity/(rock_density*core_heat_capacity);  % m2/s 
activation_energy = 5.94e4;  % J/mol
surface_temp = 40;  % should we worry about faint young sun? 
%adding constants to calculate dissipation in silicate core
C1=4.05e11; %T&S 7-192 (assumes diffusion creep and grainsize of 3 mm, flow laws from Karato and Wu 1993)
rock_activation_energy=3e5; %J/mol
rock_activation_volume=1e-6; %m^3/mol
core_pressure=triton_gravity*triton_density*(0.5*triton_radius);
core_change_est=0;
last_core_temp=0;
core_change_tol=10;

%special indicators when including Kozai oscilations
yo2=0;
Hammy=0;
visc_change=0;
Ham1_chng=1e4; 
Ham2_chng=1e4;
Ham1_tol=1e3;
Ham2_tol=1e3;

%Variables to calc melting temperature for given ammonia concentration
a1=9.6e-17;  %constants from equation 1 of Leliwa-Kopystynski et al 2002
a2=53.8;
b1=7.95e-8;
b2=650;
c2=4e-8;
max_w_NH3=32.1;
max_pressure=ice_density*water_depth*triton_gravity; %Pressure at top of core, in Pa

%IMPORTANT: If running parameters sweep, set param_sweep=1; set fig_var_num
%to 1 or appropriate value if restarting paramsweep in the middle
param_sweep=0;
num_sim=1;
if param_sweep
    fig_var_num=22;
    num_sim=6*6;
    ice_grain_size_array1=[2e-4,5e-4,1e-3,2e-3,5e-3,1e-2];
    semi_major_axis_array=[20,50,100,200,500,1000];
end

%If doing param sweep, loop through simulations
for i2=1:num_sim

%%%%%%% Initialize simulation parameters %%%%%%%
seconds_in_year = 60*60*24*365.25;
dt0 = 2e3*seconds_in_year;  % years    %%%%%%%%%% TIME STEP %%%%%%%%%%
dt = dt0;
imk2_change_est_tol = 0.05;  % recalculate im(k2) if it is likely to change by more than this multiplicative factor due to changes in layer thicknesses, tidal period, and viscosity
imk2_change_est = 1;  % want it to run a calculation of im(k2) the first time through the loop
max_time = 4.5e9*seconds_in_year;  % cut off simulation after this many years''

%%%%%%% Initialize arrays %%%%%%%
max_i = floor(max_time/dt);
time = zeros(1,max_i);
bottom_viscosity = zeros(1,max_i);
rayleigh = zeros(1,max_i);
convecting_ice_temp = zeros(1,max_i);
imk2_keep =  zeros(1,max_i);
imk2_core_keep=imk2_keep;
h2_keep =  zeros(1,max_i);
ice_thickness = zeros(1,max_i);
nusselt = zeros(1,max_i);
stagnant_lid_thickness = zeros(1,max_i);
ocean_thickness = zeros(1,max_i);
triton_eccentricity = zeros(1,max_i);
triton_inclination=zeros(1,max_i);
triton_argument_of_periapsis=zeros(1,max_i);
core_dissipation_saved=zeros(1,max_i);
radiogenic_heating=zeros(1,max_i);
semi_major_axis = zeros(1,max_i);
tidal_period = zeros(1,max_i);
core_surface_temp = zeros(1,max_i);
ocean_temp = zeros(1,max_i);
w_NH3 = zeros(1,max_i);
visc_nondimen=zeros(1,max_i);
melt_temp = zeros(1,max_i);
core_heat_flux = zeros(1,max_i);
surface_heat_flux = zeros(1,max_i);
tidal_heat_flux  = zeros(1,max_i);
ocean_obliquity_flux = zeros(1,max_i);
tidal_stress = zeros(1,max_i);
grain_size= zeros(1,max_i);
core_viscosity=zeros(1,max_i);
%Special viscosity/stress arrays
visc_array=zeros(1,1e2);
stress_array=zeros(1,1e2);

%Add parameters for kozai-oscilations
triton_inclination(1)=140*pi/180; %Inclination to neptune orbital plane
wp=15*pi/180; %argument of periapsis
triton_argument_of_periapsis(1)=wp;

%Lets save the temperature array every 100kyr
nt_s=45000; 
nodes=202+808; %202 nodes in core + 808 in ice shell/ porportional to 5.44 km per in core/ 4x resolutin in ice shell
nodes=1354; %resolution of 1 km?
Temps_saved=zeros(nodes,nt_s);

%Current orbit position
semi_major_axis_current=  3.54759e8; 
triton_eccentricity_current=1.6e-5;
L0=sqrt(grav_const*neptune_mass*semi_major_axis_current*(1-triton_eccentricity_current^2)); %Angular momentum

%Initial Orbital Conditions
%Choose initial semi-major axis here, but gets reset if doing param sweep
%to right value
semi_major_axis(1)=200*neptune_radius; 
triton_eccentricity(1)=sqrt(1-(L0^2/(grav_const*neptune_mass*semi_major_axis(1))));
orbital_frequency = sqrt(grav_const*neptune_mass/(semi_major_axis(1)^3));
a_old=semi_major_axis(1);

%%%%%%% Initialize ocean/ice shell conditions %%%%%%%
%Choose ice Grainsize Here!
ice_grain_size = 1e-3;  % m
grain_size(1)=ice_grain_size;
bulk_w_NH3 = 0; %  wt% of NH3 in ocean (water&ice) (values should be from 0 to 32.1)
%Choose initial core temperature here!
core_surface_temp(1) = 220; %  K
ocean_mass = triton_mass-rock_mass;
bulk_M_NH3  = (bulk_w_NH3/100)*(ocean_mass);
composite_rheology=1;    

    if param_sweep
    ice_grain_size=ice_grain_size_array1(floor(i2/6)+1); %1,1,1,1,2,2,2,2
   
    i3=i2-6*floor((i2-1)/6); %1,2,3,4,1,2,3,4
    semi_major_axis(1)=semi_major_axis_array(i3)*neptune_radius;
    ice_grain_size
    semi_major_axis(1,1)/neptune_radius
    pause(3)
    triton_eccentricity(1)=sqrt(1-(L0^2/(grav_const*neptune_mass*semi_major_axis(1))));
    orbital_frequency = sqrt(grav_const*neptune_mass/(semi_major_axis(1)^3));
    ocean_flag=0;
    end

pp=max_pressure;
wn=max_w_NH3/100;
if bulk_w_NH3 == 0
    wn=0;
end
%Calculate eutectic melting temp
melt_temp_min=273.1-b1*pp-a1*pp^2-a2*wn-b2*wn^2-c2*wn*pp;

%Code to determine self consistent ocean thickness given ammonia
%concentration and starting temp
%Assumes ammonia remains concentrated in the ocean (none in the ice shell)
%(thus ammonia conc. increases as ocean decreases (up to eutectic conc.)
dx_ocean=0;
if core_surface_temp(1) < melt_temp_min
    ocean_flag = 0;
    M_liq=0;
    if bulk_w_NH3>0
        w_NH3(1)=max_w_NH3;
    end
else
    ocean_flag = 1;
    %NH 6/3/20: Wrote new routine to determine initial ocean depth from NH3
    %concentration and initial core temperature
    z=linspace(0,water_depth,1e5); %Depth array, 10 m spacing
    M_ice_z=ice_density*(triton_radius^3-(triton_radius-z).^3); %Mass of ice shell as function of thickness
    M_ocean_z=ice_density.*((triton_radius-z).^3-(triton_radius-water_depth)^3); %Mass of ocean as function of ice thickness
    w_NH3_z=bulk_w_NH3+bulk_w_NH3.*(M_ice_z./M_ocean_z); % %How ammonia concentration in ocean varies as ocean depth varies
    w_NH3_z(w_NH3_z>max_w_NH3)=max_w_NH3;  %Cap ammonia concentration at eutectic composition
    wn=w_NH3_z./100;
    pp=ice_density*triton_gravity.*z;
    melt_temp_z=273.1-b1.*pp-a1.*pp.^2-a2.*wn-b2.*wn.^2-c2.*wn.*pp; %Equation from LK et al. 2002 for melt temp with NH3 and pressure
    
    [junk,mini]=min(abs(core_surface_temp(1)-melt_temp_z)); %Identify the depth where starting depth intersects melt temp w/depth
    dx_ocean=water_depth-z(mini); %Determines ocean thickness
    
    M_liq=water_density*((4*pi/3)*(dx_ocean+rock_radius)^3-rock_volume);
    %w_NH3(1) = (bulk_M_NH3/M_liq)*100;
    melt_temp(1)=melt_temp_z(mini);
    melt_temp(1)=core_surface_temp(1);
    w_NH3(1)=w_NH3_z(mini);
    w_NH3(1) = (bulk_M_NH3/M_liq)*100;
end
ocean_temp(1) = core_surface_temp(1);
ocean_thickness(1) = dx_ocean;
ice_thickness(1) = water_depth-ocean_thickness(1);
boundary_layer_temp_diff = 10;  % K, within the likely convective temp. ranges this keeps all viscosities within an order of magnitude
basal_ice_temp = core_surface_temp(1);
convecting_ice_temp(1) = basal_ice_temp - boundary_layer_temp_diff;

%%%%%%% Initialize core heat flow %%%%%%%
radtime = 4.5e9 - (time(1)/seconds_in_year);
Hcore = (0.9928*0.008e-6*9.46e-5*exp(radtime*ln2/4.47e9)) + (0.0071*0.008e-6*5.69e-4*exp(radtime*ln2/7.04e8)) + (0.029e-6*2.64e-5*exp(radtime*ln2/1.4e10)) + (0.056e-2*3.48e-9*exp(radtime*ln2/1.25e9));  % T&S eq. 4.8, assuming chondritic core
num_shells = floor(0.1*rock_radius/sqrt(4*dt*core_therm_conductivity/(rock_density*core_heat_capacity)));
shell_length = rock_radius/num_shells;  % gives a whole number of shells
shell_radius = (0:shell_length:rock_radius);
core_temps=zeros(1,num_shells+1);

% Initializing to ideal, steady-state temperatures
%ideal_temps = ocean_temp(1)+0.1*(rock_density*Hcore/(6*core_therm_conductivity))*(rock_radius^2-shell_radius.^2);  % Equation 4.42 from Turcotte 3rd edition
%core_temps = ideal_temps;
core_temps(:) = ocean_temp(1);
core_top_temp(1) = core_temps(end);

% shell_vol is the list of volumes around each shell. Every shell has a  
%   thickness shell_length(L) and shell_vol finds the volume of that shell 
%   plus or minus L/2.
% The final volume is not a full shell due to the volumes occuring at half
%   shell intervals, and the final shell being the exterior of the sphere.
% shell_area finds the surface area of the shell at these L/2 intervals. 
%   That way, when the heat flows, it flows from shell to shell across
%   the L/2 boundries between them.
shell_vol = zeros(1,num_shells+1);
shell_vol(1) = (4*pi/3)*(shell_length/2)^3;
for j = 2:num_shells
    shell_vol(j) = (4*pi/3)*((shell_radius(j)+(shell_length/2))^3-(shell_radius(j)-(shell_length/2))^3);  % m3
end
shell_vol(j+1) = (4*pi/3)*((shell_radius(j))^3-(shell_radius(j)-(shell_length/2))^3);
shell_area = 4*pi*(shell_radius(1:end-1)+(shell_length/2)).^2;  % m2
shell_mass = shell_vol*rock_density;  % kg

% Annoying thing where MATLAB looks in the wrong libraries so fortran
% code won't run.  Bookend at end of script sets things back to normal.
setenv('DYLD_LIBRARY_PATH','/usr/local/gfortran/bin')
%setenv('GFORTRAN_STDIN_UNIT', '3')
%setenv('GFORTRAN_STDOUT_UNIT', '6')
setenv('GFORTRAN_STDERR_UNIT', '0')

%More variables for Kozai Oscilations
Lz0=sqrt(1-triton_eccentricity(1)^2)*cos(triton_inclination(1));
Lz=Lz0;

% Loop until Triton reaches final position or time runs out
i = 2;
while triton_eccentricity(i-1) > 0.000016  % MAIN LOOP
    orbital_frequency = sqrt(grav_const*neptune_mass/(semi_major_axis(i-1)^3));
    tidal_period(i-1) = 2*pi/orbital_frequency;
    
    %%%%%%%%%%  1) Timestep stuff
    %Allow timestep to change depending on if there is an ocean and rapid tidal evolution 

    %Also go slower in the very beginning
    if i<100
        dt=dt0/2000;
    else
        dt=dt0;
    end
    if ocean_flag
        dt=dt/10;
    end


    %%%%%%%%% 2) Determine tidal stress and viscosity 
%     %New routine to get viscosity for compositie ice rehology
if composite_rheology
   %Shortcut, set to maximum dissipation
   %bottom_viscosity(i-1)=ice_youngs_modulus/orbital_frequency;

   if i==2
        stress=3e5; %Have to guess the tidal stress in first timestep to get first viscosity calculation
        bottom_viscosity(i-1) = get_ice_visc_composite(convecting_ice_temp(i-1),ice_grain_size,stress);
        grain_size(i)=ice_grain_size;
        tidal_stress(i)=stress;
   else
        %Need to get h2_real from the fortran code...
        tidal_bulge_perm=h2_real*triton_radius*(neptune_mass/triton_mass)*((triton_radius/semi_major_axis(i-1))^3); %permanent tidal bulge amplitude
        tidal_strain=tidal_bulge_perm*triton_eccentricity(i-1)*3/triton_radius;

        %Calculating tidal strain
        rp=semi_major_axis(i-1)*(1-triton_eccentricity(i-1));
        ra=semi_major_axis(i-1)*(1+triton_eccentricity(i-1));
        tidal_bulge_rp=h2_real*triton_radius*(neptune_mass/triton_mass)*((triton_radius/rp)^3); %permanent tidal bulge amplitude
        tidal_bulge_ra=h2_real*triton_radius*(neptune_mass/triton_mass)*((triton_radius/ra)^3); %permanent tidal bulge amplitude
        tidal_strain=(tidal_bulge_rp-tidal_bulge_ra)/triton_radius;

        %New Routine to find self-consistent stress and grainsize
        %Grainsize found through paleo-wattmeter (Behn et al. 2020/Caswell
        %and Cooper 2022)
        %
        visc_change = visc_change + abs(1 - (ocean_temp(i-2)/ocean_temp(i-1)));
        visc_change = visc_change + abs(1 - (semi_major_axis(i-2)/semi_major_axis(i-1)));
        visc_change = visc_change + abs(1-(tidal_stress(i-2)/tidal_stress(i-1)));
        visc_change = visc_change + abs(1-(ice_thickness(i-2)/ice_thickness(i-1)));
        visc_change = visc_change + abs(1-(imk2_keep(i-2)/imk2_keep(i-1)));


        %Only update basal ice viscosity if it is expected to change by ~1%
        if (visc_change > imk2_change_est_tol/3 || i<1000 || imk2_change_est==0)
            [stress,new_viscosity,ice_grain_size]=Calc_Ice_Grainsize5(ocean_temp(i-1),tidal_strain,orbital_frequency,grain_size(i-1));
            visc_change = 0;
        end
        current_grainsize=ice_grain_size; 
        
        %Damp the rate at which these parameters evolve to smooth
        %instabilities. %Grainsize evolves slower than tidal stress and
        %viscosity
        grain_size(i)=0.999*grain_size(i-1)+0.001*current_grainsize; %make grainsize evolution lag behind?
        tidal_stress(i)=0.99*tidal_stress(i-1)+0.01*stress;
        bottom_viscosity(i-1)=0.99*bottom_viscosity(i-2)+0.01*new_viscosity;
        visc_nondimen(i-1)=bottom_viscosity(i-1)/(ice_youngs_modulus/orbital_frequency);
   end

    %If assuming diffusion creep, easy peasy no stress dependence:
    else
        bottom_viscosity(i-1) = get_ice_visc(convecting_ice_temp(i-1),ice_grain_size);    
        %Calc stress
        tau = bottom_viscosity(i-1)/ice_youngs_modulus; %relaxation time
        E_real = ice_youngs_modulus*((tau.^2).*(orbital_frequency^2))./((tau.^2)*(orbital_frequency^2)+1);
        tidal_bulge_perm=h2_real*triton_radius*(neptune_mass/triton_mass)*((triton_radius/semi_major_axis(i-1))^3); %permanent tidal bulge amplitude
        tidal_strain=tidal_bulge_perm*triton_eccentricity(i-1)*3/triton_radius;
        tidal_stress(i)=E_real.*tidal_strain;    
    end
        
    %%%%%%%     3)Define convection parameters,
    %New code! For rayleigh number, recalculate viscosity. Driving stress
    %is not tides so tidal stress should not affect convective viscosity...
    convective_viscosity=bottom_viscosity(i-1);
    if i>2
    convective_stress=0.02*ice_thickness(i-2)*triton_gravity*water_density*ice_thermal_expansion*230;
    convective_viscosity = get_ice_visc_composite(convecting_ice_temp(i-1),ice_grain_size,convective_stress);
    end
    rayleigh(i-1) = water_density*triton_gravity*ice_thermal_expansion*(basal_ice_temp-surface_temp)*(ice_thickness(i-1)^3)/(ice_kappa*convective_viscosity);
    theta = activation_energy*(basal_ice_temp-surface_temp)/(gas_const*(convecting_ice_temp(i-1)^2));
    Racr=20.9*(theta^4);  %The critical Rayleigh number ~ 20.9 theta^4 (Solomatov 1995

    %Nusselt at cutoff = 
    nusselt_at_cr=0.53*(Racr^(1/3))*(theta^(-4/3));
    %create smooth transition between conduction/convection
    Ra_factor=0.1*Racr; %Ra_factor is rayleigh number range over wich nusselt number gradualy increases from 1 to Nu_at_cr
    smooth_transition=((Racr-rayleigh(i-1))/Ra_factor);

    if(rayleigh(i-1)) > Racr
        nusselt(i-1) = 0.53*(rayleigh(i-1)^(1/3))*(theta^(-4/3));
        stagnant_lid_thickness(i-1) = (ice_thickness(i-1)/nusselt(i-1))*(1 + (1/theta));
    elseif (rayleigh(i-1)<Racr && rayleigh(i-1)>0.9*Racr)
        nusselt(i-1) = nusselt_at_cr-(nusselt_at_cr-1)*smooth_transition;
        stagnant_lid_thickness(i-1) = (ice_thickness(i-1)/nusselt(i-1))*(1 + (1/theta));
    else
        nusselt(i-1)=1.0;
        stagnant_lid_thickness(i-1) = ice_thickness(i-1)*(1 - (2*boundary_layer_temp_diff/(basal_ice_temp - surface_temp)));  % assume linear temperature profile in conductive case, "convecting" ice is the warmest sublayer at the bottom
    end
    if stagnant_lid_thickness(i-1) > ice_thickness(i-1)*(1 - (2*boundary_layer_temp_diff/(basal_ice_temp - surface_temp)))
        stagnant_lid_thickness(i-1) = ice_thickness(i-1)*(1 - (2*boundary_layer_temp_diff/(basal_ice_temp - surface_temp)));  % assume linear temperature profile in conductive case, "convecting" ice is the warmest sublayer at the bottom
        nusselt(i-1)=1;
    end

    % Let's evolve stagnant lid thickness gradually to stop it from hopping
    % back and forth every timestep
    if i>2
    stagnant_lid_thickness(i-1)=0.02*stagnant_lid_thickness(i-1)+0.98*stagnant_lid_thickness(i-2);
    nusselt(i-1)=nusselt(i-1)*0.02+0.98*nusselt(i-2);
    end

    convecting_ice_top = triton_radius - stagnant_lid_thickness(i-1);
    convecting_ice_bottom = rock_radius + ocean_thickness(i-1);
    convecting_ice_mass = water_density*(4*pi/3)*((convecting_ice_top^3) - (convecting_ice_bottom^3));
    
    %%%%%%  4) Calculate love number based on internal structure and frequency
    % Call Fortran code to get appropriate k2 if other parameters have changed sufficiently to be suspicious that it may have changed -- being conservative
    if i > 2
        imk2_change_est = imk2_change_est + 0.1*abs(1 - (bottom_viscosity(i-2)/bottom_viscosity(i-1)));
        imk2_change_est = imk2_change_est + abs(1 - (tidal_period(i-2)/tidal_period(i-1)));
        imk2_change_est = imk2_change_est + abs((ocean_thickness(i-2)-ocean_thickness(i-1))/1e5);  % imk2 approx. doubles for 100 km change in ocean thickness
        imk2_change_est = imk2_change_est + abs(1 - (stagnant_lid_thickness(i-2)/stagnant_lid_thickness(i-1)));
    end
    if imk2_change_est > imk2_change_est_tol
        %if ocean_flag
            [imk2, h2_real] = run_lnc3(triton_radius/1000,triton_density/1000,water_density/1000,rock_density/1000,tidal_period(i-1)/(60*60*24),ocean_thickness(i-1)/1000,bottom_viscosity(i-1),stagnant_lid_thickness(i-1)/1000,(convecting_ice_top - convecting_ice_bottom)/1000,ocean_flag);  % need to match the weird units of Love number code
        %else

        %Previous itration used Tirade to calculate love numbers but we
        %discovered serious instabilities with Tirade, where k2/Q would
        %fluctuate by factors of 100 for small changes in ice shell
        %thickness
%             [real_k2, imk2_tirade] = run_tirade3(triton_radius,rock_radius,triton_density,water_density,rock_density,tidal_period(i-1),ocean_thickness(i-1),bottom_viscosity(i-1),very_bottom_viscosity,stagnant_lid_thickness(i-1),ice_thickness(i-1)-stagnant_lid_thickness(i-1),ocean_flag);  % need to match the weird units of Love number code
%             h2_real=real_k2*1.5;
%             imk2=imk2_tirade;
        %end
        %What if dissipation is 2x as Maxwell model?
        %imk2=imk2*2;
        imk2_change_est = 0;
    end
       
    %Let's add dissipation in the silicate core...
    %Sat Stress code can only handle 4 layers with ridid core
    %So we separately calculate k2/Q of silicate core as if it was separate
    %body with out ice shell
    %This should be a decent approximation due to Gauss's theorem
    %(symmetric mass distributed outide core radius has no gravitational influence on the
    %interior), but neglects how ice shell deformation influences core tidal response 
    T_rock_sol=1450;
    rock_viscosity=1e19*exp((-1+T_rock_sol/core_temps(1))*(rock_activation_energy+rock_activation_volume*core_pressure)/(gas_const*T_rock_sol));
    rock_viscosity(rock_viscosity>1e32)=1e32; %Max out core viscosity for stability purposes
    core_viscosity(i)=rock_viscosity;
    core_change_est=core_temps(1)-last_core_temp;
    if core_change_est > core_change_tol
        last_core_temp=core_temps(1);
        %core's core has 10% radius and density of 6,000
        [imk2_core, h2_real_core] = run_lnc5(rock_radius/1000,1.001*rock_density/1000,rock_density/1000,2*rock_density/1000,tidal_period(i-1)/(60*60*24),0,rock_viscosity,0.2*rock_radius/1000,0.7*rock_radius/1000,0);  % need to match the weird units of Love number code
        %error('look at out.love')
        imk2_core=imk2_core; %Amplify to see if it every becomes significant
    end
    imk2_core_keep(i)=imk2_core;

    %try slowly evolving imk2
    imk2_keep(i)=0.99*imk2_keep(i-1)+0.01*imk2;
    h2_keep(i)=0.99*h2_keep(i-1)+0.01*h2_real;
    
    mu=triton_mass/neptune_mass;
    %%%%%   5) Add obliquity heating
    %NEW SECTION: ADD EQUATION FOR OCEAN DISSIPATION (Nimmo and Spencer 2014, EQ 11)
    if ocean_thickness(i-1)>0
    % obliquity=0.7; %degree
    % obliquity=obliquity*pi/180; %convert to radians
    if mod(i,100)==0
    obliquity=obliquity_funk(semi_major_axis(i-1),triton_inclination(i-1),orbital_frequency);
    end
    cd=0.001; %ocean-drag coefficient
    mu=triton_mass/neptune_mass;
    ocean_heating_rate=8*pi*water_density*triton_radius^5*cd*orbital_frequency^3*obliquity^3;
    %ocean_heating_rate=ocean_heating_rate/2; %Eq 11 is approximate solution, full solution is half when obliquity 0.7
    ocean_obliquity_flux(i)=ocean_heating_rate;
    else 
        ocean_heating_rate=0;
    end
    %Set ocean heating rate to zero
    %ocean_heating_rate=0;

   
    %%%%%   6) Evolve semi-major axis using new love numbers
    %Use McDonald 1964 Method
    k2_real=0.6*h2_real;
    del=atan(imk2/(k2_real));
    term1=-6*k2_real*(orbital_frequency/mu)*(triton_radius^5/(semi_major_axis(i-1)^4));
    term2=(1-triton_eccentricity(i-1)^2)^-6;
    e1=triton_eccentricity(i-1);
    term3=(1.5*sin(del))*e1^2+0.75*(sin(del)+sin(2*del))*e1^4+(1/32)*(3*sin(del)+sin(3*del))*e1^6;
    
    %Below terms arise from dissipation in Neptune
    del2=atan(1/neptune_Q);
    term4=((3*neptune_k2*mu*orbital_frequency*neptune_radius^5)/semi_major_axis(i-1)^4)*(-3*cos(del2)^2+1); %dissipation in Neptune?
    term5=(1.5*sin(del2))*e1^2+0.75*(sin(del2)+sin(2*del2))*e1^4+(1/32)*(3*sin(del2)+sin(3*del2))*e1^6;
    
    term_k3=(2/3)*(triton_radius^2/(semi_major_axis(i-1)^2))*(1-e1^2)^-2;
    term_k32=2*e1^2*sin(del)+3*e1^4*(sin(del)+(3/4)*sin(2*del))+e1^6*(3/8)*(3*sin(del)+2*sin(2*del)+sin(3*del))+e1^8*(1/32)*(2*sin(2*del)+sin(4*del));
    dadt_m=term1*term2*term3+term1*term_k3*term_k32;%+term4*term2*term5; %neglecting neptune dissipation

    yoyo=2*semi_major_axis(i-1)^2/(grav_const*neptune_mass*triton_mass);
    tidal_dissipation_energy=-dadt_m*dt/yoyo;

     %evolve orbital inclination based on ocean obliquity heating
    %triton_inclination(i)=triton_inclination(i-1)-dt*(cos(triton_inclination(i-1))/sin(triton_inclination(i-1)))*(((ocean_heating_rate)*semi_major_axis(i-1))/(grav_const*neptune_mass*triton_mass));
    %if mod(i,10)==0
    triton_inclination(i)=inclination_funk(triton_inclination(i-1),triton_eccentricity(i-1),semi_major_axis(i-1),orbital_frequency,del,dt);

    %Now add incilation evolution due to obliquity tides 
    triton_inclination(i)=triton_inclination(i)-dt*(cos(triton_inclination(i-1))/sin(triton_inclination(i-1)))*(((ocean_heating_rate)*semi_major_axis(i-1))/(grav_const*neptune_mass*triton_mass));
    %end
    %reduce semi-major axis due to dissipation in core
    %da/dt=dE/dt*da/dE
    mu2=rock_mass/neptune_mass;
    k2_real_core=0.6*h2_real_core;
    del=atan(imk2_core/(k2_real_core));
    term1=-6*k2_real_core*(orbital_frequency/mu2)*(rock_radius^5/(semi_major_axis(i-1)^4));
    term2=(1-triton_eccentricity(i-1)^2)^-6;
    e1=triton_eccentricity(i-1);
    term3=(1.5*sin(del))*e1^2+0.75*(sin(del)+sin(2*del))*e1^4+(1/32)*(3*sin(del)+sin(3*del))*e1^6;
    
    %Below terms arise from dissipation in Neptune
    del2=atan(1/neptune_Q);
    term4=((3*neptune_k2*mu2*orbital_frequency*neptune_radius^5)/semi_major_axis(i-1)^4)*(-3*cos(del2)^2+1); %dissipation in Neptune?
    term5=(1.5*sin(del2))*e1^2+0.75*(sin(del2)+sin(2*del2))*e1^4+(1/32)*(3*sin(del2)+sin(3*del2))*e1^6;
    
    term_k3=(2/3)*(rock_radius^2/(semi_major_axis(i-1)^2))*(1-e1^2)^-2;
    term_k32=2*e1^2*sin(del)+3*e1^4*(sin(del)+(3/4)*sin(2*del))+e1^6*(3/8)*(3*sin(del)+2*sin(2*del)+sin(3*del))+e1^8*(1/32)*(2*sin(2*del)+sin(4*del));
    dadt_core=term1*term2*term3+term1*term_k3*term_k32;%+term4*term2*term5; %neglecting neptune dissipation

    yoyo=2*semi_major_axis(i-1)^2/(grav_const*neptune_mass*rock_mass);
    
    core_dissipation=-dadt_core/yoyo;
    core_dissipation_saved(i-1)=core_dissipation;
    
    %dadt_core=-core_dissipation_saved(i-1)*2*semi_major_axis(i-1)^2/(grav_const*neptune_mass*triton_mass);
    dadt_m=dadt_m+dadt_core;

    semi_major_axis(i)=semi_major_axis(i-1)+dadt_m*dt;


    %Update eccentricity using conservation of angular momemtum
    %de/dt=da/dt*dL/da*de/dL
    e1=triton_eccentricity(i-1);
    a0=semi_major_axis(i-1);
    triton_eccentricity(i)=e1+dt*(dadt_m*(1-e1^2)/(2*a0*e1));
    de_m=triton_eccentricity(i)-e1;
    if i==2
        de_m=0;
    end
    
    %7) Conduct heat flow calculations
    %%%%%%% Core heat flow %%%%%%%
    radtime = 4.5e9 - (time(i-1)/seconds_in_year);
    Hcore = (0.9928*0.008e-6*9.46e-5*exp(radtime*ln2/4.47e9)) + (0.0071*0.008e-6*5.69e-4*exp(radtime*ln2/7.04e8)) + (0.029e-6*2.64e-5*exp(radtime*ln2/1.4e10)) + (0.056e-2*3.48e-9*exp(radtime*ln2/1.25e9));  % T&S eq. 4.8, assuming chondritic core

    % Heat added by radioactive elements
    heat = shell_mass*Hcore;  % W radiogenic heat produced by each shell of core
    T_heat = (heat*dt)./(shell_mass*core_heat_capacity);  % Temperature each shell gains via radiogenic heating    core_temps(end) = ocean_temp(i-1);
    radiogenic_heating(i)=Hcore*rock_mass;

    % Heat added by heat flow through core shells
    Q = zeros(1,num_shells+1);
    Q(1) = (core_therm_conductivity/shell_length)*(core_temps(2)-core_temps(1))*shell_area(1);
    for k = 2:num_shells
        Ql = -(core_therm_conductivity/shell_length)*(core_temps(k)-core_temps(k-1))*shell_area(k-1);  % Q from inner (left) shell
        Qr = -(core_therm_conductivity/shell_length)*(core_temps(k)-core_temps(k+1))*shell_area(k);    % Q from outer (right) shell
        Q(k) = Qr+Ql;
    end
    Q(k+1) = -(core_therm_conductivity/shell_length)*(core_temps(k)-core_temps(k+1))*shell_area(k); 
    T_flow = (Q*dt)./(shell_mass*core_heat_capacity);  % Temperature each shell gains via heat flowing between shells
    
    %add head from tidal dissipation
    Q_core_dissipation=core_dissipation/rock_volume;
    T_flow=T_flow+Q_core_dissipation*dt/(rock_density*core_heat_capacity);
    % Setting core temperatures. Top of the core is set to current ocean temp
    %core_temps = core_temps+T_flow;
    core_temps = core_temps+T_heat+T_flow;    

    %Lets cap the maximum core temp at 1600 K
    %thermal energy over that we can add to core out? 
    %Assumes excess energy causes melting, which goes to surface, cools and
    %re-releases energy
    excess_temps=core_temps-1600;
    excess_temps(excess_temps<0)=0;
    core_temps(core_temps>1600)=1600; 
    excess_core_energy=excess_temps.*shell_vol*rock_density*core_heat_capacity;
    core_temps(end) = ocean_temp(i-1);

    %If timestep near 100 kyr mark, save complete temp profile
    if mod(time(i-1)/(pi*1e7*1e5),1) <0.1
        iii=floor(time(i-1)/(pi*1e7*1e5))+1;
        Temps_saved(1:202,iii)=core_temps;

        ocean_nodes=floor(ocean_thickness(i-1)/5414);
        if ocean_nodes
            Temps_saved(202:202+ocean_nodes,iii)=core_temps(end);
            ice_nodes=202+ocean_nodes:1:250;
            if isempty(ice_nodes)
                ice_nodes=250;
            end
            Temps_saved(ice_nodes,iii)=core_temps(end)-(core_temps(end)-40)*((ice_nodes-ice_nodes(1))./(250-ice_nodes(1))).^0.7;
        else
            ice_nodes=202:1:250;
            Temps_saved(ice_nodes,iii)=core_temps(end)-(core_temps(end)-40)*((ice_nodes-ice_nodes(1))./(250-ice_nodes(1))).^0.7;
        end

       

    end


    % Apply heat to appropriate layers (tidal and radiogenic)
    mu=grav_const*(neptune_mass+triton_mass);

    %Save energy arrays
    %tidal_dissipation_energy = dt*(21/2)*imk2*(grav_const*(neptune_mass)^2*(triton_radius^5)*orbital_frequency*(triton_eccentricity(i-1)^2))/(semi_major_axis(i-1)^6);
    core_heat_flow_energy = Q(k+1)*dt;

    %Modify core_heat_flow based on possible melting in core
    core_heat_flow_energy=core_heat_flow_energy-sum(excess_core_energy);

    lid_energy_loss = triton_surface_area*dt*ice_therm_conductivity*(convecting_ice_temp(i-1) - boundary_layer_temp_diff - surface_temp)/stagnant_lid_thickness(i-1);
    available_energy = tidal_dissipation_energy - core_heat_flow_energy - lid_energy_loss; %NH: reversed the sign on core heat flow here.
    %NEW: Add ocean dissipation
    if i>2
    available_energy = tidal_dissipation_energy +  (ocean_heating_rate*dt) - core_heat_flow_energy - lid_energy_loss;
    end
    
    %Save energy arrays
    core_heat_flux(i)=-core_heat_flow_energy/(1e-3*dt*4*pi*triton_radius^2); %all in mW/m2 
    surface_heat_flux(i)=lid_energy_loss/(1e-3*dt*4*pi*triton_radius^2);
    tidal_heat_flux(i)=tidal_dissipation_energy/(1e-3*dt*4*pi*triton_radius^2);
    %ocean_obliquity_flux(i)= eccentricity_ocean_heating_rate*dt/(1e-3*dt*4*pi*triton_radius^2);
    
    %%%%%%      8) Change ocean thickness, concentration.
    if ocean_flag == 1
        % If there is an ocean, the temperatures will stay constant and all
        % the energy (or loss) will go into melting (or freezing)

        %%%%%%% Variable Melting Temperature %%%%%%%
        dT = melt_temp(i-1) - ocean_temp(i-1);  % Needed temp change to reach melting temp
        melting_energy = M_liq*water_heat_capacity*dT;  % How much energy needed             %%%%% M_liq ain't Cp for solution its for H2O %%%%%
        dE_left = available_energy-melting_energy;  % energy leftover

        % when dE_left and available_energy have same sign, melting of ice or freezing of water is occuring
        if (dE_left*available_energy) > 0  
            ocean_temp(i) = ocean_temp(i-1)+dT;  % assign new ocean temp
            new_melt = dE_left/ice_heat_of_fusion;  % change in mass of water  
            ocean_thickness(i) = ocean_thickness(i-1) + new_melt/(water_density*4*pi*(convecting_ice_bottom^2));  % using thin-shell approximation
            
            %Make sure ocean doesn't get thicker than water layer
            if ocean_thickness(i)>water_depth
               ocean_thickness(i)=water_depth-1e2;
            end
            ice_thickness(i) = water_depth - ocean_thickness(i);
  
            M_liq = M_liq+new_melt;
            w_NH3(i) = (bulk_M_NH3/M_liq)*100;  % New weight percent of ammonia

            if M_liq<0 && bulk_M_NH3>0
                ocean_flag=0;
                w_NH3(i)=max_w_NH3;
            end
            
            if w_NH3(i) > max_w_NH3 
                w_NH3(i) = max_w_NH3;  % max_w_NH3 is the eutectic (max) concentration of ammonia
            end        
            
            %Calculate new melting temperature from new pressure and
            %ammonia concentration. Eq 1 from LK 2002
            P_ice = (ice_density*triton_gravity*ice_thickness(i))/1e6;  % MPa Calculates new pressure due to ice thickness change
            wn=w_NH3(i)/100;
            pp=P_ice*1e6;
            melt_temp(i)=273.1-b1*pp-a1*pp^2-a2*wn-b2*wn^2-c2*wn*pp;
            ocean_temp(i)=ocean_temp(i-1)+dT;
            core_surface_temp(i)=ocean_temp(i);

        else % heating or cooling of water is occuring (ie not enough energy to freeze or melt)
            ocean_temp(i) = ocean_temp(i-1)+(available_energy/(M_liq*water_heat_capacity)); % All energy goes to heating or cooling (depending on sign of dE_total)
            ocean_thickness(i) = ocean_thickness(i-1);
            core_surface_temp(i)=ocean_temp(i);
            core_temps(end)=core_surface_temp(i);
            w_NH3(i) = w_NH3(i-1);
            if bulk_w_NH3==0
                w_NH3(i)=0;
            end
            melt_temp(i) = melt_temp(i-1);
            ice_thickness(i) = water_depth - ocean_thickness(i);
            if ice_thickness(i)<1e2
                ice_thickness(i)=1e2;
            end
        end

        % Check if the ocean has frozen out, and if so partition the energy
        % between freezing and cooling
        if ocean_thickness(i) < 0
            excess_melt = ocean_thickness(i)*water_density*4*pi*(convecting_ice_bottom^2);  % keeping ocean_thickness as a negative number so excess energy in next line will also be negative
            excess_energy = excess_melt*ice_heat_of_fusion;
            core_heating_energy = available_energy/((ice_heat_capacity*convecting_ice_mass/(core_heat_capacity*rock_mass))+1);  % see below in other half of if statement
            delta_T = core_heating_energy/(core_heat_capacity*rock_mass);
            core_surface_temp(i) = core_surface_temp(i-1) + delta_T; 
            basal_ice_temp = core_surface_temp(i);
            convecting_ice_temp(i) = basal_ice_temp - boundary_layer_temp_diff;
            ocean_thickness(i) = 0;
            ice_thickness(i) = water_depth;
            imk2_change_est = 2*imk2_change_est_tol;  % make sure it recalculates love number
            ocean_flag = 0;
        else
            basal_ice_temp = ocean_temp(i);
            convecting_ice_temp(i) = basal_ice_temp - boundary_layer_temp_diff;
        end

    else
        % If there is no ocean, no melting or freezing will happen and the
        % energy (or loss) will go into warming (or cooling) the interior
        % Figure out how the heat would lead to temperature change of
        % the core and convecting ice, keeping in mind that the core
        % temperature has to equal the basal ice temperature (and the ocean
        % temperature, if applicable)
        core_heating_energy = available_energy/((ice_heat_capacity*convecting_ice_mass/(core_heat_capacity*rock_mass))+1);  % partition energy so the rock and ice come out to the same delta_T, could calculate ice heating energy same way and get same delta_T answer
        delta_T = core_heating_energy/(core_heat_capacity*rock_mass);
        core_surface_temp(i) = core_surface_temp(i-1) + delta_T;
        
        % Check if the melting boundary has been crossed at the base of the
        % ice, and if so partition the energy into heating and melting
        %if core_surface_temp(i) > Tm(end,w_end)
        if core_surface_temp(i) > melt_temp_min
            excess_temp = core_surface_temp(i) - melt_temp_min;
            excess_energy = excess_temp*((core_heat_capacity*rock_mass) + (ice_heat_capacity*convecting_ice_mass));
            new_melt = excess_energy/ice_heat_of_fusion;  % kg
            ocean_thickness(i) = ocean_thickness(i-1) + new_melt/(water_density*4*pi*(convecting_ice_bottom^2));  % see above in other half of if statement
            ice_thickness(i) = water_depth - ocean_thickness(i);

        %%%%%%% Variable Melting Temperature %%%%%%%
            M_liq = M_liq+new_melt;
            w_NH3(i) = (bulk_M_NH3/M_liq)*100; % New weight percent of ammonia

            if w_NH3(i) > max_w_NH3 
                w_NH3(i) = max_w_NH3; % max_w_NH3 is the eutectic concentration
            end
            if bulk_w_NH3==0
                w_NH3(i)=0;
            end

            %Calculate new melting temperature from new pressure and
            %ammonia concentration. Eq 1 from LK 2002
            P_ice = (ice_density*triton_gravity*ice_thickness(i-1))/1e6;  % MPa Calculates new pressure due to ice thickness change
            wn=w_NH3(i)/100;
            pp=P_ice*1e6;
            melt_temp(i)=273.1-b1*pp-a1*pp^2-a2*wn-b2*wn^2-c2*wn*pp;

            core_surface_temp(i) = melt_temp(i);
            imk2_change_est = 2*imk2_change_est_tol;  % make sure it recalculates love number
            ocean_flag = 1;

        else %Still no ocean
            ocean_thickness(i) = ocean_thickness(i-1);
            ice_thickness(i) = water_depth - ocean_thickness(i);
            w_NH3(i) = w_NH3(i-1);
            melt_temp(i) = melt_temp(i-1);
        end

        % Equalize the temperature where rock and water touch
        ocean_temp(i) = core_surface_temp(i);  %Update ocean temp even though there is no ocean. This is so that once an ocean forms, there isnt a big discrepency in ocean temp from i to i+1
        basal_ice_temp = core_surface_temp(i);
        convecting_ice_temp(i) = basal_ice_temp - boundary_layer_temp_diff;
    end
    
    %%%%%   9) Seek Thee Errors and You Shall Find
    if core_surface_temp(i)==0
        core_surface_temp(i)=core_surface_temp(i-1);
    end
    if isnan(core_surface_temp(i))
        error('Core_temp is NAN')
    end
    if M_liq<0
        M_liq=0;
    end

    %Now adjust orbital paramters due to Kozai effect. Adds LOTS of run
    %time and needs much smaller timestep.
    % dLz=sqrt(1-triton_eccentricity(i)^2)*cos(triton_inclination(i-1))-sqrt(1-triton_eccentricity(i-1)^2)*cos(triton_inclination(i-1));
   
    % if semi_major_axis(i)>100*neptune_radius %after this Neptune's tidal figure dampes Kozai oscilations
    % %[triton_eccentricity(i),triton_inclination(i),wp]=Kozai_func3(orbital_frequency,dt,triton_eccentricity(i),triton_inclination(i-1),wp);
    % e4=sqrt(1-(L0^2/(grav_const*neptune_mass*semi_major_axis(i))));
    % [triton_eccentricity(i),triton_inclination(i),wp]=Kozai_func5(orbital_frequency,dt,e4,triton_inclination(1),triton_eccentricity(i),triton_inclination(i-1),wp);
    % triton_argument_of_periapsis(i)=wp;
    % else
    %    triton_eccentricity(i)=sqrt(1-(L0^2/(grav_const*neptune_mass*semi_major_axis(i))));
    %    triton_inclination(i)=triton_inclination(i-1);
    %    triton_argument_of_periapsis(i)=triton_argument_of_periapsis(i-1);
    % end

    %%%% 10) Update stuff and check for finish
    ice_heat_capacity = 74.11 + (7.56*convecting_ice_temp(i));  % J/kg K; Choukroun and Grasset J. Chem. Phys. 2010 eq. 4

    time(i) = time(i-1) + dt;
    time(i)/(pi*1e7)    
    
    i = i+1;
    if time(i-1) > max_time    % See if we've reached the end of solar system history
        break
    end

end
% Clean up lengths of vectors that were calculated for i-1 originally
stagnant_lid_thickness(i-1) = stagnant_lid_thickness(i-2);
rayleigh(i-1) = rayleigh(i-2);
nusselt(i-1) = nusselt(i-2);
bottom_viscosity(i-1) = bottom_viscosity(i-2);
tidal_period(i-1) = tidal_period(i-2);

% Bookend to set environment back to normal for MATLAB
%setenv('DYLD_LIBRARY_PATH','/Applications/MATLAB_R2013a.app/sys/os/maci64/libgfortran.3.dylib')
%setenv('GFORTRAN_STDIN_UNIT', '-1')
%setenv('GFORTRAN_STDOUT_UNIT', '-1')
%setenv('GFORTRAN_STDERR_UNIT', '-1')

%%%%%%% Plotting stuff %%%%%%%
% delete(graphic_thingy)  % The most important line of the entire code!!!
f1 = figure('Position',[100,1050,1300,900]);
subplot(3,3,1)
plot(time(1:i-1)/(1e6*seconds_in_year),triton_eccentricity(1:i-1),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Triton Eccentricity')
grid on
set(gca,'Fontsize',14,'Yscale','log')
yticks([1.6e-5,1e-4,1e-3,1e-2,1e-1,1e0])
yticklabels({'1.6e-5','1e-4','1e-3','1e-2','1e-1','1'})

%
subplot(3,3,2);
plot(time(1:i-1)/(1e6*seconds_in_year),semi_major_axis(1:i-1)/neptune_radius,'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Semi-major-axis (R_N)')
grid on
set(gca,'Fontsize',14)
%
subplot(3,3,3);
plot(time(1:i-1)/(1e6*seconds_in_year),triton_inclination(1:i-1)*180/pi,'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Orbital Period (days)')
grid on
set(gca,'Fontsize',14)
%
subplot(3,3,4);
plot(time(1:i-1)/(1e6*seconds_in_year),ice_thickness(1:i-1)/1000,'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Ice thickness (km)')
grid on
set(gca,'Fontsize',14)
% subplot(3,3,5);
% plot(time(1:i-1)/seconds_in_year,ocean_thickness(1:i-1)/1000)
% xlabel('Time (yr)')
% ylabel('Ocean thickness (km)')
% grid on
% set(gca,'Fontsize',14)
%
subplot(3,3,6);
plot(time(1:i-1)/(1e6*seconds_in_year),log10(imk2_keep(1:i-1)),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Triton k_2/Q ')
grid on
set(gca,'Fontsize',14)
% subplot(3,4,6);
% % plot(time(1:i-1)/seconds_in_year,tidal_stress(1:i-1))
% % xlabel('Time (yr)')
% % ylabel('Tidal Stress (Pa)')
% % grid on
% % set(gca,'Fontsize',14)
subplot(3,3,7);
plot(time(1:i-1)/(1e6*seconds_in_year),log10(rayleigh(1:i-1)),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Rayleigh number')
grid on
set(gca,'Fontsize',14)
%
subplot(3,3,8);
plot(time(1:i-1)/(1e6*seconds_in_year),nusselt(1:i-1),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Ice shell Nusselt Number')
grid on
set(gca,'Fontsize',14)
%
subplot(3,3,5);
plot(time(1:i-1)/(1e6*seconds_in_year),log10(bottom_viscosity(1:i-1)),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Basal ice viscosity (Pa s)')
grid on
set(gca,'Fontsize',14)
%
subplot(3,3,9);
%plot(time(1:i-1)/seconds_in_year,log10(visc_nondimen(1:i-1)))
plot(time(1:i-2)/seconds_in_year,log10(core_viscosity(1:i-2)),'k-','Linewidth',1.5)
xlabel('Time (Myr)')
ylabel('Core viscosity (Pa s)')
grid on
set(gca,'Fontsize',14)

annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Initial Semi-Major-Axis = ',num2str(semi_major_axis(1)/neptune_radius),'R_N   ', ...
    'Core Temp = ',num2str(core_surface_temp(1)),'K   ','Bulk wNH3 = ',num2str(bulk_w_NH3), '   ', 'Grain Size=',num2str(1e3*ice_grain_size),'mm'], ...
    'EdgeColor', 'none','HorizontalAlignment', 'center','fontsize',20)

f2 = figure();
plot(time(1:i-1)/(1e6*seconds_in_year),surface_heat_flux(1:i-1),'k--','Linewidth',2)
hold on
xlabel('Time (Myr)')
ylabel('Heat Source/Sinks (mW/m^2)')
plot(time(1:i-1)/(1e6*seconds_in_year),core_heat_flux(1:i-1),'r--','Linewidth',2)
plot(time(30:i-1)/(1e6*seconds_in_year),tidal_heat_flux(30:i-1),'k-','Linewidth',1.5)
plot(time(1:i-1)/(1e6*seconds_in_year),1e3*radiogenic_heating(1:i-1)/(4*pi*triton_radius^2),'r-','Linewidth',1.5)
plot(time(1:i-2)/(1e6*seconds_in_year),1e3*core_dissipation_saved(1:i-2)/(4*pi*triton_radius^2),'g-','Linewidth',1.5)
plot(time(1:i-2)/(1e6*seconds_in_year),1e3*ocean_obliquity_flux(1:i-2)/(4*pi*triton_radius^2),'b-','Linewidth',2)
legend('Surface Flux','Core Flux','Tidal Heating Ice', 'Radiogenic Flux','Tidal heating core','Ocean obliquity heating','Fontsize',12)
xlabel('Time (Myr)')
grid on
set(gca,'Fontsize',18,'Yscale','log')
axis([0 time(i-1)/(1e6*seconds_in_year), 1e0 1e4])

decay_time=(2/21)*(triton_mass/neptune_mass)*((semi_major_axis(i-1)/triton_radius)^5)*(1/imk2)*(1/orbital_frequency);


if param_sweep
clear output1 output2 output3
fig_var_num = fig_var_num+1;
output1=[semi_major_axis(1)/neptune_radius,ice_grain_size*1e3,core_surface_temp(1),bulk_w_NH3,ocean_thickness(1)/1e3,time(i-1)/seconds_in_year,triton_eccentricity(i-1),ocean_thickness(i-1)/1e3,max(ocean_thickness),bottom_viscosity(i-1)];
output2=[time(1:100:i-1)/(pi*1e7);semi_major_axis(1:100:i-1)/neptune_radius;triton_eccentricity(1:100:i-1);surface_heat_flux(1:100:i-1);ice_thickness(1:100:i-1);imk2_keep(1:100:i-1)];
output3=[time(1:i-1)/(pi*1e7);semi_major_axis(1:i-1)/neptune_radius;triton_eccentricity(1:i-1);surface_heat_flux(1:i-1);ice_thickness(1:i-1);imk2_keep(1:i-1)];

dlmwrite('ParamSweep1_Triton12.csv',output1,'-append')
writematrix(output3,['Timed_output.run',num2str(fig_var_num),'.csv'])

outputdir='/Triton_ParamSweep12/';


saveas(f1,[pwd, [outputdir,num2str(fig_var_num),'.jpg']]);
saveas(f2,[pwd [outputdir,num2str(fig_var_num),'b.jpg']])
end
end
