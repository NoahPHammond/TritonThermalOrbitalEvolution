function [new_stress,new_viscosity,new_grainsize] = Calc_Ice_Grainsize5(temp,tidal_strain,orbital_frequency,grain_size)

E=9e9; %Youngs Modulus of Ice

%Ice flow law parameters
R=8.31;
A_z=[1.2*10^-10,2.2*10^-7,6.2*10^-14,4*10^-19];
n=[1,2.4,1.8,4];
p=[2,0,1.4,0];
Qe=[59400,60000,49000,60000];
d=grain_size;

%frequency time parameters
nt=400;
w=orbital_frequency;
Period=2*pi/w;
okay=0;
while okay<1
cycles=7;
dt=Period/nt;

%Forced strain
x0=tidal_strain*sin(2*pi*(1:nt)/nt);
x1=zeros(1,nt*cycles)+1; %used to calculate strain in visc layer
s2=zeros(1,nt*cycles); %stress array
visc=s2;               %visc array

for j=1:cycles
    if j==1
        st=2;
    else
        st=1;
    end

    for i=st:nt
        k=((j-1)*nt+i);
        ramp=min(k/(4*nt),1);
        s2(k)=(1-(x1(k-1)-x0(i)*ramp))*E;
        stress=sqrt(s2(k)^2);
        %calculate strain rates
        edot_disl=A_z(4)*exp(-Qe(4)/(R*temp))*stress.^(n(4));
        edot_gsb=A_z(3)*exp(-Qe(3)/(R*temp))*d.^(-p(3)).*stress.^(n(3));
        edot_diff=A_z(1)*exp(-Qe(1)/(R*temp))*d.^(-p(1)).*stress.^(n(1));
        edot_visc=edot_diff+edot_gsb+edot_disl;

        %Now use strain rates to evolve position of x1
        x1(k)=x1(k-1)+dt*edot_visc*sign(s2(k));

        visc(k)=stress/edot_visc;

    end
end

new_stress=max(s2(end-nt:end));
new_viscosity=mean(visc(end-nt:end));

if isnan(new_viscosity)
    nt=round(nt*1.5); %reduce timestep by 2 if nan errors produced: seems to be working...
elseif isinf(new_stress)
    nt=round(nt*1.5);
elseif new_stress>E
    nt=round(nt*1.5);
else
    okay=1;
end
end



%Dissipation rate by grain growth
grain_surf_energy=0.065;
grain_growth_constant=1e-15;


%Dissipation rate by dislocation creep
disl_dissipation_rate=A_z(4)*exp(-Qe(4)/(R*temp))*d^(-p(4))*new_stress^(1+n(4));

%Or is it that if internal dissipation rate from dislocation creep is
%greater than disspation rate from grain growth
gamma=0.01; % Behn, Hirth and Goldby (2020)

d_field_boundary=((12*grain_surf_energy*grain_growth_constant)/(gamma*disl_dissipation_rate))^(1/3);

new_grainsize=grain_size;
%new_grainsize=d_field_boundary;


end

