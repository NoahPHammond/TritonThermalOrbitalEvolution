function[obliquity_solution]=obliquity_funk(a,i,w)

G=6.67e-11; %gran const
triton_density=2060; %kg/m3

moment_inertia=0.33; %normalized moment of inertia estimate
decoupled_shell_factor=2;

f1=1+(5/2-(15/4)*moment_inertia)^2;
C22=0.25*(w^2/(G*4*pi*triton_density/3))*((5/f1)-1);%Eq 70 Chen and Nimmo 2014, (Darwin-Radau relation)
J2=10*C22/3;

neptune_spin_freq=2*pi/(24*3600*0.67125);
neptune_radius = 2.4633e7;
neptune_mass = 1.024e26;
neptune_k2=0.19;
J_neptune=(neptune_radius^3*neptune_spin_freq^2/(G*neptune_mass))*neptune_k2/3;

precession_rate=1.5*w*J_neptune*(neptune_radius/a)^2;
p=w/precession_rate;

obliquity=linspace(0,pi/4,2e3);





left_side=1.5*((J2+C22).*cos(obliquity)+C22).*p.*sin(obliquity);
right_side=moment_inertia*decoupled_shell_factor.*sin(i-obliquity);

% figure
% plot(obliquity,left_side-right_side)

[junk,yo]=min(abs(left_side-right_side));
obliquity_solution=obliquity(yo);
end


    