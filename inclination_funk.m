function [inclination_new]=inclination_funk(inclination_old,eccentricity,semi_major_axis,orbital_frequency,phase_lag,dt)
grav_const=6.67e-11; %gran const
triton_mass = 2.139e22;
neptune_spin_freq=2*pi/(24*3600*0.67125);
neptune_radius = 2.4633e7;
neptune_mass = 1.024e26;
neptune_k2=0.19;

A=(3/2)*grav_const*triton_mass*neptune_radius^5*neptune_k2;
mu=grav_const*(neptune_mass+triton_mass);

epsilon=inclination_old; %I think????
alpha=orbital_frequency/neptune_spin_freq;

q=sqrt(sin(epsilon)^2/(1+alpha^2-2*alpha*cos(epsilon)));
q_prime=sqrt(1-q^2);

[Fq,Eq]=ellipke(q);

Bq=(1/q^2)*(Eq-q_prime^2*Fq);

lag_term=(1+eccentricity^2*(1.5+1.5*cos(phase_lag))+(3/8)*eccentricity^4*cos(phase_lag));
dI=dt*((2*A*sin(2*phase_lag)*q*Bq)/(pi*semi_major_axis^(13/2)*sqrt(mu)*(1-eccentricity^2)^5))*lag_term;
inclination_new=inclination_old+dI;

end