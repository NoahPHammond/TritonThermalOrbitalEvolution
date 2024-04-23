%Plot some cool results?
clear all
A=load('ParamSweep1_Triton12b.txt');

%Input parameters
initial_a=A(:,1);%Semi-major axis
grain_size=A(:,2);
initial_temp=A(:,3);
bulk_NH3=A(:,4);
initial_o_thick=A(:,5);

%Output paramters
evol_time=A(:,6)./1e9;
evol_time2=A(:,11)./1e9;
viscosity_final=A(:,10);
e_final=A(:,7);

%add a missing simulation
evol_time=[4.5;evol_time];
evol_time2=[4.5;evol_time2];
grain_size=[0.2;grain_size];
initial_a=[1e3;initial_a];
e_final=[0.5;e_final];

aa=[20,50,100,200,500,1000];
cc=['k','b','r','g','c','m'];
p=zeros(1,6);
figure
subplot(1,2,1)
for i=1:6
p1=scatter(grain_size(initial_a==aa(i)),evol_time2(initial_a==aa(i)),'d','Filled',cc(i));
hold on
p(i)=plot(grain_size(initial_a==aa(i)),evol_time2(initial_a==aa(i))-(i-1)*0.01,cc(i),'DisplayName',['a_0=',num2str(aa(i)),' R_N'],'Linewidth',2);
end
scatter(grain_size(e_final>0.01),evol_time2(e_final>0.01),200,'x','k','Linewidth',3);

set(gca,'Fontsize',18,'Xscale','log')
ylabel('Eccentricity Damping Time (Gyr)')
xlabel('Grain size (mm)')
title('Time until e=10^{-2}')
legend(p)
axis([0.2 10 0 5])
xticks([0.2 0.5 1 2 5 10])

subplot(1,2,2)
for i=1:6
p1=scatter(grain_size(initial_a==aa(i)),evol_time(initial_a==aa(i)),'d','Filled',cc(i));
hold on
p(i)=plot(grain_size(initial_a==aa(i)),evol_time(initial_a==aa(i))-(i-1)*0.01,cc(i),'DisplayName',['a_0=',num2str(aa(i)),' R_N'],'Linewidth',2);
end
scatter(grain_size(e_final>0.0001),evol_time(e_final>0.0001),200,'x','k','Linewidth',3);

set(gca,'Fontsize',18,'Xscale','log')
ylabel('Circularization Time (Gyr)')
xlabel('Grain size (mm)')
title('Time until e=1.6x10^{-5}')
legend(p)
axis([0.2 10 0 5])
xticks([0.2 0.5 1 2 5 10])

tau=zeros(6,6);
tau2=zeros(6,6);
for i=1:6
    k=(i-1)*6+1;
    % tau(i,1:5)=evol_time(k+1:k+5);
    % tau(i,6)=evol_time(k);
    % tau2(i,1:5)=evol_time2(k+1:k+5);
    % tau2(i,6)=evol_time2(k);
    tau(i,:)=evol_time(k:k+5);
    tau2(i,:)=evol_time2(k:k+5);
end

tau(tau==4.5)=5;
figure
subplot(1,2,1)
ds=[0.2,0.5,1,2,5,10];
as=[20,50,100,200,500,1000];
contourf(ds,as,tau2',linspace(0.5,4.5,9),"ShowText",true,"LabelFormat","%0.1f Byr")
set(gca,'Fontsize',18,'Xscale','log','Yscale','log')
xlabel('Grain size (mm)')
ylabel('Initial Semi-Major Axis')
colormap('sky')
xticks(ds)
xticklabels({'0.2', '0.5', '1', '2', '5', '10'})
yticks(as)
yticklabels({'20', '50', '100', '200', '500', '1000'})
c = colorbar;
c.Label.String = 'Evol. Time (Byr)';
axis square
title('Circularization Time e=1e-2')



subplot(1,2,2)
ds=[0.2,0.5,1,2,5,10];
as=[20,50,100,200,500,1000];
contourf(ds,as,tau',linspace(0.5,4.5,9),"ShowText",true,"LabelFormat","%0.1f Byr")
set(gca,'Fontsize',18,'Xscale','log','Yscale','log')
xlabel('Grain size (mm)')
ylabel('Initial Semi-Major Axis')
colormap('sky')
xticks(ds)
xticklabels({'0.2', '0.5', '1', '2', '5', '10'})
yticks(as)
yticklabels({'20', '50', '100', '200', '500', '1000'})
c = colorbar;
c.Label.String = 'Evol. Time (Byr)';
axis square
title('Circularization Time e=1.6e-5')


