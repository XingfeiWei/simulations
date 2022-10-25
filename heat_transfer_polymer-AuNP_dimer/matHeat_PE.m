clear all
clc

logfile=importdata('PE_C100_L1run5_flux.txt',' ',2);
pe100=logfile.data;

t=(pe100(:,1)-pe100(1,1))*0.25/1e6;

J(:,1)=-pe100(:,3)*4186.6/6.022140857e23*1e18;  %%heat put-in 1e-18*J
J(:,2)=pe100(:,4)*4186.6/6.022140857e23*1e18;   %%heat taken-out 1e-18*J

figure;
plot(t,J,'r');
xlabel('Simulation time (ns)');
ylabel('Heat (10^{-18}J)');
legend('PE-C100 in','PE-C100 out');

%linear fitting for the heat flux in nW
ln=2000;%data for average result

x=t(size(t,1)-ln:size(t,1),1);

clear y
y=J(size(t,1)-ln:size(t,1),1);
fit1=fit(x,y,'poly1');
result(1,1)=fit1.p1;
clear y
y=J(size(t,1)-ln:size(t,1),2);
fit2=fit(x,y,'poly1');
result(1,2)=fit2.p1;
