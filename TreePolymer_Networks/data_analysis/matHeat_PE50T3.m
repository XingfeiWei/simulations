clear all
clc
logfile1=importdata('PE50T3_F9c_flux.txt',' ',2);
tree1=logfile1.data;
logfile2=importdata('PE50T3_B9c_flux.txt',' ',2);
tree2=logfile2.data;
tt = 0.25;
t1=(tree1(:,1)-tree1(1,1)+10000)*tt/1e6;
t2=(tree2(:,1)-tree2(1,1)+10000)*tt/1e6;

J1(:,1)=tree1(:,3)*4186.6/6.022140857e23*1e18;%%heat put-in 1e-18*J
J1(:,2:9)=-tree1(:,4:11)*4186.6/6.022140857e23*1e18;%%heat taken-out 1e-18*J
J2(:,1)=-tree2(:,3)*4186.6/6.022140857e23*1e18;%%heat put-int 1e-18*J
J2(:,2:9)=tree2(:,4:11)*4186.6/6.022140857e23*1e18;%%heat taken-out 1e-18*J

%for i = 1: size(J1,2)
%   J1(:,i) = J1(:,i) - J1(1,i);
%   J2(:,i) = J2(:,i) - J2(1,i);
%end

figure;
subplot(2,2,1)
plot(t1,J1(:,2:end),t1,J1(:,1),'r','LineWidth',1);
xlabel('Simulation time (ns)');
ylabel('Forward Heat (10^{-18}J)');
ylim([0 10]);
legend('#1 in','#2 in', '#3 in', '#4 in','#5 in','#6 in', '#7 in', '#8 in','#1 out');
subplot(2,2,2)
plot(t2,J2(:,2:end),t2,J2(:,1),'r','LineWidth',1);
xlabel('Simulation time (ns)');
ylabel('Backward Heat (10^{-18}J)');
ylim([0 10]);
legend('#1 out','#2 out', '#3 out', '#4 out','#5 out','#6 out', '#7 out', '#8 out','#1 in');

for i = 1:size(J1)
    J1h(i,1)=sum(J1(i,2:end));
    J1h(i,2)=J1(i,1);
end
for i = 1:size(J2)
    J2h(i,1)=sum(J2(i,2:end));
    J2h(i,2)=J2(i,1);
end

%J1h(:,1)=J1h(:,1)-J1h(1,1);J1h(:,2)=J1h(:,2)-J1h(1,2);
%J2h(:,1)=J2h(:,1)-J2h(1,1);J2h(:,2)=J2h(:,2)-J2h(1,2);

subplot(2,2,3)
plot(t1,J1h(:,1),'.r',t1,J1h(:,2),'.b',t2,J2h(:,1),'.',t2,J2h(:,2),'.g','LineWidth',1);
xlabel('Simulation time (ns)');
ylabel('Heat (10^{-18}J)');
ylim([0 10]);
legend('Forward 8 in','Forward 1 out','Backward 8 out', 'Backward 1 in');

%linear fitting for the heat flux in nW
ln=800;%data for average result
clear x y
x=t1(size(t1,1)-ln:size(t1,1),1);
y=J1h(size(t1,1)-ln:size(t1,1),1);
fit1=fit(x,y,'poly1');
result(1,1)=fit1.p1;
clear y
y=J1h(size(t1,1)-ln:size(t1,1),2);
fit2=fit(x,y,'poly1');
result(2,1)=fit2.p1;
clear x y
x=t2(size(t2,1)-ln:size(t2,1),1);
y=J2h(size(t2,1)-ln:size(t2,1),1);
fit3=fit(x,y,'poly1');
result(1,2)=fit3.p1;
clear y
y=J2h(size(t2,1)-ln:size(t2,1),2);
fit4=fit(x,y,'poly1');
result(2,2)=fit4.p1;

%% Rg plot
logfile1=importdata('PE50T3_F9c_heat.log',' ',35);
logfile2=importdata('PE50T3_B9c_heat.log',' ',35);

ac1=logfile1.data;
t1l=(ac1(:,1)-ac1(1,1))*tt/1000000; %convert to ns
ac2=logfile2.data;
t2l=(ac2(:,1)-ac2(1,1))*tt/1000000; %convert to ns

Rg1=ac1(:,12);Rg2=ac2(:,12);
llog1=size(ac1,1);
for i = 1:llog1
    Dist1(i,1)=((ac1(i,9)-ac1(i,6))^2+(ac1(i,10)-ac1(i,7))^2+(ac1(i,11)-ac1(i,8))^2)^0.5;
end
llog2=size(ac2,1);
for i = 1:llog2
    Dist2(i,1)=((ac2(i,9)-ac2(i,6))^2+(ac2(i,10)-ac2(i,7))^2+(ac2(i,11)-ac2(i,8))^2)^0.5;
end

Pair1=ac1(:,3);Pair2=ac2(:,3);
subplot(2,2,4);
plot(t1l,Pair1,'r',t2l,Pair2,'b');

result(3,1)= mean(Rg1(size(Rg1)-ln:size(Rg1),1));
result(3,2)= mean(Rg2(size(Rg2)-ln:size(Rg2),1));
result(4,1)= mean(Dist1(size(Dist1)-ln:size(Dist1),1));
result(4,2)= mean(Dist2(size(Dist2)-ln:size(Dist2),1));
result(5,1)= mean(Pair1(size(Pair1)-ln:size(Pair1),1));
result(5,2)= mean(Pair2(size(Pair2)-ln:size(Pair2),1));

%plot(t1,Rg1,'r',t2,Rg2,'b');
%plot(t1,Dist1,'r',t2,Dist2,'b');
xlabel('Simulation time (ns)');
ylabel('E pair (kcal/mol)');
legend('PE-T3-F','PE-T3-B');
