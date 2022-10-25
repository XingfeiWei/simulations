clear all
clc
tt=0.25; %time step 1fs;

logfile=importdata('PE_C100_L1run5_heat.log',' ',35);

ac=logfile.data;
t=(ac(:,1)-ac(1,1))*tt/1000000; %convert to ns

%Pair=ac(:,3); % kcal/mol
%PotE=ac(:,4); % kcal/mol

Rg=ac(:,12);
figure;
subplot(2,1,1);
plot(t,Rg);
%title('Rg vs Time');
xlabel('Simulation time (ns)');
ylabel('Rg (Å)');
legend('PE-C100');
%xticks([5:5:50])
%ylim([18 20])
%yticks([18.5 19 19.5 20])
hold on

%subplot(4,1,2);
%plot(tlog,Pair/1e5,'r');
%title('PotEnergy vs Time');
%xlabel('Time (ns)');
%set(gca,'XTickLabel',[]);
%ylabel('Pair E (kcal/mol)');
%xticks([5:5:50])
%ylim([-1.605 -1.59])
%yticks([-1.60,-1.595, -1.590])
llog=size(ac,1);
for i = 1:llog
    Dist(i,1)=((ac(i,9)-ac(i,6))^2+(ac(i,10)-ac(i,7))^2+(ac(i,11)-ac(i,8))^2)^0.5;
end

subplot(2,1,2);
plot(t,Dist);
%xlabel('Time (ns)');
ylabel('Distance (Å)');
%xticks([5:5:50])
%ylim([104 108])
%yticks([105:1:108])

%subplot(4,1,4);
%plot(tlog,PotE/1e5,'b');
xlabel('Simulation time (ns)');
legend('PE-C100');
%ylabel('Potential E (kcal/mol)');
%xticks([0:5:50])
%ylim([-1.605 -1.59])
%yticks([-1.605,-1.60,-1.595, -1.590])

ln=2000;
result(1,1)=mean(Dist(llog-ln:llog,1));
result(1,2)=mean(Rg(llog-ln:llog,1));
