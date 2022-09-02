clear all;
clc

lo=[12.7332 12.2502 15.747];% box xyz lower
hi=[142.867 143.35 139.853];% box xyz higher
S=(hi(1,1)-lo(1,1))*(hi(1,2)-lo(1,2))*2;%calculate area of xy plane
%% heat flux
filename = 'heatfluxnve1.log';%change data file name;
delimiterIn = ' ';
A = importdata(filename, ' ', 9);
B = A.data;
n=size(B,1);
ts=1000;%fs
t1=5;t2=10;%ns
for i=1:n
    t(i,1)=(i-1)*ts/1000000;%time in ns
    qo(i,1)=(B(i,4)+B(i,5))*4186.6/6.022140857e23/S*10^20; % J/m^2
    qi(i,1)=-B(i,3)*4186.6/6.022140857e23/S*10^20;% J/m^2
    dq(i,1)=qo(i,1)-qi(i,1);
end
figure;
subplot(2,1,1);
plot(t,qi,t,qo);
title('Heat flux')
xlabel('Time in ns');
ylabel('Energy in J/m^2');
hold on;
subplot(2,1,2);
plot(t,B(:,2));
title('Temperature average');
xlabel('Time in ns');
ylabel('Temperature in K');
figure
plot(t,dq);
title('Qout-Qin');
xlabel('Time in ns');
ylabel('Energy in J/m^2');
aveT = mean(B(:,2))
aveT(1,2) = std(B(:,2))
%piecewise and use the data 60 stepsize
nn=100;
mn=size(qi,1);
k=round(mn/nn);
for in=1:k
        kn=1+nn*(in-1);
        An=qi(kn:kn+nn-1,1);
        qin(in,1)=sum(An)/nn;
        Bn=qo(kn:kn+nn-1,1);
        qio(in,1)=sum(Bn)/nn;
        ttn(in,1)=t(kn,1);
end
figure;
plot(ttn,qin,'+',ttn,qio,'o');
title('100 points step average plot');%change with differen step


%find the dQ/dt use t1-t2 ns data
tq=t(t1/ts*1000000:t2/ts*1000000);
qiq=qi(t1/ts*1000000:t2/ts*1000000);
qoq=qo(t1/ts*1000000:t2/ts*1000000);

ft1=fit(tq,qiq,'poly1');
ft2=fit(tq,qoq,'poly1');
figure;
plot(ft1,tq,qiq,'b-o');
hold on;
plot(ft2,tq,qoq,'k--s');
legend('Qin','fit Qin','Qout','fit Qout');
title('dQ/dt plot from 3ns to 5ns')%change the time used to calculate dQ/dt
xlabel('Time in ns');
ylabel('Energy in J/m^2');
dqdt(1,1)=ft1.p1;
dqdt(1,2)=ft2.p1;
ft1
ft2
dQdt=mean(dqdt)/1e-9 %in W/m^2

%% temperature 
filename = 'tempnve1';%change data file name;
delimiterIn = ' ';
A = importdata(filename);
B = A.data;
clear tt;
clear T;
dt=1;%stepsize
t=10000000/dt; %input a time scale in fs
m=t/10000-1; %input a step size
n=B(1,2);
for j=1:m
    l=(j-1)*(n*2+1);
    tt(1,j)=10000/1000000*(j-1); % time scale in ns
for i=1:n
    k=2*(i-1)+1;
    X(i,j)=B(k+1+l,2);
    T(i,j)=B(k+2+l,1);
end
end
sf=X(1,1);
figure; %1
a=sf+20; %position for temperature plot in A
l=round((a-sf)/2+1); %change the initial position
subplot(3,1,1);
plot(tt,T(l,:));
title('20A');
a=sf+30; %position for temperature plot in A
l=round((a-sf)/2+1);
subplot(3,1,2);
plot(tt,T(l,:));
title('30A');
a=sf+40; %position for temperature plot in A
l=round((a-sf)/2+1);
subplot(3,1,3);
plot(tt,T(l,:));
title('40A')

nn=90;%step smooth for calcualte average Tn
mn=size(T,2);
k=mn/nn;
clear ttn;
clear Tn;
for jn=1:n
    for in=1:1:k
        kn=1+nn*(in-1);
        An=T(jn,kn:kn+nn-1);
        Tn(jn,in)=sum(An)/nn;
        ttn(1,in)=tt(1,kn);
    end
end
sf=X(1,1);
figure; %2
a=sf+20; %position for temperature plot in A
l=round((a-sf)/2+1);%XXXXXchange the shift distance!!!!!
subplot(3,1,1);
plot(ttn,Tn(l,:));
title('20A');
a=sf+30; %position for temperature plot in A
l=round((a-sf)/2+1);%XXXXXchange the shift distance!!!!!
subplot(3,1,2);
plot(ttn,Tn(l,:));
title('30A');
a=sf+40; %position for temperature plot in A
l=round((a-sf)/2+1);%XXXXXchange the shift distance!!!!!
subplot( 3,1,3);
plot(ttn,Tn(l,:));
title('40A');
Xn=size(Tn,2);

figure %3 figure
plot(X(:,1),Tn(:,Xn),'*',X(:,1),Tn(:,Xn-1),'o',X(:,1),Tn(:,Xn-2),'+',X(:,1),Tn(:,Xn-3),'o',X(:,1),Tn(:,Xn-4),'*');
title('Temperature on z-axis'); % with nn point time average
xlabel('Z axis in A');
ylabel('Temperature in K');
%delete points at the end and in the middle, and only use the last 5data
clear Tp;
clear Tp2;
clear Xp;
clear Xp2;
clear dT;
a=9;b=26; c=39;d=56;%change the fit position for temperature
for ip = a:b %change position needs plot x:y
    for jp = Xn-4:Xn%(x+1)total data points
        itp = ip-(a-1);%change the initial ip-(x-1)
        jtp = jp+5-Xn;%(x+1)total data points
        Tp(itp,jtp) = Tn(ip,jp);
    end
    Xp(itp,1) = X(ip,1);
end
xlo=X(1,1);l=size(X,1);xhi=X(l,1);
for ip =c:d%change position needs plot x:y
    for jp = Xn-4:Xn%(4+1)total data points
        itp = ip-(c-1);%change the initial ip-(x-1)
        jtp = jp+5-Xn;%(5+1)total data points
        Tp2(itp,jtp) = Tn(ip,jp);
    end
    Xp2(itp,1) =xhi-X(ip,1);
end

figure %4th figure for dT 
for i=1:5
ft=fit(Xp,Tp(:,i),'poly1');
hold on
plot(ft,Xp,Tp(:,i));
xlabel('Z axis in A');
ylabel('Temperature in K');
dT(1,i)=ft.p1;
end
%figure
for i=1:5
ft=fit(Xp2,Tp2(:,i),'poly1');
hold on
plot(ft,Xp2,Tp2(:,i));
xlabel('Z axis in A');
ylabel('Temperature in K');
dT(1,i+5)=ft.p1;
end
data(1,1)=mean(dT)/1e-10 % in K/m
data(1,2)=std(dT)/1e-10 % in K/m std

%% calculate the kappa
J=dQdt;%W/m^2
Th= dT/1e-10;%K/m
n=size(Th,2);
for i=1:n
k(1,i)=J/Th(1,i);%in W/K.m
end
kk=J/data(1,1)
k(2,1)=mean(k(1,:))
k(2,2)=std(k(1,:))

result(1,1)=dQdt;
result(1,2)=S;
result(1,3)=data(1,1);
result(1,4)=data(1,2);
result(1,5)=result(1,1)/result(1,3);
result(1,6)=k(2,1);
result(1,7)=k(2,2);
result(1,8)=aveT(1,1); %box average temperature during the nemd
result(1,9)=aveT(1,2);

