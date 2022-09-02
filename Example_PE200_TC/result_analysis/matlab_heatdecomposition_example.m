clear all;
clc;
A1 = importdata('flux1nve1.data');
A2 = importdata('flux2nve1.data');
Ak=importdata('keallnve1.data');
heat=importdata('heatfluxnve1.log',' ',9);

lo=[12.9967 12.4073 15.4892];%change y axis
hi=[142.603 143.193 140.111];%change z axis for volume match
S=(hi(1,1)-lo(1,1))*(hi(1,2)-lo(1,2))*2;%calculate area of two sides

bz=[23.4891991  69.8];%value used in the box x-axis
v=(hi(1,3)-lo(1,3))/(bz(1,2)-bz(1,1));%the vol fact in the 600K.in

B1 = A1.data;
B2 = A2.data;
Bk=Ak.data;
hl=heat.data;
ke=mean(Bk(:,2));
pe=mean(Bk(:,3));
n=size(B1,1); %maybe steps in B1 is different in B2

for i=1:n
    for j=2:7
    hf0(i,j-1)=-B1(i,j)+B2(i,j); %trans data to hf by adding flux1 and flux2
    end
    TB(i,1)=i/1000; %in ns
end
%transfer data to hf
for i=1:n
    hf(i,1)=hf0(i,1)*v/2;%multiply the vol scale for the out put dif
    hf(i,2)=hf0(i,2)*v/2;
    %hf(i,1)=((hf0(i,1))^2)^0.5;%make the data positive
    %hf(i,2)=((hf0(i,2))^2)^0.5;%make the data positive
    for j=3:6
        hf(i,j)=hf0(i,j)*v/2;
        %hf(i,j)=((hf0(i,j))^2)^0.5+((hf0(i,2)^2)^0.5);%make the data positive
    end
end
% integrate the flux in J/m^2
for k=1:6
    fi(1,k)=hf(1,k)*1000/2/(1e-20)*4186.6/6.022140857e23;
for i=2:n
    fi(i,k)=(hf(i,k)+hf(i-1,k))*1000/2/(1e-20)*4186.6/6.022140857e23+fi(i-1,k);%trapzoidal rule
end
end

%% calculate heat flux with 5 time interval 1ns each time dt
for k=1:6 % (v_Jz1 v_Jzk1 v_Jzp1 v_Jzb1 v_Jza1 v_Jzd1)
for i=1:5
    i1=(i-1)*2000+1;
    i2=i*2000;
    tt=TB(i1:i2,1);
    qq=fi(i1:i2,k);
    kfit=fit(tt,qq,'poly1');
    dqq(i,k)=kfit.p1/(1e-9);% heat flux in W/m^2
    dt=TB(i2,1)-TB(i1,1);%time difference ns
    dQ(i,k)=(fi(i2,k)-fi(i1,k))/dt/(1e-9);% heat flux in W/m^2
end
end
for i=1:6
    ans(1,i*2-1)=mean(dQ(:,i));%mean heat flux W/M^2
    ans(1,i*2)=std(dQ(:,i));%std of heat flux
    ans2(1,i*2-1)=mean(dqq(:,i));
    ans2(1,i*2)=std(dqq(:,i));
end
%put pe ke in dQ
ans(2,1)=ke;
ans(2,2)=pe;
ans2(2,1)=ke;
ans2(2,2)=pe;
%%

%auto corelation 1000fs
for k=1:6
    for j=1:n/2
        fa(j,k)=0; ta(j,1)=TB(j,1);
    end
    for i=1:n/2
    for j=1:n/2
        fa(j,k)=hf(j+i-1,k)*hf(i,k)/(hf(i,k).^2)+fa(j,k);
    end
    end
end
for k=1:6
    for i=1:n/2
        fa(i,k)=fa(i,k)/(n/2); %average to 1
    end
end

%compare the lammps heat rate J vs ns:
nhl=size(hl,1);
for i=2:nhl
    TB1(i-1,1)=i/1000; %in ns
    qo(i-1,1)=(hl(i,4)+hl(i,5))*4186.6/6.022140857e23/S*10^20;
    qi(i-1,1)=-hl(i,3)*4186.6/6.022140857e23/S*10^20;
end
qc(:,1)=fi(:,1);
figure;%plot compare heat rate
plot(TB1,qi,TB1,qo,TB,qc);
title('Heat compare')
xlabel('Time in ns');
ylabel('Energy in J/m^2');
legend('heat source','heat sink','work in z-axis');

%plot raw data
figure;
subplot(3,2,1);
plot(TB,hf(:,1));
title('z total raw in lammps units');
hold on;
subplot(3,2,2);
plot(TB,hf(:,2));
title('kinetic');
hold on;
subplot(3,2,3);
plot(TB,hf(:,3));
title('pair');
hold on;
subplot(3,2,4);
plot(TB,hf(:,4));
title('bond');
hold on;
subplot(3,2,5);
plot(TB,hf(:,5));
title('angle');
hold on;
subplot(3,2,6);
plot(TB,hf(:,6));
title('torsion');

%plot integrated data
figure;
subplot(3,2,1);
plot(TB,fi(:,1));
title('z total integrated');ylabel('heat J/m^2');xlabel('time ns');
hold on;
subplot(3,2,2);
plot(TB,fi(:,2));
title('ke heat');
hold on;
subplot(3,2,3);
plot(TB,fi(:,3));
title('pair heat');
hold on;
subplot(3,2,4);
plot(TB,fi(:,4));
title('bond heat');
hold on;
subplot(3,2,5);
plot(TB,fi(:,5));
title('angle heat');
hold on;
subplot(3,2,6);
plot(TB,fi(:,6));
title('torsion heat');


%transform to frequency domain for the flux in W/m^2
%flux unit convert W/m2
for i=1:n
    for j=1:6
        ff(i,j)=hf(i,j)/(1e-20)*4186.6/6.022140857e23/(1e-15);%unit convert
    end
end

Fs=1e12;%frequency per sample
L=length(ff);
nf=2^nextpow2(L);
Y=fft(ff,nf,1);
for i=1:nf/2+1 F(i,1)=Fs*(i-1)/nf/1e12; end %to THz
for k=1:6
P(:,k)=abs(Y(1:(nf/2+1),k)/nf); %transformed flux
end

figure; %plot flux vs THz
subplot(3,2,1);
plot(F,P(:,1));
title('FFT z-total heat');
xlabel('frequency in THz');
ylabel('heat flux in W/m^2');ylim([0,0.15e9]);
hold on;
subplot(3,2,2);
plot(F,P(:,2));
title('kinetic');ylim([0,0.15e9]);
hold on;
subplot(3,2,3);
plot(F,P(:,3));
title('pair');ylim([0,0.15e9]);
hold on;
subplot(3,2,4);
plot(F,P(:,4));
title('bond');ylim([0,0.15e9]);
hold on;
subplot(3,2,5);
plot(F,P(:,5));
title('angle');ylim([0,0.15e9]);
hold on;
subplot(3,2,6);
plot(F,P(:,6));
title('torsion');ylim([0,0.15e9]);
%inverse fourier tranform
%separate the Y to 0-0.2THZ, 0.2-0.4THZ, 0.4-0.5THZ
lF=size(F,1);
for i=1:lF
    if F(i,1)<0.02 %low limit
        lcut=i;
    end
    if F(i,1)<0.05 %high limit
        hcut=i;
    end
end
for k=1:6
    for i=1:L
        lY(i,k)=0;
        mY(i,k)=0;
        hY(i,k)=0;
        if i<lcut
            lY(i,k)=Y(i,k);
        else
            mY(i,k)=Y(i,k);
            if i>hcut
                mY(i,k)=0;
                hY(i,k)=Y(i,k);
            end
        end
    end
end
%inverse fft
lff=real(ifft(lY,nf,1));
mff=real(ifft(mY,nf,1));
hff=real(ifft(hY,nf,1));
% integrate the flux low  <0.2THz
for k=1:6
    fil(1,k)=lff(1,k)*(1e-12)/2;
for i=2:n
    fil(i,k)=(lff(i,k)+lff(i-1,k))*(1e-12)/2+fil(i-1,k);%trapzoidal rule
end
end
%plot integrated data <2.4THz
figure;
subplot(3,2,1);
plot(TB,fil(:,1));
title('z integrated flux f<0.02THz');ylabel('heat J/m^2');xlabel('time ns');
hold on;
subplot(3,2,2);
plot(TB,fil(:,2));
title('ke heat');
hold on;
subplot(3,2,3);
plot(TB,fil(:,3));
title('pair heat');
hold on;
subplot(3,2,4);
plot(TB,fil(:,4));
title('bond heat');
hold on;
subplot(3,2,5);
plot(TB,fil(:,5));
title('angle heat');
hold on;
subplot(3,2,6);
plot(TB,fil(:,6));
title('torsion heat');
% integrate the flux medium 0.2-0.4THz
for k=1:6
    fim(1,k)=mff(1,k)*(1e-12)/2;
for i=2:n
    fim(i,k)=(mff(i,k)+mff(i-1,k))*(1e-12)/2+fim(i-1,k);%trapzoidal rule
end
end
%plot integrated data 2.4<f<17THz
figure;
subplot(3,2,1);
plot(TB,fim(:,1));
title('z integrated flux 0.02<f<0.05THz');ylabel('heat J/m^2');xlabel('time ns');
hold on;
subplot(3,2,2);
plot(TB,fim(:,2));
title('ke heat');
hold on;
subplot(3,2,3);
plot(TB,fim(:,3));
title('pair heat');
hold on;
subplot(3,2,4);
plot(TB,fim(:,4));
title('bond heat');
hold on;
subplot(3,2,5);
plot(TB,fim(:,5));
title('angle heat');
hold on;
subplot(3,2,6);
plot(TB,fim(:,6));
title('torsion heat');
% integrate the flux medium 0.2-0.4THz
for k=1:6
    fih(1,k)=hff(1,k)*(1e-12)/2;
for i=2:n
    fih(i,k)=(hff(i,k)+hff(i-1,k))*(1e-12)/2+fih(i-1,k);%trapzoidal rule
end
end
%plot integrated data f>0.4THz
figure;
subplot(3,2,1);
plot(TB,fih(:,1));
title('z integrated flux f>0.05THz');ylabel('heat J/m^2');xlabel('time ns');
hold on;
subplot(3,2,2);
plot(TB,fih(:,2));
title('ke heat');
hold on;
subplot(3,2,3);
plot(TB,fih(:,3));
title('pair heat');
hold on;
subplot(3,2,4);
plot(TB,fih(:,4));
title('bond heat');
hold on;
subplot(3,2,5);
plot(TB,fih(:,5));
title('angle heat');
hold on;
subplot(3,2,6);
plot(TB,fih(:,6));
title('torsion heat');
