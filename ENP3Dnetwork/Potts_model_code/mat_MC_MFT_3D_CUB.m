%% MCMC model 3D Cubic 5*5*5 network connection: H- V| Z^
%Three states model H(0 or +1, -1), V(0 or +1, -1), Z(0 or +1, -1)
%Considering impossible links: each AuNP only has one polymer reach out
clc
clear all

nH=0;nV=0;
for i=1:10
    for j=1:10
        for k=1:3
            A(i,j,k)=255; %all white
        end
                if mod(i,2)==0 && mod(j,2)==0
                A(i,j,1)=255;A(i,j,2)=0; A(i,j,3)=0; %red
%        elseif mod(i+j,2)==0
%                A(i,j,1)=200;A(i,j,2)=200;A(i,j,3)=200;%grey
                end
    end
end
figure;
subplot(2,2,1)
image(uint8(A)); 
title(['AuNP 2D Square array']); 
%
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,2)
for i = 1:N
    for j=1:N
        for k=1:N
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2;
    surf(x, y, z, 'FaceColor','r','EdgeColor','none') 
    hold on
    patch(x', y', z', 'r','EdgeColor','none') 
            end
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['AuNP 3D Cubic array']); 
%figure;
nH1=1;nH2=1;nV1=1;nV2=1;nZ1=1;nZ2=1;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,3)
for i = 1:N
    for j=1:N
        for k=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
               col=[0.5 0.5 1.0]; %blue H
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
               col=[1.0 1.0 0.1]; %green V
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1
               col=[0.75 0.75 0.75]; %grey Z
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            end
      end 
   end
end

xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['AuNP 3D Cubic array: nH=125, nV=125, nZ=125']); 
%
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nZ1=0;nZ2=0;p0=0.00;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,4)
for i = 1:N
    for j=1:N
        for k=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
               iniN(ti,1)=0;
                if rand(1)<p0
                if rand(1)<0.5
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;
                iniN(ti,1)=1;
                else
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                iniN(ti,1)=-1;
                end
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
               %C(i,j,:)=[0.6 1.0 0.6]; %green V
               iniN(tj,2)=0;
                if rand(1)<p0
                if rand(1)<0.5
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                iniN(tj,2)=1;
                else
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                iniN(tj,2)=-1;
                end
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
               %C(i,j,:)=[0.75 0.75 0.75]; %grey U
               iniN(tk,3)=0;
                if rand(1)<p0
                if rand(1)<0.5
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ1=nZ1+1;
                iniN(tk,3)=1;
                else
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ2=nZ2+1;
                iniN(tk,3)=-1;
                end
                end
                tk=tk+1;
            end
            
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['Random link p=',num2str(p0),': nH=',num2str(nH1+nH2),', nV=',num2str(nV1+nV2),', nZ=',num2str(nZ1+nZ2)]); 
%% MC simulation neighbors
% list 2H 2V 2Z link ids for each NP 
% neighbor list 125H 125V 125Z*first layer NP 5*5=25
NPx=zeros(125,2,3);
for i=1:125
    temp = ceil(i/5);
    NPx(i,1,:)=i;
    NPx(i,2,1)=i+1;
    NPx(i,2,2)=i+5;
    NPx(i,2,3)=i+25;
    if mod(i,5)==0
    NPx(i,2,1)=i+1-5;    end
    if mod(temp,5)==0
    NPx(i,2,2)=i+5-25;   end
    if ceil(ceil(i/5)/5)==5
    NPx(i,2,3)=i+25-125;     end
end
% build neighbour list for H(2h,4v,4z) V(4h,2v,4z) Z(4h,4v,2z)
% Horizontal links
clear nbH nbV nbZ
for i=1:125 %for nbH
    k=0;
    for j=1:125
    if NPx(j,1,1)==i
        np1=j;end
    if NPx(j,2,1)==i
        np2=j;end
    end
    for i1=1:1
        for i2=1:2
            if NPx(np1,i2,i1) ~= i
                k=k+1;
                nbH(i,k)= NPx(np1,i2,i1);end
            if NPx(np2,i2,i1) ~= i
                k=k+1;
                nbH(i,k)= NPx(np2,i2,i1);end
        end
    end
    for i1=2:3
        for i2=1:2
                k=k+1;
                nbH(i,k)= NPx(np1,i2,i1);
                k=k+1;
                nbH(i,k)= NPx(np2,i2,i1);
        end
    end
end
% vertical links
for i=1:125 % for nbV
    k=0;
    for j=1:125
    if NPx(j,1,2)==i
        np1=j;end
    if NPx(j,2,2)==i
        np2=j;end
    end
        for i2=1:2
                k=k+1;
                nbV(i,k)= NPx(np1,i2,1);
                k=k+1;
                nbV(i,k)= NPx(np2,i2,1);
        end
        for i2=1:2
            if NPx(np1,i2,2) ~= i
                k=k+1;
                nbV(i,k)= NPx(np1,i2,2);end
            if NPx(np2,i2,2) ~= i
                k=k+1;
                nbV(i,k)= NPx(np2,i2,2);end
        end
        for i2=1:2
                k=k+1;
                nbV(i,k)= NPx(np1,i2,3);
                k=k+1;
                nbV(i,k)= NPx(np2,i2,3);
        end
end
% Z-links
for i=1:125 % for nbZ
    k=0;
    for j=1:125
    if NPx(j,1,3)==i
        np1=j;end
    if NPx(j,2,3)==i
        np2=j;end
    end
        for i2=1:2
                k=k+1;
                nbZ(i,k)= NPx(np1,i2,1);
                k=k+1;
                nbZ(i,k)= NPx(np2,i2,1);
        end
        for i2=1:2
                k=k+1;
                nbZ(i,k)= NPx(np1,i2,2);
                k=k+1;
                nbZ(i,k)= NPx(np2,i2,2);
        end
        for i2=1:2
            if NPx(np1,i2,3) ~= i
                k=k+1;
                nbZ(i,k)= NPx(np1,i2,3);end
            if NPx(np2,i2,3) ~= i
                k=k+1;
                nbZ(i,k)= NPx(np2,i2,3);end
        end
end

%% MC simulationTemperature effect: 
% E1 linking energy 
% E2 neighbor correlation energy
% EE E-field energy
kB = 1.380649E-23; % in J/K
NA=6.02214076E+23;  % #/mol
E1=5;%7*4186.798188/NA; % connection energy kcal/mol to J
E2=2;%3*4186.798188/NA; % neighbor correlation energy
EE=-6;%-7*4186.798188/NA; %in kcal/mol
H=iniN(:,1);V=iniN(:,2);Z=iniN(:,3);%use the initial set ratio p0
clear nH nV nZ nTot ET ET0 ET1 ET2
eh=0;ev=0;ez=0;
for i=1:125
    if H(i,1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    z1=nbH(i,7);z2=nbH(i,8);z3=nbH(i,9);z4=nbH(i,10);
    % nbH(2h,4v,4u) nbV(4h,2v,4u) nbU(4h,4v,2u)
    eh=eh+(E1+0.5*E2*(abs(H(h1,1))+abs(H(h2,1))+ ...
        abs(V(v1,1))+abs(V(v2,1))+abs(V(v3,1))+abs(V(v4,1))+ ...
        abs(Z(z1,1))+abs(Z(z2,1))+abs(Z(z3,1))+abs(Z(z4,1))))*abs(H(i,1));
    end

    if V(i,1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    z1=nbV(i,7);z2=nbV(i,8);z3=nbV(i,9);z4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V(v1,1))+abs(V(v2,1))+ ...
        abs(H(h1,1))+abs(H(h2,1))+abs(H(h3,1))+abs(H(h4,1))+ ...
        abs(Z(z1,1))+abs(Z(z2,1))+abs(Z(z3,1))+abs(Z(z4,1))))*abs(V(i,1));
    end

    if Z(i,1)~=0
    h1=nbZ(i,1);h2=nbZ(i,2);h3=nbZ(i,3);h4=nbZ(i,4);
    v1=nbZ(i,5);v2=nbZ(i,6);v3=nbZ(i,7);v4=nbZ(i,8);
    z1=nbZ(i,9);z2=nbZ(i,10);
    ez=ez+(E1+0.5*E2*(abs(Z(z1,1))+abs(Z(z2,1))+ ...
        abs(H(h1,1))+abs(H(h2,1))+abs(H(h3,1))+abs(H(h4,1))+ ...
        abs(V(v1,1))+abs(V(v2,1))+abs(V(v3,1))+abs(V(v4,1))))*abs(Z(i,1));
    end
end
et=eh+ev+ez;

it=0;
for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
    it=it+1;
    ET(it,1)=et;
kBT(it,1)=kB*T/4186.798188*NA;%kcal/mol
nit=50000;itl=0; nave=nit/5;
step=0;
while nit>0
    nit=nit-1;
    itl=itl+1;
    step(1,itl+1)=step(1,itl)+1;
H(:,itl+1)=H(:,itl);
V(:,itl+1)=V(:,itl);
Z(:,itl+1)=Z(:,itl);
%clear EH0 EV0 ET0 Eref and specify three conditions
ii=randi(125);
if H(ii,itl)==0
    if rand(1) < 0.5
        H(ii,itl+1)=1; %%0 -> +1
    else 
        H(ii,itl+1)=-1; %%0 -> -1
    end
end
if H(ii,itl)==1
    if rand(1) < 0.5
        H(ii,itl+1)=-1; %%+1 -> -1
    else 
        H(ii,itl+1)=0; %%+1 -> 0
    end
end
if H(ii,itl)==-1
    if rand(1) < 0.5
        H(ii,itl+1)=1; %%-1 -> 1
    else 
        H(ii,itl+1)=0; %%-1 -> 0
    end
end
% jh=nbH(i,1:2);jv=nbH(i,3:6);
eh=0;ev=0;ez=0;
for i=1:125
    if H(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    z1=nbH(i,7);z2=nbH(i,8);z3=nbH(i,9);z4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H(h1,itl+1))+abs(H(h2,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(H(i,itl+1));
    end
    if V(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    z1=nbV(i,7);z2=nbV(i,8);z3=nbV(i,9);z4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V(v1,itl+1))+abs(V(v2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(V(i,itl+1));
    end
    if Z(i,itl+1)~=0
    h1=nbZ(i,1);h2=nbZ(i,2);h3=nbZ(i,3);h4=nbZ(i,4);
    v1=nbZ(i,5);v2=nbZ(i,6);v3=nbZ(i,7);v4=nbZ(i,8);
    z1=nbZ(i,9);z2=nbZ(i,10);
    ez=ez+(E1+0.5*E2*(abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))))*abs(Z(i,itl+1));
    end
end

ET0(it,itl)=ev+eh+ez;
%    if ii==1
    Eref=ET(it,itl);  
    delE=ET0(it,itl)-Eref; 
    if delE <=0
       H(ii,itl+1)=H(ii,itl+1);ET(it,itl+1)=ET0(it,itl);
    else if rand(1) <= exp((-delE)/kBT(it,1))
            H(ii,itl+1)=H(ii,itl+1);ET(it,itl+1)=ET0(it,itl);
        else
            H(ii,itl+1)=H(ii,itl);ET(it,itl+1)=ET(it,itl);
        end
    end
%the same for V links 
%clear EH0 EV0 ET0 Eref 
%specify 3 different types of links
jj=randi(125);
if V(jj,itl)==0
    if rand(1) < 0.5
        V(jj,itl+1)=1; %%0 -> +1
    else 
        V(jj,itl+1)=-1; %%0 -> -1
    end
end
if V(jj,itl)==1
    if rand(1) < 0.5
        V(jj,itl+1)=-1; %%+1 -> -1
    else 
        V(jj,itl+1)=0; %%+1 -> 0
    end
end
if V(jj,itl)==-1
    if rand(1) < 0.5
        V(jj,itl+1)=1; %%-1 -> 1
    else 
        V(jj,itl+1)=0; %%-1 -> 0
    end
end
%V(jj,itl+1)=1-V(jj,itl);
eh=0;ev=0;ez=0;
for i=1:125
    if H(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    z1=nbH(i,7);z2=nbH(i,8);z3=nbH(i,9);z4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H(h1,itl+1))+abs(H(h2,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(H(i,itl+1));
    end
    if V(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    z1=nbV(i,7);z2=nbV(i,8);z3=nbV(i,9);z4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V(v1,itl+1))+abs(V(v2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(V(i,itl+1));
    end
    if Z(i,itl+1)~=0
    h1=nbZ(i,1);h2=nbZ(i,2);h3=nbZ(i,3);h4=nbZ(i,4);
    v1=nbZ(i,5);v2=nbZ(i,6);v3=nbZ(i,7);v4=nbZ(i,8);
    z1=nbZ(i,9);z2=nbZ(i,10);
    ez=ez+(E1+0.5*E2*(abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))))*abs(Z(i,itl+1));
    end
end

ET1(it,itl)=eh+ev+ez;
    Eref=ET(it,itl+1);  
    delE1=ET1(it,itl)-Eref; 
    if delE1 <=0
       V(jj,itl+1)=V(jj,itl+1);ET(it,itl+1)=ET1(it,itl);
    else if rand(1) <= exp((-delE1)/kBT(it,1))
            V(jj,itl+1)=V(jj,itl+1);ET(it,itl+1)=ET1(it,itl);
        else
            V(jj,itl+1)=V(jj,itl);ET(it,itl+1)=ET(it,itl+1);
        end
    end

%the same for U links 
j2=randi(125);
if Z(j2,itl)==0
    if rand(1) < 0.5
        Z(j2,itl+1)=1; %%0 -> +1
    else 
        Z(j2,itl+1)=-1; %%0 -> -1
    end
end
if Z(j2,itl)==1
    if rand(1) < 0.5
        Z(j2,itl+1)=-1; %%+1 -> -1
    else 
        Z(j2,itl+1)=0; %%+1 -> 0
    end
end
if Z(j2,itl)==-1
    if rand(1) < 0.5
        Z(j2,itl+1)=1; %%-1 -> 1
    else 
        Z(j2,itl+1)=0; %%-1 -> 0
    end
end

eh=0;ev=0;ez=0;
for i=1:125
    if H(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    z1=nbH(i,7);z2=nbH(i,8);z3=nbH(i,9);z4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H(h1,itl+1))+abs(H(h2,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(H(i,itl+1));
    end
    if V(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    z1=nbV(i,7);z2=nbV(i,8);z3=nbV(i,9);z4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V(v1,itl+1))+abs(V(v2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+abs(Z(z3,itl+1))+abs(Z(z4,itl+1))))*abs(V(i,itl+1));
    end
    if Z(i,itl+1)~=0
    h1=nbZ(i,1);h2=nbZ(i,2);h3=nbZ(i,3);h4=nbZ(i,4);
    v1=nbZ(i,5);v2=nbZ(i,6);v3=nbZ(i,7);v4=nbZ(i,8);
    z1=nbZ(i,9);z2=nbZ(i,10);
    ez=ez+(E1+0.5*E2*(abs(Z(z1,itl+1))+abs(Z(z2,itl+1))+ ...
        abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))+ ...
        abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))))*abs(Z(i,itl+1));
    end
end

ET2(it,itl)=eh+ev+ez;
    Eref=ET(it,itl+1);  
    delE2=ET2(it,itl)-Eref; 
    if delE2 <=0
       Z(j2,itl+1)=Z(j2,itl+1);ET(it,itl+1)=ET2(it,itl);
    else if rand(1) <= exp((-delE2)/kBT(it,1))
            Z(j2,itl+1)=Z(j2,itl+1);ET(it,itl+1)=ET2(it,itl);
        else
            Z(j2,itl+1)=Z(j2,itl);ET(it,itl+1)=ET(it,itl+1);
        end
    end

nH(it,1)=sum(abs(H(:,1)));
nV(it,1)=sum(abs(V(:,1)));
nZ(it,1)=sum(abs(Z(:,1)));
nTot(it,1)=nH(it,1)+nV(it,1)+nZ(it,1);
H_T(:,it)=H(:,itl+1);
V_T(:,it)=V(:,itl+1);
Z_T(:,it)=Z(:,itl+1);
nH(it,itl+1)=sum(abs(H(:,itl+1)));
nV(it,itl+1)=sum(abs(V(:,itl+1)));
nZ(it,itl+1)=sum(abs(Z(:,itl+1)));
nTot(it,itl+1)=nH(it,itl+1)+nV(it,itl+1)+nZ(it,itl+1);
end
meanHVZ(it,1)=mean(nH(it,itl-nave:itl+1));
meanHVZ(it,2)=mean(nV(it,itl-nave:itl+1));
meanHVZ(it,3)=mean(nZ(it,itl-nave:itl+1));
meanHVZ(it,4)=mean(nTot(it,itl-nave:itl+1));
TT(it,1)=T;
end

%plot the reuslt MC T effect
figure;
subplot(2,2,1)
plot(step,nTot); ylim([0,125]);
legend(['T=10K:H=',num2str(meanHVZ(1,1),'%0.1f'),',V=',num2str(meanHVZ(1,2),'%0.1f'),', Z=',num2str(meanHVZ(1,3),'%0.1f'),],...
['T=100K'],['T=200K'],...
['T=300K'],['T=400K'],['T=500K'],['T=600K'],['T=700K'],...
['T=800K'],...
['T=1200K'],...
['T=1600K'],...
['T=2000K'],...
['T=2500K']);
title(['MC T = 10-2500K, ',num2str(itl),' MC steps']); 

% plot map to 3D cubic network scheme at 300K last step
Tii=4;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nZ1=0;nZ2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,2)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if Z_T(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ1=nZ1+1;
                end
                if Z_T(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ2=nZ2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nZ1+nZ2)]); 


% plot map to 3D cubic network scheme at 1200K last step
Tii=10;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nZ1=0;nZ2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,3)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if Z_T(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ1=nZ1+1;
                end
                if Z_T(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ2=nZ2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nZ1+nZ2)]); 

% plot map to 3D cubic network scheme at 2000K last step
Tii=12;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nZ1=0;nZ2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,4)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if Z_T(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ1=nZ1+1;
                end
                if Z_T(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nZ2=nZ2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nZ1+nZ2)]); 
%% E-field effect: EE is the potential energy from E-field
H_E=iniN(:,1);V_E=iniN(:,2);U_E=iniN(:,3);%use the initial set ratio p0
clear nH_E nV_E nU_E nTot_E ET_E ET0_E ET1_E ET2_E
nbU=nbZ;
eh=0;ev=0;eu=0;itl=0;
for i=1:125
    if H_E(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    u1=nbH(i,7);u2=nbH(i,8);u3=nbH(i,9);u4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    end
    if V_E(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    u1=nbV(i,7);u2=nbV(i,8);u3=nbV(i,9);u4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(V_E(i,itl+1));
    end
    if U_E(i,itl+1)~=0
    h1=nbU(i,1);h2=nbU(i,2);h3=nbU(i,3);h4=nbU(i,4);
    v1=nbU(i,5);v2=nbU(i,6);v3=nbU(i,7);v4=nbU(i,8);
    u1=nbU(i,9);u2=nbU(i,10);
    eu=eu+(E1+0.5*E2*(abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(U_E(i,itl+1));
    end
end
et=eh+ev+eu;

it=0;
for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
    it=it+1;
kBT(it,1)=kB*T/4186.798188*NA;%kcal/mol
nit=50000;itl=0;nave=nit/5;
step=0;
ET_E(it,1)=et;
while nit>0
    nit=nit-1;
    itl=itl+1;
    step(1,itl+1)=step(1,itl)+1;
H_E(:,itl+1)=H_E(:,itl);
V_E(:,itl+1)=V_E(:,itl);
U_E(:,itl+1)=U_E(:,itl);
%flip links of 3 types
ii=randi(125);
if H_E(ii,itl)==0
    if rand(1) < 0.5
        H_E(ii,itl+1)=1; %%0 -> +1
    else 
        H_E(ii,itl+1)=-1; %%0 -> -1
    end
end
if H_E(ii,itl)==1
    if rand(1) < 0.5
        H_E(ii,itl+1)=-1; %%+1 -> -1
    else 
        H_E(ii,itl+1)=0; %%+1 -> 0
    end
end
if H_E(ii,itl)==-1
    if rand(1) < 0.5
        H_E(ii,itl+1)=1; %%-1 -> 1
    else 
        H_E(ii,itl+1)=0; %%-1 -> 0
    end
end
eh=0;ev=0;eu=0;
for i=1:125
    if H_E(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    u1=nbH(i,7);u2=nbH(i,8);u3=nbH(i,9);u4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    end
    if V_E(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    u1=nbV(i,7);u2=nbV(i,8);u3=nbV(i,9);u4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(V_E(i,itl+1));
    end
    if U_E(i,itl+1)~=0
    h1=nbU(i,1);h2=nbU(i,2);h3=nbU(i,3);h4=nbU(i,4);
    v1=nbU(i,5);v2=nbU(i,6);v3=nbU(i,7);v4=nbU(i,8);
    u1=nbU(i,9);u2=nbU(i,10);
    eu=eu+(E1+0.5*E2*(abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(U_E(i,itl+1));
    end
end
ET0_E(it,itl)=ev+eh+eu;
    Eref=ET_E(it,itl);  
    delE=ET0_E(it,itl)-Eref; 
    if delE <=0
       H_E(ii,itl+1)=H_E(ii,itl+1);ET_E(it,itl+1)=ET0_E(it,itl);
    else if rand(1) <= exp((-delE)/kBT(it,1))
            H_E(ii,itl+1)=H_E(ii,itl+1);ET_E(it,itl+1)=ET0_E(it,itl);
        else
            H_E(ii,itl+1)=H_E(ii,itl);ET_E(it,itl+1)=ET_E(it,itl);
        end
    end
%Vertical links 
jj=randi(125);
if V_E(jj,itl)==0
    if rand(1) < 0.5
        V_E(jj,itl+1)=1; %%0 -> +1
    else 
        V_E(jj,itl+1)=-1; %%0 -> -1
    end
end
if V_E(jj,itl)==1
    if rand(1) < 0.5
        V_E(jj,itl+1)=-1; %%+1 -> -1
    else 
        V_E(jj,itl+1)=0; %%+1 -> 0
    end
end
if V_E(jj,itl)==-1
    if rand(1) < 0.5
        V_E(jj,itl+1)=1; %%-1 -> 1
    else 
        V_E(jj,itl+1)=0; %%-1 -> 0
    end
end
%V_E(jj,itl+1)=1-V_E(jj,itl);
eh=0;ev=0; eu=0;
for i=1:125
    if H_E(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    u1=nbH(i,7);u2=nbH(i,8);u3=nbH(i,9);u4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    end
    if V_E(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    u1=nbV(i,7);u2=nbV(i,8);u3=nbV(i,9);u4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(V_E(i,itl+1));
    end
    if U_E(i,itl+1)~=0
    h1=nbU(i,1);h2=nbU(i,2);h3=nbU(i,3);h4=nbU(i,4);
    v1=nbU(i,5);v2=nbU(i,6);v3=nbU(i,7);v4=nbU(i,8);
    u1=nbU(i,9);u2=nbU(i,10);
    eu=eu+(E1+0.5*E2*(abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(U_E(i,itl+1));
    end
end
ET1_E(it,itl)=eh+ev+eu;
    Eref=ET_E(it,itl+1);  
    delE1=ET1_E(it,itl)-Eref; 
    if delE1 <=0
       V_E(jj,itl+1)=V_E(jj,itl+1);ET_E(it,itl+1)=ET1_E(it,itl);
    else if rand(1) <= exp((-delE1)/kBT(it,1))
            V_E(jj,itl+1)=V_E(jj,itl+1);ET_E(it,itl+1)=ET1_E(it,itl);
        else
            V_E(jj,itl+1)=V_E(jj,itl);ET_E(it,itl+1)=ET_E(it,itl+1);
        end
    end
%the same for U links 
j2=randi(125);
if U_E(j2,itl)==0
    if rand(1) < 0.5
        U_E(j2,itl+1)=1; %%0 -> +1
    else 
        U_E(j2,itl+1)=-1; %%0 -> -1
    end
end
if U_E(j2,itl)==1
    if rand(1) < 0.5
        U_E(j2,itl+1)=-1; %%+1 -> -1
    else 
        U_E(j2,itl+1)=0; %%+1 -> 0
    end
end
if U_E(j2,itl)==-1
    if rand(1) < 0.5
        U_E(j2,itl+1)=1; %%-1 -> 1
    else 
        U_E(j2,itl+1)=0; %%-1 -> 0
    end
end

eh=0;ev=0;eu=0;
for i=1:125
    if H_E(i,itl+1)~=0
    h1=nbH(i,1);h2=nbH(i,2);
    v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    u1=nbH(i,7);u2=nbH(i,8);u3=nbH(i,9);u4=nbH(i,10);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    end
    if V_E(i,itl+1)~=0
    h1=nbV(i,1);h2=nbV(i,2);h3=nbV(i,3);h4=nbV(i,4);
    v1=nbV(i,5);v2=nbV(i,6);
    u1=nbV(i,7);u2=nbV(i,8);u3=nbV(i,9);u4=nbV(i,10);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+abs(U_E(u3,itl+1))+abs(U_E(u4,itl+1))))*abs(V_E(i,itl+1));
    end
    if U_E(i,itl+1)~=0
    h1=nbU(i,1);h2=nbU(i,2);h3=nbU(i,3);h4=nbU(i,4);
    v1=nbU(i,5);v2=nbU(i,6);v3=nbU(i,7);v4=nbU(i,8);
    u1=nbU(i,9);u2=nbU(i,10);
    eu=eu+(E1+0.5*E2*(abs(U_E(u1,itl+1))+abs(U_E(u2,itl+1))+ ...
        abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))+ ...
        abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(U_E(i,itl+1));
    end
end

ET2_E(it,itl)=eh+ev+eu;
    Eref=ET_E(it,itl+1);  
    delE2=ET2_E(it,itl)-Eref; 
    if delE2 <=0
       U_E(j2,itl+1)=U_E(j2,itl+1);ET_E(it,itl+1)=ET2_E(it,itl);
    else if rand(1) <= exp((-delE2)/kBT(it,1))
            U_E(j2,itl+1)=U_E(j2,itl+1);ET_E(it,itl+1)=ET2_E(it,itl);
        else
            U_E(j2,itl+1)=U_E(j2,itl);ET_E(it,itl+1)=ET_E(it,itl+1);
        end
    end

H_T_E(:,it)=H_E(:,itl+1);
V_T_E(:,it)=V_E(:,itl+1);
U_T_E(:,it)=U_E(:,itl+1);
nH_E(it,1)=sum(abs(H_E(:,1)));
nV_E(it,1)=sum(abs(V_E(:,1)));
nU_E(it,1)=sum(abs(U_E(:,1)));
nTot_E(it,1)=nH_E(it,1)+nV_E(it,1)+nU_E(it,1);
nH_E(it,itl+1)=sum(abs(H_E(:,itl+1)));
nV_E(it,itl+1)=sum(abs(V_E(:,itl+1)));
nU_E(it,itl+1)=sum(abs(U_E(:,itl+1)));
nTot_E(it,itl+1)=nH_E(it,itl+1)+nV_E(it,itl+1)+nU_E(it,itl+1);
end
meanHVU_E(it,1)=mean(nH_E(it,itl-nave:itl+1));
meanHVU_E(it,2)=mean(nV_E(it,itl-nave:itl+1));
meanHVU_E(it,3)=mean(nU_E(it,itl-nave:itl+1));
meanHVU_E(it,4)=mean(nTot_E(it,itl-nave:itl+1));
end

figure;
subplot(2,2,1)
plot(step,nTot_E); ylim([0,125]);
legend(['T=10K'],...
['T=100K'],['T=200K'],...
['T=300K'],['T=400K'],['T=500K'],['T=600K'],['T=700K'],...
['T=800K'],...
['T=1200K'],...
['T=1600K'],...
['T=2000K'],...
['T=2500K']);
title(['E-field, T = 10-2500K,',num2str(itl),'MC steps total Links']); 

subplot(2,2,2)
plot(step,nH_E); ylim([0,125]);
legend(['T=10K'],...
['T=100K'],['T=200K'],...
['T=300K'],['T=400K'],['T=500K'],['T=600K'],['T=700K'],...
['T=800K'],...
['T=1200K'],...
['T=1600K'],...
['T=2000K'],...
['T=2500K']);
title(['E-field, T = 10-2500K,',num2str(itl),'MC steps H-Links']); 

subplot(2,2,3)
plot(step,nV_E); ylim([0,125]);
legend(['T=10K'],...
['T=100K'],['T=200K'],...
['T=300K'],['T=400K'],['T=500K'],['T=600K'],['T=700K'],...
['T=800K'],...
['T=1200K'],...
['T=1600K'],...
['T=2000K'],...
['T=2500K']);
title(['E-field, T = 10-2500K,',num2str(itl),'MC steps V-Links']); 

subplot(2,2,4)
plot(step,nU_E); ylim([0,125]);
legend(['T=10K'],...
['T=100K'],['T=200K'],...
['T=300K'],['T=400K'],['T=500K'],['T=600K'],['T=700K'],...
['T=800K'],...
['T=1200K'],...
['T=1600K'],...
['T=2000K'],...
['T=2500K']);
title(['E-field, T = 10-2500K,',num2str(itl),'MC steps U-Links']); 

% map to new network figure at 100-500K
% plot map to 3D cubic network scheme at 100K last step
figure;
Tii=2;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nU1=0;nU2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,1)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T_E(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T_E(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T_E(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T_E(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if U_T_E(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU1=nU1+1;
                end
                if U_T_E(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU2=nU2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nU1+nU2)]); 

% 300K
% plot map to 3D cubic network scheme at 300K last step
Tii=4;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nU1=0;nU2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,2)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T_E(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T_E(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T_E(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T_E(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if U_T_E(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU1=nU1+1;
                end
                if U_T_E(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU2=nU2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nU1+nU2)]); 

% 400K 
% plot map to 3D cubic network scheme at 400K last step
Tii=5;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nU1=0;nU2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,3)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T_E(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T_E(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T_E(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T_E(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if U_T_E(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU1=nU1+1;
                end
                if U_T_E(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU2=nU2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nU1+nU2)]); 

% 500K 
% plot map to 3D cubic network scheme at 500K last step
Tii=6;
ti=1;tj=1;tk=1;
nH1=0;nH2=0;nV1=0;nV2=0;nU1=0;nU2=0;
N=10;
a=-pi:pi/2:pi;
ph=pi/4;
X = 1:N; Y = 1:N; Z = 1:N;
subplot(2,2,4)
for k = 1:N
    for j=1:N
        for i=1:N
    clear col 
            if mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==0
               col=[1.0 0.0 0.0]; %red NP
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
            elseif mod(i,2)==1 && mod(j,2)==0 && mod(k,2)==0
                if H_T_E(ti,Tii)==1
                col=[0.75 0.75 1.0]; %blue H light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH1=nH1+1;end
                if H_T_E(ti,Tii)==-1
                col=[0.3 0.3 1.0]; %blue H dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nH2=nH2+1;
                end
                ti=ti+1;
            elseif mod(i,2)==0 && mod(j,2)==1 && mod(k,2)==0 
                if V_T_E(tj,Tii)==1
                col=[1.0 1.0 0.75]; %yellow V light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV1=nV1+1;
                end
                if V_T_E(tj,Tii)==-1
                col=[1.0 0.9 0.01]; %yellow V dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nV2=nV2+1;
                end
                tj=tj+1;
            elseif mod(i,2)==0 && mod(j,2)==0 && mod(k,2)==1 
                if U_T_E(tk,Tii)==1
                col=[0.75 0.75 0.75]; %grey U light
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU1=nU1+1;
                end
                if U_T_E(tk,Tii)==-1
                col=[0.5 0.5 0.5]; %grey U dark
    x=X(i)+[cos(a+ph); cos(a+ph)]/cos(ph)/2.5;
    y=Y(j)+[sin(a+ph); sin(a+ph)]/sin(ph)/2.5;
    z=Z(k)+[-ones(size(a)); ones(size(a))]/2.5;
    surf(x, y, z, 'FaceColor',col,'EdgeColor','none') 
    hold on
    patch(x', y', z', col,'EdgeColor','none') 
                nU2=nU2+1;
                end
                tk=tk+1;
            end 
        end
    end
end
xlim([0,11]);ylim([0,11]);zlim([0,11]);
title(['T=',num2str(TT(Tii)),'K:H=',num2str(nH1+nH2),',V=',num2str(nV1+nV2),',U=',num2str(nU1+nU2)]); 

%% MFT approximation
clear result
clc
del=0.001;
mx=0:del:1; % average network per box
kB = 1.380649E-23; % in J/K
it=0;
%for T =50:50:3000 % in K 
for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
it=it+1;
NA=6.02214076E+23;  % #/mol
E1=5*4186.798188/NA; % connection energy kcal/mol to J
E2=2*4186.798188/NA; % neighbor correlation energy
EE=-6*4186.798188/NA;%-7*4186.798188/NA; %in kcal/mol

kBT(it,1)=kB*T/4186.798188*NA;
Ti(it,1)=T;
%D=80;
lm=size(mx,2);
for i = 1:lm
    Z(1,i)=(1+6*exp(-1/(kB*T)*(E1+6*mx(1,i)*E2))+12*exp(-1/(kB*T)*2*(E1+6*mx(1,i)*E2))+8*exp(-1/(kB*T)*3*(E1+6*mx(1,i)*E2)));
    my(1,i)=(2*exp(-1/(kB*T)*(E1+6*mx(1,i)*E2))+8*exp(-1/(kB*T)*2*(E1+6*mx(1,i)*E2))+8*exp(-1/(kB*T)*3*(E1+6*mx(1,i)*E2)))/Z(1,i); %network in Horizontal
end

for i = 2:lm
   if (my(1,i-1)-mx(1,i-1))*(my(1,i)-mx(1,i))<=0
       result(it,1)=(mx(1,i)+mx(1,i-1)+my(1,i)+my(1,i-1))/4; %Horizontal
       result(it,2)=3*result(it,1);%Total
   end
end

% with E-field
%EE=-Efield*e*200*L*1E-10*1e6;e=1.602176634E-19;Efield=0.05;L=100; % E-field energy
 % in C
 % in V/m
 %linker length in A
for i = 1:lm
    mH(1,i)=mx(1,i);
    for j = 1:lm
    mV(1,j)=mx(1,j);
    ZE(j,i)=1+exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2+EE))+exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2-EE))+...
        exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+...
        exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+ exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE));

    mx1(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2+EE))+exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)))/ZE(j,i);

    mx2(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2-EE))+exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)))/ZE(j,i);

    my1(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)))/ZE(j,i);
    
    my2(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+4*mV(1,i))*E2))+...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+8*mV(1,i))*E2)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2+EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)) + ...
        exp(-1/(kB*T)*(3*E1+(6*mH(1,i)+12*mV(1,i))*E2-EE)))/ZE(j,i);

    mHs(j,i)=(mx1(j,i)+mx2(j,i));
    mVs(j,i)=(my1(j,i)+my2(j,i));
    del(j,i)=abs(mHs(j,i)-mH(1,i))+2*abs(mVs(j,i)-mV(1,j));
    end
end
del0(it,1)=min(min(del));
for j = 1:lm
    for i = 1:lm
       if del(j,i)==del0(it,1)
result(it,3)=mHs(j,i);
result(it,4)=mHs(j,i)+mVs(j,i)*2;break;
       end
   end
end
end
%
figure1 = figure;% plot minimize mx my MFT surface
axes1 = axes('Parent',figure1);
hold(axes1,'on');
surf(mx,my,del,'Parent',axes1,'EdgeColor','none');
view(axes1,[-190.424929686292 90]);
hold(axes1,'off');
set(axes1,'Colormap',[0 1 1;0.0416666666666667 1 0.958333333333333;0.0833333333333333 1 0.916666666666667;0.125 1 0.875;0.166666666666667 1 0.833333333333333;0.208333333333333 1 0.791666666666667;0.25 1 0.75;0.291666666666667 1 0.708333333333333;0.333333333333333 1 0.666666666666667;0.375 1 0.625;0.416666666666667 1 0.583333333333333;0.458333333333333 1 0.541666666666667;0.5 1 0.5;0.541666666666667 1 0.458333333333333;0.583333333333333 1 0.416666666666667;0.625 1 0.375;0.666666666666667 1 0.333333333333333;0.708333333333333 1 0.291666666666667;0.75 1 0.25;0.791666666666667 1 0.208333333333333;0.833333333333333 1 0.166666666666667;0.875 1 0.125;0.916666666666667 1 0.0833333333333334;0.958333333333333 1 0.0416666666666666;1 1 0;0.988492134099517 0.957328662567213 0.00432800365526536;0.976984268199034 0.914657325134426 0.00865600731053072;0.965476402298552 0.871985987701638 0.0129840109657961;0.953968536398069 0.829314650268851 0.0173120146210614;0.942460670497586 0.786643312836064 0.0216400182763268;0.930952804597103 0.743971975403277 0.0259680219315921;0.91944493869662 0.70130063797049 0.0302960255868575;0.907937072796138 0.658629300537703 0.0346240292421229;0.896429206895655 0.615957963104915 0.0389520328973882;0.884921340995172 0.573286625672128 0.0432800365526536;0.873413475094689 0.530615288239341 0.0476080402079189;0.861905609194206 0.487943950806554 0.0519360438631843;0.850397743293724 0.445272613373767 0.0562640475184497;0.838889877393241 0.40260127594098 0.060592051173715;0.827382011492758 0.359929938508192 0.0649200548289804;0.815874145592275 0.317258601075405 0.0692480584842457;0.804366279691792 0.274587263642618 0.0735760621395111;0.79285841379131 0.231915926209831 0.0779040657947764;0.781350547890827 0.189244588777044 0.0822320694500418;0.768294506943918 0.189020458627642 0.103806191070405;0.755238465997008 0.188796328478239 0.125380312690769;0.742182425050099 0.188572198328837 0.146954434311132;0.72912638410319 0.188348068179435 0.168528555931496;0.716070343156281 0.188123938030033 0.190102677551859;0.703014302209371 0.187899807880631 0.211676799172223;0.689958261262462 0.187675677731229 0.233250920792586;0.676902220315553 0.187451547581826 0.25482504241295;0.663846179368644 0.187227417432424 0.276399164033313;0.650790138421734 0.187003287283022 0.297973285653677;0.637734097474825 0.18677915713362 0.319547407274041;0.624678056527916 0.186555026984218 0.341121528894404;0.611622015581007 0.186330896834816 0.362695650514768;0.598565974634097 0.186106766685413 0.384269772135131;0.585509933687188 0.185882636536011 0.405843893755495;0.572453892740279 0.185658506386609 0.427418015375858;0.55939785179337 0.185434376237207 0.448992136996222;0.54634181084646 0.185210246087805 0.470566258616585;0.533285769899551 0.184986115938403 0.492140380236949;0.520229728952642 0.184761985789 0.513714501857312;0.507173688005733 0.184537855639598 0.535288623477676;0.494117647058824 0.184313725490196 0.556862745098039;0.488434023836375 0.198599273929414 0.547227864836937;0.482750400613926 0.212884822368631 0.537592984575834;0.477066777391477 0.227170370807849 0.527958104314732;0.471383154169029 0.241455919247066 0.518323224053629;0.46569953094658 0.255741467686284 0.508688343792527;0.460015907724131 0.270027016125501 0.499053463531424;0.454332284501682 0.284312564564719 0.489418583270322;0.448648661279234 0.298598113003936 0.479783703009219;0.442965038056785 0.312883661443154 0.470148822748116;0.437281414834336 0.327169209882371 0.460513942487014;0.431597791611888 0.341454758321589 0.450879062225911;0.425914168389439 0.355740306760806 0.441244181964809;0.42023054516699 0.370025855200024 0.431609301703706;0.414546921944541 0.384311403639242 0.421974421442604;0.408863298722093 0.398596952078459 0.412339541181501;0.403179675499644 0.412882500517677 0.402704660920399;0.397496052277195 0.427168048956894 0.393069780659296;0.391812429054747 0.441453597396112 0.383434900398194;0.386128805832298 0.455739145835329 0.373800020137091;0.380445182609849 0.470024694274547 0.364165139875989;0.3747615593874 0.484310242713764 0.354530259614886;0.369077936164952 0.498595791152982 0.344895379353784;0.363394312942503 0.512881339592199 0.335260499092681;0.357710689720054 0.527166888031417 0.325625618831579;0.352027066497605 0.541452436470634 0.315990738570476;0.346343443275157 0.555737984909852 0.306355858309374;0.340659820052708 0.57002353334907 0.296720978048271;0.334976196830259 0.584309081788287 0.287086097787168;0.329292573607811 0.598594630227504 0.277451217526066;0.323608950385362 0.612880178666722 0.267816337264963;0.317925327162913 0.62716572710594 0.258181457003861;0.312241703940464 0.641451275545157 0.248546576742758;0.306558080718016 0.655736823984375 0.238911696481656;0.300874457495567 0.670022372423592 0.229276816220553;0.295190834273118 0.68430792086281 0.219641935959451;0.289507211050669 0.698593469302027 0.210007055698348;0.283823587828221 0.712879017741245 0.200372175437246;0.278139964605772 0.727164566180462 0.190737295176143;0.272456341383323 0.74145011461968 0.181102414915041;0.266772718160875 0.755735663058897 0.171467534653938;0.261089094938426 0.770021211498115 0.161832654392836;0.256794409922144 0.77765514122189 0.160293281553758;0.252499724905862 0.785289070945666 0.15875390871468;0.24820503988958 0.792923000669441 0.157214535875602;0.243910354873299 0.800556930393216 0.155675163036524;0.239615669857017 0.808190860116992 0.154135790197446;0.235320984840735 0.815824789840767 0.152596417358368;0.231026299824453 0.823458719564542 0.15105704451929;0.226731614808172 0.831092649288318 0.149517671680211;0.22243692979189 0.838726579012093 0.147978298841133;0.218142244775608 0.846360508735869 0.146438926002055;0.217934850334018 0.852360168293936 0.142009089157246;0.217727455892429 0.858359827852004 0.137579252312436;0.217520061450839 0.864359487410072 0.133149415467626;0.21731266700925 0.87035914696814 0.128719578622816;0.21710527256766 0.876358806526208 0.124289741778006;0.216897878126071 0.882358466084275 0.119859904933196;0.216690483684481 0.888358125642343 0.115430068088386;0.216483089242891 0.894357785200411 0.111000231243577;0.216275694801302 0.900357444758479 0.106570394398767;0.216068300359712 0.906357104316547 0.102140557553957;0.215860905918123 0.912356763874615 0.0977107207091469;0.215653511476533 0.918356423432683 0.0932808838643371;0.215446117034943 0.92435608299075 0.0888510470195272;0.215238722593354 0.930355742548818 0.0844212101747173;0.215031328151764 0.936355402106886 0.0799913733299074;0.214823933710175 0.942355061664954 0.0755615364850976;0.214616539268585 0.948354721223022 0.0711316996402877;0.214409144826996 0.954354380781089 0.0667018627954779;0.214201750385406 0.960354040339157 0.0622720259506681;0.213994355943816 0.966353699897225 0.0578421891058581;0.213786961502227 0.972353359455293 0.0534123522610482;0.213579567060637 0.978353019013361 0.0489825154162383;0.213372172619048 0.984352678571429 0.0445526785714285;0.222222222222222 1 0;0.259259259259259 1 0;0.296296296296296 1 0;0.333333333333333 1 0;0.37037037037037 1 0;0.407407407407407 1 0;0.444444444444444 1 0;0.481481481481481 1 0;0.518518518518518 1 0;0.555555555555556 1 0;0.592592592592593 1 0;0.62962962962963 1 0;0.666666666666667 1 0;0.703703703703704 1 0;0.740740740740741 1 0;0.777777777777778 1 0;0.814814814814815 1 0;0.851851851851852 1 0;0.888888888888889 1 0;0.925925925925926 1 0;0.962962962962963 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;0.981908333333333 0.979825 0.0134;0.963816666666667 0.95965 0.0268;0.945725 0.939475 0.0402;0.927633333333333 0.9193 0.0536;0.909541666666667 0.899125 0.067;0.89145 0.87895 0.0804;0.873358333333333 0.858775 0.0938;0.855266666666667 0.8386 0.1072;0.837175 0.818425 0.1206;0.819083333333333 0.79825 0.134;0.800991666666667 0.778075 0.1474;0.7829 0.7579 0.1608;0.7945 0.7554 0.157;0.806 0.7529 0.1546;0.743623529411765 0.639182745098039 0.235052549019608;0.681247058823529 0.525465490196078 0.315505098039216;0.618870588235294 0.411748235294118 0.395957647058824;0.556494117647059 0.298030980392157 0.476410196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.53652513368984 0.233785204991087 0.527257040998217;0.578932620320856 0.283256684491979 0.497651336898396;0.621340106951872 0.33272816399287 0.468045632798574;0.663747593582888 0.382199643493761 0.438439928698752;0.706155080213904 0.431671122994652 0.40883422459893;0.74856256684492 0.481142602495544 0.379228520499109;0.790970053475936 0.530614081996435 0.349622816399287;0.833377540106952 0.580085561497326 0.320017112299465;0.875785026737968 0.629557040998218 0.290411408199643;0.918192513368984 0.679028520499109 0.260805704099822;0.9606 0.7285 0.2312;0.9689 0.7292 0.2373;0.977 0.7304 0.2418;0.9842 0.733 0.2446;0.99 0.7365 0.2429;0.9946 0.7407 0.2394;0.9966 0.7458 0.2351;0.9971 0.7513 0.2309;0.9972 0.7569 0.2267;0.9971 0.7626 0.2224;0.9969 0.7683 0.2181;0.9966 0.774 0.2138;0.9962 0.7798 0.2095;0.9957 0.7856 0.2053;0.9949 0.7915 0.2012;0.9938 0.7974 0.1974;0.9923 0.8034 0.1939;0.9906 0.8095 0.1906;0.9885 0.8156 0.1875;0.9861 0.8218 0.1846;0.9835 0.828 0.1817;0.9807 0.8342 0.1787;0.9778 0.8404 0.1757;0.9748 0.8467 0.1726;0.972 0.8529 0.1695;0.9694 0.8591 0.1665;0.9671 0.8654 0.1636;0.9651 0.8716 0.1608;0.9634 0.8778 0.1582;0.9619 0.884 0.1557;0.9608 0.8902 0.1532;0.9601 0.8963 0.1507;0.9596 0.9023 0.148;0.9595 0.9084 0.145;0.9597 0.9143 0.1418;0.9601 0.9203 0.1382;0.9608 0.9262 0.1344;0.9618 0.932 0.1304;0.9629 0.9379 0.1261;0.9642 0.9437 0.1216;0.9657 0.9494 0.1168;0.9674 0.9552 0.1116;0.9692 0.9609 0.1061;0.9711 0.9667 0.1001;0.973 0.9724 0.0938;0.9749 0.9782 0.0872;0.9769 0.9839 0.0805]);
colorbar(axes1);
%% plot and compare Networks
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plotMDT1 = plot([300 800	1200	1600	2000	2500], [4.1	7.1	44.0	78.4	99.2	102.7],'p-k','MarkerSize',18,'Parent',axes1);
plotMDT2 = plot([300	800	1200	1600	2000	2500], [1.0	4.0	16.5	19.0	17.6	11.2],'p-k','MarkerSize',2,'Parent',axes1);
plotMDE1 = plot ([10 100 200 300 400 500], [3.0	11.0	14.8	23.2	32.6	57.4],'p', 'Color',[1 0 1],'MarkerSize',18,'Parent',axes1);
plotMDE2 = plot ([10 100 200 300 400 500], [2.0	8.0	14.8	19.2	30.6	55.4],'p', 'Color',[1 0 1],'MarkerSize',2,'Parent',axes1);
plotMFT1 = plot(Ti,result(:,2)*125,'o','MarkerSize',18,'Parent',axes1);
plotMFT2 = plot(Ti,result(:,1)*125,'o','MarkerSize',2,'Parent',axes1);
plotMFE1 = plot(Ti,result(:,4)*125,'o','MarkerSize',18,'Parent',axes1);
plotMFE2 = plot(Ti,result(:,3)*125,'o','MarkerSize',2,'Parent',axes1);
%plot3 = plot ([10 200	300	400], [0 5 11 14],'p', 'Color',[1 0 1],'MarkerSize',18,'Parent',axes1);
plotMCT1 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHVZ(:,4),'s','MarkerSize',18,'Parent',axes1);
plotMCT2 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHVZ(:,1),'s','MarkerSize',2,'Parent',axes1);
plotMCE1 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHVU_E(:,4),'^','MarkerSize',18,'Parent',axes1);
plotMCE2 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHVU_E(:,1),'^','MarkerSize',2,'Parent',axes1);
%plot1 = plot(Ti,result*125,'Parent',axes1);

set(plotMDT1(1),'DisplayName','DPD T-effect, Total','LineWidth',1.5,'LineStyle','-');
set(plotMDT2(1),'DisplayName','DPD T-effect, x-link','LineWidth',1.5,'LineStyle','--');
set(plotMDE1(1),'DisplayName','DPD E-field, Total','LineWidth',1.5,'LineStyle','-');
set(plotMDE2(1),'DisplayName','DPD E-field, x-link','LineWidth',1.5,'LineStyle','--');
set(plotMFT1(1),'DisplayName','MFT T-effect, Total','LineWidth',1.5,'LineStyle','-','Color',[0 0 1]);
set(plotMFT2(1),'DisplayName','MFT T-effect, x-link','LineWidth',1.5,'LineStyle','--','Color',[0 0 1]);
set(plotMFE1(1),'DisplayName','MFT E-field, Total','LineWidth',1.5,'LineStyle','-','Color',[1 0 0]);
set(plotMFE2(1),'DisplayName','MFT E-field, x-link','LineWidth',1.5,'LineStyle','--','Color',[1 0 0]);

%set(plot3(1),'DisplayName','MD, E-field 0.001 V/Å, Total','LineWidth',3,'LineStyle','none');
%set(plot4(1),'DisplayName','MC, No E-field, X','LineWidth',1.5,'LineStyle','--','Color',[0 1 1]);
%set(plot4(2),'DisplayName','MC, No E-field, Y','LineWidth',1.5,'LineStyle','-.','Color',[0 1 1]);
set(plotMCT1(1),'DisplayName','MC T-effect, Total','LineWidth',1.5,'LineStyle','-','Color',[0 1 1]);
set(plotMCT2(1),'DisplayName','MC T-effect, x-link','LineWidth',1.5,'LineStyle','--','Color',[0 1 1]);
%set(plot5(1),'DisplayName','MC, E-field, X','LineWidth',1.5,'LineStyle','--','Color',[1 0 0.7]);
%set(plot5(2),'DisplayName','MC, E-field, Y','LineWidth',1.5,'LineStyle','-.','Color',[1 0 0.7]);
set(plotMCE1(1),'DisplayName','MC E-field, Total','LineWidth',1.5,'LineStyle','-','Color',[1 0 0.7]);
set(plotMCE2(1),'DisplayName','MC E-field, x-link','LineWidth',1.5,'LineStyle','--','Color',[1 0 0.7]);
box(axes1,'on');
hold(axes1,'off');
legend1 = legend(axes1,'show');
ylabel('Likers (#)');
%ylim([0,120]);
xlabel('Temperature (K)');
set(axes1,'FontSize',20);
title(['3D Cubic network: Ising model MC & MFT vs DPD']); 
result*125;
