%MCMC model 2D Square 10*10 network connection H-> V|^
%Three states model H(0 or +1, -1), V(0 or +1, -1)
%Considering impossible links: each AuNP only has one polymer reach out
clc
clear all
nH1=0;nH2=0;nV1=0;nV2=0;
for i=1:20
    for j=1:20
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
title(['AuNP 2D array: nH=0, nV=0']); 
ti=0;tj=0;
for i=1:20
    for j=1:20
        if mod(i,2)==0 && mod(j,2)==1 
            ti=i/2; tj=ceil(j/2);
            iniH(ti,tj)=0;
                if rand(1)<1
                    if rand(1)<0.5
                A(i,j,1)=128;A(i,j,2)=128;A(i,j,3)=255;%Horizontal+1  blue
                nH1=nH1+1;
                iniH(ti,tj)=1;
                else
                A(i,j,1)=128;A(i,j,2)=128;A(i,j,3)=255;%Horizontal-1  blue
                nH2=nH2+1;
                iniH(ti,tj)=-1;
                    end
                end
        elseif mod(i,2)==1 && mod(j,2)==0 
            ti=ceil(i/2); tj=j/2;
            iniV(ti,tj)=0;
            r1=rand(1);r2=rand(1);
            if r1<1
                if r2<0.5
                A(i,j,1)=255;A(i,j,2)=255;A(i,j,3)=25.5;% Vertical+1  yellow
                nV1=nV1+1;
                iniV(ti,tj)=1;
            else
                A(i,j,1)=255;A(i,j,2)=255;A(i,j,3)=25.5;%Vertical-1  yellow
                nV2=nV2+1;
                iniV(ti,tj)=-1;
                end
            end
        end
    end
end
%A(11,5,:)=0;
subplot(2,2,2)
image(uint8(A));
title(['Random initial state: nH=',num2str(nH1+nH2),', nV=',num2str(nV1+nV2)]); 

% MC simulation 
ti=0;tj=0;p0=0.0;
for i=1:20
    for j=1:20
        if mod(i,2)==0 && mod(j,2)==1 
            ti=i/2; tj=ceil(j/2);
            iniH(ti,tj)=0;
                if rand(1)<p0
                    if rand(1)<0.5
                A(i,j,1)=192;A(i,j,2)=192;A(i,j,3)=255;%Horizontal+1 light blue
                nH1=nH1+1;
                iniH(ti,tj)=1;
                else
                A(i,j,1)=76.8;A(i,j,2)=76.8;A(i,j,3)=255;%Horizontal-1 dark blue
                nH2=nH2+1;
                iniH(ti,tj)=-1;
                    end
                end
        elseif mod(i,2)==1 && mod(j,2)==0 
            ti=ceil(i/2); tj=j/2;
            iniV(ti,tj)=0;
            r1=rand(1);r2=rand(1);
            if r1<p0
                if r2<0.5
                A(i,j,1)=255;A(i,j,2)=255;A(i,j,3)=192;% Vertical+1 light yellow
                nV1=nV1+1;
                iniV(ti,tj)=1;
            else
                A(i,j,1)=255;A(i,j,2)=230.4;A(i,j,3)=2.55;%Vertical-1 dark yellow
                nV2=nV2+1;
                iniV(ti,tj)=-1;
                end
            end
        end
    end
end
%H(:,1) and V(:,1) initialized random link
for i=1:10
    for j=1:10
        k=(i-1)*10+j;
        H(k,1)=iniH(i,j); %links-horizontal
        V(k,1)=iniV(j,i); %links-vertical
        H_E(k,1)=iniH(i,j); %links-horizontal E-field
        V_E(k,1)=iniV(j,i); %links-vertical E-field
    end
end
% neighbor list for 100 H and 100 V ** first row: -> H1-H10 and first column: |^ V1-V91
clear nbH nbV
for i=1:100
    nbH(i,1)=i-1;nbH(i,2)=i+1; %2 H neighbor column 1,2
    if mod(i,10)==1
        nbH(i,1)=i-1+10;end    %2 H neighbor 1,2
    if mod(i,10)==0
        nbH(i,2)=i+1-10;end    %2 H neighbor 1,2
    nbH(i,3)=(mod(i,10)-1)*10+floor(i/10);nbH(i,4)=(mod(i,10)-1)*10+floor(i/10)+1; %4 V neighbor column 3-6
    nbH(i,5)=(mod(i,10)-0)*10+floor(i/10);nbH(i,6)=(mod(i,10)-0)*10+floor(i/10)+1; 
    if mod(i,10)==0
    nbH(i,3)=(mod(i,10)-1)*10+100+floor(i/10)-1;nbH(i,4)=(mod(i,10)-1)*10+100+floor(i/10); 
    nbH(i,5)=(mod(i,10)-0)*10+floor(i/10)-1;nbH(i,6)=(mod(i,10)-0)*10+floor(i/10); 
        if i==10
            nbH(i,3)=(mod(i,10))*10+100+floor(i/10)-1; nbH(i,5)=(mod(i,10)+1)*10+floor(i/10)-1;end
    end  %4 V neighbor 3-6
    if i < 10
    nbH(i,3)=(mod(i,10)-1)*10+floor(i/10)+10;nbH(i,4)=(mod(i,10)-1)*10+floor(i/10)+1; %4 V neighbor column 3-6
    nbH(i,5)=(mod(i,10)-0)*10+floor(i/10)+10;nbH(i,6)=(mod(i,10)-0)*10+floor(i/10)+1;  %4 V neighbor 3&5        
    end 
    %if i>90
    %nbH(i,4)=(mod(i,10)-2)*10+ceil(i/10)+1-10; 
    %nbH(i,6)=(mod(i,10)-1)*10+ceil(i/10)+1-10; end         %4 V neighbor 3-6
   % if i==100
   % nbH(i,4)=(mod(i,10)-2)*10+100+ceil(i/10)+1-10; 
   % nbH(i,6)=(mod(i,10)-1)*10+100+ceil(i/10)+1-10; end     %4 V neighbor 3-6
   % if i==91
   %     nbH(i,4)=(mod(i,10)-2)*10+100+ceil(i/10)+1-10; end %4 V neighbor 3-6
    nbV(i,:)=nbH(i,:); %  the same as V neighbors
    % nbVT(i,:)=nbH(i,:);H neighbors are the same as V neighbors, 2 V neighbor column 1-2, and 4 H neighbor column 3-6
end
% transpose V neighbor ids
%for i=1:10
%    for j=1:10
%        idV1(i,j)=(i-1)*10+j; %new V linker ids
%        idV0(j,i)=(i-1)*10+j; %old V linker ids 
%    end
%end
%for i=1:10
%    for j=1:10
%        k1=idV1(i,j);
%        k0=idV0(i,j);
%        nbV(k1,:)=nbVT(k0,:);
%    end
%end

%% Temperature effect: E1 linking energy and E2 correlation energy
kB = 1.380649E-23; % in J/K
NA=6.02214076E+23;  % #/mol
E1=5;%7*4186.798188/NA; % connection energy kcal/mol to J
E2=2;%3*4186.798188/NA; % neighbor correlation energy
EE=-10;%-4.5, -7*4186.798188/NA; %in kcal/mol

eh=0;ev=0;
for i=1:100
    h1=nbH(i,1);h2=nbH(i,2);v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6); % 
    %    V2   V4
    % H1-O-H?-O-H2
    %    V1   V3
    eh=eh+(E1+0.5*E2*(abs(H(h1,1))+abs(H(h2,1))+abs(V(v1,1))+abs(V(v2,1))+abs(V(v3,1))+abs(V(v4,1))))*abs(H(i,1));
    v1=nbV(i,1);v2=nbV(i,2);h1=nbV(i,3);h2=nbV(i,4);h3=nbV(i,5);h4=nbV(i,6);
    ev=ev+(E1+0.5*E2*(abs(V(v1,1))+abs(V(v2,1))+abs(H(h1,1))+abs(H(h2,1))+abs(H(h3,1))+abs(H(h4,1))))*abs(V(i,1));
    %Z(i,itl)=1+2*exp(-EH(i,itl)/kBT(it,1))+2*exp(-EV(i,itl)/kBT(it,1))+4*exp(-(EH(i,itl)+EV(i,itl))/kBT(it,1));
    %Z(1,i)=(1+4*exp(-1/(kB*T)*(E1+6*mx(1,i)*E2))+4*exp(-1/(kB*T)*2*(E1+6*mx(1,i)*E2)));
end
et=eh+ev;

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
%clear EH0 EV0 ET0 Eref and specify three conditions
ii=randi(100);
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
eh=0;ev=0;
    for i=1:100
    h1=nbH(i,1);h2=nbH(i,2);v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    eh=eh+(E1+0.5*E2*(abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))))*abs(H(i,itl+1));
    v1=nbV(i,1);v2=nbV(i,2);h1=nbV(i,3);h2=nbV(i,4);h3=nbV(i,5);h4=nbV(i,6);
    ev=ev+(E1+0.5*E2*(abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))))*abs(V(i,itl+1));
    end
ET0(it,itl)=ev+eh;
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
%the same for Vertical links 
%clear EH0 EV0 ET0 Eref 
%specify 3 different types of links
jj=randi(100);
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
eh=0;ev=0; 
    for i=1:100
    h1=nbH(i,1);h2=nbH(i,2);v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    %    V2
    % H3-O-H4
    %    V?
    % H1-O-H2 
    %    V1
    eh=eh+(E1+0.5*E2*(abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(V(v3,itl+1))+abs(V(v4,itl+1))))*abs(H(i,itl+1));
    v1=nbV(i,1);v2=nbV(i,2);h1=nbV(i,3);h2=nbV(i,4);h3=nbV(i,5);h4=nbV(i,6);
    ev=ev+(E1+0.5*E2*(abs(V(v1,itl+1))+abs(V(v2,itl+1))+abs(H(h1,itl+1))+abs(H(h2,itl+1))+abs(H(h3,itl+1))+abs(H(h4,itl+1))))*abs(V(i,itl+1));
    end
ET1(it,itl)=eh+ev;
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

%    if R1 <= 4*exp(-(EH(i,itl)+EV(i,itl))/kBT(it,1))/Z(i,itl)
%        H(i,itl+1)=1;V(i,itl+1)=1;
%    else if R1<=4*exp(-(EH(i,itl)+EV(i,itl))/kBT(it,1))/Z(i,itl)+2*exp(-EH(i,itl)/kBT(it,1))/Z(i,itl)
%            H(i,itl+1)=1;V(i,itl+1)=0;
%        else if R1<=4*exp(-(EH(i,itl)+EV(i,itl))/kBT(it,1))/Z(i,itl)+2*exp(-EH(i,itl)/kBT(it,1))/Z(i,itl)+2*exp(-EV(i,itl)/kBT(it,1))/Z(i,itl)
%                H(i,itl+1)=0;V(i,itl+1)=1;
%            else
%                H(i,itl+1)=0;V(i,itl+1)=0;
%            end
%        end
%end
nH(it,1)=sum(abs(H(:,1)));
nV(it,1)=sum(abs(V(:,1)));
nTot(it,1)=nH(it,1)+nV(it,1);
H_T(:,it)=H(:,itl+1);
V_T(:,it)=V(:,itl+1);
nH(it,itl+1)=sum(abs(H(:,itl+1)));
nV(it,itl+1)=sum(abs(V(:,itl+1)));
nTot(it,itl+1)=nH(it,itl+1)+nV(it,itl+1);
end
meanHV(it,1)=mean(nH(it,itl-nave:itl+1));
meanHV(it,2)=mean(nV(it,itl-nave:itl+1));
meanHV(it,3)=mean(nTot(it,itl-nave:itl+1));
end
subplot(2,2,3)
plot(step,nTot); ylim([0,50]);
legend(['T=  10: nH=',num2str(meanHV(1,1),'%0.1f'),', nV=',num2str(meanHV(1,2),'%0.1f'),', nT=',num2str(meanHV(1,3),'%0.1f'),],...
['T=100: nH=',num2str(meanHV(2,1),'%0.1f'),', nV=',num2str(meanHV(2,2),'%0.1f'),', nT=',num2str(meanHV(2,3),'%0.1f'),],...
['T=200: nH=',num2str(meanHV(3,1),'%0.1f'),', nV=',num2str(meanHV(3,2),'%0.1f'),', nT=',num2str(meanHV(3,3),'%0.1f'),],...
['T=300: nH=',num2str(meanHV(4,1),'%0.1f'),', nV=',num2str(meanHV(4,2),'%0.1f'),', nT=',num2str(meanHV(4,3),'%0.1f'),],...
['T=400: nH=',num2str(meanHV(5,1),'%0.1f'),', nV=',num2str(meanHV(5,2),'%0.1f'),', nT=',num2str(meanHV(5,3),'%0.1f'),],...
['T=500: nH=',num2str(meanHV(6,1),'%0.1f'),', nV=',num2str(meanHV(6,2),'%0.1f'),', nT=',num2str(meanHV(6,3),'%0.1f'),],...
['T=600: nH=',num2str(meanHV(7,1),'%0.1f'),', nV=',num2str(meanHV(7,2),'%0.1f'),', nT=',num2str(meanHV(7,3),'%0.1f'),],...
['T=700: nH=',num2str(meanHV(8,1),'%0.1f'),', nV=',num2str(meanHV(8,2),'%0.1f'),', nT=',num2str(meanHV(8,3),'%0.1f'),],...
['T=800: nH=',num2str(meanHV(9,1),'%0.1f'),', nV=',num2str(meanHV(9,2),'%0.1f'),', nT=',num2str(meanHV(9,3),'%0.1f'),],...
['T=1200: nH=',num2str(meanHV(10,1),'%0.1f'),', nV=',num2str(meanHV(10,2),'%0.1f'),', nT=',num2str(meanHV(10,3),'%0.1f'),],...
['T=1600: nH=',num2str(meanHV(11,1),'%0.1f'),', nV=',num2str(meanHV(11,2),'%0.1f'),', nT=',num2str(meanHV(11,3),'%0.1f'),],...
['T=2000: nH=',num2str(meanHV(12,1),'%0.1f'),', nV=',num2str(meanHV(12,2),'%0.1f'),', nT=',num2str(meanHV(12,3),'%0.1f'),],...
['T=2500: nH=',num2str(meanHV(13,1),'%0.1f'),', nV=',num2str(meanHV(13,2),'%0.1f'),', nT=',num2str(meanHV(13,3),'%0.1f'),]);
title(['MC T = 10-2500K, ',num2str(itl),' MC steps']); 

% map to new network figure at 2000K last step
Tii=12;%2000K
for i = 1:10
    for j=1:10
        Ht(i,j)=H_T((i-1)*10+j,Tii);
        ai=(11-i)*2;aj=j*2-1;
        if Ht(i,j)==1
            A(ai,aj,1)=192;A(ai,aj,2)=192;A(ai,aj,3)=255;%Horizontal+1 light blue
        else if Ht(i,j)==-1
                A(ai,aj,1)=76.8;A(ai,aj,2)=76.8;A(ai,aj,3)=255;%Horizontal-1 dark blue
            else
            A(ai,aj,1)=255;A(ai,aj,2)=255;A(ai,aj,3)=255;%blank
            end
        end
        
        Vt(j,i)=V_T((i-1)*10+j,Tii);
        bi=(11-j)*2-1;bj=i*2;
        if Vt(j,i)==1
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=192;% Vertical+1 light
        else if Vt(j,i)== -1
                A(bi,bj,1)=255;A(bi,bj,2)=230.4;A(bi,bj,3)=2.55;%Vertical-1 dark
            else 
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=255;%blank
            end
        end
    end
end

subplot(2,2,4)
image(uint8(A));
title(['Last step: T=2000K, nH=',num2str(nH(Tii,itl+1)),', nV=',num2str(nV(Tii,itl+1))]); 
%% E-field effect: EE is the potential energy from E-field
it=0;
for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
    it=it+1;
kBT(it,1)=kB*T/4186.798188*NA;%kcal/mol
nit=50000;itl=0;nave/nit/5;
step=0;
ET_E(it,1)=ET(it,1);
while nit>0
    nit=nit-1;
    itl=itl+1;
    step(1,itl+1)=step(1,itl)+1;
H_E(:,itl+1)=H_E(:,itl);
V_E(:,itl+1)=V_E(:,itl);
%flip links of 3 types
ii=randi(100);
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
%H_E(ii,itl+1)=1-H_E(ii,itl);
% jh=nbH(i,1:2);jv=nbH(i,3:6);
eh=0;ev=0;
    for i=1:100
    h1=nbH(i,1);h2=nbH(i,2);v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    v1=nbV(i,1);v2=nbV(i,2);h1=nbV(i,3);h2=nbV(i,4);h3=nbV(i,5);h4=nbV(i,6);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))))*abs(V_E(i,itl+1));
    end
ET0_E(it,itl)=ev+eh;
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
%Vertical links are similar 3 types
    jj=randi(100);
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
eh=0;ev=0; 
    for i=1:100
    h1=nbH(i,1);h2=nbH(i,2);v1=nbH(i,3);v2=nbH(i,4);v3=nbH(i,5);v4=nbH(i,6);
    eh=eh+(E1+0.5*E2*(abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(V_E(v3,itl+1))+abs(V_E(v4,itl+1))))*abs(H_E(i,itl+1))+EE*H_E(i,itl+1);
    v1=nbV(i,1);v2=nbV(i,2);h1=nbV(i,3);h2=nbV(i,4);h3=nbV(i,5);h4=nbV(i,6);
    ev=ev+(E1+0.5*E2*(abs(V_E(v1,itl+1))+abs(V_E(v2,itl+1))+abs(H_E(h1,itl+1))+abs(H_E(h2,itl+1))+abs(H_E(h3,itl+1))+abs(H_E(h4,itl+1))))*abs(V_E(i,itl+1));
    end
ET1_E(it,itl)=eh+ev;
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
H_E_T(:,it)=H_E(:,itl+1);
V_E_T(:,it)=V_E(:,itl+1);
nH_E(it,1)=sum(abs(H_E(:,1)));
nV_E(it,1)=sum(abs(V_E(:,1)));
nTot_E(it,1)=nH_E(it,1)+nV_E(it,1);
nH_E(it,itl+1)=sum(abs(H_E(:,itl+1)));
nV_E(it,itl+1)=sum(abs(V_E(:,itl+1)));
nTot_E(it,itl+1)=nH_E(it,itl+1)+nV_E(it,itl+1);
end
meanHV_E(it,1)=mean(nH_E(it,itl-nave:itl+1));
meanHV_E(it,2)=mean(nV_E(it,itl-nave:itl+1));
meanHV_E(it,3)=mean(nTot_E(it,itl-nave:itl+1));
end
figure;
subplot(2,2,1)
plot(step,nTot_E); ylim([0,50]);
legend(['T=   10: nT=',num2str(meanHV_E(1,3),'%0.1f'),],...
['T= 100: nT=',num2str(meanHV_E(2,3),'%0.1f'),],...
['T= 200: nT=',num2str(meanHV_E(3,3),'%0.1f'),],...
['T= 300: nT=',num2str(meanHV_E(4,3),'%0.1f'),],...
['T= 400: nT=',num2str(meanHV_E(5,3),'%0.1f'),],...
['T= 500: nT=',num2str(meanHV_E(6,3),'%0.1f'),],...
['T= 600: nT=',num2str(meanHV_E(7,3),'%0.1f'),],...
['T= 700: nT=',num2str(meanHV_E(8,3),'%0.1f'),],...
['T= 800: nT=',num2str(meanHV_E(9,3),'%0.1f'),],...
['T=1200: nT=',num2str(meanHV_E(10,3),'%0.1f'),],...
['T=1600: nT=',num2str(meanHV_E(11,3),'%0.1f'),],...
['T=2000: nT=',num2str(meanHV_E(12,3),'%0.1f'),],...
['T=2500: nT=',num2str(meanHV_E(13,3),'%0.1f'),]);
title(['E-field, T = 10-2500K, MC Total Links']); 

subplot(2,2,2)
plot(step,nH_E); ylim([0,50]);
legend(['T=   10: nH=',num2str(meanHV_E(1,1),'%0.1f')],...
['T= 100: nT=',num2str(meanHV_E(2,3),'%0.1f'),],...
['T= 200: nT=',num2str(meanHV_E(3,3),'%0.1f'),],...
['T= 300: nT=',num2str(meanHV_E(4,3),'%0.1f'),],...
['T= 400: nT=',num2str(meanHV_E(5,3),'%0.1f'),],...
['T= 500: nT=',num2str(meanHV_E(6,3),'%0.1f'),],...
['T= 600: nT=',num2str(meanHV_E(7,3),'%0.1f'),],...
['T= 700: nT=',num2str(meanHV_E(8,3),'%0.1f'),],...
['T= 800: nT=',num2str(meanHV_E(9,3),'%0.1f'),],...
['T=1200: nT=',num2str(meanHV_E(10,3),'%0.1f'),],...
['T=1600: nT=',num2str(meanHV_E(11,3),'%0.1f'),],...
['T=2000: nT=',num2str(meanHV_E(12,3),'%0.1f'),],...
['T=2500: nT=',num2str(meanHV_E(13,3),'%0.1f'),]);
title(['E-field, T = 10-2500K, MC steps Horizontal links']); 

% map to new network figure at 300K
Tii=3; %300K
for i = 1:10
    for j=1:10
        Ht_E(i,j)=H_E_T((i-1)*10+j,Tii);
        ai=(11-i)*2;aj=j*2-1;
        if Ht_E(i,j)==1
           A(ai,aj,1)=192;A(ai,aj,2)=192;A(ai,aj,3)=255;%Horizontal+1 light
        else if Ht_E(i,j)==-1
                A(ai,aj,1)=76.8;A(ai,aj,2)=76.8;A(ai,aj,3)=255;%Horizontal-1 dark
            else
            A(ai,aj,1)=255;A(ai,aj,2)=255;A(ai,aj,3)=255;%blank
            end
        end
  
        Vt_E(j,i)=V_E_T((i-1)*10+j,Tii);
        bi=(11-j)*2-1;bj=i*2;
        if Vt_E(j,i)==1
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=192;% Vertical+1 light
        else if Vt_E(j,i)== -1
                A(bi,bj,1)=255;A(bi,bj,2)=230.4;A(bi,bj,3)=2.55;%Vertical-1 dark
            else 
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=255;%blank
            end
        end

    end
end

subplot(2,2,3)
image(uint8(A));
title(['Last step: E-field, T=300K, nH=',num2str(nH_E(Tii,itl+1)),', nV=',num2str(nV_E(Tii,itl+1))]); 

% map to new network figure at 2000K
Tii=7; %at 2000K
for i = 1:10
    for j=1:10
        Ht_E(i,j)=H_E_T((i-1)*10+j,Tii);
        ai=i*2;aj=j*2-1;
        if Ht_E(i,j)==1
           A(ai,aj,1)=192;A(ai,aj,2)=192;A(ai,aj,3)=255;%Horizontal+1 light
        else if Ht_E(i,j)==-1
                A(ai,aj,1)=76.8;A(ai,aj,2)=76.8;A(ai,aj,3)=255;%Horizontal-1 dark
            else
            A(ai,aj,1)=255;A(ai,aj,2)=255;A(ai,aj,3)=255;%blank
            end
        end
        
        Vt_E(j,i)=V_E_T((i-1)*10+j,Tii);
        bi=j*2-1;bj=i*2;
        if Vt_E(j,i)==1
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=192;% Vertical+1 light
        else if Vt_E(j,i)== -1
                A(bi,bj,1)=255;A(bi,bj,2)=230.4;A(bi,bj,3)=2.55;%Vertical-1 dark
            else 
            A(bi,bj,1)=255;A(bi,bj,2)=255;A(bi,bj,3)=255;%blank
            end
        end
    end
end

subplot(2,2,4)
image(uint8(A));
title(['Last step: E-field, T=2000K, nH=',num2str(nH_E(Tii,itl+1)),', nV=',num2str(nV_E(Tii,itl+1))]); 
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
EE=-10*4186.798188/NA;%-4.5*4186.798188/NA; %in kcal/mol

kBT(it,1)=kB*T/4186.798188*NA;
Ti(it,1)=T;
%D=80;
lm=size(mx,2);
for i = 1:lm
    Z(1,i)=(1+4*exp(-1/(kB*T)*(E1+4*mx(1,i)*E2))+4*exp(-1/(kB*T)*2*(E1+4*mx(1,i)*E2)));
    my(1,i)=(2*exp(-1/(kB*T)*(E1+4*mx(1,i)*E2))+4*exp(-1/(kB*T)*2*(E1+4*mx(1,i)*E2)))/Z(1,i); %network in Horizontal
end

for i = 2:lm
   if (my(1,i-1)-mx(1,i-1))*(my(1,i)-mx(1,i))<=0
       result(it,1)=my(1,i); %Horizontal
       result(it,2)=2*my(1,i);%Total
   end
end

%% with E-field

%-Efield*e*200*L*1E-10*1e6;e=1.602176634E-19;Efield=0.05;L=100; % E-field energy
 % in C
 % in V/m
 %linker length in A

for i = 1:lm
    mH(1,i)=mx(1,i);
    for j = 1:lm
    mV(1,j)=mx(1,j);
    ZE(j,i)=1+2*exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,j))*E2))+exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,j))*E2+EE))+ ... 
        exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,i))*E2-EE))+2*exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,i))*E2+EE))+ ... 
        2*exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,i))*E2-EE));
    mx1(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,j))*E2+EE))+2*exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,j))*E2+EE)))/ZE(j,i);
    mx2(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,j))*E2-EE))+2*exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,j))*E2-EE)))/ZE(j,i);
    my1(j,i)=(exp(-1/(kB*T)*(E1+(2*mH(1,i)+2*mV(1,j))*E2))+exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,j))*E2+EE))+exp(-1/(kB*T)*(2*E1+(4*mH(1,i)+4*mV(1,i))*E2-EE)))/ZE(j,i);
    
    mHs(j,i)=(mx1(j,i)+mx2(j,i));
    mVs(j,i)=(my1(j,i)+my1(j,i));
    del(j,i)=abs(mHs(j,i)-mH(1,i))+abs(mVs(j,i)-mV(1,j));
    end
end
del0(it,1)=min(min(del));
for j = 1:lm
    for i = 1:lm
       if del(j,i)==del0(it,1)
result(it,3)=mHs(j,i);
result(it,4)=mHs(j,i)+mVs(j,i);break;
       end
   end
end
end
figure;
mesh(mH(:),mV(:),del)
xlabel('H')
ylabel('V')
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plotMDT1 = plot([300	800	1200	1600	2000	2500], [4.0 8.9 36.5 43.9 58.2  68.9],'p-k','MarkerSize',18,'Parent',axes1);
plotMDT2 = plot([300	800	1200	1600	2000	2500], [3.0	3.7	17.9	16.3	24.0	26.4],'p-k','MarkerSize',2,'Parent',axes1);
plotMDE1 = plot ([10 100 200 300 400 500 600 700 800], [2.0	7.0	6.0	16.8 23.6 39.0 69.6 114.4 139.9],'p', 'Color',[1 0 1],'MarkerSize',18,'Parent',axes1);
plotMDE2 = plot ([10 100 200 300 400 500 600 700 800], [1.0 6.0 6.0 15.8 23.6 38.9 68.9  88.9  51.8],'p', 'Color',[1 0 1],'MarkerSize',2,'Parent',axes1);
plotMFT1 = plot(Ti,result(:,2)*100,'o','MarkerSize',12,'Parent',axes1);
plotMFT2 = plot(Ti,result(:,1)*100,'o','MarkerSize',12,'Parent',axes1);
plotMFE1 = plot(Ti,result(:,4)*100,'o','MarkerSize',12,'Parent',axes1);
plotMFE2 = plot(Ti,result(:,3)*100,'o','MarkerSize',12,'Parent',axes1);
plotMCT1 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHV(:,3),'s','MarkerSize',18,'Parent',axes1);
plotMCT2 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHV(:,1),'s','MarkerSize',2,'Parent',axes1);
plotMCE1 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHV_E(:,3),'^','MarkerSize',18,'Parent',axes1);
plotMCE2 = plot([10,100,200,300,400,500,600,700,800,1200,1600,2000,2500], meanHV_E(:,1),'^','MarkerSize',2,'Parent',axes1);
%plot1 = plot(Ti,result*100,'Parent',axes1);
set(plotMDT1(1),'DisplayName','DPD T-effect, Total','LineWidth',1.5,'LineStyle','-','Color',[0 0 0]);
set(plotMDT2(1),'DisplayName','DPD T-effect, x-link','LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
set(plotMDE1(1),'DisplayName','DPD E-field 0.0005 V/Å, Total','LineWidth',1.5,'LineStyle','-','Color',[1 0 1]);
set(plotMDE2(1),'DisplayName','DPD E-field 0.0005 V/Å, x-link','LineWidth',1.5,'LineStyle','--','Color',[1 0 1]);
set(plotMFT1(1),'DisplayName','MFT T-effect, Total','LineWidth',3,'LineStyle','-','Color',[0 0 1]);
set(plotMFT2(1),'DisplayName','MFT T-effect, x-link','LineWidth',2,'LineStyle','--','Color',[0 0 1]);
set(plotMFE1(1),'DisplayName','MFT E-field 0.0005 V/Å, Total','LineWidth',1.5,'LineStyle','-','Color',[1 0 0]);
set(plotMFE2(1),'DisplayName','MFT E-field 0.0005 V/Å, x-link','LineWidth',2,'LineStyle','--','Color',[1 0 0]);
set(plotMCT1(1),'DisplayName','MC T-effect, Total','LineWidth',1.5,'LineStyle','-','Color',[0 1 1]);
set(plotMCT2(1),'DisplayName','MC T-effect, x-link','LineWidth',1.5,'LineStyle','--','Color',[0 1 1]);
%set(plot4(2),'DisplayName','MC, No E-field, Vertical','LineWidth',1.5,'LineStyle','none','Marker','none','Color',[0 1 1]);
set(plotMCE1(1),'DisplayName','MC E-field, Total','LineWidth',1.5,'LineStyle','-','Color',[1 0 0.7]);
set(plotMCE2(1),'DisplayName','MC E-field, x-link','LineWidth',1.5,'LineStyle','--','Color',[1 0 0.7]);
%set(plot5(2),'DisplayName','MC, E-field, Vertical','LineWidth',1.5,'LineStyle','none','Marker','none','Color',[1 0 0.7]);
box(axes1,'on');
hold(axes1,'off');
legend1 = legend(axes1,'show');
ylabel('Likers (#)');
%ylim([0,120]);
xlabel('Temperature (K)');
set(axes1,'FontSize',20);
title(['MFT: Horizontal and Total links']); 

result*100;
