clear all;
clc
%cd ../V3
logfile     =  'T2500KHEXR2au.txt';
logfile2     = 'T2500KHEXR2pah.txt';
logfile3 = 'log.T2500KHEXR2';
Edata=importdata(logfile3,' ',130);
Ea=Edata.data; PE=Ea(:,4);tlog=Ea(:,1)/1000000;
figure;plot(tlog,PE,tlog,smooth(PE,100));title('Potential Energy vs Time');

%
st=41;  %number of frames
np=100; %number of AuNP
npp=100; %number of PAH
rc=27; %cut off 27 A radius 

for ii = 1:st
    clear Ab Bb Cb R netS net Link;
    t(ii,1)=(ii-1)*0.5;
    
n=(ii-1)*40009+9; 

data=importdata(logfile,' ',n);
Aa=data.data;
Ba=sortrows(Aa,1);
for i = 1:size(Ba,1)
    Ba(i,3)=Ba(i,6)*1000+Ba(i,3);
    Ba(i,4)=Ba(i,7)*1000+Ba(i,4);
    Ba(i,5)=Ba(i,8)*1000+Ba(i,5);
end

for i=1:np
   mi=size(Ba,1)/np*(i-1)+1; ni=size(Ba,1)/np*i;
com(i,1)=mean(Ba(mi:ni,3));
com(i,2)=mean(Ba(mi:ni,4));
com(i,3)=mean(Ba(mi:ni,5));
end

data2=importdata(logfile2,' ',n);
Ab=data2.data;
Bb=sortrows(Ab,1);
for i = 1:size(Bb,1)
    Bb(i,3)=Bb(i,6)*1000+Bb(i,3);
    Bb(i,4)=Bb(i,7)*1000+Bb(i,4);
    Bb(i,5)=Bb(i,8)*1000+Bb(i,5);
end
for i=1:npp
    mi=size(Bb,1)/npp*(i-1)+1; ni=size(Bb,1)/npp*i;
    Cb(:,:,i)=Bb(mi:ni,3:5);
end
for i=1:np
    for j=1:npp
        for k=1:size(Bb,1)/npp
            Link(i,j,k)=0;
            R(i,j,k)=((Cb(k,1,j)-com(i,1))^2+(Cb(k,2,j)-com(i,2))^2+(Cb(k,3,j)-com(i,3))^2)^0.5;
            if R(i,j,k)<rc
            Link(i,j,k)=1;
            end
        end
    end
end

%% PBC 5 AuNPs on the boundary
comPBC(1,:)=com(2,:); 
comPBC(2,:)=com(22,:); 
comPBC(3,:)=com(42,:); 
comPBC(4,:)=com(62,:); 
comPBC(5,:)=com(82,:); comPBC(:,1)=comPBC(:,1)+1000;
for i=[2 22 42 62 82]
    for j=1:npp
        for k=1:size(Bb,1)/npp
            R(i,j,k)=((Cb(k,1,j)-comPBC((i-2)/20+1,1))^2+(Cb(k,2,j)-comPBC((i-2)/20+1,2))^2+(Cb(k,3,j)-comPBC((i-2)/20+1,3))^2)^0.5;
            if R(i,j,k)<rc
            Link(i,j,k)=1;
            end
        end
    end
end
%%
count=0;
for i=1:np
    for j=1:npp
        netS(i,j)=0;
        net(i,j)=sum(Link(i,j,:));
        if net(i,j)>0
            netS(i,j)=1;
            if i~=j
                count=count+1;
            end
        end
    end
end

%if t(ii,1) == 0 || t(ii,1) == 5|| t(ii,1) == 10
%figure;contour(netS);
%end
ncon(ii,1)=count;
links(1,:)=[0,0,0,0,0,0]; ncon(ii,2:7)= links(1,:);

%% 6 types of links: the link vector r_ij 
if ncon(ii,1)~=0
kl=0;
r_ij=0;
for i=1:npp
    if netS(i,i)==1
        for j=1:np
            if netS(j,i)==1 && i~=j % only ENP with initial polymer
                kl=kl+1;
                r_ij(kl,1)=i;
                r_ij(kl,2)=j;
                r_ij(kl,3:5)=com(j,:)-com(i,:);
                    if r_ij(kl,3)>890
                        r_ij(kl,3)= 1000-r_ij(kl,3);
                    end
                    if r_ij(kl,3)< -890
                        r_ij(kl,3)= 1000+r_ij(kl,3);
                    end
                    if r_ij(kl,4)>890
                        r_ij(kl,4)= 1000-r_ij(kl,4);
                    end
                    if r_ij(kl,4)< -890
                        r_ij(kl,4)= 1000+r_ij(kl,4);
                    end
            end
        end
    end
end
%6 types of links, 1:+x0y,2:-x0y,3:+x+y,4:-x-y,5:+x-y,6:-x+y.
lr=size(r_ij, 1);link1=0;link2=0;link3=0;link4=0;link5=0;link6=0;
for i = 1:lr
    if r_ij(i,3)/100>0.9 && r_ij(i,3)/100 < 1.1 && abs(r_ij(i,4)/100) < 0.1
       r_ij(i,6)=1;
       link1=link1+1; 
    end
    if r_ij(i,3)/100<-0.9 && r_ij(i,3)/100>-1.1 && abs(r_ij(i,4)/100) < 0.1
       r_ij(i,6)=2;
       link2=link2+1; 
    end
    if r_ij(i,3)/100>0.4 && r_ij(i,3)/100 < 0.6 && r_ij(i,4)/100 > 0.5
       r_ij(i,6)=3;link3=link3+1;end
    if r_ij(i,3)/100<-0.4 && r_ij(i,3)/100 >-0.6 && r_ij(i,4)/100 <-0.5
       r_ij(i,6)=4;link4=link4+1;end
    if r_ij(i,3)/100>0.4 && r_ij(i,3)/100 < 0.6 && r_ij(i,4)/100 < -0.5
       r_ij(i,6)=5;link5=link5+1;end
    if r_ij(i,3)/100<-0.4 && r_ij(i,3)/100 >-0.6 && r_ij(i,4)/100 > 0.5
       r_ij(i,6)=6;link6=link6+1;end
end
links(1,:)=[link1,link2,link3,link4,link5,link6];
ncon(ii,2:7)= links(1,:); %1st nearest neighbors 6 types of links, 1:+x0y,2:-x0y,3:+x+y,4:-x-y,5:+x-y,6:-x+y.
ncon(ii,8)= ncon(ii,1)-sum(links(1,:)); %unusual links all others
%figure; hist(r_ij(:,6),6);
end
end
figure; plot(t,ncon);
[X,Y] = meshgrid(1:100);
figure;surf(X,Y,net);
%%
npah= sum(netS,1);
nenp= sum(netS,2);
figure; subplot(1,2,1); hist(npah);subplot(1,2,2); hist(nenp);