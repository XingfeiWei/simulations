clear all;
clc
%cd ../V3
logfile =      'T300KHCPR2au.txt';
logfile2 =     'T300KHCPR2pah.txt';
logfile3 = 'log.T300KHCPR2';
Edata=importdata(logfile3,' ',127); % the starting line in log
Ea=Edata.data; PE=Ea(:,4);tlog=Ea(:,1)/1000000;
figure;plot(tlog,PE,tlog,smooth(PE,100));title('Potential Energy vs Time');
xlbox=500;ylbox=519.61524227066320;zlbox=489.89794855663558;

%%
st=41;  %number of frames
np=180; %number of AuNP
npp=180; %number of PAH
rc=27; %cut off 27 A radius 

for ii = 1:st
    clear Ab Bb Cb R netS net Link;
    t(ii,1)=(ii-1)*0.5;
    
n=(ii-1)*72009+9; % the line for each step

data=importdata(logfile,' ',n);
Aa=data.data;
Ba=sortrows(Aa,2);
for i=1:size(Ba,1)
    if Ba(i,6)~=0
        Ba(i,3) = Ba(i,3)+Ba(i,6)*xlbox;
    end
    if Ba(i,7)~=0
        Ba(i,4) = Ba(i,4)+Ba(i,7)*ylbox;
    end
    if Ba(i,8)~=0
        Ba(i,5) = Ba(i,5)+Ba(i,8)*zlbox;
    end
end

for i=1:np
   mi=size(Ba,1)/np*(i-1)+1; ni=size(Ba,1)/np*i;
com(i,1)=mean(Ba(mi:ni,3));
com(i,2)=mean(Ba(mi:ni,4));
com(i,3)=mean(Ba(mi:ni,5));
end

data2=importdata(logfile2,' ',n);
Ab=data2.data;
Bb=sortrows(Ab,2);
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

count=0;
for i=1:np
    for j=1:npp
        netS(i,j)=0;
        net(i,j)=sum(Link(i,j,:));
        if net(i,j)>0 %&& net(i,i)~=0
            netS(i,j)=1;
            if i~=j
                count=count+1;
            end
        end
    end
end

if t(ii,1) == 0 || t(ii,1) == 5|| t(ii,1) == 10
figure;contour(netS);
end
%add
%% 12 types of links: the link vector r_ij 
%clear r_ij links ncon
ncon(ii,1)=count; 
links(1,:)=[0,0,0,0,0,0,0,0,0,0,0,0]; ncon(ii,2:13)= links(1,:);
if ncon(ii,1)~=0
kl=0;
r_ij=0;
for ir=1:npp
    if netS(ir,ir)==1
        for jr=1:np
            if netS(jr,ir)==1 && ir~=jr % only ENP with initial polymer
                kl=kl+1;
                r_ij(kl,1)=ir;
                r_ij(kl,2)=jr;
                r_ij(kl,3:5)=com(jr,:)-com(ir,:);
                    if r_ij(kl,3)>350
                        r_ij(kl,3)= xlbox-r_ij(kl,3);
                    end
                    if r_ij(kl,3)< -350
                        r_ij(kl,3)= xlbox+r_ij(kl,3);
                    end
                    if r_ij(kl,4)>350
                        r_ij(kl,4)= ylbox-r_ij(kl,4);
                    end
                    if r_ij(kl,4)< -350
                        r_ij(kl,4)= ylbox+r_ij(kl,4);
                    end
                    if r_ij(kl,5)>350
                        r_ij(kl,5)= zlbox-r_ij(kl,5);
                    end
                    if r_ij(kl,5)< -350
                        r_ij(kl,5)= zlbox+r_ij(kl,5);
                    end
            end
        end
    end
end
%% 12 types of links, 1:+x,2:-x,3:+y,4:-y,5:+z,6:-z.
%6 types of in-plane HEX links, 1:+x0y,2:-x0y,3:+x+y,4:-x-y,5:+x-y,6:-x+y.
lr=size(r_ij, 1);link1=0;link2=0;link3=0;link4=0;link5=0;link6=0;link7=0;link8=0;link9=0;link10=0;link11=0;link12=0;
for i = 1:lr
    if abs(r_ij(i,5)/100) < 0.1
    if r_ij(i,3)/100>0.9 && r_ij(i,3)/100 < 1.1 && abs(r_ij(i,4)/100) < 0.1
       r_ij(i,6)=1; % store link types in r_ij(i,6)
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
 %6 types of out-of-plane links, 7:+x+y+z,8:-x+y+z,9:0x-y+z,10:+x+y-z,11:-x+y-z,12:0x-y-z.
  
 if r_ij(i,5)/100 > 0.1  %%  top 3 type links 
    if r_ij(i,3)/100>0.3 && r_ij(i,3)/100 < 0.7 && r_ij(i,4)/100 > 0.1 && r_ij(i,4)/100 < 0.5
       r_ij(i,6)=7;link7=link7+1;end
    if r_ij(i,3)/100 < -0.3 && r_ij(i,3)/100 > -0.7 && r_ij(i,4)/100 > 0.1 && r_ij(i,4)/100 < 0.5
       r_ij(i,6)=8;link8=link8+1;end
    if abs(r_ij(i,3))/100<0.1 && r_ij(i,4)/100 < -0.4 && r_ij(i,4)/100 > -0.8
       r_ij(i,6)=9;link9=link9+1;end   
    if r_ij(i,3)/100<-0.3 && r_ij(i,3)/100>-0.7 && r_ij(i,4)/100 <-0.1 && r_ij(i,4)/100 >-0.5
       r_ij(i,6)=10;link10=link10+1;end
    if r_ij(i,3)/100 > 0.3 && r_ij(i,3)/100 < 0.7 && r_ij(i,4)/100 <-0.1 && r_ij(i,4)/100 >-0.5
       r_ij(i,6)=11;link11=link11+1;end
    if abs(r_ij(i,3))/100<0.1 && r_ij(i,4)/100 > 0.4 && r_ij(i,4)/100 < 0.8
       r_ij(i,6)=12;link12=link12+1;end 
 end
 
 if r_ij(i,5)/100 < -0.1 %% bottom 3 type links 
    if r_ij(i,3)/100<-0.3 && r_ij(i,3)/100 >-0.7 && r_ij(i,4)/100 <-0.1 && r_ij(i,4)/100 >-0.5
       r_ij(i,6)=7;link7=link7+1;end
    if r_ij(i,3)/100 > 0.3 && r_ij(i,3)/100 < 0.7 && r_ij(i,4)/100 <-0.1 && r_ij(i,4)/100 >-0.5
       r_ij(i,6)=8;link8=link8+1;end
    if abs(r_ij(i,3))/100<0.1 && r_ij(i,4)/100 > 0.4 && r_ij(i,4)/100 < 0.8
       r_ij(i,6)=9;link9=link9+1;end 
    if r_ij(i,3)/100>0.3 && r_ij(i,3)/100 < 0.7 && r_ij(i,4)/100 > 0.1 && r_ij(i,4)/100 < 0.5
       r_ij(i,6)=10;link10=link10+1;end
    if r_ij(i,3)/100 < -0.3 && r_ij(i,3)/100 > -0.7 && r_ij(i,4)/100 > 0.1 && r_ij(i,4)/100 < 0.5
       r_ij(i,6)=11;link11=link11+1;end
    if abs(r_ij(i,3))/100<0.1 && r_ij(i,4)/100 < -0.4 && r_ij(i,4)/100 > -0.8
       r_ij(i,6)=12;link12=link12+1;end                 
 end
end
links(1,:)=[link1,link2,link3,link4,link5,link6,link7,link8,link9,link10,link11,link12];
ncon(ii,2:13)= links(1,:); %1st nearest neighbors 6 types of links, 1:+x0y,2:-x0y,3:+x+y,4:-x-y,5:+x-y,6:-x+y.
ncon(ii,14)= ncon(ii,1)-sum(links(1,:)); %unusual links all others
end
end
figure; plot(t,ncon);
[X,Y] = meshgrid(1:180);
figure;surf(X,Y,net);
%%
npah= sum(netS,1);
nenp= sum(netS,2);
figure; subplot(1,2,1); hist(npah);subplot(1,2,2); hist(nenp);