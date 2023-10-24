clear all;
clc
%cd ./T/T1600K
logfile     =  'T2500KSQR2au.txt';
logfile2    =  'T2500KSQR2pah.txt';
logfile3 = 'log.T2500KSQR2';
Edata=importdata(logfile3,' ',130);
Ea=Edata.data; PE=Ea(:,4);tlog=Ea(:,1)/1000000;
figure;plot(tlog,PE,tlog,smooth(PE,100));title('Potential Energy vs Time');

%%
st=41;np=100;npp=100;rc=27

for ii = 1:st
    clear Ab Bb Cb R netS net Link;
    t(ii,1)=(ii-1)*0.5;
    
n=(ii-1)*40009+9; 

data=importdata(logfile,' ',n);
Aa=data.data;
Ba=sortrows(Aa,2);

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
        if net(i,j)>0
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
ncon(ii,1)=count;

clear netS1 netS2 netS3 netS4 netS4v netS5
for i=1:np
    for j=1:npp
        netS1(i,j)=0; netS2(i,j)=0;netS3(i,j)=0;netS4(i,j)=0;netS4v(i,j)=0; netS5(i,j)=0;
        if netS(i,j)==1
            netS5(i,j)=1; % set all connection special
            if i==j
                netS1(i,j)=1; netS5(i,j)=0;end% on the diagonal
            if i-j==1 && netS(j,j)==1
                netS2(i,j)=1; netS5(i,j)=0;% forward horizontal
                if  mod(i,10)==0 && j==i+1
                    netS2(i,j)=0;netS5(i,j)=1; end; end% special
            if j-i==1 && netS(j,j)==1
                netS2(i,j)=1; netS5(i,j)=0; % backward horizontal
                if mod(j,10)==0 && i==j+1
                   netS2(i,j)=0; netS5(i,j)=1; end; end;
            if mod(i,10)==0 && j == i-9 && netS(j,j)==1
                netS2(i,j)=1;netS5(i,j)=0;end
            if mod(j,10)==0 && i == j-9 && netS(j,j)==1
                netS2(i,j)=1;netS5(i,j)=0;end

            if abs(i-j)==10 && netS(j,j)==1
                netS3(i,j)=1; netS5(i,j)=0;end% horizontal 
            if abs(i-j)==90 && netS(j,j)==1
                netS3(i,j)=1;netS5(i,j)=0;end

            if abs(i-j)==2 && netS(j,j)==1
                netS4(i,j)=1; netS5(i,j)=0;end % 2nd neighbor 
            if abs(i-j)==2 && mod(i,10)==0 && j==i+2 && netS(j,j)==1
                netS5(i,j)=1;end
            if abs(i-j)==2 && mod(j,10)==0 && i==j+2 && netS(j,j)==1
                netS5(i,j)=1; end
            if  mod(i,10)==1 && j-i == 8 && netS(j,j)==1
                netS4(i,j)=1;netS5(i,j)=0;end
            if  mod(i,10)==9 && i == j-8 && netS(j,j)==1
                netS4(i,j)=1;netS5(i,j)=0;end 
            if  mod(i,10)==2 && j-i == 8 && netS(j,j)==1
                netS4(i,j)=1;netS5(i,j)=0;end
            if  mod(i,10)==0 && i == j-8 && netS(j,j)==1
                netS4(i,j)=1;netS5(i,j)=0;end 
            if abs(i-j)==20 || abs(i-j)==80 && netS(j,j)==1
                netS4v(i,j)=1;netS5(i,j)=0;end
            end
        end
end
result(ii,1)= sum(sum(netS,2)); % totol
result(ii,2)= sum(sum(netS1,2));% self net
result(ii,3)= sum(sum(netS2,2));% 1st Horizontal
result(ii,4)= sum(sum(netS3,2));% 1st Vertical
result(ii,5)= sum(sum(netS4,2));% 2nd Horizontal
result(ii,6)= sum(sum(netS4v,2));% 2nd Vertical
result(ii,7)= sum(sum(netS5,2)); % special

links(1,:)=[0,0,0,0]; ncon(ii,2:5)= links(1,:);
%% 4 types of links: +x,-x,+y,-y,others
if ncon(ii,1)~=0
kl=0;
r_ij=0;
for ir=1:npp
    if netS(ir,ir)==1
        for jr=1:np
            if netS(jr,ir)==1 && ir~=jr
                kl=kl+1;
                r_ij(kl,1)=ir;
                r_ij(kl,2)=jr;
                r_ij(kl,3:5)=com(jr,:)-com(ir,:);
                    if r_ij(kl,3)>880
                        r_ij(kl,3)= 1000-r_ij(kl,3);
                    end
                    if r_ij(kl,3)< -880
                        r_ij(kl,3)= 1000+r_ij(kl,3);
                    end
                    if r_ij(kl,4)>880
                        r_ij(kl,4)= 1000-r_ij(kl,4);
                    end
                    if r_ij(kl,4)< -880
                        r_ij(kl,4)= 1000+r_ij(kl,4);
                    end
            end
        end        
    end
end
%6 types of links, 1:+x,2:-x,3:+y,4:-y
lr=size(r_ij, 1);link1=0;link2=0;link3=0;link4=0;
    if r_ij~=0
        for ir = 1:lr 
    if r_ij(ir,3)/100>0.8 && r_ij(ir,3)/100 < 1.2 && abs(r_ij(ir,4)/100) < 0.2 
       r_ij(ir,6)=1;
       link1=link1+1; 
    end
    if r_ij(ir,3)/100<-0.8 && r_ij(ir,3)/100>-1.2 && abs(r_ij(ir,4)/100) < 0.2 
       r_ij(ir,6)=2;
       link2=link2+1; 
    end
    if abs(r_ij(ir,3)/100) < 0.2 &&r_ij(ir,4)/100>0.8 && r_ij(ir,4)/100 < 1.2 
       r_ij(ir,6)=3;link3=link3+1;end
    if abs(r_ij(ir,3)/100) < 0.2 &&r_ij(ir,4)/100<-0.8 && r_ij(ir,4)/100 > -1.2 
       r_ij(ir,6)=4;link4=link4+1;end
    end
end
links(1,:)=[link1,link2,link3,link4];
ncon(ii,2:5)= links(1,:); %1st nearest neighbors 4 types of links, 1:+x,2:-x,3:+y,4:-y
ncon(ii,6)= ncon(ii,1)-sum(links(1,:)); %unusual links all others
end
end
figure; plot(t,ncon);
[X,Y] = meshgrid(1:100);
figure;surf(X,Y,net);
%%
figure; axes1 = axes;
hold(axes1,'on');
contour(netS1,'DisplayName','Self-connection','LineColor',[0 0 0]);
contour(netS2,'DisplayName','H-connection','LineColor',[1 0 0]);
contour(netS3,'DisplayName','V-connection','LineColor',[0 1 0]);
contour(netS4,'DisplayName','H-2nd connection','LineColor',[0 0 1]);
contour(netS4v,'DisplayName','V-2nd connection','LineColor',[0 1 1]);
contour(netS5,'DisplayName','Special-connection','LineColor',[1 0 1]);
xlabel('PAH ID');
ylabel('AuNP ID');
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'BoxStyle','full','FontSize',20,'Layer','top','LineWidth',2);
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.22 0.64 0.33 0.226]);
%%
npah= sum(netS,1);
nenp= sum(netS,2);
figure; subplot(1,2,1); hist(npah);subplot(1,2,2); hist(nenp);