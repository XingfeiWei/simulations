%% MFT approximation 3D HCP
clear all
clc
delx=0.1;
mx=0:delx:1; % average network per box
kB = 1.380649E-23; % in J/K

NA=6.02214076E+23;  % #/mol
E1=5*4186.798188/NA; % connection energy kcal/mol to J
E2=2*4186.798188/NA; % neighbor correlation energy
EE=-6*4186.798188/NA;%-7*4186.798188/NA; %in kcal/mol

% no E-field
it=0;
%for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
for T =[10] % in K 
it=it+1;
kBT(it,1)=kB*T/4186.798188*NA;
Ti(it,1)=T;
lm=size(mx,2);
for i = 1:lm
    Z(1,i)=(1+12*exp(-1/(kB*T)*(E1+22*mx(1,i)*E2))+ ...
        60*exp(-1/(kB*T)*2*(E1+22*mx(1,i)*E2))+ ...
        160*exp(-1/(kB*T)*3*(E1+22*mx(1,i)*E2))+ ...
        240*exp(-1/(kB*T)*4*(E1+22*mx(1,i)*E2))+ ...
        192*exp(-1/(kB*T)*5*(E1+22*mx(1,i)*E2))+ ...
        64*exp(-1/(kB*T)*6*(E1+22*mx(1,i)*E2)));
    my(1,i)=(2*exp(-1/(kB*T)*(E1+22*mx(1,i)*E2))+ ...
        20*exp(-1/(kB*T)*2*(E1+22*mx(1,i)*E2))+ ...
        80*exp(-1/(kB*T)*3*(E1+22*mx(1,i)*E2))+ ...
        160*exp(-1/(kB*T)*4*(E1+22*mx(1,i)*E2))+ ...
        160*exp(-1/(kB*T)*5*(E1+22*mx(1,i)*E2))+ ...
        64*exp(-1/(kB*T)*6*(E1+22*mx(1,i)*E2)))/Z(1,i); %network in Horizontal
end

for i = 2:lm
   if (my(1,i-1)-mx(1,i-1))*(my(1,i)-mx(1,i))<=0
       result(it,1)=(mx(1,i)+mx(1,i-1)+my(1,i)+my(1,i-1))/4; %Horizontal
       result(it,2)=6*result(it,1);%Total
   end
end
end

%% with E-field
%EE=-Efield*e*200*L*1E-10*1e6;e=1.602176634E-19;Efield=0.05;L=100; % E-field energy
 % in C
 % in V/m
 %linker length in A
it = 0;
%parfor (n=1:nmax, 12)
%for T =[10,100,200,300,400,500,600,700,800,1200,1600,2000,2500] % in K 
for T =[10] % in K 
it=it+1;
kBT(it,1)=kB*T/4186.798188*NA;
Ti(it,1)=T;
lm=size(mx,2);
nnx = 2;
nny = 4;

for i = 1:lm
    mH(1,i)=mx(1,i);
    for j = 1:lm
    mV(1,j)=mx(1,j);
    for ja = 1:lm
        mA(1,ja)=mx(1,ja);
        for jb = 1:lm
            mB(1,jb)=mx(1,jb);
    Elink(1,i,j,ja,jb) = E1+(nnx*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2+EE;
    Elink(2,i,j,ja,jb) = E1+(nnx*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2-EE;
    Elink(3,i,j,ja,jb) = E1+(nny*mH(1,i)+nnx*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2+EE*cos(pi()/3);
    Elink(4,i,j,ja,jb) = E1+(nny*mH(1,i)+nnx*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2-EE*cos(-pi()/3);
    Elink(5,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nnx*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2+EE*cos(pi()/3);
    Elink(6,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nnx*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2-EE*cos(-pi()/3);
    Elink(7,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nnx*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2;
    Elink(8,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nnx*mA(1,ja)+nny*mB(1,jb)+nny*mB(1,jb))*E2;
    Elink(9,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nnx*mB(1,jb)+nny*mB(1,jb))*E2+2/3*EE*cos(pi()/6);
    Elink(10,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nnx*mB(1,jb)+nny*mB(1,jb))*E2-2/3*EE*cos(-pi()/6);
    Elink(11,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nnx*mB(1,jb))*E2+2/3*EE*cos(pi()/6);
    Elink(12,i,j,ja,jb) = E1+(nny*mH(1,i)+nny*mV(1,j)+nny*mV(1,j)+nny*mA(1,ja)+nny*mB(1,jb)+nnx*mB(1,jb))*E2-2/3*EE*cos(-pi()/6);

ZE1(i,j,ja,jb) = 0;  % for the total partition function 
for kk1 = 1:6 
    for kk1b = 1:2
    kk = (kk1-1)*2+kk1b;
    ZE1(i,j,ja,jb) = ZE1(i,j,ja,jb) + (exp(-1/(kB*T)*Elink(kk,i,j,ja,jb)));
    end
end

ZE2(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        ZE2(i,j,ja,jb) = ZE2(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
            end
    end
    end
end

ZE3(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
        ZE3(i,j,ja,jb) = ZE3(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                    end
                end
            end
    end
    end
end

ZE4(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
        ZE4(i,j,ja,jb) = ZE4(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
        ZE5(i,j,ja,jb) = ZE5(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
        ZE6(i,j,ja,jb) = ZE6(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
ZE(i,j,ja,jb)=1+ZE1(i,j,ja,jb)+ZE2(i,j,ja,jb)+ZE3(i,j,ja,jb)+ZE4(i,j,ja,jb)+ZE5(i,j,ja,jb)+ZE6(i,j,ja,jb);


ZE1xp(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(1,i,j,ja,jb)));  % for x+ direction
ZE2xp(i,j,ja,jb) = 0;
for kk1 = 1:1
    for kk1b = 1:1
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        ZE2xp(i,j,ja,jb) = ZE2xp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
            end
    end
    end
end

ZE3xp(i,j,ja,jb) = 0;
for kk1 = 1:1
    for kk1b = 1:1
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
        ZE3xp(i,j,ja,jb) = ZE3xp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                    end
                end
            end
    end
    end
end

ZE4xp(i,j,ja,jb) = 0;
for kk1 = 1:1
for kk1b = 1:1
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
        ZE4xp(i,j,ja,jb) = ZE4xp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5xp(i,j,ja,jb) = 0;
for kk1 = 1:1
for kk1b = 1:1
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
        ZE5xp(i,j,ja,jb) = ZE5xp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6xp(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:1
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
        ZE6xp(i,j,ja,jb) = ZE6xp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
mx1(i,j,ja,jb)=(ZE1xp(i,j,ja,jb)+ZE2xp(i,j,ja,jb)+ZE3xp(i,j,ja,jb)+ZE4xp(i,j,ja,jb)+ZE5xp(i,j,ja,jb)+ZE6xp(i,j,ja,jb))/ZE(i,j,ja,jb);


ZE1xm(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(2,i,j,ja,jb))); % for x- direction

ZE2xm(i,j,ja,jb) = 0;
for kk1 = 1:1
    for kk1b = 2:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        ZE2xm(i,j,ja,jb) = ZE2xm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
            end
    end
    end
end

ZE3xm(i,j,ja,jb) = 0;
for kk1 = 1:1
    for kk1b = 2:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
        ZE3xm(i,j,ja,jb) = ZE3xm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                    end
                end
            end
    end
    end
end

ZE4xm(i,j,ja,jb) = 0;
for kk1 = 1:1
for kk1b = 2:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
        ZE4xm(i,j,ja,jb) = ZE4xm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5xm(i,j,ja,jb) = 0;
for kk1 = 1:1
for kk1b = 2:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
        ZE5xm(i,j,ja,jb) = ZE5xm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6xm(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 2:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
        ZE6xm(i,j,ja,jb) = ZE6xm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
mx2(i,j,ja,jb)=(ZE1xm(i,j,ja,jb)+ZE2xm(i,j,ja,jb)+ZE3xm(i,j,ja,jb)+ZE4xm(i,j,ja,jb)+ZE5xm(i,j,ja,jb)+ZE6xm(i,j,ja,jb))/ZE(i,j,ja,jb);

ZE1yp(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(3,i,j,ja,jb))); % for U+ or V+ direction

ZE2yp(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        if kk_1 == 3 || kk_2 == 3
        ZE2yp(i,j,ja,jb) = ZE2yp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
        end
            end
    end
    end
end

ZE3yp(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        if kk_1 == 3 || kk_2 == 3 || kk_3 == 3
        ZE3yp(i,j,ja,jb) = ZE3yp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                        end
                    end
                end
            end
    end
    end
end

ZE4yp(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                        if kk_1 == 3 || kk_2 == 3 || kk_3 == 3 || kk_4 == 3
        ZE4yp(i,j,ja,jb) = ZE4yp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5yp(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                            if kk_1 == 3 || kk_2 == 3 || kk_3 == 3 || kk_4 == 3 || kk_5 == 3
        ZE5yp(i,j,ja,jb) = ZE5yp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6yp(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
                                 if kk_1 == 3 || kk_2 == 3 || kk_3 == 3 || kk_4 == 3 || kk_5 == 3 || kk_6 == 3
        ZE6yp(i,j,ja,jb) = ZE6yp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
my1(i,j,ja,jb)=(ZE1yp(i,j,ja,jb)+ZE2yp(i,j,ja,jb)+ZE3yp(i,j,ja,jb)+ZE4yp(i,j,ja,jb)+ZE5yp(i,j,ja,jb)+ZE6yp(i,j,ja,jb))/ZE(i,j,ja,jb);


ZE1ym(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(4,i,j,ja,jb))); % for U- or V- direction

ZE2ym(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        if kk_1 == 4 || kk_2 == 4
        ZE2ym(i,j,ja,jb) = ZE2ym(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
        end
            end
    end
    end
end

ZE3ym(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        if kk_1 == 4 || kk_2 == 4 || kk_3 == 4
        ZE3ym(i,j,ja,jb) = ZE3ym(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                        end
                    end
                end
            end
    end
    end
end

ZE4ym(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                        if kk_1 == 4 || kk_2 == 4 || kk_3 == 4 || kk_4 == 4
        ZE4ym(i,j,ja,jb) = ZE4ym(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5ym(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                            if kk_1 == 4 || kk_2 == 4 || kk_3 == 4 || kk_4 == 4 || kk_5 == 4
        ZE5ym(i,j,ja,jb) = ZE5ym(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6ym(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
                                 if kk_1 == 4 || kk_2 == 4 || kk_3 == 4 || kk_4 == 4 || kk_5 == 4 || kk_6 == 4
        ZE6ym(i,j,ja,jb) = ZE6ym(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
my2(i,j,ja,jb)=(ZE1ym(i,j,ja,jb)+ZE2ym(i,j,ja,jb)+ZE3ym(i,j,ja,jb)+ZE4ym(i,j,ja,jb)+ZE5ym(i,j,ja,jb)+ZE6ym(i,j,ja,jb))/ZE(i,j,ja,jb);


ZE1A(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(7,i,j,ja,jb))); % for alpha direction
ZE2A(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        if kk_1 == 7 || kk_2 == 7
        ZE2A(i,j,ja,jb) = ZE2A(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
        end
            end
    end
    end
end

ZE3A(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        if kk_1 == 7 || kk_2 == 7 || kk_3 == 7
        ZE3A(i,j,ja,jb) = ZE3A(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                        end
                    end
                end
            end
    end
    end
end

ZE4A(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                        if kk_1 == 7 || kk_2 == 7 || kk_3 == 7 || kk_4 == 7
        ZE4A(i,j,ja,jb) = ZE4A(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5A(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                            if kk_1 == 7 || kk_2 == 7 || kk_3 == 7 || kk_4 == 7 || kk_5 == 7
        ZE5A(i,j,ja,jb) = ZE5A(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6A(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
                                 if kk_1 == 7 || kk_2 == 7 || kk_3 == 7 || kk_4 == 7 || kk_5 == 7 || kk_6 == 7
        ZE6A(i,j,ja,jb) = ZE6A(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
mA1(i,j,ja,jb)=(ZE1A(i,j,ja,jb)+ZE2A(i,j,ja,jb)+ZE3A(i,j,ja,jb)+ZE4A(i,j,ja,jb)+ZE5A(i,j,ja,jb)+ZE6A(i,j,ja,jb))/ZE(i,j,ja,jb);


ZE1Bp(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(9,i,j,ja,jb))); % for beta or gama direction
ZE2Bp(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        if kk_1 == 9 || kk_2 == 9
        ZE2Bp(i,j,ja,jb) = ZE2Bp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
        end
            end
    end
    end
end

ZE3Bp(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        if kk_1 == 9 || kk_2 == 9 || kk_3 == 9
        ZE3Bp(i,j,ja,jb) = ZE3Bp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                        end
                    end
                end
            end
    end
    end
end

ZE4Bp(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                        if kk_1 == 9 || kk_2 == 9 || kk_3 == 9 || kk_4 == 9
        ZE4Bp(i,j,ja,jb) = ZE4Bp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5Bp(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                            if kk_1 == 9 || kk_2 == 9 || kk_3 == 9 || kk_4 == 9 || kk_5 == 9
        ZE5Bp(i,j,ja,jb) = ZE5Bp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6Bp(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
                                 if kk_1 == 9 || kk_2 == 9 || kk_3 == 9 || kk_4 == 9 || kk_5 == 9 || kk_6 == 9
        ZE6Bp(i,j,ja,jb) = ZE6Bp(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
mB1(i,j,ja,jb)=(ZE1Bp(i,j,ja,jb)+ZE2Bp(i,j,ja,jb)+ZE3Bp(i,j,ja,jb)+ZE4Bp(i,j,ja,jb)+ZE5Bp(i,j,ja,jb)+ZE6Bp(i,j,ja,jb))/ZE(i,j,ja,jb);

ZE1Bm(i,j,ja,jb) = (exp(-1/(kB*T)*Elink(10,i,j,ja,jb))); % for beta or gama direction
ZE2Bm(i,j,ja,jb) = 0;
for kk1 = 1:5
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:6
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
        if kk_1 == 10 || kk_2 == 10
        ZE2Bm(i,j,ja,jb) = ZE2Bm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb))));
        end
            end
    end
    end
end

ZE3Bm(i,j,ja,jb) = 0;
for kk1 = 1:4
    for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:5
            for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:6
                    for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        if kk_1 == 10 || kk_2 == 10 || kk_3 == 10
        ZE3Bm(i,j,ja,jb) = ZE3Bm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb))));
                        end
                    end
                end
            end
    end
    end
end

ZE4Bm(i,j,ja,jb) = 0;
for kk1 = 1:3
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:4
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:5
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:6
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                        if kk_1 == 10 || kk_2 == 10 || kk_3 == 10 || kk_4 == 10
        ZE4Bm(i,j,ja,jb) = ZE4Bm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb))));
                        end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE5Bm(i,j,ja,jb) = 0;
for kk1 = 1:2
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:3
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:4
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:5
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:6
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                            if kk_1 == 10 || kk_2 == 10 || kk_3 == 10 || kk_4 == 10 || kk_5 == 10
        ZE5Bm(i,j,ja,jb) = ZE5Bm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb))));
                            end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end

ZE6Bm(i,j,ja,jb) = 0; 
for kk1 = 1:1
for kk1b = 1:2
        kk_1 = (kk1-1)*2+kk1b;
    for kk2 = kk1+1:2
    for kk2b = 1:2
        kk_2 = (kk2-1)*2+kk2b;
                for kk3 = kk2+1:3
                for kk3b = 1:2
                        kk_3 = (kk3-1)*2+kk3b;
                        for kk4 = kk3+1:4
                        for kk4b = 1:2
                        kk_4 = (kk4-1)*2+kk4b;
                            for kk5 = kk4+1:5
                            for kk5b = 1:2
                            kk_5 = (kk5-1)*2+kk5b;
                                 for kk6 = kk5+1:6
                                 for kk6b = 1:2
                                 kk_6 = (kk6-1)*2+kk6b;
                                 if kk_1 == 10 || kk_2 == 10 || kk_3 == 10 || kk_4 == 10 || kk_5 == 10 || kk_6 == 10
        ZE6Bm(i,j,ja,jb) = ZE6Bm(i,j,ja,jb) + (exp(-1/(kB*T)*(Elink(kk_1,i,j,ja,jb) + Elink(kk_2,i,j,ja,jb) + Elink(kk_3,i,j,ja,jb) + Elink(kk_4,i,j,ja,jb) + Elink(kk_5,i,j,ja,jb)+ Elink(kk_6,i,j,ja,jb))));
                                 end
                                 end
                                 end
                            end
                            end
                        end
                        end
                    end
                end
            end
    end
    end
end
mB2(i,j,ja,jb)=(ZE1Bm(i,j,ja,jb)+ZE2Bm(i,j,ja,jb)+ZE3Bm(i,j,ja,jb)+ZE4Bm(i,j,ja,jb)+ZE5Bm(i,j,ja,jb)+ZE6Bm(i,j,ja,jb))/ZE(i,j,ja,jb);

    mHs(i,j,ja,jb)=(mx1(i,j,ja,jb)+mx2(i,j,ja,jb));
    mVs(i,j,ja,jb)=(my1(i,j,ja,jb)+my2(i,j,ja,jb));
    mAs(i,j,ja,jb)=2*mA1(i,j,ja,jb);
    mBs(i,j,ja,jb)=(mB1(i,j,ja,jb)+mB2(i,j,ja,jb));
    del(j,i,ja,jb)=abs(mHs(i,j,ja,jb)-mH(1,i))+2*abs(mVs(i,j,ja,jb)-mV(1,j))+abs(mAs(i,j,ja,jb)-mA(1,ja))+2*abs(mBs(i,j,ja,jb)-mB(1,jb));
        end
    end
    end
end
del0(it,1)=min(min(min(min(del))));
for j = 1:lm
    for i = 1:lm
        for ja = 1:lm
            for jb = 1:lm
       if del(i,j,ja,jb)==del0(it,1)
result(it,3)=mHs(i,j,ja,jb);
result(it,4)=mHs(i,j,ja,jb)+mVs(i,j,ja,jb)*2+mAs(i,j,ja,jb)+mBs(i,j,ja,jb)*2;
result(it,5)=mVs(i,j,ja,jb);
result(it,6)=mAs(i,j,ja,jb);
result(it,7)=mBs(i,j,ja,jb);
break;
       end
            end
        end
   end
end
%save del.mat del 
%save mHs.mat mHs 
%save mVs.mat mVs 
%save mAs.mat mAs 
%save mBs.mat mBs 
end
%
fileID = fopen('HCP3D_MFT_result.txt','w');
fprintf(fileID,'MFT results: nx, ntot, nx-E, ntot-E, nV-E, nA-E, nB-E \n');
for i = 1:size(result,1)
fprintf(fileID,'%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f \n',result(i,1:7));
end
fclose(fileID);
