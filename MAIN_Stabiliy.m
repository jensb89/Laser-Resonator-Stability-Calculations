%MAIN RESONATOR STABILITY
%03.2014, Jens Brauer

clc
clear all
close all

%% Setup
nK = 1.76;
res = Resonator(3.2E-3,...                   %TiSa-crystal thickness / m
                 nK,...                      %refractive index of crystal
                 atan(nK),...                %angle of incidence on crystal (Brewster) / rad
                 0.1,...                     %focusing mirrors radius of curvature / m
                 1.84,...                    % Length of Resonator (L_ges = L1+L2+R+dK) / m
                 800e-9,...                  %Wavelength / m
                 's');                       %Polarisation ('s' or 'p')

%% Optimal folding angle Theta_ges = Theta_1 + Theta_2
theta_ges=res.getTheta();
fprintf('Theta_ges = %1.2fdeg\n',theta_ges/pi*180)

%% Set Distribution of Arm Lengths + folding angle
res.setArmLengthDist(1/3); %Arm1 = val*(L1+L2); Arm2 = (1-val)*(L1+L2)
res.setThetaDist(2/3); %Theta1 = val*theta_ges, theta2 = (1-val)*theta_ges
res.setRelCrystalPos(0.5); %Crystal in the middle of the cavity

%% Calc Stability for 's' and 'p' for different cavity mirror distances 
delta = linspace(0,5E-3); 

for i=1:length(delta)
    res.setDelta(delta(i)); %deviation from Cavity distance 2*f+dK
    
    [~,~,ws_arm1(i),ws_arm2(i),Ss_arm1(i),Ss_arm2(i)] = res.calcRoundTrip();
    res.setPolarisation('p');
    [~,~,wp_arm1(i),wp_arm2(i),Sp_arm1(i),Sp_arm2(i)] = res.calcRoundTrip();
    res.setPolarisation('s');
    
end
res.estimation_Delta(); %Prints estimated stability boundaries in console window
% Here delta = 0 == delta + (dK_eff - dK)!!
%% PLotting
%x = (res.R+delta)/1E-3;
x=delta/1E-3;

figure;
plot(x,ws_arm1/1E-3,x,wp_arm1/1E-3);
xlabel('Delta / mm');
ylabel('Strahldurchmesser / mm');
legend('S','P');
title(['Arm1 = ',num2str(res.L1),'mm']);

figure;
plot(x,ws_arm2/1E-3,x,wp_arm2/1E-3);
xlabel('Delta / mm');
ylabel('Strahldurchmesser / mm');
legend('S','P');
title(['Arm2 = ',num2str(res.L2),'mm']);

figure;
plot(x,Sp_arm1,x,Ss_arm1);
xlabel('Delta / mm');
ylabel('Stability factor S');
hold on
plot(x,ones(length(x)),'r');
xlabel('Delta / mm');
ylabel('Stability factor S');
ylim([0 3]);

%% Beam Inside Crystal
res.setDelta(2E-3);
res.calcRoundTrip();
res.setPolarisation('s')
[zs,ws] = res.beamInsideCrystal();
res.setPolarisation('p')
res.calcRoundTrip();
[zp,wp] = res.beamInsideCrystal();
figure;
plot(zs,ws,zp,wp);
xlabel('Distance z inside Crystal / m');
ylabel('Beam diameter \omega(z) / m');
legend('S','P');
crystal_end = res.dK/(cos(asin(sin(res.thetaK)/res.nK)));
hold on 
plot(crystal_end*ones(1,length(zs)),linspace(0,max(wp),length(zs)),'r');