%% Nd:Yag Pump - beam
%---()----()-----|(----/ /
%L1 f1 L2 f2 L3  f3 L4 K_in

nK = 1.76;
thetaK = atan(nK);
lambda = 532E-9;
pol = 's';

L1 = free(10E-2);
L2 = free(11E-2);
L3 = free(3.25E-2);
L4 = free(5E-2-1.6E-3);
f1 = lens(-1);
f2 = lens(10E-2);
theta = 13/180*pi;
n = 1.5195; %N-BK7 at 532nm
theta_in = asin(sin(theta)/n);
f3 = curved_mirror_transmission(0.1,10E-3,theta,n,pol);


if strcmp(pol,'p')
    nVak = 1;
    theta_out = thetaK;
    theta_in = asin(sin(theta_out)*nVak/nK);
    K_in = [cos(theta_out)/cos(theta_in) 0;0 cos(theta_in)/cos(theta_out)]; %Curved_Interface, arbitrary incidence, with R=inf
else
    K_in = eye(2);
end

%M = f2 * L2 * f1 * L1
%M = f3 * L3 * f2 * L2 * f1 * L1
M = K_in * L4 * f3 * L3 * f2 * L2 * f1 * L1;

w0 = 2.3E-3; %Millennia Pro (w0 = 2.3E-3)

z_pump=linspace(0,6E-3,300);

A=M(1,1); B=M(1,2); C=M(2,1); D=M(2,2);
q1 =  1/(-1i*lambda/pi/w0^2);
q2 = (A*q1+B)/(C*q1+D);

n=1.76;
w_pump  = sqrt(-lambda/pi./imag(1./(q2+z_pump/n))); 

figure;
plot(z_pump,w_pump,'r');
xlabel('Distance z inside Crystal / m');
ylabel('Beam diameter \omega(z) / m');