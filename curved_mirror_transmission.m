function abcd = curved_mirror_transmission(R,d,varargin)
%Returns ABCD Matrrix for Transmission through curved/spherical mirror with radius R

if nargin >=2, theta=varargin{1}; else theta=13/180*pi; end

if nargin >=3, n=varargin{2}; else n=1.5; end

if nargin >=4, pol=varargin{3}; else pol='s'; end

nVak = 1;
theta_out = theta;
theta_in = asin(sin(theta)/n);

switch pol
    case 'p'
        delta_n = (nVak * cos(theta_in) -n*cos(theta_out))/(cos(theta_in)*cos(theta_out));
        abcd = [cos(theta_out)/cos(theta_in),0;0,cos(theta_in)/cos(theta_out)]...
               *[1,d/(n*cos(theta_in));0,1]...
               *[cos(theta_in)/cos(theta_out),0;delta_n/R,cos(theta_out)/cos(theta_in)];
        %[ 1 d*n^2*(1-sin(theta)^2)/(n^2-sin(theta)^2)^1.5
        %        0                   1];
    case 's'
        delta_n = (nVak * cos(theta_in) -n*cos(theta_out));
        abcd = [1,0;0,1]*[1,d/(n*cos(theta_in));0,1]*[1,0;delta_n/R,1];
                %[1 d/(n^2-sin(theta)^2)^0.5
               % 0,1];
    otherwise
        error('Only two polarization states allowed (p or s)!');
end
end