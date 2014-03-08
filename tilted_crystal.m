function abcd = tilted_crystal(d, varargin)
%Returns ABCD matrix for a tilted crystal
%
% Syntax = tilted_crystal(d, <theta>, <n>)
%
% d=Thickness of the crystal / planparallel plate (the real thickness, not the effective!)
% Theta and n are opitonal. Default = Theta = atan(1.76), n=1.76 (TiSA)
%
% 

if nargin>=2, theta=varargin{1}; else theta=atan(1.76); end

if nargin >=3, n=varargin{2}; else n=1.76; end

if nargin >=4, pol=varargin{3}; else pol='s'; end

theta_out = theta;
theta_in = asin(sin(theta)/n);

switch pol
    case 'p'
        abcd = [cos(theta_out)/cos(theta_in),0;0,cos(theta_in)/cos(theta_out)]...
               *[1,d/(n*cos(theta_in));0,1]...
               *[cos(theta_in)/cos(theta_out),0;0,cos(theta_out)/cos(theta_in)];
        %[ 1 d*n^2*(1-sin(theta)^2)/(n^2-sin(theta)^2)^1.5
        %        0                   1];
    case 's'
        abcd = [1,0;0,1]*[1,d/(n*cos(theta_in));0,1]*[1,0;0,1];
                %[1 d/(n^2-sin(theta)^2)^0.5
               % 0,1];
    otherwise
        error('Only two polarization states allowed (p or s)!');
end
end

