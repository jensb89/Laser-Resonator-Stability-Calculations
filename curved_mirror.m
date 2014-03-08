function abcd = curved_mirror(R,varargin)
%Returns ABCD Matrrix for reflection on curved/spherical mirror with radius R

if nargin >= 2, theta = varargin{1}; else theta = 0; end

if nargin >=3, pol = varargin{2}; else pol = 's'; end

switch pol
    case 'p'
        R_eff = R*cos(theta);
    case 's'
        R_eff = R/cos(theta);
    otherwise
        error('Only two polarization states allowed (p or s)!');
end

abcd=[
    1       0
    -2/R_eff    1
];

end