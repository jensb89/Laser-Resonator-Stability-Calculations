function abcd=free(L,varargin)
% Returns the ABCD matrix for free space propagation.
%
% SYNTAX: abcd=free(L <,n>)

if nargin>=2, n=varargin{1}; else n=1; end

abcd=[
    1   L/n
    0   1
];
