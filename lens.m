function abcd=lens(f)
%Returns ABCD matrix for a lens with focal length f

abcd=[
    1       0
    -1/f    1
];