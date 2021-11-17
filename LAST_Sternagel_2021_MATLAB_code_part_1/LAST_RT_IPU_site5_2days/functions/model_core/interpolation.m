%% Function for the interpolation of v and D

function [v, D] = interpolation(vup, vlow, Dup, Dlow, zup, zlow, zstar)

zint = (abs(zstar)-abs(zup)) / (abs(zlow)-abs(zup));
v = vup - (zint*(vup-vlow));
D = Dup - (zint*(Dup-Dlow));

end