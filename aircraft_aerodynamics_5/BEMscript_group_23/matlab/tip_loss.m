function F = tip_loss(z, D, r_mean, phi)
% inputs (z, D, r, phi) calculating the tip loss factor, F, with Prandtl's
% model. phi in rads. if D/2 - r is very small, the exponent becomes really
% small and acos is complex. This is avoided by setting the exponent as 0
% if this happens

exponent = -(z*(D/2- r_mean)/(2*r_mean *sin(phi)));
F = 2/pi*acos(exp(exponent));
end
