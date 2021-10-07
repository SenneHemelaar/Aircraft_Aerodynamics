function [vx,vz] = VOR2D(Gamma,xi,zi,xj,zj)
%This function calculates the influence of a lumped vortex element.
%Therefore, vx and vz are the velocity components in the x direction and z
%direction respectively in a collocation point i induced by a vortex placed
%in a vortex point j.
rx=xi-xj;
rz=zi-zj;
r=sqrt(rx^2+rz^2);

vx=(Gamma/(2*pi*r^2))*rz;
vz=-(Gamma/(2*pi*r^2))*rx;
end

