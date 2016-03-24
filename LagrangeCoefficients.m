function [ r_vector,v_vector ] = LagrangeCoefficients( r0,v0,dtheta,muo )
% this function get r and v after true anomaly (dtheta) by Lagrange Coefficients
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% r0 : initial position vector (1x3) (km)
% v0 : initial velocity vector (1x3) (km/s)
% dtheta : true anomaly at which we calculate r and v (degree)
% muo : gravitational parameter (km^3/s^2)
%% outputs:
% r_vector : position vector at dtheta (1x3) (km)
% v_vector : velocity vector at dtheta (1x3) (km/s)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
mag_r0=norm(r0) % km
mag_v0=norm(v0) % km/s
v_r0=dot(r0,v0)/mag_r0 % km/s
h=mag_r0*sqrt(mag_v0^2-v_r0^2) %km^2/s
r=h^2/muo/(1+(h^2/muo/mag_r0-1)*cosd(dtheta)-h*v_r0/muo*sind(dtheta)) % km
f=1-muo*r/h^2*(1-cosd(dtheta))
g=r*mag_r0/h*sind(dtheta) % s
f_dot=muo/h*(1-cosd(dtheta))/sind(dtheta)*(muo/h^2*(1-cosd(dtheta))-1/mag_r0-1/r)
g_dot=1-muo*mag_r0/h^2*(1-cosd(dtheta))
r_vector=f*r0+g*v0 % km
v_vector=f_dot*r0+g_dot*v0 % km/s
end