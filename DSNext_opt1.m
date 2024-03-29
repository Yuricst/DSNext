%% DSNext back of the envelope study
% Version 02
% Trying to optimize for good orbit
% Yuri Shimane, 2019/12/07
% ============================================================= %

% house keeping
clear; close all; clc;

% parameters
mu = 132712440018; % sun gravitational parameter
au2km = 149597871; % conversion from AU to km

% Earth orbit (assume circular)
Earth.e  = 0;
Earth.a  = au2km; % [km]
Earth.h  = sqrt(Earth.a*mu*(1-Earth.e^2));
Earth.vp = sqrt(mu/Earth.a);
Earth.P  = 2*pi*sqrt(Earth.a^3/mu); % orbital period [s]
Earth.r0 = [Earth.a; 0; 0];           % initial position vector [km]
Earth.v0 = [0; Earth.vp; 0];          % initial velocity vector [km]

% Mars orbit (assume circular)
Mars.e  = 0;
Mars.a  = 1.524*au2km; % [km]
Mars.h  = sqrt(Mars.a*mu*(1-Mars.e^2));
Mars.vp = sqrt(mu/Mars.a);
Mars.P  = 2*pi*sqrt(Mars.a^3/mu);  % orbital period [s]
Mars.r0 = [Mars.a; 0; 0];            % initial position vector [km]
Mars.v0 = [0; Mars.vp; 0];           % initial velocity vector [km]


%% Relay satellite orbit
Rel01.e  = 0;
inc = 10;

a_range_AU = linspace(0.1,1.4,inc);

% iterate over allowable range of semi-major axis
for i = 1:inc
    Rel01.a  = a_range_AU(i)*au2km;
    Rel01.h  = sqrt(Rel01.a*mu*(1-Rel01.e^2));
    Rel01.theta0 = 180;                          % INPUT PARAMETER
    Rel01.vp = sqrt(mu/Rel01.a);
    Rel01.P  = 2*pi*sqrt(Rel01.a^3/mu);  % orbital period [s]
    Rel01.r0 = Rel01.h^2/(mu*(1 + Rel01.e*cosd(Rel01.theta0)))*...
        [cosd(Rel01.theta0); sind(Rel01.theta0); 0];            % initial position vector [km]
    Rel01.v0 = mu/Rel01.h *...
        [-sind(Rel01.theta0); Rel01.e + cosd(Rel01.theta0); 0];
    
    % for loop
    
    
end

% communication band
f = 34.2 * 10^9; % [Hz]
c = 2.99792458 * 10^8; % [m/s]
lambda = c/f;  % [m]
