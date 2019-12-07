%% DSNext back of the envelope study
% 
% 
% Yuri Shimane, listening to Goldberg Variation, 2019/12/07
% ============================================================= %

% house keeping
clear; close all; clc;

% parameters
mu = 132712000000; % sun gravitational parameter
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
Rel01.a  = 0.8*au2km;
Rel01.h  = sqrt(Rel01.a*mu*(1-Rel01.e^2));
Rel01.theta0 = 100;                          % INPUT PARAMETER
Rel01.vp = sqrt(mu/Rel01.a);
Rel01.P  = 2*pi*sqrt(Rel01.a^3/mu);  % orbital period [s]
Rel01.r0 = Rel01.h^2/(mu*(1 + Rel01.e*cosd(Rel01.theta0)))*...
    [cosd(Rel01.theta0); sind(Rel01.theta0); 0];            % initial position vector [km]
Rel01.v0 = mu/Rel01.h *...
    [-sind(Rel01.theta0); Rel01.e + cosd(Rel01.theta0); 0];

% communication band
f = 34.2 * 10^9; % [Hz]
c = 2.99792458 * 10^8; % [m/s]
lambda = c/f;  % [m]

%% Propagate prbit with Kepler's method
nsteps = 1000;
tmax = Mars.P*5;
dt = tmax/nsteps;
time = linspace(0,tmax,nsteps);
for i = 1:nsteps
    % compute angular position
    Earth.theta(1,i)  = kepler_forward(time(1,i),Earth.P,Earth.e);
    Mars.theta(1,i)   = kepler_forward(time(1,i),Mars.P,Mars.e);
    Rel01.theta(1,i)  = kepler_forward(time(1,i),Rel01.P,Rel01.e) + Rel01.theta0;
    % compute position vector in PF frame
    Earth.rPF(:,i) = Earth.h^2/(mu*(1+Earth.e*cosd(Earth.theta(1,i))))*...
        [cosd(Earth.theta(1,i)); sind(Earth.theta(1,i)); 0];
    Mars.rPF(:,i)  = Mars.h^2/(mu*(1+Mars.e*cosd(Mars.theta(1,i))))*...
        [cosd(Mars.theta(1,i)); sind(Mars.theta(1,i)); 0];
    Rel01.rPF(:,i)  = Rel01.h^2/(mu*(1+Rel01.e*cosd(Rel01.theta(1,i))))*...
        [cosd(Rel01.theta(1,i)); sind(Rel01.theta(1,i)); 0];
    % compute relative distance at each time step
    dr_EM(1,i) = norm(Mars.rPF(:,i) - Earth.rPF(:,i));
    dr_Rel01E(1,i) = norm(Rel01.rPF(:,i) - Earth.rPF(:,i));
    dr_Rel01M(1,i) = norm(Rel01.rPF(:,i) - Mars.rPF(:,i));
    
    % Communication free space loss
    Lfs_EM(1,i) = 4*pi*dr_EM(1,i)^2/lambda;
    Lfs_EM_dB(1,i) = -10*log10(Lfs_EM(1,i));
    
    Lfs_Rel01E(1,i) = 4*pi*dr_Rel01E(1,i)^2/lambda;
    Lfs_Rel01E_dB(1,i) = -10*log10(Lfs_Rel01E(1,i));
    
    Lfs_Rel01M(1,i) = 4*pi*dr_Rel01M(1,i)^2/lambda;
    Lfs_Rel01M_dB(1,i) = -10*log10(Lfs_Rel01M(1,i));
end

%% static plot
figure(11)
stp_interest = 1;
plot(Earth.rPF(1,:)/au2km,Earth.rPF(2,:)/au2km,'--b')
hold on
plot(Earth.rPF(1,stp_interest)/au2km,Earth.rPF(2,stp_interest)/au2km,'ob','MarkerSize',8)
hold on
plot(Mars.rPF(1,:)/au2km,Mars.rPF(2,:)/au2km,'--r')
hold on
plot(Mars.rPF(1,stp_interest)/au2km,Mars.rPF(2,stp_interest)/au2km,'or','MarkerSize',8)
hold on
plot(0,0,'*','MarkerSize',8,'color','#EDB120')
hold on
plot(Rel01.rPF(1,:)/au2km,Rel01.rPF(2,:)/au2km,'--k')
hold on
plot(Rel01.rPF(1,stp_interest)/au2km,Rel01.rPF(2,stp_interest)/au2km,'^k','MarkerSize',6)
axis equal
grid on; xlabel('x [AU]'); ylabel('[AU]')

figure(12)
subplot(2,1,1)
plot(time,dr_EM/au2km,'--k')
hold on
plot(time,dr_Rel01M/au2km,'-.m')
hold on
plot(time,dr_Rel01E/au2km,'-.c')
grid on;
xlabel('time'); ylabel('Relative distancd [AU]')
legend('Mars-Earth','SC-Mars','SC-Earth')

subplot(2,1,2)
plot(time,Lfs_EM_dB,'--k')
hold on
plot(time,Lfs_Rel01M_dB,'-.m')
hold on
plot(time,Lfs_Rel01E_dB,'-.c')
grid on;
xlabel('time'); ylabel('Free space loss [dB]')
legend('Mars-Earth','SC-Mars','SC-Earth')


%% dynamic plot
dynamicplot = input('Dynamic plot? [y/n]: ','s');
if dynamicplot == 'y'
for i = 1:10:nsteps
    figure(15)
    % entire orbit
    plot(Earth.rPF(1,:)/au2km,Earth.rPF(2,:)/au2km,':b')
    hold on
    plot(Mars.rPF(1,:)/au2km,Mars.rPF(2,:)/au2km,':r')
    hold on
    plot(Rel01.rPF(1,:)/au2km,Rel01.rPF(2,:)/au2km,':k')
    hold on
    plot(0,0,'*','MarkerSize',8,'color','#EDB120')
    % moving position
    plot(Earth.rPF(1,i)/au2km,Earth.rPF(2,i)/au2km,'ob')
    hold on
    plot(Mars.rPF(1,i)/au2km,Mars.rPF(2,i)/au2km,'or')
    hold on
    plot(Rel01.rPF(1,i)/au2km,Rel01.rPF(2,i)/au2km,'^k')
    hold on
    plot([Earth.rPF(1,i)/au2km Mars.rPF(1,i)/au2km],...
         [Earth.rPF(2,i)/au2km Mars.rPF(2,i)/au2km],'--k')
    hold on
    plot([Rel01.rPF(1,i)/au2km Mars.rPF(1,i)/au2km],...
         [Rel01.rPF(2,i)/au2km Mars.rPF(2,i)/au2km],'-.m')
    hold on
    plot([Earth.rPF(1,i)/au2km Rel01.rPF(1,i)/au2km],...
         [Earth.rPF(2,i)/au2km Rel01.rPF(2,i)/au2km],'-.c')
    axis equal
    xlim([-2 2]); ylim([-2 2]);
    grid on; 
    xlabel('x [AU]'); ylabel('[AU]')
    %pause(0.0001)
    hold off
end
end


%% Comm. budget 
% signal source (some Mars mission / architecture?)

% antenna parameters (assume identical antenna for down/up link
Rel01.antenna.diam = 3;  % antenna diameter
Rel01.antenna.eta = 0.7; % antenna efficiency
Rel01.antenna.gain = Rel01.antenna.eta*(pi*Rel01.antenna.diam/lambda)^2;
Rel01.antenna.gain_dB = 10*log10(Rel01.antenna.gain);
% downlink components
Rel01.downlink.cable_loss_dB = 3;
Rel01.downlink.LNA_gain_dB = 10;  % low-noise amplifier gain

Rel01.downlink.dPower_dB = Rel01.antenna.gain_dB - Rel01.downlink.cable_loss_dB +...
    Rel01.downlink.LNA_gain_dB;
Rel01.downlink.dPower = 10^(Rel01.downlink.dPower_dB/10);

% uplink components
Rel01.uplink.dPin_dB = Rel01.downlink.dPower_dB;




