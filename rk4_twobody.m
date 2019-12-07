function [rr,vv,time] = rk4_twobody(r0,v0,mu,dt,nsteps,a_thrust)
% function integrates two-body by runge-kutta 4th order
% INPUT
%   r0 : 3x1 initial position vector [km]
%   v0 : 3x1 initial velocity vector [km/s]
%   mu : gravitatonal parameter [km^3/s^2]
%   dt : time step size [s]
%   nsteps : number of steps to integrate
%   a_thrust : 3x1 constant acceleration vector (e.g. thrust)
% OUTPUT
%   rr : 3 x nsteps propagated position [km]
%   vv : 3 x nsteps propagated velocity [km/s]
%   time : 1 x nsteps propagated time array [s]
% Yuri Shimane, 2019.12.07
% ==================================================== %

 if ~exist('a_thrust','var')
     % if a_thrust is not given, set to 0
     a_thrust = 0;
 end

% initialize
rr = zeros(3,nsteps);
vv = zeros(3,nsteps);
time = zeros(1,nsteps);
% populate initial condition to array
rr(:,1) = r0;
vv(:,1) = v0;
sv_next = vertcat(r0,v0);

for i = 1:nsteps-1
    % Runge-Kutta 4th order coefficient k1
    dr = vv(:,i);
    dv = (-mu/norm(rr(:,i))^3 + a_thrust/norm(rr(:,i)))...
        .* rr(:,i);
    k1 = vertcat(dr,dv);
    
    % Runge-Kutta 4th order coefficient k2
    dr = vv(:,i) + (dt/2)*k1(4:6,1);
    dv = (-mu/norm(rr(:,i) + (dt/2)*k1(1:3,1))^3 ...
        + a_thrust/norm(rr(:,i) + (dt/2)*k1(1:3,1))) .* ...
        (rr(:,i) + (dt/2)*k1(1:3,1));
    k2 = vertcat(dr,dv);
    
    % Runge-Kutta 4th order coefficient k3
    dr = vv(:,i) + (dt/2)*k2(4:6,1);
    dv = (-mu/norm(rr(:,i) + (dt/2)*k1(1:3,1))^3 ...
        + a_thrust/norm(rr(:,i) + (dt/2)*k1(1:3,1))).* ...
        (rr(:,i) + (dt/2)*k2(1:3,1));
    k3 = vertcat(dr,dv);
    
    % Runge-Kutta 4th order coefficient k4
    dr = vv(:,i) + (dt)*k3(4:6,1);
    dv = (-mu/norm(rr(:,i) + (dt).*k1(1:3,1))^3 ...
        + a_thrust/norm(rr(:,i) + (dt).*k1(1:3,1))).* ...
        (rr(:,i) + (dt)*k3(1:3,1));
    k4 = vertcat(dr,dv);
    
    % propagate state-vector
    sv_next = sv_next + (dt/6)*(k1+2*k2+2*k3+k4);
    
    % populate position, velocity, time
    rr(:,i+1) = sv_next(1:3,1);
    vv(:,i+1) = sv_next(4:6,1);
    time(1,i+1) = time(1,i) + dt;
end

end