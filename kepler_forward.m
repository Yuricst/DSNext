function theta = kepler_forward(time,P,e)
% Function to propagate orbit forward with Kepler equation for elliptical
% orbit
% INPUT
%   time (float): time in the future to determine [sec]
%   P (float): orbital period [sec]
%   e (float): eccentricity
% OUTPUT
%   theta (float): true anomaly [deg]
% Yuri Shimane, 2019/12/07
% ============================================================ %

M = (2*pi/P)*time;

%...Set an error tolerance:
error = 1.e-8;

%...Select a starting value for E:
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

%...Iterate on Equation 3.17 until E is determined to within
%...the error tolerance:
ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

theta = 2*atand(sqrt((1+e)/(1-e)) * tan(E/2));

end