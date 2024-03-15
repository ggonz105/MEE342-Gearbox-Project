% Example Parameters %
P = 1500; % Applied Load (N) %
L = 1; % Length of Shaft (m)%
D = .05; % Diameter of Shaft (m)
E = 210*10^9; % Young's Modulus (Pa)
G = 79.3*10^9; % Shear Modulus (Pa)
T = 1000; % Applied Torque (N-m)
x = 0:.01:1; % Any given distance on shaft (m)

% Intial Calculations %
I = (pi/64) * D^4; % Moment of Inertia
J = (pi/32) * D^4; % Polar Moment of Inertia

% Stress Calculations%
M = P .* (L-x); % Moment Calculation for different points on the shaft
sigma = (M .* (D/2)) ./ I; % Normal Stress at different points on shaft (Pa)
tau = T * (D/2) / J; % Shear Stress (Pa) %
sigma_p1 = sigma./2 + sqrt((sigma./2).^2 + tau^2); % First Principal Stress (Pa)
sigma_p2 = sigma./2 - sqrt((sigma./2).^2 + tau^2); % Second Principal Stress (Pa)
sigma_prime = sqrt((sigma_p1.^2 - sigma_p1.*sigma_p2 + sigma_p2.^2)); % Simga prime used in fos calcualtion (Pa)

% Deflection of Shaft Calcualtion %
delta_bending = P .* (L-x).^3./(3*E*I); % Deflection due to Bending (m)

% Factor of Safety Calculation %
Sy = 370*10^6; % Typcial Yield Strength of AISI 1018 Steel
n = Sy./sigma_prime; % Factor of Safety along shaft