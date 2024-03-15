% Example Parameters %
P = 1500; 
L = 1; 
D = .05; 
E = 210*10^9; 
G = 79.3*10^9;
T = 1000; 
x = 0:.01:1; % Any given distance on the shaft (m) %


I = (pi/64) * D^4; 
J = (pi/32) * D^4; 

% Stress Calculations %
M = P .* (L-x); 
sigma = (M .* (D/2)) ./ I; % Stress due to bending %
tau = T * (D/2) / J; % Stress due to Torque %
sigma_p1 = sigma./2 + sqrt((sigma./2).^2 + tau^2); 
sigma_p2 = sigma./2 - sqrt((sigma./2).^2 + tau^2);
sigma_prime = sqrt((sigma_p1.^2 - sigma_p1.*sigma_p2 + sigma_p2.^2)); % Principal Stress used for fos Calculation %

% Deflection Calculation %
delta_bending = P .* (L-x).^3./(3*E*I); 

% Factor of Safety Calculation
Sy = 370*10^6; 
n = Sy./sigma_prime; 