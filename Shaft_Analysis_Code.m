% Example Parameters %
P = 1500; 
L = 1; 
D = .05; 
E = 210*10^9; 
G = 79.3*10^9;
T = 1000; 
x = 0:.01:1; 


I = (pi/64) * D^4; 
J = (pi/32) * D^4; 


M = P .* (L-x); 
sigma = (M .* (D/2)) ./ I; 
tau = T * (D/2) / J; 
sigma_p1 = sigma./2 + sqrt((sigma./2).^2 + tau^2); 
sigma_p2 = sigma./2 - sqrt((sigma./2).^2 + tau^2);
sigma_prime = sqrt((sigma_p1.^2 - sigma_p1.*sigma_p2 + sigma_p2.^2));


delta_bending = P .* (L-x).^3./(3*E*I); 


Sy = 370*10^6; 
n = Sy./sigma_prime; 