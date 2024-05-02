% Fatigue and Yielding Analysis Intermediate Shaft
n = 2; % Minimum Factor of Safety (ny and nf)
kc = 1; % Loading Modification Factor
ke = 0.814; % 99 Percent Reliabilty Factor
kd = 1; % Assumed Temperature Factor
kb = 1; % Used to get Intial Diameter Guesses
kf = 1; % Assumed Misc Modification Factor
Se_p = 0.5*Sut; % Sut < 200kpsi

sqrt_a = 0.246 - (3.08 * 10^-3) * (Sut * 10^-3) + (1.51 * 10^-5) * (Sut * 10^-3)^2 - (2.67 * 10^-8) * (Sut * 10^-3)^3;
sqrt_as = 0.19 - (2.51 * 10^-3) * (Sut * 10^-3) + (1.35 * 10^-5) * (Sut * 10^-3)^2 - (2.67 * 10^-8) * (Sut * 10^-3)^3;

% Surface Condition Modification Factor (CD Steel)
a = 2.7;
b = -0.265;
ka = a*(Sut*10^-3)^b;

% Diameter at Bearing (D1)
Ma_0 = 3.1158*10^3;
Mm_0 = 0;
Ta_0 = 0;
Tm_0 = Ta;
Kt_0 = 2.7; % Assumed for Sharp Fillet
Kts_0 = 2.2; % Assumed for Sharp Fillet
Kf_0 = Kt_0;
Kfs_0 = Kts_0;

Se_0 = ka*kb*kc*ke*kf*Se_p;
d1 = ((16*n/pi)*((2*Kf_0*Ma_0)/Se_0) + sqrt(3*(Kfs_0*Tm_0)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter at Keyway of Gear 4 (D2)
Ma_a = 3.411*10^3;
Mm_a = 0;
Ta_a = 0;
Tm_a = Ta;
Kt_a = 2.7; % Assumed for Sharp Fillet
Kts_a = 2.2; % Assumed for SSharp Fillet
Kf_a = Kt_a;
Kfs_a = Kts_a;

Se_a = ka*kb*kc*ke*kf*Se_p;
d2 = ((16*n/pi)*((2*Kf_a*Ma_0)/Se_a) + sqrt(3*(Kfs_a*Tm_a)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Max Moment on Middle Shaft (D3)
Ma_b = sqrt((212.25^2)^2 + (2.7102*10^3)^2);
Mm_b = 0;
Ta_b = 0;
Tm_b = 0;
Kt_b = 2.14; % Assumed for Keyseat
Kts_b = 3; % Assumed for Keyseat
Kf_b = Kt_b;
Kfs_b = Kts_b;

Se_b = ka*kb*kc*ke*kf*Se_p;
d3 = ((16*n/pi)*((2*Kf_b*Ma_b)/Se_b) + sqrt(3*(Kfs_b*Tm_b)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Shoulder of Gear 3 (D4)
Ma_c = sqrt((215.1877)^2 + (851.9034)^2);
Mm_c = 0;
Ta_c = 0;
Tm_c = Ta;
Ktc = 2.7; % Assumed for Sharp Fillet
Ktsc = 2.2; % Assumed for Sharp Fillet
Kfc = Ktc;
Kfsc = Ktsc;

Se_c = ka*kb*kc*ke*kf*Se_p;
d4 = ((16*n/pi)*((2*Kfc*Ma_c)/Se_c) + sqrt(3*(Kfsc*Tm_c)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

%%
D1 = d1.*ones(1,10);
kb_1 = ones(1,10);
qd_1 = ones(1,10);
qs_1 = ones(1,10);
kf_1 = Kf_0*ones(1,10);
kfs_1 = Kfs_0*ones(1,10);
Se_1 = Se_0*ones(1,10);

for i = 2:10
    kb_1(i) = 0.879.*(D1(i-1)).^-0.107;
    qd_1(i) = (1+sqrt_a./0.02.*D1(i-1)).^-1;
    qs_1(i) = (1+sqrt_as./0.02.*D1(i-1)).^-1;
    kf_1(i) = 1 + qd_1(i-1).*(Kt_0-1);
    kfs_1(i) = 1 + qs_1(i-1).*(Kts_0-1);
    Se_1(i) = ka.*kb_1(i-1).*kc.*ke.*kf.*Se_p;
    D1(i) = ((16*n/pi).*(((2*kf_1(i-1).*Ma_0)./Se_1(i-1)) + sqrt(3.*(kfs_1(i-1).*Tm_0).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_1(4) * Ma_0)/(pi*D1(4)^3);
sigma_mp = sqrt(3)*((16*kfs_1(4)*Tm_0)/(pi*D1(4)^3));

ny_D1 = Sy/(sigma_ap + sigma_mp);
nf_D1 = ((sigma_ap/Se_1(4) + (sigma_mp/Sut)))^-1;

fprintf('Input Shaft Yielding FOS @ End of 4in Shaft = %f\n',ny_D1)
fprintf('Input Shaft Fatigue FOS @ End of 4in Shaft = %f\n', nf_D1)
%%
D2 = d2.*ones(1,10);
kb_2 = ones(1,10);
qd_2 = ones(1,10);
qs_2 = ones(1,10);
kf_2 = Kf_a*ones(1,10);
kfs_2 = Kfs_a*ones(1,10);
Se_2 = Se_a*ones(1,10);

for i = 2:10
    kb_2(i) = 0.879.*(D2(i-1)).^-0.107;
    qd_2(i) = (1+sqrt_a./0.02.*D2(i-1)).^-1;
    qs_2(i) = (1+sqrt_as./0.02.*D2(i-1)).^-1;
    kf_2(i) = 1 + qd_2(i-1).*(Kt_a-1);
    kfs_2(i) = 1 + qs_2(i-1).*(Kts_a-1);
    Se_2(i) = ka.*kb_2(i-1).*kc.*ke.*kf.*Se_p;
    D2(i) = ((16*n/pi).*(((2*kf_2(i-1).*Ma_a)./Se_2(i-1)) + sqrt(3.*(kfs_2(i-1).*Tm_a).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_2(4) * Ma_a)/(pi*D2(4)^3);
sigma_mp = sqrt(3)*((16*kfs_2(4)*Tm_a)/(pi*D2(4)^3));

ny_D2 = Sy/(sigma_ap + sigma_mp);
nf_D2 = ((sigma_ap/Se_2(4) + (sigma_mp/Sut)))^-1;

fprintf('Input Shaft Yielding FOS @ Bearing A = %f\n',ny_D2)
fprintf('Input Shaft Fatigue FOS @ Bearing A = %f\n', nf_D2)
%%
D3 = d3.*ones(1,10);
kb_3 = ones(1,10);
qd_3 = ones(1,10);
qs_3 = ones(1,10);
kf_3 = Kf_b*ones(1,10);
kfs_3 = Kfs_b*ones(1,10);
Se_3 = Se_b*ones(1,10);

for i = 2:10
    kb_3(i) = (0.91.*D3(i-1)).^-0.157;
    qd_3(i) = (1+sqrt_a./0.02.*D3(i-1)).^-1;
    qs_3(i) = (1+sqrt_as./0.02.*D3(i-1)).^-1;
    kf_3(i) = 1 + qd_3(i-1).*(Kt_b-1);
    kfs_3(i) = 1 + qs_3(i-1).*(Kts_b-1);
    Se_3(i) = ka.*kb_3(i-1).*kc.*ke.*kf.*Se_p;
    D3(i) = ((16*n/pi).*(((2*kf_3(i-1).*Ma_b)./Se_3(i-1)) + sqrt(3.*(kfs_3(i-1).*Tm_b).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_3(6) * Ma_b)/(pi*D3(6)^3);
sigma_mp = sqrt(3)*((16*kfs_3(6)*Tm_b)/(pi*D3(6)^3));

ny_D3 = Sy/(sigma_ap + sigma_mp);
nf_D3 = ((sigma_ap/Se_3(6) + (sigma_mp/Sut)))^-1;

fprintf('Input Shaft Yielding FOS @ Gear 2 = %f\n',ny_D3)
fprintf('Input Shaft Fatigue FOS @ Gear 2 = %f\n', nf_D3)
%%
D4 = d4.*ones(1,10);
kb_4 = ones(1,10);
qd_4 = ones(1,10);
qs_4 = ones(1,10);
kf_4 = Kfc*ones(1,10);
kfs_4 = Kfsc*ones(1,10);
Sed_4 = Se_c*ones(1,10);

for i = 2:10
    kb_4(i) = 0.879.*D4(i-1).^-0.107;
    qd_4(i) = (1+sqrt_a./0.02.*D4(i-1)).^-1;
    qs_4(i) = (1+sqrt_as./0.02.*D4(i-1)).^-1;
    kf_4(i) = 1 + qd_4(i-1).*(Kfc-1);
    kfs_4(i) = 1 + qs_4(i-1).*(Kfsc-1);
    Sed_4(i) = ka.*kb_4(i-1).*kc.*ke.*kf.*Se_p;
    D4(i) = ((16*n/pi).*((2*kf_4(i-1).*Ma_c)./Sed_4(i-1)) + sqrt(3.*(kfs_4(i-1).*Tm_c).^2)./Sut).^(1/3);
end
% 
sigma_ap = (32*kf_4(3)*Ma_c)/(pi*D4(3)^3);
sigma_mp = sqrt(3)*((16*kfs_4(3)*Tm_c)/(pi*D4(3)^3));

ny_D4 = Sy/(sigma_ap + sigma_mp);
nf_D4 = ((sigma_ap/Sed_4(3) + (sigma_mp/Sut)))^-1;

fprintf('Input Shaft Yielding FOS @ Bearing C = %f\n',ny_D4)
fprintf('Input Shaft Fatigue FOS @ Bearing C = %f\n', nf_D4)