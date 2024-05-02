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
Ma_d = 0;
Mm_d = 0;
Ta_d = 0;
Tm_d = Ta;
Kt_d = 2.7; % Assumed for Sharp Fillet
Kts_d = 2.2; % Assumed for Sharp Fillet
Kf_d = Kt_d;
Kfs_d = Kts_d;

Se_d = ka*kb*kc*ke*kf*Se_p;
d9 = ((16*n/pi)*((2*Kf_d*Ma_d)/Se_d) + sqrt(3*(Kfs_d*Tm_d)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter at Keyway of Gear 4 (D2)
Ma_e = sqrt((Wr_45 * (L5+L6))^2 + (Wt_45 * (L5+L6))^2);
Mm_e = 0;
Ta_e = 0;
Tm_e = 0;
Kt_e = 2.14; % Assumed for Keyseat
Kts_e = 3; % Assumed for Keyseat
Kf_e = Kt_e;
Kfs_e = Kts_e;

Se_e = ka*kb*kc*ke*kf*Se_p;
d8 = ((16*n/pi)*((2*Kf_e*Ma_d)/Se_e) + sqrt(3*(Kfs_e*Tm_e)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Max Moment on Middle Shaft (D3)
Ma_3 = sqrt((413.24^2) + (223.32^2));
Mm_3 = 0;
Ta_Int = 0;
Tm_3 = Ta;
Kt_3 = 2.14; % Assumed for Keyseat
Kts_3 = 3; % Assumed for Keyseat
Kf_3 = Kt_3;
Kfs_3 = Kts_3;

Se_3 = ka*kb*kc*ke*kf*Se_p;
d7 = ((16*n/pi)*((2*Kf_3*Ma_3)/Se_3) + sqrt(3*(Kfs_3*Tm_3)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Shoulder of Gear 3 (D4)
Ma_f = sqrt((627.28^2) + (275.90)^2);
Mm_f = 0;
Ta_f = 0;
Tm_f = Ta;
Ktf = 2.7; % Assumed for Sharp Fillet
Ktsf = 2.2; % Assumed for Sharp Fillet
Kff = Ktf;
Kfsf = Ktsf;

Sef = ka*kb*kc*ke*kf*Se_p;
d6 = ((16*n/pi)*((2*Kff*Ma_f)/Sef) + sqrt(3*(Kfsf*Tm_f)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Shoulder to Right of Gear 3 (D5)
Ma_g = sqrt((420.21^2) + (171.24)^2);
Mm_g = 0;
Ta_g = 0;
Tm_g = 0;
Ktg = 2.7; % Assumed for Sharp Fillet
Ktsg = 2.2; % Assumed for Sharp Fillet
Kfg = Ktg;
Kfsg = Ktsg;

Seg = ka*kb*kc*ke*kf*Se_p;
d5 = ((16*n/pi)*((2*Kfg*Ma_g)/Seg) + sqrt(3*(Kfsg*Tm_g)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing
%%
D9 = d9.*ones(1,10);
kb_9 = ones(1,10);
qd_9 = ones(1,10);
qs_9 = ones(1,10);
kf_9 = Kfg*ones(1,10);
kfs_9 = Kfsg*ones(1,10);
Se_9 = Seg*ones(1,10);

for i = 2:10
    kb_9(i) = 0.879.*(D9(i-1)).^-0.107;
    qd_9(i) = (1+sqrt_a./0.02.*D9(i-1)).^-1;
    qs_9(i) = (1+sqrt_as./0.02.*D9(i-1)).^-1;
    kf_9(i) = 1 + qd_9(i-1).*(Kt_d-1);
    kfs_9(i) = 1 + qs_9(i-1).*(Kts_d-1);
    Se_9(i) = ka.*kb_9(i-1).*kc.*ke.*kf.*Se_p;
    D9(i) = ((16*n/pi).*(((2*kf_9(i-1).*Ma_d)./Se_9(i-1)) + sqrt(3.*(kfs_9(i-1).*Tm_d).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_9(4) * Ma_d)/(pi*D9(4)^3);
sigma_mp = sqrt(3)*((16*kfs_9(4)*Tm_d)/(pi*D9(4)^3));

ny_D9 = Sy/(sigma_ap + sigma_mp);
nf_D9 = ((sigma_ap/Se_9(4) + (sigma_mp/Sut)))^-1;
fprintf('Intermediate Shaft Yielding FOS @ Left Bearing = %f\n',ny_D9)
fprintf('Intermediate Shaft Fatigue FOS @ Left Bearing = %f\n', nf_D9)
%%
D8 = d8.*ones(1,10);
kb_8 = ones(1,10);
qd_8 = ones(1,10);
qs_8 = ones(1,10);
kf_8 = Kff*ones(1,10);
kfs_8 = Kfsf*ones(1,10);
See_8 = Sef*ones(1,10);

for i = 2:10
    kb_8(i) = 0.879.*(D8(i-1)).^-0.107;
    qd_8(i) = (1+sqrt_a./0.02.*D8(i-1)).^-1;
    qs_8(i) = (1+sqrt_as./0.02.*D8(i-1)).^-1;
    kf_8(i) = 1 + qd_8(i-1).*(Kt_e-1);
    kfs_8(i) = 1 + qs_8(i-1).*(Kts_e-1);
    See_8(i) = ka.*kb_8(i-1).*kc.*ke.*kf.*Se_p;
    D8(i) = ((16*n/pi).*(((2*kf_8(i-1).*Ma_e)./See_8(i-1)) + sqrt(3.*(kfs_8(i-1).*Tm_e).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_8(2) * Ma_e)/(pi*D8(2)^3);
sigma_mp = sqrt(3)*((16*kfs_8(2)*Tm_e)/(pi*D8(2)^3));

ny_D8 = Sy/(sigma_ap + sigma_mp);
nf_D8 = ((sigma_ap/See_8(2) + (sigma_mp/Sut)))^-1;

fprintf('Intermediate Shaft Yielding FOS @ Gear 4 = %f\n',ny_D8)
fprintf('Intermediate Shaft Fatigue FOS @ Gear 4 = %f\n', nf_D8)
%%
D7 = d7.*ones(1,10);
kb_7 = ones(1,10);
qd_7 = ones(1,10);
qs_7 = ones(1,10);
kf_7 = Kf_3*ones(1,10);
kfs_7 = Kfs_3*ones(1,10);
Sed_7 = Seg*ones(1,10);

for i = 2:10
    kb_7(i) = 0.879.*(D7(i-1)).^-0.107;
    qd_7(i) = (1+sqrt_a./0.02.*D7(i-1)).^-1;
    qs_7(i) = (1+sqrt_as./0.02.*D7(i-1)).^-1;
    kf_7(i) = 1 + qd_7(i-1).*(Kt_3-1);
    kfs_7(i) = 1 + qs_7(i-1).*(Kts_3-1);
    Sed_7(i) = ka.*kb_7(i-1).*kc.*ke.*kf_7(i-1).*Se_p;
    D7(i) = ((16*n/pi).*(((2*kf_7(i-1).*Ma_3)./Sed_7(i-1)) + sqrt(3.*(kfs_7(i-1).*Tm_3).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_7(7) * Ma_3)/(pi*D7(7)^3);
sigma_mp = sqrt(3)*((16*kfs_7(7)*Tm_3)/(pi*D7(7)^3));

ny_D7 = Sy/(sigma_ap + sigma_mp);
nf_D7 = ((sigma_ap/Sed_7(7) + (sigma_mp/Sut)))^-1;

fprintf('Intermediate Shaft Yielding FOS @ Center Shaft = %f\n',ny_D7)
fprintf('Intermediate Shaft Fatigue FOS @ Center Shaft = %f\n', nf_D7)
%%
D6 = d6.*ones(1,10);
kb_6 = ones(1,10);
qd_6 = ones(1,10);
qs_6 = ones(1,10);
kf_6 = Kf_e*ones(1,10);
kfs_6 = Kfs_e*ones(1,10);
Sed_6 = Se_e*ones(1,10);

for i = 2:10
    kb_6(i) = 0.879.*D6(i-1).^-0.107;
    qd_6(i) = (1+sqrt_a./0.02.*D6(i-1)).^-1;
    qs_6(i) = (1+sqrt_as./0.02.*D6(i-1)).^-1;
    kf_6(i) = 1 + qd_6(i-1).*(Kff-1);
    kfs_6(i) = 1 + qs_6(i-1).*(Kfsf-1);
    Sed_6(i) = ka.*kb_6(i-1).*kc.*ke.*kf.*Se_p;
    D6(i) = ((16*n/pi).*((2*kf_6(i-1).*Ma_f)./Sed_6(i-1)) + sqrt(3.*(kfs_6(i-1).*Tm_f).^2)./Sut).^(1/3);
end
% 
sigma_ap = (32*kf_6(3)*Ma_f)/(pi*D6(3)^3);
sigma_mp = sqrt(3)*((16*kfs_6(3)*Tm_f)/(pi*D6(3)^3));

ny_D6 = Sy/(sigma_ap + sigma_mp);
nf_D6 = ((sigma_ap/Sed_6(3) + (sigma_mp/Sut)))^-1;

fprintf('Intermediate Shaft Yielding FOS @ Gear 2 = %f\n',ny_D6)
fprintf('Intermediate Shaft Fatigue FOS @ Gear 2 = %f\n', nf_D6)
%%
% Iterations to lower D5
D5 = d5.*ones(1,10);
kb_5 = ones(1,10);
qd_5 = ones(1,10);
qs_5 = ones(1,10);
kf_5 = Kf_d*ones(1,10);
kfs_5 = Kfs_d*ones(1,10);
Sed_5 = Se_d*ones(1,10);

for i = 2:10
    kb_5(i) = 0.879.*(D5(i-1)).^-0.107;
    qd_5(i) = (1+sqrt_a./0.02.*D5(i-1)).^-1;
    qs_5(i) = (1+sqrt_as./0.02.*D5(i-1)).^-1;
    kf_5(i) = 1 + qd_5(i-1).*(Ktg-1);
    kfs_5(i) = 1 + qs_5(i-1).*(Ktsg-1);
    Sed_5(i) = ka.*kb_5(i-1).*kc.*ke.*kf.*Se_p;
    D5(i) = ((16*n/pi).*(((2*kf_5(i-1).*Ma_d)./Sed_5(i-1)) + sqrt(3.*(kfs_5(i-1).*Tm_d).^2)./Sut)).^(1/3);
end
sigma_ap = (32*kf_5(1) * Ma_g)/(pi*D5(1)^3);
sigma_mp = sqrt(3)*((16*kfs_5(1)*Tm_g)/(pi*D5(1)^3));

ny_D5 = Sy/(sigma_ap + sigma_mp);
nf_D5 = ((sigma_ap/Sed_5(1) + (sigma_mp/Sut)))^-1;

fprintf('Intermediate Shaft Yielding FOS @ Bearing G = %f\n',ny_D5)
fprintf('Intermediate Shaft Fatigue FOS @ Bearing G = %f\n', nf_D5)