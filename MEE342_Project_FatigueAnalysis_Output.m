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
Ma_10 = 0;
Mm_10 = 0;
Ta_10 = 0;
Tm_10 = Ta;
Kt_10 = 2.7; % Assumed for Sharp Fillet
Kts_10 = 2.2; % Assumed for Sharp Fillet
Kf_10 = Kt_10;
Kfs_10 = Kts_10;

Se_10 = ka*kb*kc*ke*kf*Se_p;
d10 = ((16*n/pi)*((2*Kf_10*Ma_10)/Se_10) + sqrt(3*(Kfs_10*Tm_10)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter at Keyway of Gear 4 (D2)
Ma_i = sqrt(22.9906^2 + 8.3679^2);
Mm_i = 0;
Ta_i = 0;
Tm_i = Ta;
Kt_i = 2.7; % Assumed for Sharp Fillet
Kts_i = 2.2; % Assumed for SSharp Fillet
Kf_i = Kt_i;
Kfs_i = Kts_i;

Se_i = ka*kb*kc*ke*kf*Se_p;
d11 = ((16*n/pi)*((2*Kf_i*Ma_10)/Se_i) + sqrt(3*(Kfs_i*Tm_i)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Max Moment on Middle Shaft (D3)
Ma_j = sqrt(94.3171^2 + 34.3286^2);
Mm_j = 0;
Ta_j = 0;
Tm_j = 0;
Kt_j = 2.14; % Assumed for Keyseat
Kts_j = 3; % Assumed for Keyseat
Kf_j = Kt_j;
Kfs_j = Kts_j;

Se_j = ka*kb*kc*ke*kf*Se_p;
d12 = ((16*n/pi)*((2*Kf_j*Ma_j)/Se_j) + sqrt(3*(Kfj*Tm_j)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

% Diameter Design for Shoulder of Gear 3 (D4)
Ma_k = sqrt(47.4175^2 + 17.2586^2);
Mm_k = 0;
Ta_k = 0;
Tm_k = 0;
Ktk = 2.7; % Assumed for Sharp Fillet
Ktsk = 2.2; % Assumed for Sharp Fillet
Kfk = Ktk;
Kfsk = Ktsk;

Se_k = ka*kb*kc*ke*kf*Se_p;
d13 = ((16*n/pi)*((2*Kfk*Ma_k)/Se_k) + sqrt(3*(Kfsk*Tm_k)^2)/Sut)^(1/3); % Intial guess of D1 based on critical point on bearing

%%
D10 = d10.*ones(1,10);
kb_10 = ones(1,10);
qd_10 = ones(1,10);
qs_10 = ones(1,10);
kf_10 = Kf_10.*ones(1,10);
kfs_10 = Kfs_10.*ones(1,10);
Se_10 = Se_10.*ones(1,10);

for i = 2:10
    kb_10(i) = 0.879.*(D10(i-1)).^-0.107;
    qd_10(i) = (1+sqrt_a./0.02.*D10(i-1)).^-1;
    qs_10(i) = (1+sqrt_as./0.02.*D10(i-1)).^-1;
    kf_10(i) = 1 + qd_10(i-1).*(Kt_10-1);
    kfs_10(i) = 1 + qs_10(i-1).*(Kts_10-1);
    Se_10(i) = ka.*kb_10(i-1).*kc.*ke.*kf.*Se_p;
    D10(i) = ((16*n/pi).*(((2*kf_10(i-1).*Ma_10)./Se_10(i-1)) + sqrt(3.*(kfs_10(i-1).*Tm_10).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_10(3) * Ma_10)/(pi*D10(3)^3);
sigma_mp = sqrt(3)*((16*kfs_10(3)*Tm_10)/(pi*D10(3)^3));

ny_D10 = Sy/(sigma_ap + sigma_mp);
nf_D10 = ((sigma_ap/Se_10(3) + (sigma_mp/Sut)))^-1;

fprintf('Output Shaft Yielding FOS @ End of 2in Shaft = %f\n',ny_D10)
fprintf('Output Shaft Fatigue FOS @ End of 2in Shaft = %f\n', nf_D10)
%%
D11 = d11.*ones(1,10);
kb_11 = ones(1,10);
qd_11 = ones(1,10);
qs_11 = ones(1,10);
kf_11 = Kf_i.*ones(1,10);
kfs_11 = Kfs_i.*ones(1,10);
Se_11 = Se_i.*ones(1,10);

for i = 2:10
    kb_11(i) = 0.879.*(D11(i-1)).^-0.107;
    qd_11(i) = (1+sqrt_a./0.02.*D11(i-1)).^-1;
    qs_11(i) = (1+sqrt_as./0.02.*D11(i-1)).^-1;
    kf_11(i) = 1 + qd_11(i-1).*(Kt_i-1);
    kfs_11(i) = 1 + qs_11(i-1).*(Kts_i-1);
    Se_11(i) = ka.*kb_11(i-1).*kc.*ke.*kf.*Se_p;
    D11(i) = ((16*n/pi).*(((2*kf_11(i-1).*Ma_i)./Se_11(i-1)) + sqrt(3.*(kfs_11(i-1).*Tm_i).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_11(8) * Ma_i)/(pi*D11(8)^3);
sigma_mp = sqrt(3)*((16*kfs_11(8)*Tm_i)/(pi*D11(8)^3));

ny_D11 = Sy/(sigma_ap + sigma_mp);
nf_D11 = ((sigma_ap/Se_11(8) + (sigma_mp/Sut)))^-1;

fprintf('Output Shaft Yielding FOS @ Bearing I = %f\n',ny_D11)
fprintf('Output Shaft Fatigue FOS @ Bearing I = %f\n', nf_D11)
%%
D12 = d12.*ones(1,10);
kb_12 = ones(1,10);
qd_12 = ones(1,10);
qs_12 = ones(1,10);
kf_12 = Kf_j.*ones(1,10);
kfs_12 = Kfs_j.*ones(1,10);
Se_12 = Se_j.*ones(1,10);

for i = 2:10
    kb_12(i) = 0.879.*D12(i-1).^-0.107;
    qd_12(i) = (1+sqrt_a./0.02.*D12(i-1)).^-1;
    qs_12(i) = (1+sqrt_as./0.02.*D12(i-1)).^-1;
    kf_12(i) = 1 + qd_12(i-1).*(Kt_j-1);
    kfs_12(i) = 1 + qs_12(i-1).*(Kts_j-1);
    Se_12(i) = ka.*kb_12(i-1).*kc.*ke.*kf.*Se_p;
    D12(i) = ((16*n/pi).*(((2*kf_12(i-1).*Ma_j)./Se_12(i-1)) + sqrt(3.*(kfs_12(i-1).*Tm_j).^2)./Sut)).^(1/3);
end

sigma_ap = (32*kf_12(8) * Ma_j)/(pi*D12(8)^3);
sigma_mp = sqrt(3)*((16*kfs_12(8)*Tm_j)/(pi*D12(8)^3));

ny_D12 = Sy/(sigma_ap + sigma_mp);
nf_D12 = ((sigma_ap/Se_12(8) + (sigma_mp/Sut)))^-1;

fprintf('Output Shaft Yielding FOS @ Gear 5 = %f\n',ny_D12)
fprintf('Output Shaft Fatigue FOS @ Gear 5 = %f\n', nf_D12)
%%
D13 = d13.*ones(1,10);
kb_13 = ones(1,10);
qd_13 = ones(1,10);
qs_13 = ones(1,10);
kf_13 = Kfk.*ones(1,10);
kfs_13 = Kfsk.*ones(1,10);
Sed_13 = Se_k.*ones(1,10);

for i = 2:10
    kb_13(i) = 0.879.*D13(i-1).^-0.107;
    qd_13(i) = (1+sqrt_a./0.02.*D13(i-1)).^-1;
    qs_13(i) = (1+sqrt_as./0.02.*D13(i-1)).^-1;
    kf_13(i) = 1 + qd_13(i-1).*(Kfk-1);
    kfs_13(i) = 1 + qs_13(i-1).*(Kfsk-1);
    Sed_13(i) = ka.*kb_13(i-1).*kc.*ke.*kf.*Se_p;
    D13(i) = ((16*n/pi).*((2*kf_13(i-1).*Ma_k)./Sed_13(i-1)) + sqrt(3.*(kfs_13(i-1).*Tm_k).^2)./Sut).^(1/3);
end
% 
sigma_ap = (32*kf_13(4)*Ma_k)/(pi*D13(4)^3);
sigma_mp = sqrt(3)*((16*kfs_13(4)*Tm_k)/(pi*D13(4)^3));

ny_D13 = Sy/(sigma_ap + sigma_mp);
nf_D13 = ((sigma_ap/Sed_13(4) + (sigma_mp/Sut)))^-1;

fprintf('Output Shaft Yielding FOS @ Bearing K = %f\n',ny_D13)
fprintf('Output Shaft Fatigue FOS @ Bearing K = %f\n', nf_D13)