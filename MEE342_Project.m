% To run code for plots, run top section, then run section whichever plots
% are needed

% Inputs
Sut = 85*10^3; % [psi] Ultimate Tensile Strength 1040 CD Steel
Sy = 71*10^3; % [psi] Yield Strength 1040 CD Steel
E = 29*10^3; % [psi] Young's Modulus of 1040 CD Steel

% Guesses/Assumptions
Dc = 6.08; Df = 1.33; Dg = 6.08; Dj = 1.33; % Diameters of gears (Guessed)
GR_BF = 16/73; % gear ratio = (Number of teeth of Gear 2)/(Number of teeth of Gear 3)
GR_EJ = 16/73; % gear ratio = (Number of teeth of Gear 4)/(Number of teeth of Gear 5)
theta = pi/9;
theta2 = pi/9;
Px = 50.25*3; % Axial Load From Thrust of Blades [lbs]
Vy = 260*3; % Lateral Load from Weight of Blades [lbs]
Ta = 38.6*3; % Torque Applied by Blades [ft-lbs]

% Torque Calculations
Wt_23 = (Ta*12)/(Dc/2); % Tangential force from gear 2 on 3
Wr_23 = Wt_23 * tan(theta); % Radial force from gear 2 on 3
Wt_45 = (Wt_23 * (Df/2))/(Dg/2); % Tangential force from gear 4 on 5
Wr_45 = Wt_45 * tan(theta2); % Radial force from gear 4 on 5

% Reaction Force Calculation for Input Shaft

L0 = 4; % Length of Shaft out of Gearbox
L1 = 0.75; % Length Between Bearing 1 and Gear 2
L2 = 1.875; % Length of Gear 2
L3 = 1.875; % Length Between Bearing 2 and Gear 2
L4 = 0.625; % Length From Bearing 2 Reaction to end of Shaft
L_Input = L0 + L1 + L2 + L3 + L4;

Tb = Ta; % Torque at Gear B
Ay = (Vy*(L_Input-L4) + Wr_23*(L3))/(L2 + L3); % Reaction at 1st Bearing
Cy = Vy + Wr_23 - Ay; % Reaction at 2nd Bearing
Az = (Wt_23*(L3))/(L2+L3);
Cz = Wt_23 - Az;
Cx = Px;

x1 = linspace(0,L_Input,1075);
In_Shaft = [0, L0, L0+L1, L0+L1+L2, L0+L1+L2+L3/2, L_Input+L4];
Torque1 = [Ta,Ta,Ta-Tb,Ta-Tb];
T1 = zeros(size(x2));
T1(1:(In_Shaft(5)*100)) = Torque1(1);

F1_x = [Px,Px,Px,Cx];
Fx_InShaft1 = zeros(size(x1));
Fx_InShaft1(1:(In_Shaft(6)*100)) = F1_x(1);

Vy_1 = -Vy*(x1>0) + Ay*(x1>(L0+L1)) - Wr_23*(x1>(L0+L1+L2)) + Cy*(x1>(L_Input-L4));
Vz_1 = -Az*(x1>(L0+L1)) + Wt_23*(x1>(L0+L1+L2)) - Cz*(x1>(L_Input-L4));
My_1 = -Az*(x1-(L0+L1)).*(x1>(L0+L1)) + Wt_23*(x1-(L0+L1+L2)).*(x1>(L0+L1+L2)) - Cz*(x1-(L_Input-L4)).*(x1>(L_Input-L4));
Mz_1 = -Vy*(x1-0).*(x1>0) + Ay*(x1-(L0+L1)).*(x1>(L0+L1)) - Wr_23*(x1-(L0+L1+L2)).*(x1>(L0+L1+L2)) + Cy*(x1-(L_Input-L4)).*(x1>(L_Input-L4));;

% Reaction Force Calculations for Intermediate Shaft
L5 = 0.75; % Length Bearing 4 Reaction 
L6 = 2; % Length to Gear 3 Reaction
L7 = 5.75; % Length to Gear 2 Reaction
L8 = 2.25; % Length to Bearing 3 Reaction
L9 = 0.75; % Length from Bearing 3 Reaction to End of Shaft
L_Int = L5 + L6 + L7 + L8 + L9;
Tf = Ta * GR_BF;
Te = Tf;
Gy = (Wr_45 * L6 + Wr_23 * (L6 + L7)) / (L6 + L7 + L8);
Dy = -Gy + Wr_23 + Wr_45;
Gz = (-Wt_45 * L6 + Wt_23 * (L6 + L7)) / (L6 + L7 + L8);
Dz = -Gz + Wt_23 - Wt_45;

% Setting Singularity Functions and Torque as function of x
x2 = linspace(0,L_Int,1075);
Int_Shaft = [0, L5, L5+L6, L5+L6+L7, L_Int-L9, L_Int];
Torque2 = [Tf-Te,Tf-Te,Tf,Tf-Te];
T2 = zeros(size(x2));
T2((Int_Shaft(1)*100+1):(Int_Shaft(3)*100)) = Torque2(2);
T2((Int_Shaft(3)*100+1):(Int_Shaft(4)*100)) = Torque2(3);
T2((Int_Shaft(4)*100+1):(Int_Shaft(5)*100)) = Torque2(4);

Vy_2 = Dy*(x2>L5) - Wr_45*(x2>(L5+L6)) - Wr_23*(x2>(L5+L6+L7)) + Gy*(x2>(L_Int - L9));
Vz_2 = -Dz*(x2>L5) - Wt_45*(x2>(L5+L6)) + Wt_23*(x2>(L5+L6+L7)) - Gz*(x2>(L_Int - L9));
My_2 = -Dz*(x2-(L5)).*(x2>(L5)) - Wt_45*(x2-(L5+L6)).*(x2>(L5+L6)) + Wt_23*(x2-(L5+L6+L7)).*(x2>(L5+L6+L7)) - Gz*(x2-(L_Int-L9)).*(x2>(L_Int-L9));
Mz_2 = Dy*(x2-(L5)).*(x2>(L5)) - Wr_45*(x2-(L5+L6)).*(x2>(L5+L6)) - Wr_23*(x2-(L5+L6+L7)).*(x2>(L5+L6+L7)) + Gy*(x2-(L_Int-L9)).*(x2>(L_Int-L9));

% Recation Force Calculations for Output Shaft
L10 = 2; % Length of Shaft Sticking out of Gearbox
L11 = 0.75; % Length from end of Output Shaft to Bearing 6 Reaction
L12 = 2.25; % Length from Bearing 6 to Gear 5 Reaction
L13 = 1.625; % Length from Gear 5 to Bearing 5 Reaction
L14 = 0.625; % Length from Bearing 5 Reaction to end of Output Shaft
L_Out = L10 + L11 + L12 + L13 + L14;
Tf = Ta * GR_BF;
Te = Tf;
Tj = Te * GR_EJ;
Th = Tj;
Ky = (Wr_45 * L12) / (L12 + L13);
Iy = -Ky + Wr_45;
Kz = (Wt_45 * L12) / (L12 + L13);
Iz = -Kz + Wt_45;

% Setting Singularity Functions and Torque as Functions of x
x3 = linspace(0,L_Out,700);
Out_Shaft = [0, L10, L10+L11, L10+L11+L12, L10+L11+L12+L13, L_Out];
Torque3 = [Th,Th,Th-Tj,Th-Tj];
T3 = zeros(size(x3));
T3((Out_Shaft(1)*100+1):(Out_Shaft(3)*100)) = Torque3(2);
T3((Out_Shaft(4)*100):(Out_Shaft(5)*100)) = Torque3(4);

Vy = Iy*(x3>(L10 + L11)) - Wr_45*(x3>(L10 + L11 + L12)) + Ky*(x3>(L10 + L11 + L12 +L13));
Vz = -Iz*(x3>(L10 + L11)) + Wt_45*(x3>(L10 + L11 + L12)) - Kz*(x3>(L10 + L11 + L12 +L13));
My = -Iz*(x3-(L10+L11)).*(x3>(L10+L11)) + Wt_45*(x3-(L10+L11+L12)).*(x3>(L10+L11+L12)) - Kz*(x3-(L10+L11+L12+L13)).*(x3>(L10+L11+L12+L13));
Mz = Iy*(x3-(L10+L11)).*(x3>(L10+L11)) - Wr_45*(x3-(L10+L11+L12)).*(x3>(L10+L11+L12)) + Ky*(x3-(L10+L11+L12+L13)).*(x3>(L10+L11+L12+L13));

%% Plots for Shaft 1
close all;

figure(1)
plot(x1,T1,"k-","LineWidth",2)
xlabel('Shaft Distance(in)')
ylabel("Torque (lbs-in)")
title('Torque Diagram (Input Shaft)')

figure(2)
plot(x1,Fx_Shaft1,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Axial Force (lbs)')
title('Axial Force Diagram (Input Shaft)')

figure(3)
plot(x1,Vy_1,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Y Diagram (Input Shaft)')

figure(4)
plot(x1,Vz_1,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Z Diagram (Input Shaft)')

figure(5)
plot(x1,My_1,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Y Diagram (Input Shaft)')

figure(6)
plot(x1,Mz_1,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Z Diagram (Input Shaft)')
%% Plots for Shaft 2
close all;

figure(1)
plot(x2,T2,"k-","LineWidth",2)
xlabel('Shaft Distance(in)')
ylabel("Torque (lbs-in)")
title('Torque Diagram (Intermediate Shaft)')

figure(2)
plot(x2,Vy_2,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Y Diagram (Intermediate Shaft)')

figure(3)
plot(x2,Vz_2,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Z Diagram (Intermediate Shaft)')

figure(4)
plot(x2,My_2,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Y Diagram (Intermediate Shaft)')
ylim([-800,100])

figure(5)
plot(x2,Mz_2,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Z Diagram (Intermediate Shaft)')
%% Plots for Shaft 3
close all;
figure(1)
plot(x3,T3,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel("Torque (lbs-in)")
title('Torque Diagram (Output Shaft)')

figure(2)
plot(x3,Vy,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Y Diagram (Output Shaft)')

figure(3)
plot(x3,Vz,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Shear Force (lbs)')
title('Shear Force Z Diagram (Output Shaft)')

figure(4)
plot(x3,My,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Y Diagram (Output Shaft)')

figure(5)
plot(x3,Mz,'k-','LineWidth',2)
xlabel('Distance on Shaft (in)')
ylabel('Moment (lb-in)')
title('Moment in Z Diagram (Output Shaft)')