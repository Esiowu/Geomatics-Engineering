clear all
close all
clc

% Load GCPs Geodetic Coords (Known) model coordinates in mm

% Load Results from Lab3

GCPs_Pt_Model = [
    102, -2.274491853, -5.934950174, -151.6811648;
    105, 87.43594397, -88.14744037, -148.4980106;
    200, 18.2144714, 109.654938, -153.5805529;
    100, -9.474484647, 96.32036807, -153.5457957;
    104, 18.38865314, -79.46286962, -148.9211622;
    201, 43.99032779, 7.374288007, -150.9729891;
    202, -7.476281135, -48.38480639, -151.2056239;
    203, 50.98013036, -90.06114869, -148.2546889
];

% Load Known GCPS

GCPs_Known = [
    102, 109.70, -642.35, 1086.43;
    105, 517.62, -194.43, 1090.65;
    200, -466.39, -542.31, 1091.55
    100, -399.28, -679.72, 1090.96;
    104, 475.55, -538.18, 1090.50;
    201, 42.73, -412.19, 1090.82;
    202, 321.09, -667.45, 1083.49;
    203, 527.78, -375.72, 1092.00
];




% Absolute orientation
% Absolutely orient your model using the three control points:
% Calculate approximate parameter values and estimate the absolute orientation parameters of your model space from the three control points

% Extract Values GCPs
% Image Space
Xm_full = GCPs_Pt_Model (:, 2);
Ym_full = GCPs_Pt_Model (:, 3);
Zm_full = GCPs_Pt_Model(:, 4);

% Object Space
Xo_full = GCPs_Known (:, 2);
Yo_full = GCPs_Known (:, 3);
Zo_full = GCPs_Known (:, 4);


% Now Lets Move to residual calculation (ANGLES IN RAD ,others in mm)

Omega_Final = -0.0257602254186736;
Phi_Final = -0.00256615368895036;
Kappa_Final = -1.57415400701183;
tx_Final = 99.4422567961494;
ty_Final =-629.195935451632 ;
tz_Final = 1842.12314947580;
Lambda_Final =4.97685643192124 ;


% Get DHat Values From RO Least Squares
DHat_Omega = -1.17904282893057e-12;
DHat_Phi = 1.82150287596531e-12;
DHat_Kappa = -1.10572041008367e-14;
DHat_tx = 8.93292106520172e-10;
DHat_ty = -1.42664496658397e-09;
DHat_tz =-3.61610278233854e-10 ;
DHat_Lambda = 3.88386992261676e-15;

						

% Get Object Space Rotation Matrix
M_Object = AO_M3_Matrix(Omega_Final, Phi_Final, Kappa_Final);

    % Access each element of the matrix M
    M11 = M_Object(1,1); % Element in the first row, first column
    M12 = M_Object(1,2); % Element in the first row, second column
    M13 = M_Object(1,3); % Element in the first row, third column
    M21 = M_Object(2,1); % Element in the second row, first column
    M22 = M_Object(2,2); % Element in the second row, second column
    M23 = M_Object(2,3); % Element in the second row, third column
    M31 = M_Object(3,1); % Element in the third row, first column
    M32 = M_Object(3,2); % Element in the third row, second column
    M33 = M_Object(3,3); % Element in the third row, third column



    % For W - Misclosure
    % Initialize W_Final as a 18x1 vector
    W_Final = zeros(3*8, 1); % Preallocate W_Final for 6 observations, each contributing 3 elements

    % Calculate the misclosure vector W_Final for each set of coordinates
    for i = 1:8
        W_Final(3*i-2) = Lambda_Final * (M11 * Xm_full(i) + M12 * Ym_full(i) + M13 * Zm_full(i)) + tx_Final - Xo_full(i);
        W_Final(3*i-1) = Lambda_Final * (M21 * Xm_full(i) + M22 * Ym_full(i) + M23 * Zm_full(i)) + ty_Final - Yo_full(i);
        W_Final(3*i)   = Lambda_Final * (M31 * Xm_full(i) + M32 * Ym_full(i) + M33 * Zm_full(i)) + tz_Final - Zo_full(i);
    end


    % For Derivatives from A matrix_final is a 24 by 7 (eight points)
    % For A- Matrix (Partial Derivatives expressed as determinants)
    % Initialize A_Matrix
    A_Matrix_Final = zeros(24, 7);

    % Initialize the partial derivative output arrays outside the loop
    partial_Xo_Omega_Final = zeros(8, 1);
    partial_Xo_Phi_Final = zeros(8, 1);
    partial_Xo_K_Final = zeros(8, 1);
    partial_Xo_Lambda_Final = zeros(8, 1);
    partial_Xo_tx_Final = zeros(8, 1);
    partial_Xo_ty_Final = zeros(8, 1);
    partial_Xo_tz_Final = zeros(8, 1);

    partial_Yo_Omega_Final = zeros(8, 1);
    partial_Yo_Phi_Final = zeros(8, 1);
    partial_Yo_K_Final = zeros(8, 1);
    partial_Yo_Lambda_Final = zeros(8, 1);
    partial_Yo_tx_Final = zeros(8, 1);
    partial_Yo_ty_Final = zeros(8, 1);
    partial_Yo_tz_Final = zeros(8, 1);

    partial_Zo_Omega_Final = zeros(8, 1);
    partial_Zo_Phi_Final = zeros(8, 1);
    partial_Zo_K_Final = zeros(8, 1);
    partial_Zo_Lambda_Final = zeros(8, 1);
    partial_Zo_tx_Final = zeros(8, 1);
    partial_Zo_ty_Final = zeros(8, 1);
    partial_Zo_tz_Final = zeros(8, 1);


% Initialize A_Matrix_Final as a 24x7 matrix for the A matrix
A_Matrix_Final = zeros(24, 7);

% Fill in the A matrix with the partial derivatives for each point
for i = 1:8
    % Calculate the partial derivatives for X°
    partial_Xo_Omega_Final = Lambda_Final * Ym_full * (-sin(Omega_Final) * sin(Kappa_Final) + cos(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final* Zm_full * (cos(Omega_Final) * sin(Kappa_Final) + sin(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final));   
    partial_Xo_Phi_Final = -Lambda_Final * Xm_full * (sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Ym_full * (sin(Omega_Final) * cos(Phi_Final) * cos(Kappa_Final)) - Lambda_Final * Zm_full * (cos(Omega_Final) * cos(Phi_Final) * cos(Kappa_Final));
    partial_Xo_K_Final = -Lambda_Final * Xm_full * (cos(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Ym_full * ((cos(Omega_Final) * cos(Kappa_Final)) - sin(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * ((sin(Omega_Final) * cos(Kappa_Final) + cos(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)));
    partial_Xo_Lambda_Final = Xm_full * M11 + Ym_full * M12 + Zm_full * M13;
    partial_Xo_tx_Final = [1;1;1;1;1;1;1;1];
    partial_Xo_ty_Final = [0;0;0;0;0;0;0;0];
    partial_Xo_tz_Final = [0;0;0;0;0;0;0;0];

    % Calculate the partial derivatives for Y°
    partial_Yo_Omega_Final = Lambda_Final * Ym_full* (-sin(Omega_Final) * cos(Kappa_Final) - cos(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * (cos(Omega_Final) * cos(Kappa_Final) - sin(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final));
    partial_Yo_Phi_Final = Lambda_Final * Xm_full * (sin(Phi_Final) * sin(Kappa_Final)) - Lambda_Final * Ym_full * (sin(Omega_Final) * cos(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * ((cos(Omega_Final) * cos(Phi_Final) * sin(Kappa_Final)));
    partial_Yo_K_Final = -Lambda_Final * Xm_full* (cos(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Ym_full * (-cos(Omega_Final) * sin(Kappa_Final) - sin(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Zm_full * (-sin(Omega_Final) * sin(Kappa_Final)) + (cos(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final));
    partial_Yo_Lambda_Final = Xm_full * M21 + Ym_full * M22 + Zm_full * M23;
    partial_Yo_tx_Final = [0;0;0;0;0;0;0;0];
    partial_Yo_ty_Final = [1;1;1;1;1;1;1;1];
    partial_Yo_tz_Final = [0;0;0;0;0;0;0;0];


    % Calculate the partial derivatives for Z°
    partial_Zo_Omega_Final = -Lambda_Final * Ym_full * (cos(Omega_Final) * cos(Phi_Final)) - Lambda_Final * Zm_full * (sin(Omega_Final) * cos(Phi_Final));
    partial_Zo_Phi_Final = Lambda_Final * Xm_full * cos(Phi_Final) + Lambda_Final * Ym_full * sin(Omega_Final) * sin(Phi_Final) - Lambda_Final * Zm_full * cos(Omega_Final) * sin(Phi_Final);
    partial_Zo_K_Final =  [0;0;0;0;0;0;0;0];
    partial_Zo_Lambda_Final = Xm_full * M31 + Ym_full * M32 + Zm_full * M33;
    partial_Zo_tx_Final = [0;0;0;0;0;0;0;0];
    partial_Zo_ty_Final = [0;0;0;0;0;0;0;0];
    partial_Zo_tz_Final = [1;1;1;1;1;1;1;1];

    % Construct the 24 x 7 matrix A with each partial derivative occupying three rows per point
    % Partial Derivatives for each point is within the A matrix
    A_Matrix_Final = [partial_Xo_Omega_Final(1),partial_Xo_Phi_Final(1),partial_Xo_K_Final(1),partial_Xo_tx_Final(1), partial_Xo_ty_Final(1),partial_Xo_tz_Final(1),partial_Xo_Lambda_Final(1);
                    partial_Yo_Omega_Final(1),partial_Yo_Phi_Final(1),partial_Yo_K_Final(1), partial_Yo_tx_Final(1),partial_Yo_ty_Final(1),partial_Yo_tz_Final(1), partial_Yo_Lambda_Final(1);
                    partial_Zo_Omega_Final(1),partial_Zo_Phi_Final(1),partial_Zo_K_Final(1),partial_Zo_tx_Final(1), partial_Zo_ty_Final(1), partial_Zo_tz_Final(1), partial_Zo_Lambda_Final(1);
                    partial_Xo_Omega_Final(2),partial_Xo_Phi_Final(2), partial_Xo_K_Final(2),partial_Xo_tx_Final(2), partial_Xo_ty_Final(2),partial_Xo_tz_Final(2),partial_Xo_Lambda_Final(2);
                    partial_Yo_Omega_Final(2),partial_Yo_Phi_Final(2), partial_Yo_K_Final(2), partial_Yo_tx_Final(2),partial_Yo_ty_Final(2),partial_Yo_tz_Final(2), partial_Yo_Lambda_Final(2);
                    partial_Zo_Omega_Final(2),partial_Zo_Phi_Final(2), partial_Zo_K_Final(2),partial_Zo_tx_Final(2), partial_Zo_ty_Final(2), partial_Zo_tz_Final(2), partial_Zo_Lambda_Final(2);
                    partial_Xo_Omega_Final(3),partial_Xo_Phi_Final(3), partial_Xo_K_Final(3),partial_Xo_tx_Final(3),partial_Xo_ty_Final(3),partial_Xo_tz_Final(3),partial_Xo_Lambda_Final(3);
                    partial_Yo_Omega_Final(3),partial_Yo_Phi_Final(3), partial_Yo_K_Final(3),partial_Yo_tx_Final(3),partial_Yo_ty_Final(3),partial_Yo_tz_Final(3),partial_Yo_Lambda_Final(3);
                    partial_Zo_Omega_Final(3),partial_Zo_Phi_Final(3), partial_Zo_K_Final(3),partial_Zo_tx_Final(3),partial_Zo_ty_Final(3),partial_Zo_tz_Final(3),partial_Zo_Lambda_Final(3);
                    partial_Xo_Omega_Final(4),partial_Xo_Phi_Final(4), partial_Xo_K_Final(4),partial_Xo_tx_Final(4),partial_Xo_ty_Final(4),partial_Xo_tz_Final(4),partial_Xo_Lambda_Final(4);
                    partial_Yo_Omega_Final(4),partial_Yo_Phi_Final(4), partial_Yo_K_Final(4),partial_Yo_tx_Final(4),partial_Yo_ty_Final(4),partial_Yo_tz_Final(4),partial_Yo_Lambda_Final(4);
                    partial_Zo_Omega_Final(4),partial_Zo_Phi_Final(4), partial_Zo_K_Final(4),partial_Zo_tx_Final(4),partial_Zo_ty_Final(4),partial_Zo_tz_Final(4),partial_Zo_Lambda_Final(4);
                    partial_Xo_Omega_Final(5),partial_Xo_Phi_Final(5), partial_Xo_K_Final(5),partial_Xo_tx_Final(5),partial_Xo_ty_Final(5),partial_Xo_tz_Final(5),partial_Xo_Lambda_Final(5);
                    partial_Yo_Omega_Final(5),partial_Yo_Phi_Final(5), partial_Yo_K_Final(5),partial_Yo_tx_Final(5),partial_Yo_ty_Final(5),partial_Yo_tz_Final(5),partial_Yo_Lambda_Final(5);
                    partial_Zo_Omega_Final(5),partial_Zo_Phi_Final(5), partial_Zo_K_Final(5),partial_Zo_tx_Final(5),partial_Zo_ty_Final(5),partial_Zo_tz_Final(5),partial_Zo_Lambda_Final(5);
                    partial_Xo_Omega_Final(6),partial_Xo_Phi_Final(6), partial_Xo_K_Final(6),partial_Xo_tx_Final(6),partial_Xo_ty_Final(6),partial_Xo_tz_Final(6),partial_Xo_Lambda_Final(6);
                    partial_Yo_Omega_Final(6),partial_Yo_Phi_Final(6), partial_Yo_K_Final(6),partial_Yo_tx_Final(6),partial_Yo_ty_Final(6),partial_Yo_tz_Final(6),partial_Yo_Lambda_Final(6);
                    partial_Zo_Omega_Final(6),partial_Zo_Phi_Final(6), partial_Zo_K_Final(6),partial_Zo_tx_Final(6),partial_Zo_ty_Final(6),partial_Zo_tz_Final(6),partial_Zo_Lambda_Final(6);
                    partial_Xo_Omega_Final(7),partial_Xo_Phi_Final(7), partial_Xo_K_Final(7),partial_Xo_tx_Final(7),partial_Xo_ty_Final(7),partial_Xo_tz_Final(7),partial_Xo_Lambda_Final(7);
                    partial_Yo_Omega_Final(7),partial_Yo_Phi_Final(7), partial_Yo_K_Final(7),partial_Yo_tx_Final(7),partial_Yo_ty_Final(7),partial_Yo_tz_Final(7),partial_Yo_Lambda_Final(7);
                    partial_Zo_Omega_Final(7),partial_Zo_Phi_Final(7), partial_Zo_K_Final(7),partial_Zo_tx_Final(7),partial_Zo_ty_Final(7),partial_Zo_tz_Final(7),partial_Zo_Lambda_Final(7);
                    partial_Xo_Omega_Final(8),partial_Xo_Phi_Final(8), partial_Xo_K_Final(8),partial_Xo_tx_Final(8),partial_Xo_ty_Final(8),partial_Xo_tz_Final(8),partial_Xo_Lambda_Final(8);
                    partial_Yo_Omega_Final(8),partial_Yo_Phi_Final(8), partial_Yo_K_Final(8),partial_Yo_tx_Final(8),partial_Yo_ty_Final(8),partial_Yo_tz_Final(8),partial_Yo_Lambda_Final(8);
                    partial_Zo_Omega_Final(8),partial_Zo_Phi_Final(8), partial_Zo_K_Final(8),partial_Zo_tx_Final(8),partial_Zo_ty_Final(8),partial_Zo_tz_Final(8),partial_Zo_Lambda_Final(8)];
    end




    % Now lets calculate residuals
    % VHat = ( VxHat, VyHat ,VzHat)'

    % Initiate partials For r vector 3x8 matrices for points
    
    % Define the row ranges for concatenation
    row_ranges = {[1:3], [4:6], [7:9], [10:12], [13:15], [16:18], [19:21], [22:24]};
    
    % Initialize empty matrices for partial derivatives
    partial_ro_partial_omega = [];
    partial_ro_partial_Phi = [];
    partial_ro_partial_kappa = [];
    partial_ro_partial_Lambda = [];
    partial_ro_partial_tx = [];
    partial_ro_partial_ty = [];
    partial_ro_partial_tz = [];
    
    % Concatenate submatrices for each partial derivative
    for i = 1:numel(row_ranges)
        partial_ro_partial_omega = [partial_ro_partial_omega, A_Matrix_Final(row_ranges{i}, 1)];
        partial_ro_partial_Phi = [partial_ro_partial_Phi, A_Matrix_Final(row_ranges{i}, 2)];
        partial_ro_partial_kappa = [partial_ro_partial_kappa, A_Matrix_Final(row_ranges{i}, 3)];
        partial_ro_partial_Lambda = [partial_ro_partial_Lambda, A_Matrix_Final(row_ranges{i}, 7)];
        partial_ro_partial_tx = [partial_ro_partial_tx, A_Matrix_Final(row_ranges{i}, 4)];
        partial_ro_partial_ty = [partial_ro_partial_ty, A_Matrix_Final(row_ranges{i}, 5)];
        partial_ro_partial_tz = [partial_ro_partial_tz, A_Matrix_Final(row_ranges{i}, 6)];
    end


    VHat1 = (W_Final(1:3,1) + (partial_ro_partial_omega(1:3,1) * DHat_Omega) + (partial_ro_partial_Phi(1:3,1) * DHat_Phi) + (partial_ro_partial_kappa(1:3,1) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,1) * DHat_Lambda) + (partial_ro_partial_tx(1:3,1) * DHat_tx) + (partial_ro_partial_ty(1:3,1) * DHat_ty) + (partial_ro_partial_tz(1:3,1) * DHat_tz))'/1000;
    VHat2 = (W_Final(4:6,1) + (partial_ro_partial_omega(1:3,2) * DHat_Omega) + (partial_ro_partial_Phi(1:3,2) * DHat_Phi) + (partial_ro_partial_kappa(1:3,2) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,2) * DHat_Lambda) + (partial_ro_partial_tx(1:3,2) * DHat_tx) + (partial_ro_partial_ty(1:3,2) * DHat_ty) + (partial_ro_partial_tz(1:3,2) * DHat_tz))'/1000;
    VHat3 = (W_Final(7:9,1) + (partial_ro_partial_omega(1:3,3) * DHat_Omega) + (partial_ro_partial_Phi(1:3,3) * DHat_Phi) + (partial_ro_partial_kappa(1:3,3) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,3) * DHat_Lambda) + (partial_ro_partial_tx(1:3,3) * DHat_tx) + (partial_ro_partial_ty(1:3,3) * DHat_ty) + (partial_ro_partial_tz(1:3,3) * DHat_tz))'/1000;
    VHat4 = (W_Final(10:12,1) + (partial_ro_partial_omega(1:3,4) * DHat_Omega) + (partial_ro_partial_Phi(1:3,4) * DHat_Phi) + (partial_ro_partial_kappa(1:3,4) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,4) * DHat_Lambda) + (partial_ro_partial_tx(1:3,4) * DHat_tx) + (partial_ro_partial_ty(1:3,4) * DHat_ty) + (partial_ro_partial_tz(1:3,4) * DHat_tz))'/1000;
    VHat5 = (W_Final(13:15,1) + (partial_ro_partial_omega(1:3,5) * DHat_Omega) + (partial_ro_partial_Phi(1:3,5) * DHat_Phi) + (partial_ro_partial_kappa(1:3,5) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,5) * DHat_Lambda) + (partial_ro_partial_tx(1:3,5) * DHat_tx) + (partial_ro_partial_ty(1:3,5) * DHat_ty) + (partial_ro_partial_tz(1:3,5) * DHat_tz))'/1000;
    VHat6 = (W_Final(16:18,1) + (partial_ro_partial_omega(1:3,6) * DHat_Omega) + (partial_ro_partial_Phi(1:3,6) * DHat_Phi) + (partial_ro_partial_kappa(1:3,6) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,6) * DHat_Lambda) + (partial_ro_partial_tx(1:3,6) * DHat_tx) + (partial_ro_partial_ty(1:3,6) * DHat_ty) + (partial_ro_partial_tz(1:3,6) * DHat_tz))'/1000;
    VHat7 = (W_Final(19:21,1) + (partial_ro_partial_omega(1:3,7) * DHat_Omega) + (partial_ro_partial_Phi(1:3,7) * DHat_Phi) + (partial_ro_partial_kappa(1:3,7) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,7) * DHat_Lambda) + (partial_ro_partial_tx(1:3,7) * DHat_tx) + (partial_ro_partial_ty(1:3,7) * DHat_ty) + (partial_ro_partial_tz(1:3,7) * DHat_tz))'/1000;
    VHat8 = (W_Final(22:24,1) + (partial_ro_partial_omega(1:3,8) * DHat_Omega) + (partial_ro_partial_Phi(1:3,8) * DHat_Phi) + (partial_ro_partial_kappa(1:3,8) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,8) * DHat_Lambda) + (partial_ro_partial_tx(1:3,8) * DHat_tx) + (partial_ro_partial_ty(1:3,8) * DHat_ty) + (partial_ro_partial_tz(1:3,8) * DHat_tz))'/1000;


    VHat = [VHat1;VHat2;VHat3;VHat4;VHat5;VHat6;VHat7;VHat8];
    residuals_x = VHat(:,1);
    residuals_y = VHat(:,2);
    residuals_z = VHat(:,3);

    % Calculate RMSE for x ,y and z
    rmse_x = sqrt(mean(residuals_x.^2));
    rmse_y = sqrt(mean(residuals_y.^2));
    rmse_z = sqrt(mean(residuals_z.^2));
    
    % Display RMSE values
    disp(['RMSE for x: ', num2str(rmse_x)]);
    disp(['RMSE for y: ', num2str(rmse_y)]);
    disp(['RMSE for z: ', num2str(rmse_z)]);
    disp(' '); % Empty line for separation
    disp(' '); % Empty line for separation
