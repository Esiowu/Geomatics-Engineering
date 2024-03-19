clear all
close all
clc

% Load GCPs Geodetic Coords (Known) model coordinates in mm

% Read the contents of the text file into a table
dataTable1 = readtable('TESTDATA1.txt');

% Convert the table to a numeric matrix
ModelMatrix = table2array(dataTable1);

% Display the matrix
disp(ModelMatrix);

% Load GCPs Model space  Coords (Test Data) object coordinates in m

% Read the contents of the text file into a table
dataTable2 = readtable('TESTDATA2.txt');

% Convert the table to a numeric matrix
ObjectMatrix = table2array(dataTable2);

% Display the matrix
disp(ObjectMatrix);

% Convert values from meters to millimeters
ObjectMatrix(:, 2:4) = ObjectMatrix(:, 2:4) * 1000;

% Display the modified matrix
disp(ObjectMatrix);

% Extract Values GCPs
% Image Space
Xm_full = ModelMatrix (:, 2);
Ym_full = ModelMatrix (:, 3);
Zm_full = ModelMatrix (:, 4);

% Object Space
Xo_full = ObjectMatrix (:, 2);
Yo_full = ObjectMatrix (:, 3);
Zo_full = ObjectMatrix (:, 4);


% Now Lets Move to residual calculation (ANGLES IN RAD ,others in mm)

DHat_Omega = -1.72240441697721e-15;
DHat_Phi = -4.10967731402799e-15;
DHat_Kappa = -3.74256122379315e-16;
DHat_tx = 3.99152853001091e-09;
DHat_ty = 1.83678533530704e-08;
DHat_tz =8.39027022511067e-10 ;
DHat_Lambda = 3.33795785089366e-13;


% Calculate misclosure for each point (ANGLES IN RAD ,others in mm)

Omega_Final = -0.014410375447377;
Phi_Final = -0.012846372769068;
Kappa_Final = 0.329698792362283;
tx_Final = 6.349935469653413e+06;
ty_Final =3.964457991333090e+06 ;
tz_Final = 1.458257928987210e+06;
Lambda_Final =7.585354580619940e+03 ;



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
    W_Final = zeros(3*6, 1); % Preallocate W_Final for 6 observations, each contributing 3 elements

    % Calculate the misclosure vector W_Final for each set of coordinates
    for i = 1:6
        W_Final(3*i-2) = Lambda_Final * (M11 * Xm_full(i) + M12 * Ym_full(i) + M13 * Zm_full(i)) + tx_Final - Xo_full(i);
        W_Final(3*i-1) = Lambda_Final * (M21 * Xm_full(i) + M22 * Ym_full(i) + M23 * Zm_full(i)) + ty_Final - Yo_full(i);
        W_Final(3*i)   = Lambda_Final * (M31 * Xm_full(i) + M32 * Ym_full(i) + M33 * Zm_full(i)) + tz_Final - Zo_full(i);
    end


    % For Derivatives from A matrix_final is a 24 by 7 (six points)
    % For A- Matrix (Partial Derivatives expressed as determinants)
    % Initialize A_Matrix
    A_Matrix_Final = zeros(18, 7);

    % Initialize the partial derivative output arrays outside the loop
    partial_Xo_Omega_Final = zeros(6, 1);
    partial_Xo_Phi_Final = zeros(6, 1);
    partial_Xo_K_Final = zeros(6, 1);
    partial_Xo_Lambda_Final = zeros(6, 1);
    partial_Xo_tx_Final = zeros(6, 1);
    partial_Xo_ty_Final = zeros(6, 1);
    partial_Xo_tz_Final = zeros(6, 1);

    partial_Yo_Omega_Final = zeros(6, 1);
    partial_Yo_Phi_Final = zeros(6, 1);
    partial_Yo_K_Final = zeros(6, 1);
    partial_Yo_Lambda_Final = zeros(6, 1);
    partial_Yo_tx_Final = zeros(6, 1);
    partial_Yo_ty_Final = zeros(6, 1);
    partial_Yo_tz_Final = zeros(6, 1);

    partial_Zo_Omega_Final = zeros(6, 1);
    partial_Zo_Phi_Final = zeros(6, 1);
    partial_Zo_K_Final = zeros(6, 1);
    partial_Zo_Lambda_Final = zeros(6, 1);
    partial_Zo_tx_Final = zeros(6, 1);
    partial_Zo_ty_Final = zeros(6, 1);
    partial_Zo_tz_Final = zeros(6, 1);


% Initialize A_Matrix_Final as a 18x7 matrix for the A matrix
A_Matrix_Final = zeros(18, 7);

% Fill in the A matrix with the partial derivatives for each point
for i = 1:6
    % Calculate the partial derivatives for X°
    partial_Xo_Omega_Final = Lambda_Final * Ym_full * (-sin(Omega_Final) * sin(Kappa_Final) + cos(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final* Zm_full * (cos(Omega_Final) * sin(Kappa_Final) + sin(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final));   
    partial_Xo_Phi_Final = -Lambda_Final * Xm_full * (sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Ym_full * (sin(Omega_Final) * cos(Phi_Final) * cos(Kappa_Final)) - Lambda_Final * Zm_full * (cos(Omega_Final) * cos(Phi_Final) * cos(Kappa_Final));
    partial_Xo_K_Final = -Lambda_Final * Xm_full * (cos(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Ym_full * ((cos(Omega_Final) * cos(Kappa_Final)) - sin(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * ((sin(Omega_Final) * cos(Kappa_Final) + cos(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)));
    partial_Xo_Lambda_Final = Xm_full * M11 + Ym_full * M12 + Zm_full * M13;
    partial_Xo_tx_Final = [1;1;1;1;1;1];
    partial_Xo_ty_Final = [0;0;0;0;0;0];
    partial_Xo_tz_Final = [0;0;0;0;0;0];

    % Calculate the partial derivatives for Y°
    partial_Yo_Omega_Final = Lambda_Final * Ym_full* (-sin(Omega_Final) * cos(Kappa_Final) - cos(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * (cos(Omega_Final) * cos(Kappa_Final) - sin(Omega_Final) * sin(Phi_Final) * sin(Kappa_Final));
    partial_Yo_Phi_Final = Lambda_Final * Xm_full * (sin(Phi_Final) * sin(Kappa_Final)) - Lambda_Final * Ym_full * (sin(Omega_Final) * cos(Phi_Final) * sin(Kappa_Final)) + Lambda_Final * Zm_full * ((cos(Omega_Final) * cos(Phi_Final) * sin(Kappa_Final)));
    partial_Yo_K_Final = -Lambda_Final * Xm_full* (cos(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Ym_full * (-cos(Omega_Final) * sin(Kappa_Final) - sin(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final)) + Lambda_Final * Zm_full * (-sin(Omega_Final) * sin(Kappa_Final)) + (cos(Omega_Final) * sin(Phi_Final) * cos(Kappa_Final));
    partial_Yo_Lambda_Final = Xm_full * M21 + Ym_full * M22 + Zm_full * M23;
    partial_Yo_tx_Final = [0;0;0;0;0;0];
    partial_Yo_ty_Final = [1;1;1;1;1;1];
    partial_Yo_tz_Final = [0;0;0;0;0;0];


    % Calculate the partial derivatives for Z°
    partial_Zo_Omega_Final = -Lambda_Final * Ym_full * (cos(Omega_Final) * cos(Phi_Final)) - Lambda_Final * Zm_full * (sin(Omega_Final) * cos(Phi_Final));
    partial_Zo_Phi_Final = Lambda_Final * Xm_full * cos(Phi_Final) + Lambda_Final * Ym_full * sin(Omega_Final) * sin(Phi_Final) - Lambda_Final * Zm_full * cos(Omega_Final) * sin(Phi_Final);
    partial_Zo_K_Final =  [0;0;0;0;0;0];
    partial_Zo_Lambda_Final = Xm_full * M31 + Ym_full * M32 + Zm_full * M33;
    partial_Zo_tx_Final = [0;0;0;0;0;0];
    partial_Zo_ty_Final = [0;0;0;0;0;0];
    partial_Zo_tz_Final = [1;1;1;1;1;1];
 
    % Construct the 18 x 7 matrix A with each partial derivative occupying three rows per point
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
                    partial_Zo_Omega_Final(6),partial_Zo_Phi_Final(6), partial_Zo_K_Final(6),partial_Zo_tx_Final(6),partial_Zo_ty_Final(6),partial_Zo_tz_Final(6),partial_Zo_Lambda_Final(6)];
    end




    % Now lets calculate residuals
    % VHat = ( VxHat, VyHat ,VzHat)'

    % Initiate partials For r vector 6x3 matrices for point 1
    partial_ro_partial_omega =[A_Matrix_Final(1:3,1),A_Matrix_Final(4:6,1),A_Matrix_Final(7:9,1),A_Matrix_Final(10:12,1),A_Matrix_Final(13:15,1),A_Matrix_Final(16:18,1)];
    partial_ro_partial_Phi =[A_Matrix_Final(1:3,2),A_Matrix_Final(4:6,2),A_Matrix_Final(7:9,2),A_Matrix_Final(10:12,2),A_Matrix_Final(13:15,2),A_Matrix_Final(16:18,2)];
    partial_ro_partial_kappa =[A_Matrix_Final(1:3,3),A_Matrix_Final(4:6,3),A_Matrix_Final(7:9,3),A_Matrix_Final(10:12,3),A_Matrix_Final(13:15,3),A_Matrix_Final(16:18,3)];
    partial_ro_partial_Lambda =[A_Matrix_Final(1:3,7),A_Matrix_Final(4:6,7),A_Matrix_Final(7:9,7),A_Matrix_Final(10:12,7),A_Matrix_Final(13:15,7),A_Matrix_Final(16:18,7)];
    partial_ro_partial_tx =[A_Matrix_Final(1:3,4),A_Matrix_Final(4:6,4),A_Matrix_Final(7:9,4),A_Matrix_Final(10:12,4),A_Matrix_Final(13:15,4),A_Matrix_Final(16:18,4)];
    partial_ro_partial_ty =[A_Matrix_Final(1:3,5),A_Matrix_Final(4:6,5),A_Matrix_Final(7:9,5),A_Matrix_Final(10:12,5),A_Matrix_Final(13:15,5),A_Matrix_Final(16:18,5)];
    partial_ro_partial_tz =[A_Matrix_Final(1:3,6),A_Matrix_Final(4:6,6),A_Matrix_Final(7:9,6),A_Matrix_Final(10:12,6),A_Matrix_Final(13:15,6),A_Matrix_Final(16:18,6)];


    VHat1 = (W_Final(1:3,1) + (partial_ro_partial_omega(1:3,1) * DHat_Omega) + (partial_ro_partial_Phi(1:3,1) * DHat_Phi) + (partial_ro_partial_kappa(1:3,1) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,1) * DHat_Lambda) + (partial_ro_partial_tx(1:3,1) * DHat_tx) + (partial_ro_partial_ty(1:3,1) * DHat_ty) + (partial_ro_partial_tz(1:3,1) * DHat_tz))'/1000;
    VHat2 = (W_Final(4:6,1) + (partial_ro_partial_omega(1:3,2) * DHat_Omega) + (partial_ro_partial_Phi(1:3,2) * DHat_Phi) + (partial_ro_partial_kappa(1:3,2) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,2) * DHat_Lambda) + (partial_ro_partial_tx(1:3,2) * DHat_tx) + (partial_ro_partial_ty(1:3,2) * DHat_ty) + (partial_ro_partial_tz(1:3,2) * DHat_tz))'/1000;
    VHat3 = (W_Final(7:9,1) + (partial_ro_partial_omega(1:3,3) * DHat_Omega) + (partial_ro_partial_Phi(1:3,3) * DHat_Phi) + (partial_ro_partial_kappa(1:3,3) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,3) * DHat_Lambda) + (partial_ro_partial_tx(1:3,3) * DHat_tx) + (partial_ro_partial_ty(1:3,3) * DHat_ty) + (partial_ro_partial_tz(1:3,3) * DHat_tz))'/1000;
    VHat4 = (W_Final(10:12,1) + (partial_ro_partial_omega(1:3,4) * DHat_Omega) + (partial_ro_partial_Phi(1:3,4) * DHat_Phi) + (partial_ro_partial_kappa(1:3,4) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,4) * DHat_Lambda) + (partial_ro_partial_tx(1:3,4) * DHat_tx) + (partial_ro_partial_ty(1:3,4) * DHat_ty) + (partial_ro_partial_tz(1:3,4) * DHat_tz))'/1000;
    VHat5 = (W_Final(13:15,1) + (partial_ro_partial_omega(1:3,5) * DHat_Omega) + (partial_ro_partial_Phi(1:3,5) * DHat_Phi) + (partial_ro_partial_kappa(1:3,5) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,5) * DHat_Lambda) + (partial_ro_partial_tx(1:3,5) * DHat_tx) + (partial_ro_partial_ty(1:3,5) * DHat_ty) + (partial_ro_partial_tz(1:3,5) * DHat_tz))'/1000;
    VHat6 = (W_Final(16:18,1) + (partial_ro_partial_omega(1:3,6) * DHat_Omega) + (partial_ro_partial_Phi(1:3,6) * DHat_Phi) + (partial_ro_partial_kappa(1:3,6) * DHat_Kappa) + (partial_ro_partial_Lambda(1:3,6) * DHat_Lambda) + (partial_ro_partial_tx(1:3,6) * DHat_tx) + (partial_ro_partial_ty(1:3,6) * DHat_ty) + (partial_ro_partial_tz(1:3,6) * DHat_tz))'/1000;

    VHat = [VHat1;VHat2;VHat3;VHat4;VHat5;VHat6];
