
% ENGO 431
% Principles of Photogrammetry
% Laboratory Assignment 4
% Absolute Orientation TEST DATA VALIDATION

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
Xm = ModelMatrix ([1, 2, 5], 2);
Ym = ModelMatrix ([1, 2, 5], 3);
Zm = ModelMatrix ([1, 2, 5], 4);

% Object Space
Xo = ObjectMatrix ([1, 2, 5], 2);
Yo = ObjectMatrix ([1, 2, 5], 3);
Zo = ObjectMatrix ([1, 2, 5], 4);

% Extract Values Check Points
% Image Space
Xmc = ModelMatrix ([3, 4, 6], 2);
Ymc = ModelMatrix ([3, 4, 6], 3);
Zmc = ModelMatrix ([3, 4, 6], 4);

% Object Space
Xoc = ObjectMatrix ([3, 4, 6], 2);
Yoc = ObjectMatrix ([3, 4, 6], 3);
Zoc = ObjectMatrix ([3, 4, 6], 4);


% To get Initial values of Parameters

% Choose Furthest Points
% Calculate pairwise distances
numPoints = size(ObjectMatrix, 1);
distances = zeros(numPoints, numPoints);
for i = 1:numPoints
    for j = i+1:numPoints
        distances(i, j) = sqrt(sum((ObjectMatrix(i, 2:4) - ObjectMatrix(j, 2:4)).^2));
    end
end

% Find the maximum distance
[maxDistance, maxIndices] = max(distances(:));
[rowIdx, colIdx] = ind2sub(size(distances), maxIndices);
furthestPoints = [ObjectMatrix(rowIdx, 1); ObjectMatrix(colIdx, 1)];

% Display the furthest points and the distance between them
disp('Furthest Points:');
disp(furthestPoints);
disp(['Maximum Distance: ', num2str(maxDistance), ' mm']);



% Absolute orientation
% Q1 - Calculate approximate parameter values

% Step 1: Estimate Initial Values of Omega,Phi,Kappa,Lambda,tx,ty,tz
% Initial Values of Omega and Phi = 0 Assuming straight Line Flight
OmegaO = 0;
PhiO = 0;

% Initial Values of Kappa
% Choosing two GCPs Furthest Apart
% KO = AlphaOij - Alphami

% AlphaOij 
Numerator1 = Xo(3,:) - Xo(1,:);
Denominator1 = Yo(3,:) - Yo(1,:);
AlphaOij = atan2(Numerator1, Denominator1);

% Alphamij
Numerator2 = Xm(3,:) - Xm(1,:);
Denominator2 = Ym(3,:) - Ym(1,:);
Alphamij = atan2(Numerator2, Denominator2);

KappaO = (AlphaOij - Alphamij);

% KappaOdeg = KappaO*(180/pi);

% Initial Value of LambdaO = Doij/Dmij

% For Doij
ao = Xo(1) - Xo(3);
bo = Yo(1) - Yo(3);
co = Zo(1) - Zo(3);
Doij = sqrt(ao^2 + bo^2 + co^2);

% For Dmij
am = Xm(1) - Xm(3);
bm = Ym(1) - Ym(3);
cm = Zm(1) - Zm(3);
Dmij = sqrt(am^2 + bm^2 + cm^2);

% Now solve for initial LambdaO
LambdaO = Doij / Dmij;


% Initial Value of Translation (txo,tyo,tzo)
% Using one GCPs to determine approximate translation values

% Lets get initial Mo = R3(KappaO) Recall Initial Values of Omega and Phi = 0 Assuming straight Line Flight

Mo = AO_M_Matrix(OmegaO, PhiO, KappaO);

rio =[Xo(1,:);Yo(1,:);Zo(1,:)];
rim =[Xm(1,:);Ym(1,:);Zm(1,:)];

to = rio - LambdaO * Mo *rim;
txo = to(1,:);
tyo = to(2,:);
tzo = to(3,:);


% Initiate Xo
X_o = [OmegaO,PhiO,KappaO,txo,tyo,tzo,LambdaO];

% Display the initial values parameters (OmegaO,PhiO,KappaO,txo,tyo,tzo,LambdaO)
disp('Initial values parameters (OmegaO,PhiO,KappaO,txo,tyo,tzo,LambdaO)');
disp(X_o);


% One iteration least squares
% Set convergence thresholds
angle_convergence_threshold = 1e-6;
translation_convergence_threshold = 1e-6;
scale_convergence_threshold = 1e-6;
max_iterations = 10;

% Initialize iteration counter
iteration = 3;

% Perform iterative least squares adjustment
while iteration < max_iterations

    % Store the current parameter vector for comparison
    Xo_prev = X_o;
    Omega = Xo_prev(1);
    Phi = Xo_prev(2);
    Kappa = Xo_prev(3);
    tx = Xo_prev(4);
    ty = Xo_prev(5);
    tz = Xo_prev(6);
    Lambda = Xo_prev(7);
    
    % Perform one iteration of the least squares algorithm to update Dhat
    % Get Least Squares Parameters Required
    
    % For P ( Weight Matrix)
    % For all points CL (Variance/Covariance Matrix)
    CL = eye(9);
    
    % Calculating for P, the observation weight matrix
    P = inv(CL); 
    
    % For M Rotation matrix
    M = AO_M2_Matrix(Omega, Phi, Kappa);
    
    
    % Access each element of the matrix M
    m11 = M(1,1); % Element in the first row, first column
    m12 = M(1,2); % Element in the first row, second column
    m13 = M(1,3); % Element in the first row, third column
    m21 = M(2,1); % Element in the second row, first column
    m22 = M(2,2); % Element in the second row, second column
    m23 = M(2,3); % Element in the second row, third column
    m31 = M(3,1); % Element in the third row, first column
    m32 = M(3,2); % Element in the third row, second column
    m33 = M(3,3); % Element in the third row, third column
    
    % For W - Misclosure
    % Initialize a cell array to store the misclosure for each tie point
    Misclosure = cell(3, 1);

    % Initialize the W matrix
    W = zeros(9, 1);

    % Fill in the W matrix with the given formula
    for i = 1:3
        W(3*i-2) = Lambda * (m11 * Xm(i) + m12 * Ym(i) + m13 * Zm(i)) + tx - Xo(i);
        W(3*i-1) = Lambda * (m21 * Xm(i) + m22 * Ym(i) + m23 * Zm(i)) + ty - Yo(i);
        W(3*i)   = Lambda * (m31 * Xm(i) + m32 * Ym(i) + m33 * Zm(i)) + tz - Zo(i);
    end
    

    
    % For A- Matrix (Partial Derivatives expressed as determinants)
    % Initialize A_Matrix
    A_Matrix = zeros(9, 7);
    
    for i = 1:3
        % Initialize the partial_Xo output arrays
        partial_Xo_Omega = zeros(3, 1);
        partial_Xo_Phi = zeros(3, 1);
        partial_Xo_K = zeros(3, 1);
        partial_Xo_Lambda = zeros(3, 1);
        partial_Xo_tx = zeros(3, 1);
        partial_Xo_ty = zeros(3, 1);
        partial_Xo_tz = zeros(3, 1);
    
        % Initialize the partial_Yo output arrays
        partial_Yo_Omega = zeros(3, 1);
        partial_Yo_Phi = zeros(3, 1);
        partial_Yo_K = zeros(3, 1);
        partial_Yo_Lambda = zeros(3, 1);
        partial_Yo_tx = zeros(3, 1);
        partial_Yo_ty = zeros(3, 1);
        partial_Yo_tz = zeros(3, 1);
    
        % Initialize the partial_Zo output arrays
        partial_Zo_Omega = zeros(3, 1);
        partial_Zo_Phi = zeros(3, 1);
        partial_Zo_K = zeros(3, 1);
        partial_Zo_Lambda = zeros(3, 1);
        partial_Zo_tx = zeros(3, 1);
        partial_Zo_ty = zeros(3, 1);
        partial_Zo_tz = zeros(3, 1);
    
    
        % Define the partial derivatives for X° as given in page 27 of Lecture Notes
        partial_Xo_Omega = Lambda * Ym * (-sin(Omega)*sin(Kappa) + cos(Omega)*sin(Phi)*cos(Kappa)) + Lambda * Zm * (cos(Omega)*sin(Kappa) + sin(Omega)*sin(Phi)*cos(Kappa));
        partial_Xo_Phi = -Lambda* Xm *(sin(Phi)*cos(Kappa)) + Lambda* Ym *(sin(Omega)*cos(Phi)*cos(Kappa)) - Lambda*Zm *(cos(Omega)*cos(Phi)*cos(Kappa));
        partial_Xo_K = -Lambda* Xm *(cos(Phi)*sin(Kappa)) + Lambda*Ym *((cos(Omega)*cos(Kappa)) - sin(Omega)*sin(Phi)*sin(Kappa)) + Lambda* Zm *((sin(Omega)*cos(Kappa)+ cos(Omega)*sin(Phi)*sin(Kappa)));
        partial_Xo_Lambda = Xm * m11 + Ym*m12 + Zm*m13;
        partial_Xo_tx = [1;1;1];
        partial_Xo_ty = [0;0;0];
        partial_Xo_tz = [0;0;0];
    
    
        % Define the partial derivatives for Y° as given in page 28 of Lecture Notes
        partial_Yo_Omega = (Lambda*Ym*(-sin(Omega)*cos(Kappa) - cos(Omega)*sin(Phi)*sin(Kappa))) + (Lambda*Zm*(cos(Omega)*cos(Kappa) - sin(Omega)*sin(Phi)*sin(Kappa)));
        partial_Yo_Phi = (Lambda*Xm*sin(Phi)*sin(Kappa)) - (Lambda*Ym*(sin(Omega)*cos(Phi)*sin(Kappa))) + (Lambda*Zm*(cos(Omega)*cos(Phi)*sin(Kappa)));
        partial_Yo_K = (-Lambda*Xm*(cos(Phi)*cos(Kappa))) + (Lambda*Ym*(-cos(Omega)*sin(Kappa) - sin(Omega)*sin(Phi)*cos(Kappa))) + (Lambda*Zm*(-sin(Omega)*sin(Kappa)) + (cos(Omega)*sin(Phi)*cos(Kappa)));
        partial_Yo_Lambda = Xm*m21 + Ym*m22 + Zm*m23;
        partial_Yo_tx = [0;0;0];
        partial_Yo_ty = [1;1;1];
        partial_Yo_tz = [0;0;0];
    
    
        % Define the partial derivatives for Z° as given in page 29 of Lecture Notes
        partial_Zo_Omega = (-Lambda*Ym*(cos(Omega)*cos(Phi)) - Lambda*Zm*(sin(Omega)*cos(Phi)));
        partial_Zo_Phi = (Lambda*Xm*cos(Phi))+ (Lambda*Ym*sin(Omega)*sin(Phi))- (Lambda*Zm*cos(Omega)*sin(Phi));
        partial_Zo_K = [0;0;0];
        partial_Zo_Lambda = Xm*m31 + Ym*m32 + Zm*m33;
        partial_Zo_tx = [0;0;0];
        partial_Zo_ty = [0;0;0];
        partial_Zo_tz = [1;1;1];
    
    
        % Construct the 9 x 7 matrix A with each partial derivative occupying three rows per point
        A_Matrix = [partial_Xo_Omega(1),partial_Xo_Phi(1),partial_Xo_K(1),partial_Xo_tx(1), partial_Xo_ty(1),partial_Xo_tz(1),partial_Xo_Lambda(1);
                    partial_Yo_Omega(1),partial_Yo_Phi(1),partial_Yo_K(1),partial_Yo_tx(1),partial_Yo_ty(1),partial_Yo_tz(1),partial_Yo_Lambda(1);
                    partial_Zo_Omega(1),partial_Zo_Phi(1),partial_Zo_K(1),partial_Zo_tx(1), partial_Zo_ty(1),partial_Zo_tz(1),partial_Zo_Lambda(1);
                    partial_Xo_Omega(2),partial_Xo_Phi(2),partial_Xo_K(2),partial_Xo_tx(2), partial_Xo_ty(2),partial_Xo_tz(2),partial_Xo_Lambda(2);
                    partial_Yo_Omega(2),partial_Yo_Phi(2),partial_Yo_K(2),partial_Yo_tx(2),partial_Yo_ty(2),partial_Yo_tz(2),partial_Yo_Lambda(2);
                    partial_Zo_Omega(2),partial_Zo_Phi(2),partial_Zo_K(2),partial_Zo_tx(2), partial_Zo_ty(2),partial_Zo_tz(2),partial_Zo_Lambda(2);
                    partial_Xo_Omega(3),partial_Xo_Phi(3),partial_Xo_K(3),partial_Xo_tx(3),partial_Xo_ty(3),partial_Xo_tz(3),partial_Xo_Lambda(3);
                    partial_Yo_Omega(3),partial_Yo_Phi(3),partial_Yo_K(3),partial_Yo_tx(3),partial_Yo_ty(3),partial_Yo_tz(3),partial_Yo_Lambda(3);
                    partial_Zo_Omega(3),partial_Zo_Phi(3),partial_Zo_K(3),partial_Zo_tx(3),partial_Zo_ty(3),partial_Zo_tz(3),partial_Zo_Lambda(3)];
    end
    
    % Calculate for DeltaHat (solve for Dhat = [DhatOmega,DhatPhi,DhatKappa,Dhattx,Dhatty,Dhattz,DhatLambda)')
    N = -inv(A_Matrix' * P * A_Matrix);
    U = A_Matrix' * P * W;
    DHat = (N * U)';
    
  
    % Update the parameter vector
    X_o = Xo_prev + DHat;


    % Calculate change in parameters
    delta_angles = norm(X_o(1:3) - Xo_prev(1:3));
    delta_translations = norm(X_o(4:6) - Xo_prev(4:6));
    delta_scale = abs(X_o(7) - Xo_prev(7));
    
    % Check for convergence
    if delta_angles < angle_convergence_threshold && ...
            delta_translations < translation_convergence_threshold && ...
            delta_scale < scale_convergence_threshold
        disp('Converged');
        break;
    end
        % Increment iteration counter
        iteration = iteration + 1;
    end

% Check if maximum iterations reached without convergence
if iteration == max_iterations
    disp('Maximum iterations reached without convergence');
end

% Display the final estimated parameters
X_final = X_o;
XHat = [X_final(1)*180/pi, X_final(2)*180/pi, X_final(3)*180/pi, X_final(4)/1000, X_final(5)/1000, X_final(6)/1000, X_final(7)/1000];
disp('Final estimated parameters (Omega, Phi, Kappa, tx, ty, tz, Lambda):');
disp(XHat);













% Now Lets Move to Similarity Transformation of Points ( Using the parameters in Xhat from LS above)
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

DHat_Omega = -1.559232253770598e-15;
DHat_Phi = -3.802807221619803e-15;
DHat_Kappa = 4.796392450088256e-15;
DHat_tx = -9.932125616708791e-07;
DHat_ty =-7.977939767209832e-08;
DHat_tz =6.927680425736572e-09 ;
DHat_Lambda = 3.656734557831265e-11;


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




% X_Transformed = Lambda_Final * M_Object * ([Xm_full(1,:);Ym_full(1,:);Zm_full(1,:)]) + [tx_Final;ty_Final;tz_Final],

% Preallocate the transformed coordinates matrix
Coords_Transformed = zeros(3, 6); % 3 rows for X, Y, Z and 6 columns for 6 points

% Apply the transformation to each of the 6 points
for i = 1:6
    % Extract the i-th point's coordinates
    point_coordinates = [Xm_full(i, :); Ym_full(i, :); Zm_full(i, :)];

    % Apply the transformation
    Coords_Transformed(:, i) = Lambda_Final * M_Object * point_coordinates + [tx_Final; ty_Final; tz_Final];
end

Transformed_ObjectPoints = Coords_Transformed'/1000 ;

% Write final parameters (Transformed_ObjectPoints) to an Excel file
writematrix(Transformed_ObjectPoints, 'Transformed_TestObjectPoints.xlsx', 'Sheet', 1);


% TEST RESIDUALS
Residuals = Transformed_ObjectPoints - (ObjectMatrix(:,2:4))/1000;

% Calculate Vx
Vx = Transformed_ObjectPoints(:, 1) - ObjectMatrix(:, 2) / 1000;
Vy = Transformed_ObjectPoints(:, 2) - ObjectMatrix(:, 3) / 1000;
Vz = Transformed_ObjectPoints(:, 3) - ObjectMatrix(:, 4) / 1000;

% Now lets Calculate RMSE

% Calculate RMSE for Vx and Vy
% rmse_x = sqrt(mean(Vx^2));
% rmse_y = sqrt(mean(Vy^2));
% rmse_Z = sqrt(mean(Vz^2));
% 
% % test
% ZZ = rmse_y - rmse_x;



% Transform Left and Right PC

Pc_Left_Xm = 0;
Pc_Left_Ym = 0;
Pc_Left_Zm = 0;

PC_left_Transformed = Lambda_Final * M_Object * ([Pc_Left_Xm;Pc_Left_Ym;Pc_Left_Zm]) + [tx_Final;ty_Final;tz_Final];
PC_left_Transformed_m = PC_left_Transformed /1000;



Pc_Right_Xm = 92;
Pc_Right_Ym = 5.0455;
Pc_Right_Zm = 2.1725;

PC_Right_Transformed = Lambda_Final * M_Object * ([Pc_Right_Xm;Pc_Right_Ym;Pc_Right_Zm]) + [tx_Final;ty_Final;tz_Final];

PC_Right_Transformed_m = PC_Right_Transformed/1000;




% Lets Get Resultant Matrix Products


M_LeftRO = [1, 0, 0; 
            0, 1, 0; 
            0, 0, 1];


M_RightRO = [0.9981, 0.0553, -0.0259; 
             -0.0551, 0.9984, 0.0091; 
             0.0263, -0.0077, 0.9996];



% For Left MIO
% Moi = Mim * Mom
Mo_i_LEFT = M';


% For Right MIO
% Moi = Mim * Mom
Mo_i_RIGHT =M_RightRO * M';


% Extract Angles

% LEFT
% Access each element of the matrix Mo_i_LEFT
    m11_left = Mo_i_LEFT(1,1); % Element in the first row, first column
    m12_left = Mo_i_LEFT(1,2); % Element in the first row, second column
    m13_left = Mo_i_LEFT(1,3); % Element in the first row, third column
    m21_left = Mo_i_LEFT(2,1); % Element in the second row, first column
    m22_left = Mo_i_LEFT(2,2); % Element in the second row, second column
    m23_left = Mo_i_LEFT(2,3); % Element in the second row, third column
    m31_left = Mo_i_LEFT(3,1); % Element in the third row, first column
    m32_left = Mo_i_LEFT(3,2); % Element in the third row, second column
    m33_left = Mo_i_LEFT(3,3); % Element in the third row, third column

% Omega Left (angles in dd)
w_Left = (atan (-m32_left/m33_left))*(180/pi);

% Phi Left (angles in dd)
Phi_Left = (asin(m31_left))*(180/pi);

% Kappa Left (angles in dd)
k_left = (atan(-m21_left/m11_left)) *(180/pi);



% RIGHT

% Access each element of the matrix Mo_i_RIGHT
    m11_Right = Mo_i_RIGHT(1,1); % Element in the first row, first column
    m12_Right = Mo_i_RIGHT(1,2); % Element in the first row, second column
    m13_Right = Mo_i_RIGHT(1,3); % Element in the first row, third column
    m21_Right = Mo_i_RIGHT(2,1); % Element in the second row, first column
    m22_Right = Mo_i_RIGHT(2,2); % Element in the second row, second column
    m23_Right = Mo_i_RIGHT(2,3); % Element in the second row, third column
    m31_Right = Mo_i_RIGHT(3,1); % Element in the third row, first column
    m32_Right = Mo_i_RIGHT(3,2); % Element in the third row, second column
    m33_Right = Mo_i_RIGHT(3,3); % Element in the third row, third column


% Omega Right (angles in dd)
w_Right = (atan (-m32_Right/m33_Right))*(180/pi);

% Phi Right (angles in dd)
Phi_Right= (asin(m31_Right))*(180/pi);

% Kappa Right (angles in dd)
k_Right = (atan(-m21_Right/m11_Right)) *(180/pi);





