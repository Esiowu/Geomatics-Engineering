
% ENGO 431
% Principles of Photogrammetry
% Laboratory Assignment 4
% Absolute Orientation

clear all
close all
clc
format long
% Load Results from Lab3

GCPs_Pt_Model = [
    102, -2.274491853, -5.934950174, -151.6811648;
    105, 87.43594397, -88.14744037, -148.4980106;
    200, 18.2144714, 109.654938, -153.5805529
];

Check_Pt_Model = [
    100, -9.474484647, 96.32036807, -153.5457957;
    104, 18.38865314, -79.46286962, -148.9211622;
    201, 43.99032779, 7.374288007, -150.9729891;
    202, -7.476281135, -48.38480639, -151.2056239;
    203, 50.98013036, -90.06114869, -148.2546889
];

Tie_Pt_Model = [
    1, 8.95457365032259,85.640745706275, -152.885158235898;
    2, 90.5053940751437, 75.4705347580149, -151.494508969955;
    3, 5.16111207234995, 10.421240189128, -150.552935570933;
    4, 86.3843785451478, 4.28528767281028, -153.66614677389;
    5, 12.5158036392059, -48.7386130699065, -149.914506458617;
    6, 89.0118176296829, -49.9323655074747, -148.784309390745
];

% Load Known GCPS

GCPs_Known = [
    102, 109.70, -642.35, 1086.43;
    105, 517.62, -194.43, 1090.65;
    200, -466.39, -542.31, 1091.55
    
];

Check_Pt_Known = [
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
Xm = GCPs_Pt_Model (:, 2);
Ym = GCPs_Pt_Model (:, 3);
Zm = GCPs_Pt_Model(:, 4);

% Object Space
Xo = GCPs_Known (:, 2);
Yo = GCPs_Known (:, 3);
Zo = GCPs_Known (:, 4);

% Extract Values Check Points
% Image Space
Xmc = Check_Pt_Model (:, 2);
Ymc = Check_Pt_Model (:, 3);
Zmc = Check_Pt_Model (:, 4);

% Object Space
Xoc = Check_Pt_Known (:, 2);
Yoc = Check_Pt_Known (:, 3);
Zoc = Check_Pt_Known (:, 4);



% To get Initial values of Parameters

% Choose Furthest Points
% Calculate pairwise distances
numPoints = size(GCPs_Pt_Model, 1);
distances = zeros(numPoints, numPoints);
for i = 1:numPoints
    for j = i+1:numPoints
        distances(i, j) = sqrt(sum((GCPs_Pt_Model(i, 2:4) - GCPs_Pt_Model(j, 2:4)).^2));
    end
end

% Find the maximum distance
[maxDistance, maxIndices] = max(distances(:));
[rowIdx, colIdx] = ind2sub(size(distances), maxIndices);
furthestPoints = [GCPs_Pt_Model(rowIdx, 1); GCPs_Pt_Model(colIdx, 1)];

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

KappaOdeg = KappaO*(180/pi);

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
XHat = [X_final(1)*180/pi, X_final(2)*180/pi, X_final(3)*180/pi, X_final(4), X_final(5), X_final(6), X_final(7)];
disp('Final estimated parameters (Omega, Phi, Kappa, tx, ty, tz, Lambda):');
disp(XHat);




% After the iterative process is completed, compute the covariance matrix
Cx = inv(A_Matrix' * A_Matrix);

% Compute standard deviations from the diagonal elements of the covariance matrix
std_devs = sqrt(diag(Cx));

% Initialize the correlation coefficient matrix
correlation_coefficient_matrix = zeros(length(std_devs));

% Fill in the correlation coefficient matrix
for i = 1:length(std_devs)
    for j = 1:length(std_devs)
        correlation_coefficient_matrix(i, j) = Cx(i, j) / (std_devs(i) * std_devs(j));
    end
end

% Display the correlation coefficient matrix
disp('Correlation Coefficient Matrix of the Parameters:');
disp(correlation_coefficient_matrix);

% Write final parameters (correlation_coefficient_matrix) to an Excel file
writematrix(correlation_coefficient_matrix, 'AOcorrelation_coefficient_matrix.xlsx', 'Sheet', 1);

% Calculate redundancy number for x, y, and z of GCP points

% A, P, CL are already defined matrices
% Calculate A'
A = A_Matrix;
A_transpose = transpose(A_Matrix);
I = eye(size(A, 1));
% Calculate the product A' * P * A
product_APA = A_transpose * P * A;

% Calculate the inverse of the matrix APA
inverse_APA = pinv(product_APA);

% % Calculate the product A * (inverse_APA) * A'
% product_APA_A = A * inverse_APA * A_transpose;

% Calculate the inverse of matrix C
inverse_CL = inv(CL);

% Calculate the final result R
R = I - A * inverse_APA * A_transpose * inverse_CL;

% Calculate redundancy numbers
redundancy_numbers = diag(R);

% Display the redundancy numbers
disp('Redundancy numbers:');
disp(redundancy_numbers);

% Now Lets Move to Similarity Transformation of Points ( Using the parameters in Xhat from LS above)
% Extract Values GCPs
% Image Space
Xm_full = Xm;
Ym_full = Ym;
Zm_full = Zm;

% Object Space
Xo_full = Xo;
Yo_full = Yo;
Zo_full = Zo;



% Now lets Transform all points (ANGLES IN RAD ,others in mm)
% X_Transformed = Lambda_Final * M_Object * ([Xm_full(1,:);Ym_full(1,:);Zm_full(1,:)]) + [tx_Final;ty_Final;tz_Final],
Omega_Final = X_o(:,1);
Phi_Final = X_o(:,2);
Kappa_Final = X_o(:,3);
tx_Final = X_o(:,4);
ty_Final = X_o(:,5);
tz_Final = X_o(:,6);
Lambda_Final =X_o(:,7);



% Get Object Space Rotation Matrix
M_Object = AO_M3_Matrix(Omega_Final, Phi_Final, Kappa_Final);

% Preallocate the transformed coordinates matrix
Coords_Transformed = zeros(3, 3); % 3 rows for X, Y, Z and 6 columns for 6 points

% Apply the transformation to each of the 3 GCP points
for i = 1:3
    % Extract the i-th point's coordinates
    point_coordinates = [Xm_full(i, :); Ym_full(i, :); Zm_full(i, :)];

    % Apply the transformation
    Coords_Transformed(:, i) = Lambda_Final * M_Object * point_coordinates + [tx_Final; ty_Final; tz_Final];
end

Transformed_GCPsPoints = Coords_Transformed';

% Write final parameters (correlation_coefficient_matrix) to an Excel file
writematrix(Transformed_GCPsPoints, 'AOTransformed_GCPsPoints.xlsx', 'Sheet', 1);

% Apply the transformation to each of the 3 Check Points 

for i = 1:5
    % Extract the i-th point's coordinates
    Check_Coords = [Xmc(i, :); Ymc(i, :); Zmc(i,:)];

    % Apply the transformation
    Check_Coords_Transformed(:, i) = Lambda_Final * M_Object * Check_Coords + [tx_Final; ty_Final; tz_Final];

end

Transformed_CheckPoints = Check_Coords_Transformed';

% Write final parameters (correlation_coefficient_matrix) to an Excel file
writematrix(Transformed_CheckPoints, 'AOTransformed_CheckPoints.xlsx', 'Sheet', 1);


% Apply the transformation to each of the 6 Tie Points

XmT = Tie_Pt_Model (:, 2);
YmT = Tie_Pt_Model (:, 3);
ZmT = Tie_Pt_Model(:, 4);

% Preallocate the transformed coordinates matrix
Tie_Coords_Transformed= zeros(3, 6); % 6 rows for X, Y, Z and 6 columns for 6 points

% Apply the transformation to each of the 3 Check Points 

for i = 1:6
    % Extract the i-th point's coordinates
    Tie_Coords = [XmT(i, :); YmT(i, :); ZmT(i,:)];

    % Apply the transformation
    Tie_Coords_Transformed(:, i) = Lambda_Final * M_Object * Tie_Coords + [tx_Final; ty_Final; tz_Final];

end

Transformed_TiePoints = Tie_Coords_Transformed';
% Write final parameters (correlation_coefficient_matrix) to an Excel file
writematrix(Transformed_TiePoints, 'AOTransformed_TiePoints.xlsx', 'Sheet', 1);



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
% From Ro in lab 2
bx = 92;
by = -1.10554050327635;
bz = -1.28785487369518;
Omega_Ro = -0.0187144073515171;
Phi_Ro = 0.00540798410465166;
Kappa_Ro = -0.0293520731606641;

M_LeftRO = [1, 0, 0; 
            0, 1, 0; 
            0, 0, 1];

M_RightRO = [0.999554642015508, -0.0294438767965824, -0.00485548601454187;
             0.0293474295005143, 0.999391254816230, -0.0188639385747797;
             0.00540795774410976, 0.0187130413356140, 0.999810270039776];

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





% Calculate Accuracies

% From Ro in lab 2

% Given data
baseVector = [92; -1.10554050327635;-1.28785487369518];
c = 153.358; % mm

% Call function for  Accuracy 
[differences, accuracies] = LAB4AO_Accuracy(GCPs_Known, Transformed_GCPsPoints, Check_Pt_Known, Transformed_CheckPoints, baseVector, c);

% Display the results
disp('Differences between real and calculated checkpoints:');
disp(differences);
disp('Accuracies of the transformed points (in m):');
disp(accuracies);



