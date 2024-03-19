
% ENGO 431
% Principles of Photogrammetry
% Laboratory Assignment 4
% Absolute Orientation
% Function to define the rotation matrix MR based on initial estimates of rotation angles Omega, Phi, Kappa

function M = AO_M2_Matrix(Omega, Phi, Kappa)
    % Calculate the sine and cosine of the angles
    c_phi = cos(Phi);
    s_phi = sin(Phi);
    c_omega = cos(Omega);
    s_omega = sin(Omega);
    c_k = cos(Kappa);
    s_k = sin(Kappa);

    % Construct the rotation matrix
    M = [c_phi * c_k, c_omega * s_k + s_omega * s_phi * c_k, s_omega * s_k - c_omega * s_phi * c_k;
         -c_phi * s_k, c_omega * c_k - s_omega * s_phi * s_k, s_omega * c_k + c_omega * s_phi * s_k;
         s_phi, -s_omega * c_phi, c_omega * c_phi];
end