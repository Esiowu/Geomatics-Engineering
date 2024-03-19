
% ENGO 431
% Principles of Photogrammetry
% Laboratory Assignment 4
% Absolute Orientation
% Function to define the rotation matrix MR based on initial estimates of rotation angles Omega, Phi, Kappa

function M = AO_M3_Matrix(Omega_Final, Phi_Final, Kappa_Final)
    % Calculate the sine and cosine of the angles
    c_phi = cos(Phi_Final);
    s_phi = sin(Phi_Final);
    c_omega = cos(Omega_Final);
    s_omega = sin(Omega_Final);
    c_k = cos(Kappa_Final);
    s_k = sin(Kappa_Final);

    % Construct the rotation matrix
    M = [c_phi * c_k, c_omega * s_k + s_omega * s_phi * c_k, s_omega * s_k - c_omega * s_phi * c_k;
         -c_phi * s_k, c_omega * c_k - s_omega * s_phi * s_k, s_omega * c_k + c_omega * s_phi * s_k;
         s_phi, -s_omega * c_phi, c_omega * c_phi];
end