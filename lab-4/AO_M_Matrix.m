
% ENGO 431
% Principles of Photogrammetry
% Laboratory Assignment 4
% Absolute Orientation
% Function to define the rotation matrix MR based on initial estimates of rotation angles Omega, Phi, Kappa

function M = AO_M_Matrix(OmegaO, PhiO, KappaO)
    % Calculate the sine and cosine of the angles
    c_phi = cos(PhiO);
    s_phi = sin(PhiO);
    c_omega = cos(OmegaO);
    s_omega = sin(OmegaO);
    c_k = cos(KappaO);
    s_k = sin(KappaO);

    % Construct the rotation matrix
    M = [c_phi * c_k, c_omega * s_k + s_omega * s_phi * c_k, s_omega * s_k - c_omega * s_phi * c_k;
         -c_phi * s_k, c_omega * c_k - s_omega * s_phi * s_k, s_omega * c_k + c_omega * s_phi * s_k;
         s_phi, -s_omega * c_phi, c_omega * c_phi];
end