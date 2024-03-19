

%  Inputs:
%   - GCPs_Known:       Object coordinates of ground control points (3 x N matrix)
%   - Transformed_GCPsPoints: Transformed coordinates of ground control points (K x 3 matrix)
%   - Check_Pt_Known:   Real coordinates of check points (K x 3 matrix)
%   - Transformed_CheckPoints: Transformed coordinates of check points (K x 3 matrix)
%   - baseVector:       Base vector (3 x 1 column vector)
%   - c:                Constant value
%
%   Outputs:
%   - differences:      Differences between real and calculated checkpoints (K x 3 matrix)
%   - accuracies:       Accuracy values for X, Y, and Z coordinates (3 x 1 column vector)
%
%   Explanation:
%   This function calculates the differences between real and calculated checkpoints and
%   computes the accuracies of X, Y, and Z based on the provided parameters.
%
%   Author: [Simisola Oyetunde]
%   Date: [2024-03-17]
%   Version: [04]


function [differences, accuracies] = LAB4AO_Accuracy(GCPs_Known, Transformed_GCPsPoints, Check_Pt_Known, Transformed_CheckPoints, baseVector, c)
    % Calculate differences between real and calculated checkpoints
    differences = Check_Pt_Known(:, 2:4) - Transformed_CheckPoints(:, 1:3);

    % Define image point precision
    imagePointPrecision = 0.004E-3;

    % Calculate mean heights
    H = mean(GCPs_Known (3, :));
    h = mean(Transformed_GCPsPoints(:, 3));

    % Calculate scale factor
    S = (H - h) / (c * 10^-3);

    % Calculate base vector magnitude and BH
    baseMagnitude = sqrt(baseVector(1, 1)^2 + baseVector(2, 1)^2 + baseVector(3, 1)^2);
    BH = baseMagnitude / (H - h);

    % Calculate accuracies
    accuracyX = S * imagePointPrecision;
    accuracyY = accuracyX;
    accuracyZ = ((sqrt(2) * S) / BH) * imagePointPrecision;

    % Store accuracies in an array
    accuracies = [accuracyX; accuracyY; accuracyZ];
end
