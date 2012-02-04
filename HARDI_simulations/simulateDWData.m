function S = simulateDWData (bVal, gradientDirections, orientations, weights, isAnisotropic)

if isAnisotropic
    L = diag([1.7e-3, 3e-4, 3e-4]);
else
    L = diag([7e-4, 7e-4, 7e-4]);
end

S = zeros(length(gradientDirections), 1);

for i=1:length(gradientDirections)
    for j=1:length(orientations)
        angle = orientations(j);
        R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
        D = R'*L*R;
        S(i) = S(i) + weights(j)*exp(-bVal*(gradientDirections(i,:)*D*gradientDirections(i,:)'));
    end
end

end