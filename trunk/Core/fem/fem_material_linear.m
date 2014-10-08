function [ D ] = fem_material_linear( E, nu )
% FEM_MATERIAL_LINEAR creates the linear elasticity
%   matrix given a young's modulus E and poisson ration nu
    
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));
    
    D = mu*eye(6);
    D(1:3,1:3) = lambda+2*mu*eye(3);
    
end

