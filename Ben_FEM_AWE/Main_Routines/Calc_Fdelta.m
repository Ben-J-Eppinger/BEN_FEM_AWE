function [F_delta] = Calc_Fdelta(V,phi,f)
%Calc_Fdelta Cacluates the local load vector on a particular element
%  
%   Inputs: 
%   V:  the array of vertices for the element in question. should be 2x3
%       with x cords in row 1 and z cords in row 2 with each collumn being a
%       node
%   phi: A Data structure with the basis fucntions
%   f:  The forcing function in terms of space f(x,z)
%
%   Outputs: 
%   F_delta: The local load vector in the elemtn in question

% define the jacobiain

B = [V(:,2) - V(:,1) V(:,3) - V(:,1)];
det_B = det(B);

% numerical quadrature

% get the gaussian quadriture points
prec=3;
G=GaussTriangle(prec,0);

% initialize the local stiffness matrix
F_delta = zeros(3,1);

% compute each element of F_delta
for n = 1:3
    
    I=0; % initialize the integral
    
    for i=1:size(G,1) % loop through the gaussian quadrature
        
        cord =  B*(G(i,1:2)') + V(:,1); % find the physical cordinates at each quadrature point
        
        I = I + G(i,3)*... % weight
            f(cord(1),cord(2))*... % value
            phi{n}(G(i,1), G(i,2)); % product
        
    end
    
    I = I*det_B; % scalte by the determinant
    F_delta(n) = I; % put into local stifness matrix
    
end

end

