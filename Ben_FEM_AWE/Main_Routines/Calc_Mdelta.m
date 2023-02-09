function [M_delta] = Calc_Mdelta(V,phi)
%Calc_Mdelta Cacluates the local mass matrix on a particular element
%  
%   Inputs: 
%   V:  the array of vertices for the element in question. should be 2x3
%       with x cords in row 1 and z cords in row 2 with each collumn being a
%       node
%   phi: A Data structure with the basis fucntions
%
%   Outputs: 
%   M_delta: The local load vector in the elemtn in question

% define the jacobiain

B = [V(:,2) - V(:,1) V(:,3) - V(:,1)];
det_B = det(B);

% numerical quadrature

% get the gaussian quadriture points
prec=3;
G=GaussTriangle(prec,0);

% initialize the local mass matrix
M_delta = zeros(3,3);

% compute each elemtn of m_delta
for n = 1:3
    for m = 1:3
        
        I=0; % initialize the integral
        
        for i=1:size(G,1) % loop through the gaussian quadrature
            
            I = I + G(i,3)*... % weight
                phi{n}(G(i,1), G(i,2)) * phi{m}(G(i,1), G(i,2)); % product
       
        end
        
        I = I*det_B; % scalte by the determinant 
        M_delta(n,m) = I; % put into local stifness matrix 
        
    end
end

end

