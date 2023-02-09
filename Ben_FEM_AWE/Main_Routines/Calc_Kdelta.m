function [K_delta] = Calc_Kdelta(V,c,dphi)
%Calc_Kdelta Cacluates the local Stiffness matrix on a particular element
%  
%   Inputs: 
%   V:  the array of vertices for the element in question. should be 2x3
%       with x cords in row 1 and z cords in row 2 with each collumn being a
%       node
%   c:  The squared wave speed function in terms of space c(x,z)
%   dphi: an array where each collumn in the gradient of a shape function
%
%   Outputs: 
%   K_delta: The local load stiffness matrix on the element in question

% define the jacobiain

B = [V(:,2) - V(:,1) V(:,3) - V(:,1)];
det_B = det(B);
inv_B = inv(B);

% numerical quadrature

% get the gaussian quadriture points
prec=3;
G=GaussTriangle(prec,0);

% initialize the local stiffness matrix
K_delta = zeros(3,3);

% compute each elemtn of k_delta
for n = 1:3
    for m = 1:3
        
        I=0; % initialize the integral
        
        for i=1:size(G,1) % loop through the gaussian quadrature
            
            cord =  B*(G(i,1:2)') + V(:,1); % find the physical cordinates at each quadrature point
            
            I = I + G(i,3)*... % weight
                c(cord(1), cord(2))*... % c
                dot((inv_B'*dphi(:,n)),(inv_B'*dphi(:,m))); % product
       
        end
        
        I = I*det_B; % scalte by the determinant 
        K_delta(n,m) = I; % put into local stifness matrix 
        
    end
end

end

