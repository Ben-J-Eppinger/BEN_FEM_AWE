function [M_tilda_delta] = Calc_M_tilda_delta(V,c,h,phi,G)
%Calc_M_tilda_delta Cacluates the local absorbing matrix on a particular 
%element
%  
%   Inputs: 
%   V:  the array of vertices for the element in question. should be 2x3
%       with x cords in row 1 and z cords in row 2 with each collumn being a
%       node
%   c:  squared wave speed function c(x,z)
%   h: length of physical edge
%   phi: A Data structure with the basis fucntions on the reference edge
%   G: array with Gauss points and weights for 1D quadrature
%
%   Outputs: 
%   M_tilda_delta: The local load absorbing matrix on the elemtn in question

% initialize the local mass matrix
M_tilda_delta = zeros(2,2);

% compute each elemtn of m_tilda_delta
for n = 1:2
    for m = 1:2
        
        I=0; % initialize the integral
        
        for i=1:size(G,1) % loop through the gaussian quadrature
            
            cord = V(:,1) + h*G(i,1)*V(:,2); 
            
            I = I + G(i,2)*... % weight
                sqrt(c(cord(1), cord(2)))*... % c
                phi{n}(G(i,1)) * phi{m}(G(i,1)); % product
       
        end
        
        I = I*h; % scalte by the determinant 
        M_tilda_delta(n,m) = I; % put into local stifness matrix 
        
    end
end



end

