function [F,M,K] = Assembly(Elements,Nodes,f,c)
%Assembly Assembles the F, M, and F array needed to solve the system
%
%   Inputs:
%   Elements: the elements array returned by GenerateLinearMesh
%   Nodes: the nodes array returned by GenerateLinearMesh
%   f: the force function interms of space f(x,z)
%   c: the squared wave speed as a function of space c(x,z)
%
%   Outputs:
%   F: the load vector
%   M: the mass matrix
%   K: the stiffness matrix

tic

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);

%define the shape funcitons and their derivatives

phi{1} = @(xi,eta) 1 - xi - eta;
phi{2} = @(xi,eta) xi;
phi{3} = @(xi,eta) eta;

dphi = [-1 1 0;
    -1 0 1];

% initialize M, K and F

M = sparse(N,N); 

K = sparse(N,N);

F = zeros(N,1);

% loop over all elements

p_done = 0;
fprintf('Asembling Matrices: %6.0f%%',p_done)

for k = 1:N_el
    
    IV = Elements(:,k);
    V = Nodes(:,IV); 
    F_delta = Calc_Fdelta(V,phi,f);
    K_delta = Calc_Kdelta(V,c,dphi);
    M_delta = Calc_Mdelta(V,phi);
    
    for i = 1:n
        
        F(IV(i)) = F(IV(i)) + F_delta(i);
        
        for j = 1:n
            
            K(IV(i),IV(j)) = K(IV(i),IV(j)) + K_delta(i,j);            
            M(IV(i),IV(j)) = M(IV(i),IV(j)) + M_delta(i,j);
            
        end
    end
    
     % display progress
    %if mod(k,round(N_el/10)) == 0 
    p_done = 100*k/N_el;
    fprintf('\b\b\b\b\b\b\b\b%3.0f%% || ',p_done)
    %end
        
end

toc

end


