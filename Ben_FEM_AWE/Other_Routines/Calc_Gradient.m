function [d_x_u, d_z_u] = Calc_Gradient(Elements, Nodes, elements_per_node, C)
%Calc_Gradient Summary of this function goes here
%   Detailed explanation goes here

n = size(Elements,1); 
N_el = size(Elements,2); 
N = size(Nodes,2); 

dphi = [-1 1 0;
    -1 0 1];

% initialize components of the the gradient
d_x_u = zeros(N,1);
d_z_u = zeros(N,1); 

% loop over elements
for k = 1:N_el
    
    % grab the node locations
    IV = Elements(:,k);
    V = Nodes(:,IV); 
    % compute B and its inverse
    B = [V(:,2) - V(:,1) V(:,3) - V(:,1)];
    inv_B = inv(B);
    % compute the gradient on the element
    grad_u_delta = sum((inv_B')*dphi.*C(IV)',2);
    
    for i = 1:n
        
        % populate the components by averaging togethere the gradients
        % calculaed on all the elements touching a given node
        fact = elements_per_node(IV(i)); 
        d_x_u(IV(i)) = d_x_u(IV(i)) + grad_u_delta(1)/fact;
        d_z_u(IV(i)) = d_z_u(IV(i)) + grad_u_delta(2)/fact;
        
    end
   
end

end

