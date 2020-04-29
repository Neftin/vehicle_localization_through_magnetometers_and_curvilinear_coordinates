function [ c , ceq ] = constraints(z,auxdata,f)
 
%   c = ...   % Compute nonlinear inequalities at x.
% ceq = ...   % Compute nonlinear equalities at x.
    
    c   = [];
    ceq = discrete_diff_constraints(z,auxdata,f);
    
    
    
end