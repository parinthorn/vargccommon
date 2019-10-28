
function [x, x2,history] = grouplasso(G, b, P,lambda, PARAMETER,rho,Num,varargin)
% group_lassooff  Solve group lasso problem via ADMM
%
% [x, history] = grouplasso(G, b, P,lambda, PARAMETER, rho);
%
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Gx - b ||_2^2 + \lambda sum(norm(P*x))
%
%
% P is a weight matrix that maps x to the penalized entries weighted by
% p_ij > 1 
% 
% PARAMETER = [n p K] where n is the number of variables and p is the order
% of AR model and K is the number of models 
% 
% if K > 1 then this program estimates K models jointly 
% (they all have common Granger graph)
% 
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% BLOCK_SIZE_SUM2NORM is an integer indicating the block size when computing the sum of norm
%
% varargin is 'initial condition' for x (optional)


t_start = tic;

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

PRINT_RESULT = 1;
FREQ_PRINT = 100;
MAXITERS = 1e9; 
ABSTOL = 1e-6;
RELTOL = 1e-4;

% store variables
[ng,nn] = size(G);
np = size(P,1);

% nn = n^2*p; np = (n^2-n)*p;
Gtb = G'*b;
Pt = P';

L = chol(sparse(G'*G+rho*Pt*P) ,'lower'); 
L = sparse(L); U = L';

%% ADMM solver

optargin = size(varargin,2);

if optargin == 0,
    x1 = zeros(nn,1); 
else
    x1 = varargin{1};
end
x2 = P*x1;
z = zeros(np,1);

if PRINT_RESULT
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for k = 1:MAXITERS
ta = tic;
    % x1-update
    q = Gtb + Pt*(rho*x2-z);    % temporary value   
    x1 = U \ (L \ q); % x1 is not generally sparse 

    % x2-update
    x2old = x2;
    Px1 = P*x1;
%     x2 = prox_sumof2norm(Px1+z/rho,p,lambda/rho);
      x2 = prox_sumof2norm(Px1+z/rho,p*K,lambda/rho); % joint estimation of K models
      
    % z-udpate
    z = z + rho*(Px1 - x2);

    % stopping criterion
    
%     obj = 0.5*norm(G*x1-b)^2+lambda*norm21(x2,p);
    obj = 0.5*norm(G*x1-b)^2+lambda*norm21(x2,p*K);

    history.objval(k) = obj;  
    history.r_norm(k) = norm(Px1-x2);
    history.s_norm(k)  = norm(rho*Pt*(x2-x2old));    
      
    history.eps_pri(k) = sqrt(np)*ABSTOL + RELTOL*max(norm(Px1), norm(x2));
    history.eps_dual(k)= sqrt(nn)*ABSTOL + RELTOL*norm(Pt*z);
%     history.x{k} = x2;
history.tpi(k,1) = toc(ta);

    if (PRINT_RESULT && mod(k,FREQ_PRINT) == 0)
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
   history.fit = 0.5*norm(G*x1-b)^2;
         break;
    end

end

toc(t_start);


% return sparse result
% since we penalize the l1 norm on Px then x2 = Px is sparse 
% (where P maps x to off-diagonal entries) but x1 is not
% generally sparse. We can extract the off-diagonal entries from x2.
% 

if nn ~= n^2*p*K,  
    K = nn/(n^2*p); % this is actual K
end
X = reshape(x1,p*K,n^2); % x1 is not sparse
X2 = reshape(x2,p*K,n^2-n); % x2 is offdiagonal entries only and sparse
IND_offdiag = setdiff((1:n^2)',(1:n+1:n^2)','rows');
X(:,IND_offdiag) = X2;
x = X(:);

end



