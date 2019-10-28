function [lambda_max] = lambdamax_grouplasso(G,b,PARAMETER,varargin)
% [lambda_max] = lambdamax_grouplass(G,b,[n p K],P)
% 
% lambda_max calculates the maximum value of penalty parameter in the
% problem
%
%   minimize 1/2*|| Gx - b ||_2^2 + \lambda sum(norm(P*x))
%
%
% P = kron(P1,I_{pK}) is a weight matrix that maps x to the penalized entries weighted by
% p_ij > 1 
% 
% PARAMETER = [n p K] where n is the number of variables and p is the order
% of AR model and K is the number of models 
% 
% if K > 1 then the program gives the lambda1 max for the estimation of common GC graph
% 
% By default, if the option of giving P is neglected then P is the projection matrix
% maps all entries to off-diagonal entries (zero rows of P are removed).
% 

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

IND_DIAG = (1:n+1:n^2)';% indices of diagonal elements
IND_OFFDIAG = setdiff((1:n^2)',IND_DIAG,'rows'); % indices of off-diagonal

optargin = size(varargin,2);

if optargin == 0,
    IND_DIAG = 1:n+1:n^2; 
    P1 = speye(n^2); 
    P1(IND_DIAG,:) = []; % select only off-diagonal entries
else
    P1 = varargin{1};
end


Q1 = speye(n^2,n^2); 
Q1(IND_OFFDIAG,:) = 0; % Q is the projection matrix on the diagonal entries
Q = kron(Q1,speye(p*K));

[IND3_OFFDIAG,J] = find(Q);
Gc = G(:,IND3_OFFDIAG');
diagz = zeros(n^2*p*K,1); diagz(IND3_OFFDIAG) = Gc\b;

vecP = max((P1),[],1);
ind_z_vecP = find(any((P1),1) == 0); % extract the zero entries in P
vecP(ind_z_vecP) = 1; % replace zero by 1 (to be the divider)

c = G'*b-G'*G*diagz;
C = reshape(c,p*K,n^2);

lambda_max = max( abs(sqrt(sum(C.^2,1)) )./vecP);





