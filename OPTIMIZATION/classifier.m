function L = classifier(Y,G,y,P,template_GC,p,idx,indPx)
% G, y is data from single entity that has K=1
% template_vect is off-diagonal VAR coefficient that learned from the estimation process
[n,~,K] = size(Y);
Atmp = repmat(template_GC,1,1,p,K);
% Atmp = Atmp(indPx);
template_vect = Atmp(idx);
template_vect = template_vect(indPx);
% p = length(template_vect)/(n^2-n)/K;
x = refit(G,y,P,template_vect,1);
A = devect(x,n,p,K);
L = FITTING(Y,A,n,p,K);

end
function LogLikeLihood = FITTING(data,A,n,p,K)
T = size(data,2);
LogLikeLihood = 0;
tmpA = reshape(A,[n,n*p,K]);
for kk=1:K
    [H,Y] = H_gen(data(:,:,kk),p);
    Ek = Y - tmpA(:,:,kk)*H; % Error term
    Sigma = Ek*Ek'/(T-p);
%     disp(size(Sigma))
    LogLikeLihood = LogLikeLihood+(-1/2)*(T-p)*log(det(Sigma))+(-1/2)*n*(T-p);
end

end