function indPx = efficient_Px(P,n,p,K)
ind = (1:1:n^2*p*K)';
indPx = P*ind;
end