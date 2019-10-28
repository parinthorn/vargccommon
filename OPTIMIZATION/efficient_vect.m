function idx = efficient_vect(PARAMETER)
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
A = zeros(n,n,p,K);

A(:) = 1:1:n^2*p*K;


ind = linindex(n,n,p,'row');
x = zeros(p*n^2*K,1);
for kk=1:K
    A0 = A(:,:,:,kk);
    Ek = zeros(K,1); Ek(kk) = 1;
    x0 = A0(ind); % size pn^2 x 1
    x0 = reshape(x0,p,n^2); % size p x n^2
    x0 = kron(x0,Ek'); % insert zero columns
    x0 = reshape(x0,p*n^2*K,1);
    x = x+x0;
end
idx = x;

end

