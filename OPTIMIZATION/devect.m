function A = devect(x,n,p,k)
tmp = reshape(x,[p*k,n*n]);
tmpA(k).A = zeros(n,n);
idk = 0;
for ii = 1:k
    tmpA2 = [];
    for jj=1:p
        idk = idk+1;
        tmpA2 = [tmpA2 (reshape(tmp(idk,:)',[n,n]))'];
%         tmpA(ii).A = [tmpA(ii).A (reshape(tmp(idk,:)',[n,n]))'];
    end
    tmpA(ii).A = tmpA2;
end
A = reshape([tmpA.A],[n,n,p,k]);
end