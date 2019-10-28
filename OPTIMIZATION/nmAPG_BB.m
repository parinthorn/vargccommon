function [x,history]= nmAPG_BB(G,b,PInd,reg,pp,qq,PARAMETERS,gLen,d,eta,rho,varargin)
IT_MAX = 200000;
TOL = 1e-6;
objTOL = 1e-6;
ALLPRINT = 0;FREQ = 100;
% L = max(eig(G'*G));
n = PARAMETERS(1);
p = PARAMETERS(2);
K = PARAMETERS(3);
% gLen = p*K;
% ax = (1/L);
% ay = ax;
t_start = tic;
% r = prox_pq_Px(z,PInd,v,reg,pp,qq,PARAMETER)
if isempty(varargin)
    x0 = zeros(n^2*p*K,1);
else
    x0 = varargin{1};
end
obj = @(u) 0.5*(norm(G*u-b))^2+reg*normpq(u(PInd),pp,qq,gLen);
% APG INITIALIZATION
k=1;
xkm1 = x0;
zk = xkm1;
xk = xkm1;
tk = 1;
tkm1 = 0;
ck = obj(xk);
qk = 1;
ykm1 = zeros(n^2*p*K,1);
GtG = G'*G;
axITR = 0;
ayITR = 0;
if ALLPRINT
    fprintf('%3s\t%10s\t%10s\n','iter', 'objective','step size');
end
while k<IT_MAX
    ta = tic;
    yk = xk+(tkm1/tk)*(zk-xk)+((tkm1-1)/tk)*(xk-xkm1);
    
    sk = yk-ykm1;
    rk = GtG*(sk);
    ay = norm(sk,2)^2/(sk'*rk);
%     zkp1 = prox_pq_Px(yk-ay*(GtG*yk-G'*b),PInd,ay,reg,pp,qq,PARAMETERS,gLen);
    while 1
        ayITR = ayITR+1;
        zkp1 = prox_pq_Px(yk-ay*(GtG*yk-G'*b),PInd,ay,reg,pp,qq,PARAMETERS,gLen);
        ay = ay*rho;
        if (obj(zkp1)<=obj(yk)-d*norm(zkp1-yk)^2 )|| (obj(zkp1)<=ck-d*norm(zkp1-yk)^2 )
            break
        end
    end
    if (obj(zkp1)<=ck-d*norm(zkp1-yk)^2 )
        xkp1 = zkp1;
    else
        sk = xk-ykm1;
        rk = GtG*(xk-ykm1);
        ax = norm(sk,2)^2/(sk'*rk);
        while 1
            axITR = axITR+1;
            vkp1 = prox_pq_Px(xk-ax*(GtG*xk-G'*b),PInd,ax,reg,pp,qq,PARAMETERS,gLen);
            ax=ax*rho;
            if obj(vkp1) <= ck-d*norm(vkp1-xk,2)^2
                if obj(zkp1)<=obj(vkp1)
                    xkp1 = zkp1;
                else
                    xkp1 = vkp1;
                end
                break
            end
        end
    end
    
    
    objval = obj(xkp1);
    tkp1 = 0.5*(sqrt(4*tk^2+1)+1);
    qkp1 = eta*qk+1;
    ckp1 = (eta*qk*ck+objval)/qkp1;
    ykm1 = yk;
    history.tpi(k,1) = toc(ta);
    history.obj(k,1) = objval;
    if (norm(xkp1-xk)<(TOL*norm(xk))) || (norm(xkp1)==0) || ((k>1) && (abs(objval-history.obj(k-1))<objTOL*objval))
        x = sparse(xkp1);
        break
    else
        x = nan;
    end
    
    
    k = k+1;
    if ((mod(k,FREQ)==0) && ALLPRINT)
        fprintf('%3d\t%10.4f\t%10.4f\n',k,objval,norm(xkp1-xk)/norm(xk))
%         toc(ta);
    end
    tk = tkp1;
    tkm1 = tk;
    zk = zkp1;
%     vk = vkp1;
    xkm1 = xk;
    xk = xkp1;
    qk = qkp1;
    ck = ckp1;
    
end
toc(t_start);
end
function z = normpq(x,p,q,gLen)
n = length(x);
gNo = n/gLen;
tmp = reshape(x,gLen,gNo);
z = norm((sum(abs(tmp).^(p),1)).^(1/p),q)^q;
end