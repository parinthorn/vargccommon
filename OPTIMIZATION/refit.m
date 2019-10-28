function x = refit(G,y,P,Px,varargin)
if size(varargin)~=0
    options = optimset('Display','off');
else
    options = optimset('Display','notify');
end
tmpP = P;
% disp(size(tmpP))
% disp(tmpP(Px~=0,:))
tmpP(Px~=0,:) = [];
b=zeros(size(tmpP,1),1);
x = lsqlin(G,y,[],[],tmpP,b,[],[],[],options);
% tmp_G = G;
% xc = zeros(length(x),1);
% Zindex = (x==0);
% nz_ind = (x~=0);
% tmp_G(:,Zindex) = [];
% % disp(size(x))
% % disp(size(tmp_G))
% % disp(size(y))
% tmp = tmp_G\y;
% xc(nz_ind) = tmp;

end