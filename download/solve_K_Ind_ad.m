N = 31;
K = 2;
L0 = ones(N,1);
% A = ones(N,K);
% Am = repmat(linspace(1,3,N)',1,K);
% A = [a(1:15,:);fliplr(a(1:16,:)')'];
A = ones(N,K);
A(:,2) = A(:,2)+0.1;
s = 6.4;
theta = 7;
t = 0.01;
T = tau(t,N);
M = ones(N);
w0 = ones(N,1); 
options=optimset('Maxfuneval',6.2e+04,'Algorithm','Levenberg-Marquardt');
w = fsolve(@(w) f(w,A,M,T,L0),w0,options);
w1 = w;
L1 = Labor(w1,A,M,T,L0);
P1 = ds_price(w1,A,T);
plot(L1);
% hold on
% plot(w1);
% plot(P1);
legend('labor','wage','price');
function F = f(w,A,M,T,L0)
    F = eq_condition(w,A,M,T,L0)'*eq_condition(w,A,M,T,L0);
end

function eq = eq_condition(w,A,M,T,L0)
    s = 6.4;
    eta = 2.5;
    wm = repmat(w,1,length(w));
    L = Labor(w,A,M,T,L0);
    P = ds_price(w,A,T);
    Pm = (((wm).^(1-s).*T.^(1-s)*(A.^(-1)).^(1-s)).^(1/(1-s)));
    omega0 = w.*L.*(P'.^(eta-1));
    omega1 = repmat(omega0,1,size(Pm,2));
    omega2 = omega1 .* Pm.^(s-eta);
    omega3 = T.^(1-s)*omega2;
%         omega3 = w'.*L'.*(P.^(eta-1))*T.^(1-s)*Pm.^(s-eta);
%     omega = (A).^(1-s).*L;
%     omega = omega.^(-1);
%     psi = w.*L;
%     psi = psi./(P'.^(1-s));
%     psi = T.^(1-s)*psi;
    wn = (A.^(s-1).*omega3*ones(size(Pm,2),1)).*(L.^(-1));
    eq = w - wn.^(1/s);
end

function L = Labor(w,A,M,T,L0)
    P = Pi(w,A,M,T);
    L = P * L0;
end

function Pi = Pi(w,A,M,T)
    theta = 1.1;
    p = ds_price(w,A,T);
    phi = (w./p').^(theta);
    Phi = M.^(-theta);
    Phi = Phi*phi;
    Phi = Phi.^(-1);
    Pi = phi*Phi';
    Pi = Pi.*M.^(-theta);
end

function p = ds_price(w,A,T)
    s = 6.4;
    eta = 2.5;
    wm = repmat(w,1,length(w));
    p = (((((wm).^(1-s).*T.^(1-s)*(A.^(-1)).^(1-s)).^(1/(1-s))).^(1-eta)*ones(size(A,2),1)).^(1/(1-eta)))';
end

function T = tau(t,N)
    T = ones(N);
    for i = 1:N
        for j = 1:N
            T(i,j) = exp(t*abs(i-j));
        end
    end
end
