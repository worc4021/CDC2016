function cand = approxAlphaStar(Lambda,A,D,G,kMax,tol)

calc = zeros(size(Lambda,1),kMax);
k = 1;

model = struct('obj',ones(1,size(Lambda,2)),...
            'A',G,...
            'rhs', ones(size(G,1),1),...
            'sense', char(ones(size(G,1),1)*'<'),...
            'lb',ones(size(G,2),1)*-inf,...
            'modelsense','max');
param = struct('OutputFlag', 0);
res = gurobi(model,param);
wMax = res.x;

Ak = D;

for i = 1:size(Lambda,1)
    calc(i,1) = Lambda(i,:)*Ak*diag(sign(Lambda(i,:)*Ak))*wMax;
end


while (k<kMax && max(calc(:,k))>tol)
    Ak = A*Ak;
    k = k+1;
    for i = 1:size(Lambda,1)
        calc(i,k) = Lambda(i,:)*Ak*diag(sign(Lambda(i,:)*Ak))*wMax;
    end
end

cand = max(sum(calc,2));