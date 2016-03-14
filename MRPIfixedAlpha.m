function [Aout,bout,iter] = MRPIfixedAlpha(Ain,bin,A,D,G,iMax,alpha)

Lambda = Ain;
lambda = bin;

LambdaNext = [];
lambdaNext = [];

model = struct('A',G,...
            'rhs', ones(size(G,1),1),...
            'sense', char(ones(size(G,1),1)*'<'),...
            'lb',ones(size(G,2),1)*-inf,...
            'modelsense','max');
param = struct('OutputFlag', 0);

iter = 1;
while and(~isContained(Lambda,lambda,LambdaNext,lambdaNext),iter<iMax)
    
    if mod(iter,10) == 0
        clc
        fprintf('Currently at iteration %d, there are currently %d constraints\n',iter,length(lambda))
    end
    
    if iter~=1
        Lambda = LambdaNext;
        lambda = lambdaNext;
    end
    
    LambdaNext = zeros(size(Lambda));
    lambdaNext=zeros(size(lambda));
    for i = 1:length(lambda)
        model.obj = (Lambda(i,:)*D)';
        res = gurobi(model,param);
        LambdaNext(i,:) = Lambda(i,:)*A;
        lambdaNext(i) = lambda(i) - res.objval*alpha;
    end
    
    [LambdaNext,lambdaNext] = bigReduce([Lambda;LambdaNext],...
        [lambda;lambdaNext]);
    iter = iter+1;

end

Aout = LambdaNext;
bout = lambdaNext;