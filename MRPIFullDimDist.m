function [Aout,bout,i] = MRPIFullDimDist(Ain,bin,G,A,D,iMax,chars)

if nargin<6
    chars = ' ';
end

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

i = 1;
while and(~isContained(Lambda,lambda,LambdaNext,lambdaNext),i<iMax)
    
    if mod(i,10) == 0
        clc
        fprintf([chars,'\n']);
        fprintf('Currently at iteration %d, there are currently %d constraints\n',i,length(lambda))
    end
    
    if i~=1
        Lambda = LambdaNext;
        lambda = lambdaNext;
        
    end
    
    LambdaNext = zeros(size(Lambda));
    lambdaNext=lambda;
    for i = 1:length(lambda)
        model.obj = (Lambda(i,1:2)*D)';
        res = gurobi(model,param);
        LambdaNext(i,:) = [Lambda(i,1:2)*A,res.objval];
    end
    
    
    [LambdaNext,lambdaNext] = bigReduce([Lambda;LambdaNext],...
        [lambda;lambdaNext]);
    i = i+1;

end

Aout = LambdaNext;
bout = lambdaNext;