clear all
pF = [0.05:0.001:0.95];
pS = [0.05:0.001:0.95];
nu = 2.5;
lambda = 0.1115;
beta = 0.0065;
k = 1; % optima léptetõ
for f = 1:901
    for s = 1:901
        qopt = 1 - pS(s);
        betaopt = 1 - (1 - beta)*pF(f)*nu/qopt;
        lambdaopt = beta*pF(f)*nu*lambda/(qopt*betaopt);
        if (1 - beta)*nu*pF(f) + pS(s) < 1 && nu*pF(f) + pS(s) >= 1 && betaopt <= 0.95 && betaopt >= 0.001 && lambdaopt <= 0.95 && lambdaopt >= 0.05
            optima(:,k) = [pF(f),pS(s),qopt,betaopt,lambdaopt];
            k = k + 1;
        end
    end
end
