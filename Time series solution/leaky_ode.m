function [sol] = leaky_ode(t, init, theta)
    R_0 = theta(1);
    gamma = theta(2);
    mu = theta(3);
    epsilon = theta(4);
    beta = R_0*(gamma+mu);
    q = gamma/(gamma+mu);
    S = init(1);
    I = init(2);
    R = init(3);
    %dydt = zeros(size(y))
    lambda = beta*I;
    sol = [mu-(mu+lambda)*S,
        lambda*S+epsilon*lambda*R-mu*I-gamma*I,
        gamma*I-mu*R-epsilon*lambda*R];
    
end