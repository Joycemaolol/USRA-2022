clear;
N = 100000;
gamma = 365/8;
mu = 1/60;
q = gamma/(gamma+mu);
eps_l = linspace(0, 1, N);
beta = 400;
R0=beta/(gamma+mu);


realparts = zeros(N,2);
imagparts = zeros(N,2);
check = zeros(N,2);

for n=1:N
    epsilon = eps_l(n);
    q=gamma/(gamma+mu);
    lambda = (-mu*(epsilon*(R0 + q - 1) - 1) - sqrt(-4*(q - 1)*epsilon*mu^2*(R0 - 1) + mu^2*(epsilon*(R0 + q - 1) - 1)^2))/(2*(q - 1)*epsilon);
    S_leaky = mu/(lambda + mu);
    I_leaky = lambda/beta;
    J_leaky = [-(lambda+mu), -beta*S_leaky; 
        lambda-epsilon*lambda, beta*S_leaky+epsilon*beta*(1-S_leaky)-2*epsilon*lambda-(mu+gamma)];    
    e_leaky = eig(J_leaky);
    %e = cplxpair(e);
    realparts_leaky(n,:) = real(e_leaky);
    imagparts_leaky(n,:) = imag(e_leaky);



    I_sir = mu/(mu+gamma)-mu/beta;
    lambda_sir = beta*I_sir;
    S_sir = mu/(beta*I_sir + mu);
    J_sir = [-lambda_sir-mu, -beta*S_sir;
        lambda_sir, beta*S_sir-mu-gamma];    
    e_sir = eig(J_sir);
    %e = cplxpair(e);
    realparts_sir(n,:) = real(e_sir);
    imagparts_sir(n,:) = imag(e_sir);


    I_sis = (beta*mu-mu^2-gamma*mu)/(mu*beta);
    lambda_sis = beta*I_sis;
    S_sis = (mu+gamma*I_sis)/(mu+lambda_sis);
    J_sis = [-lambda_sis-mu, gamma-beta*S_sis;
        lambda_sis, beta*S_sis-mu-gamma];    
    e_sis = eig(J_sis);
    %e = cplxpair(e);
    realparts_sis(n,:) = real(e_sis);
    imagparts_sis(n,:) = imag(e_sis);
end
%%
%realparts_leaky(215:end,[1,2,3,4]) = realparts_leaky(215:end,[2,3,1,4]);
%imagparts_leaky(215:end,[1,2,3,4]) = imagparts_leaky(215:end,[2,3,1,4]);
subplot(1, 2, 1)

hold all
p1 = plot(eps_l,realparts_leaky(:,1),'.',Color='red', LineWidth=2);
plot(eps_l,realparts_leaky(:,2),'.',Color='red', LineWidth=2)


p2 = plot(eps_l,realparts_sir(:,1),'--',Color='blue', LineWidth=5);
plot(eps_l,realparts_sir(:,2),'--',Color='blue', LineWidth=5)


p3 = plot(eps_l,realparts_sis(:,1),':',Color='green', LineWidth=2);
plot(eps_l,realparts_sis(:,2),':',Color='green', LineWidth=2)

plot([1/(R0+q-1),1/(R0+q-1)],[-500,500],'k:','LineWidth',2)
legend([p1, p2, p3], {"Leaky", "SIR", "SIS"}, 'Fontsize',30, 'Position',[0.2 0.6 0.1 0.2])

grid

text(1/R0+0.05,-270,'$\varepsilon_{L}=\frac{1}{R_{0}}$','Interpreter','latex','Fontsize',30)
ylim([-360 10])
set(gca,"FontSize",20)

xlabel('$\varepsilon_{L}$','Interpreter','latex','Fontsize',30)
ylabel('Real Parts of Eigenvalues', 'Fontsize',30)
%leg=legend('Eigenvalue 1','Eigenvalue 2','Eigenvalue 3','Eigenvalue 4','Fontsize',20);
%set(leg,'Position',[0.3,0.4,0,0]);
% x0=10;
% y0=10;
% width=440;
% height=400;
% set(gcf,'position',[x0,y0,width,height])

subplot(1, 2, 2)

hold all
p4 = plot(eps_l,imagparts_leaky(:,1),'.',color='red', LineWidth=2);
plot(eps_l,imagparts_leaky(:,2),'.',color='red', LineWidth=2)


p5 = plot(eps_l,imagparts_sir(:,1),'--',color='blue', LineWidth=2);
plot(eps_l,imagparts_sir(:,2),'--',color='blue', LineWidth=2)


p6 = plot(eps_l,imagparts_sis(:,1),':',Color='green', LineWidth=2);
plot(eps_l,imagparts_sis(:,2),':',Color='green', LineWidth=2)

plot([1/(R0+q-1),1/(R0+q-1)],[-150,16],'k:','LineWidth',2)
legend([p4, p5, p6], {"Leaky", "SIR", "SIS"}, 'Fontsize',30, 'Position',[0.2 0.6 0.1 0.2])
grid
ylim([-3 3])
set(gca,"FontSize",20)

xlabel('$\varepsilon_{L}$','Interpreter','latex', 'Fontsize',30)
ylabel('Imaginary Parts of Eigenvalues', 'Fontsize',30)
text(1/R0+0.05,-1,'$\varepsilon_{L}=\frac{1}{R_{0}}$','Interpreter','latex','Fontsize',30)
%leg=legend('Eigenvalue 1','Eigenvalue 2','Eigenvalue 3','Eigenvalue 4','Fontsize',20);
% x0=10;
% y0=10;
% width=440;
% height=400;
% set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'eigen_compared.jpg')