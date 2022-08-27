clear
epsilon = 0.25;
gamma = 365/8;
mu = 1/60;
q = gamma/(gamma+mu);


R_0 = 10;
beta = R_0*(gamma+mu);

theta = [R_0, gamma, mu, epsilon];

lambda = mu * (-(epsilon*(1-R_0-q)+1)+sqrt((epsilon*(1-R_0-q)-1)^2 - 4*epsilon*q*R_0))/(2*epsilon*(1-q));%at endemic equilibrium
S0_temp = mu/(lambda + mu);
I0_temp = lambda/beta*100;
R0_temp = gamma*I0_temp/(epsilon*lambda + mu);
S0_temp = 0.9
I0_temp = 0.1
R0_temp = 0


N = S0_temp+I0_temp+R0_temp;
S0 = S0_temp/N;
I0 = I0_temp/N;
R0 = R0_temp/N;


t_start = 0;
t_stop = 15;
%% system of ode solver
%solve leaky
 odeoptions = odeset('RelTol',1e-13, ...
     'AbsTol',1e-16, ...
     'Refine', 1);
 [t, y] = ode45(@(t, y) leaky_ode(t, y, theta), [t_start t_stop],[S0, I0, R0], odeoptions);
 %t=transpose(linspace(0,100,1000));
 %y =transpose(deval(leakymodel,t));

% solve SEIR
odeoptions = odeset('RelTol',1e-13, ...
     'AbsTol',1e-16, ...
     'Refine', 1);
 [t_r, y_r] = ode45(@(t, y) SIR(t, y, theta), [t_start t_stop],[S0, I0, R0], odeoptions);

% solve SEIS
 odeoptions = odeset('RelTol',1e-13, ...
     'AbsTol',1e-16, ...
     'Refine', 1);
 [t_s, y_s] = ode45(@(t, y) SIS(t, y, theta), [t_start t_stop],[S0, I0, R0], odeoptions);

%% The effective S compartment
figure(3)
plot(t, y(:,1)+epsilon*y(:,3),'--', LineWidth=2)
hold on
plot(t_r, y_r(:,1), '-', LineWidth=2)
hold on
plot(t_s, y_s(:,1)+y_s(:,3), ':', LineWidth=2)
legend("Leaky", "SIR", "SIS", 'FontSize',30)
set(gca,"FontSize",30)
ylabel('Effective S compartment','Interpreter','LaTeX','FontSize',30);
xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',30);
ylim([-0.01 0.5])

%% The effective R compartment 
figure(2)
plot(t, y(:,3)-epsilon*y(:,3),'-.', LineWidth=2)
hold on
plot(t_r, y_r(:,3), '-b', LineWidth=2)
hold on
plot(t_s, y_s(:,3)-y_s(:,3), '-r', LineWidth=2)
legend("Leaky", "SIR", "SIS", 'FontSize',20)
ylabel('Solution','Interpreter','LaTeX','FontSize',20);
xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',20);
title('Effective R compartment', 'FontSize',20)
%% Plot the S, E, I, R compartment solution of the leaky, SIR, and SIS model
 
subplot(2, 2, 1)
plot(t, y(:,1),'-.', LineWidth=2)
hold on
plot(t_r, y_r(:,1), '-b', LineWidth=2)
hold on
plot(t_s, y_s(:,1), '-r', LineWidth=2)
legend("leaky", "seir", "seis")
ylabel('Solution','Interpreter','LaTeX','FontSize',14);
xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',14);
title('S compartment')

subplot(2, 2, 2)
plot(t, y(:,2),'-.', LineWidth=2)
hold on
plot(t_r, y_r(:,2), '-b', LineWidth=2)
hold on
plot(t_s, y_s(:,2), '-r', LineWidth=2)
legend("leaky", "seir", "seis")
ylabel('Solution','Interpreter','LaTeX','FontSize',14);
xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',14);
title('I compartment')

subplot(2, 2, 3)
plot(t, y(:,3),'-.', LineWidth=2)
hold on
plot(t_r, y_r(:,3), '-b', LineWidth=2)
hold on
plot(t_s, zeros(size(t_s)), '-r',LineWidth=2)
legend("leaky", "seir", "seis")
ylabel('Solution','Interpreter','LaTeX','FontSize',14);
xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',14);
title('R compartment')


%% plot the reinfection probability

% I = exp(y(:,3));
% 
% rp = q*epsilon*beta*I./(epsilon*beta*I+repelem(mu,length(I)));
% 
% figure(4)
% plot(t, rp,'-o')
% ylabel('Reinfection probability','Interpreter','LaTeX','FontSize',14);
% xlabel('Time (years), $t$','Interpreter','LaTeX','FontSize',14);
% title('Reinfection Probability')


