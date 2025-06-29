clc; clear; close all;

T = readtable('solar.csv');
simtime = height(T);
delt = 1;

Varmin = [0, 0, 0, 0, 0000000];
Varmax = [10, 10000, 1000, 1, 3063777.09];

CostFunction = @(x) netpc(x, Varmax, T, simtime);

nVar = 5;
Varsize = [1, nVar];
%% Parameters of PSO
Simtime = 60;      % Max number of iterations
npop = 30;          % population size
kappa = 1.0;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
w = chi;
wdamp = 1;
c1 = chi*phi1;
c2 = chi*phi2;

%% Initialization
empty_particle.Position = [];
empty_particle.velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
GlobalBest.Cost = inf;
particle = repmat(empty_particle, npop, 1);

for i = 1:npop
    particle(i).Position = unifrnd(Varmin, Varmax, Varsize);
    
    try
        particle(i).Cost = CostFunction(particle(i).Position);
    catch err
        disp('Error in CostFunction:');
        disp(err.message);
        rethrow(err);  % optional: remove in production
    end
    
    particle(i).velocity = zeros(Varsize);
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    if particle(i).Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end

BestCosts = zeros(Simtime,1);

%% Main loop of PSO
for it = 1:Simtime    
    for i = 1:npop
        % Update velocity
        particle(i).velocity = w * particle(i).velocity ...
            + c1 * rand(Varsize) .* (particle(i).Best.Position - particle(i).Position) ...
            + c2 * rand(Varsize) .* (GlobalBest.Position - particle(i).Position);

        % Update position
        particle(i).Position = particle(i).Position + particle(i).velocity;
        particle(i).Position = max(particle(i).Position, Varmin);
        particle(i).Position = min(particle(i).Position, Varmax);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);

        % Update the personal best
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;

            % Update the global best
            if particle(i).Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end

    BestCosts(it) = GlobalBest.Cost;
end

%% Result Section
[cost,extra] = netpc(GlobalBest.Position,Varmax,T,simtime,true);

figure;

Bestcost = BestCosts-(sum(extra.penaltyg)+sum(extra.penalty));

plot(1:Simtime, Bestcost, 'b-', 'LineWidth', 2);
ylim([1e7,1e10]);
xlabel('Iteration');
ylabel('Net Present Cost (NPC)($)');
title('Convergence of PSO');
grid on;
display(['Final cost = ', num2str(Bestcost(Simtime))]);
time = 1:height(T);  % or time = 1:simtime if simtime = 288

figure;
plot(time, extra.P_totalf, 'r', 'LineWidth', 2); hold on;
plot(time, T.P_el, 'b--', 'LineWidth', 2);hold on;
plot(time, extra.Diff,'g','LineWidth', 2);
legend('P_{total}', 'Load demand');
xlabel('Time (hr)');
ylabel('Power (W)');
title('Optimized Power Output vs Load Demand');
grid on;
xlim([1, max(time)]);  % Ensure full x-axis coverage

time2 = 1:height(T);
figure;
plot(time2, extra.BESSstate, 'b', 'LineWidth', 2);
% plot(time, T.P_el, 'b--', 'LineWidth', 2);
legend('State of charge');

xlabel('Time (hr)');
ylabel('SOC(%)');
title('BESS variation');
grid on;
xlim([1, max(time2)]);  % Ensure full x-axis coverage
