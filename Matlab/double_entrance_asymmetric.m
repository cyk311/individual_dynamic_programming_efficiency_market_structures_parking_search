clear;

N     = 100;
x     = 0:N;
mu    = 50;
sigma = 30;
m     = 0:100;
E     = 0:100;

%% k_1 is a parameter that guarantees the integeration of f_M being 1.
% M is a random variable representing the total number of vacant spaces in the parking lot.
sum_1 = 0;
sum = zeros(1,N+1);
for i = 1:N+1
    sum(i) = normpdf(i-1,mu,sigma);
    sum_1 = sum_1 + sum(i);
end
k_1 = 1 / sum_1;

%% Parameter Scanning.
equilibrium = zeros(1, 101);
for e = E
    B     = 0:e;
    peak  = max(e, 100-e);
    left  = max(0, 100-2*e);
    right = max(0, 2*e-100);
    t = [left:peak-1, peak, peak-1:-1:right];

    %% k_2 is a parameter that guarantees the integeration of (p0-p1*x) on [0,100] being 1.
    p0_0 = 3/202;   
    p1_0 = 1/10100; % p0_0 and p1_0 matches with p0 and p1 in the single entrance scenario.
    p_dx = p0_0 - p1_0 * t;
    k_2 = 1 / dot(ones(size(t)), p_dx);
    
    p0 = k_2 * p0_0;
    p1 = k_2 * p1_0;
    
    %% This is how p_x is generated
    f_M   = k_1 * normpdf(m, mu, sigma);
    
    p_x   = zeros(size(t));
    
    for i = 1:length(t)
        ti = t(i);
        sum_term_0 = 0;
        for j = 1:length(m)
            mj = m(j);
            term_0 = f_M(j) * (1 - (1 - p0 + p1 * ti)^mj);
            sum_term_0 = sum_term_0 + term_0;
        end
        p_x(i) = sum_term_0;
    end

    x_0 = zeros(size(B));
    for b = B
        du  = 5;
        u   = du * (100 - abs(e-x));
        c   = x;
        
        % Stop function.
        stop = u - c - (u(b+1) - c(b+1) - 2*(c - c(b+1)));
        S = N;
        for i =1:N
            if (stop(i) >= 0) && (stop(i+1) < 0)
                S = i-1;
            end
        end
    
        % E_net1
        E_net1 = zeros(1,S+1);
        for i = 1:S
            sum_tobe = (u(i+1:S+1)-c(i+1:S+1)) .* p_x(i+1:S+1);
            E_net1(i) = ones(1,S+1-i) * sum_tobe';
        end
        
        % E_net2
        E_net2 = zeros(1,S+1);
        pxt = 1. - p_x;
        A = u(b+1)-c(b+1) - 2*(c(S+1)-c(b+1));
        
        for i = 1:S
            prod_E_net2 = prod(pxt(i+1:S+1));
            E_net2(i) = A * prod_E_net2;
        end
        
        % E_net3
        E_net3_0 = (u(b+1)-c(b+1))-2*(c(x+1)-c(b+1));
        E_net3   = E_net3_0(1:S+1);
        
        E_net = E_net1 - E_net2 - E_net3;
    
        % Zero of E_net
        x_0(b+1) = S; % Initialize x_0 in case that no solution exists in the domain.
        for i = 1:S
            if E_net(i) > 0 && E_net(i+1) < 0
                x_0(b+1) = i;
            end
        end
    end
    
    for i = 1:(length(B)-1)
        if (x_0(i)-B(i))^2 < (x_0(i+1)-B(i+1))^2
            equilibrium(e+1) = x_0(i);
        end
    end
end

distance = E - equilibrium;
total_distance = distance + flip(distance);

figure;
plot(E, total_distance);
grid on;
title('Total Distance From the Equilibrium Point to the Destination')

figure;
plot(E, equilibrium, E, E);
legend('equilibrium point','e', Location='northwest')
xlabel('e')

grid on;
title('Equilibrium Points for Different Destinations')