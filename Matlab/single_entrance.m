clear;

N     = 100;
B     = 0:100;
u     = 5;
c     = 1;
x     = 0:N;
mu    = 50;
sigma = 30;
m     = 0:100;
p0    = 3/202; 
p1    = 1/10100; 

%% k is a parameter that guarantees the integeration of f_M being 1.
sum = 0;
s = zeros(1,N+1);
for i = 1:N+1
    s(i) = normpdf(i-1,50,30);
    sum = sum + s(i);
end
k = 1 / sum;

%% "p(x)" (i.e. "p_x") is the possibility of a vacant position existing at point x.
% This is how p_x is generated.
f_M = k * normpdf(m, mu, sigma);

p_x = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    sum_term_0 = 0;
    for j = 1:length(m)
        mj = m(j);
        term_0 = f_M(j) * (1 - (1 - p0 + p1 * xi)^mj);
        sum_term_0 = sum_term_0 + term_0;
    end
    p_x(i) = sum_term_0;
end

%% The latex expression of the net expectation of utility (E_net) is as follow:
%{

E_{net} = [(u-c)\sum_{x+1}^N x p(x)] -
          [b(u-c)-2c(N-b)]\prod_{x+1}^N [1-p(x)] - 
          [b(u-c) - 2c(x-b)]
Which is devided into 3 parts for convenience of calculation. (i.e. E_net1, E_net2 and E_net3)

%}

x_0 = zeros(size(B));
for b = B
    E_net1 = zeros(1,N+1);
    for i = 1:N
        sum_E_net1 = x(i+1:N+1) * p_x(i+1:N+1)';
        E_net1(i) = (u-c) * sum_E_net1;
    end
    
    E_net2 = zeros(1,N+1);
    pxt = 1. - p_x;
    A = b*(u-c) - 2*c*(N-b);
    
    for i = 1:N
        prod_E_net2 = prod(pxt(i+1:N+1));
        E_net2(i) = A * prod_E_net2;
    end
    
    E_net3 = b*(u-c)-2*c*(x-b);
    
    E_net = E_net1 - E_net2 - E_net3;

    % Zero of E_net
    x_0(b+1) = N; % Initialize x_0 in case that no solution exists in the domain.
    for i = 1:N
        if E_net(i) > 0 && E_net(i+1) < 0
            x_0(b+1) = i-1;
        end
    end
    
    %% Plot
    if b == 94
        figure;
        plot(x,E_net);
        ylabel('$E_{net}$',Interpreter='latex');
        grid on;
    end

end

figure;
plot(B,x_0);
xlabel('$b$',Interpreter='latex');
ylabel('$x_0(b)$',Interpreter='latex')

hold on;


plot(B,B);