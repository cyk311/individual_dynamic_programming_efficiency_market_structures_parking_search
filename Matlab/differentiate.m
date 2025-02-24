clear;

N = 100;

p0 = 3/202; 
p1 = 1/10100; 

x = 0:100;

m = 0:100;

mu = 50;
sigma = 17;

%% k is a parameter that guarantees the integeration of f_M being 1.
sum = 0;
s = zeros(N+1);
for i = 1:N+1
    s(i) = normpdf(i-1,50,30);
    sum = sum + s(i);
end
k = 1 / sum;

% 计算正态分布的概率密度函数 f_M(m)
f_M = k * normpdf(m, mu, sigma);

% This creates a symmetric space.
half_1 = cat(2, 0:N/2-1, zeros(1,N/2+1)); % (0 1 2...N/2-1 0...0)
mid    = cat(2, cat(2, zeros(1,N/2), N/2), zeros(1,N/2)); % (0...0 N/2 0...0)
half_2 = cat(2, zeros(1,N/2+1), (N/2-1:-1:0)); % (0...0 N/2-1 N/2-2 ... 1 0)

t = half_1 + mid + half_2; % (0 1 2...N/2-1 N/2 N/2-1...2 1 0)

%% Define p_dx
p_dx(1,:) = p0 - p1 * x;
p_dx(2,:) = 1 / dot(ones(size(t)), p0 - p1 * t) * (p0 - p1 * t);


%% Calculate partial derivatives and the original function
dEp_dx = zeros(size(x));
p_x   = zeros(size(p_dx));

for i = 1:length(x)
    sum_term_1 = 0;
    sum_term_0 = 0;
    for j = 1:length(m)
        term_1 = f_M(j) * (- m(j) * p1 * (1 - p_dx(1,i))^(m(j) - 1));
        term_0 = f_M(j) * (1 - (1 - p_dx(:,i)).^m(j));
        sum_term_1 = sum_term_1 + term_1;
        sum_term_0 = sum_term_0 + term_0;
    end
    dEp_dx(i) = sum_term_1;
    p_x(:,i)  = sum_term_0;
end

p_x_1 = p_x(1,:);
p_x_2 = p_x(2,:);

%% Plot
% Single Entrance
figure;   % 创建绘图框
subplot(1,2,1);   % 绘图框共有1行2列，副图绘制在第1个位置
plot(x, dEp_dx);
xlabel('x');
title('$\frac{dE[p(M|x)]}{dx}$','Interpreter','latex','FontSize',14);

subplot(1,2,2);
plot(x, p_x_1);     % 绘图框共有1行2列，副图绘制在第2个位置
xlabel('x');
title('$E[p(M|x)]$','Interpreter','latex','FontSize',14);

% Double Entrance
figure;
plot(x, p_x_2);
xlabel('x');
title('$E[p(M|x)]$','Interpreter','latex','FontSize',14);
