clear;

%% 设置需要涂色的坐标
A = [(93:-1:90), 95, (88:-1:75), 96, (74:-1:62), 97, (61:-1:51), 98, (50:-1:42), 99, (41:-1:35), 100];
B = [49, 51, 50];
C = [2, 4, 5, 12, 29, 36, 44, 47, 52, 54, 62, 63, 65, 68, 77, 79, 86, 87, 89, 94];

% 使用 ismember 找到 A 中与 C 重复的元素
toRemove = ismember(A, C);

% 删除重复的元素
A(toRemove) = [];

% 创建一个 的零矩阵
matrix_1 = zeros(length(A)+1, 101); % Single Entrance
matrix_2 = zeros(length(B)+1, 101); % Double Entrance

%% 停车者停靠位置
% Single Entrance
for i = 1:length(A)
    matrix_1(i:length(A),A(i)) = 2;
end

% Double Entrance
for i = 1:length(B)
    matrix_2(i:length(B),B(i)) = 2;
end

% 预设已占位置
for i = C
    matrix_1(:,i) = 1;
    matrix_2(:,i) = 1;
end

%% 使用 pcolor 绘制矩阵
figure;

subplot(2,1,1);
pcolor(matrix_1);

% 设置颜色映射（矩阵每一行是颜色的RGB）
colormap([1 1 1; 1 0 0; 0 0 1]); % 第一行是白色，第二行是红色，第三行是蓝色

grid on;
set(gca, 'GridColor', [0 0 0], 'LineWidth', 0.3); % 设置网格线颜色和宽度

% 设置坐标轴
axis equal;
axis tight;
set(gca, 'XAxisLocation', 'bottom'); % 将 X 轴刻度放在底部

% 添加标题
ylabel('Time')
title('Single Entrance Occupation');

subplot(2,1,2);
pcolor(matrix_2);
colormap([1 1 1; 1 0 0; 0 0 1]); 
grid on;
set(gca, 'GridColor', [0 0 0], 'LineWidth', 0.3);
axis equal;
axis tight;
set(gca, 'XAxisLocation', 'bottom'); 
ylabel('Time')
title(sprintf('Double Entrance Occupation \n (Parkers enter from both sides alternately)')); 

