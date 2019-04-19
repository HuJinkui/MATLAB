%%  清空环境变量
clear all
clc
%% 导入数据
load faciesdatawithoutxy45.mat           %导入数据，无坐标，只有特征值（一列为一个特征）
attributes = faciesdatawithoutxy45';     %数据旋转，一行为一个特征
%% 数据归一化
attributes = mapminmax(attributes);
%%
% 60个样本
P_train = attributes(:,1:end);
%%  SOM神经网络创建、训练及仿真测试
%%
% 1. 创建网络
net = newsom(P_train,[3 4]);
%plotsom(hextop(3,4));
%%
% 2. 设置训练参数
net.trainParam.epochs = 300;
%%
% 3. 训练网络
net = train(net,P_train);

%%
t_sim_sofm_1 = sim(net,P_train);
T_sim_sofm_1 = vec2ind(t_sim_sofm_1);  %向量和稀疏矩阵之间转换
result_som = T_sim_sofm_1' ;

% save seismicfacies89 result_som;     %将分类结果输出为mat格式（MATLAB数据）
