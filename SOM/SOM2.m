%%  ��ջ�������
clear all
clc
%% ��������
load faciesdatawithoutxy45.mat           %�������ݣ������ֻ꣬������ֵ��һ��Ϊһ��������
attributes = faciesdatawithoutxy45';     %������ת��һ��Ϊһ������
%% ���ݹ�һ��
attributes = mapminmax(attributes);
%%
% 60������
P_train = attributes(:,1:end);
%%  SOM�����紴����ѵ�����������
%%
% 1. ��������
net = newsom(P_train,[3 4]);
%plotsom(hextop(3,4));
%%
% 2. ����ѵ������
net.trainParam.epochs = 300;
%%
% 3. ѵ������
net = train(net,P_train);

%%
t_sim_sofm_1 = sim(net,P_train);
T_sim_sofm_1 = vec2ind(t_sim_sofm_1);  %������ϡ�����֮��ת��
result_som = T_sim_sofm_1' ;

% save seismicfacies89 result_som;     %�����������Ϊmat��ʽ��MATLAB���ݣ�
