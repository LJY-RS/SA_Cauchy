
%%%% 基于同名点的三维点集配准示例 %%%%
clc;clear;close all;warning off

addpath(genpath('RigidEstimation'));
N_in = 100;                    % 点集中正确匹配点个数
noise=0.1;                     % 点集测量误差，假设0.1m
nTest=100;                     % 100次独立实验
O_rate = 0.9;                  % 粗差比率，假设90%

success=0;                     %成功率
for i=1:nTest
    i
    N_ou=round(N_in/(1-O_rate))-N_in;         % 粗差观测数量
    X=randn(3,N_in+N_ou)*100;                 % 随机生成高斯分布点集X ~ N(0,100^2)
    r=xrand(3,1,[-pi/2 pi/2]);                % 随机生成旋转向量
    R= rodrigues(r);                          % 构成旋转矩阵
    t=xrand(3,1,[-100 100]);                  % 随机生成平移向量
    
    Ygt=R*X+repmat(t,1,N_in+N_ou);         % 生成真值匹配点Ygt
    randomvals=randn(3,N_in+N_ou)*noise;   % 生成高斯分布噪声 randomvals ~ N(0,0.1^2)
    Y=Ygt+randomvals;                         % 生成带噪声匹配点
    
    maxYn = max(Y,[],2);
    minYn = min(Y,[],2);
    
    if N_ou~=0                                % 随机生成粗差点，X和Yo之间没有任何关系，不满足任何模型
        Yo=randn(3,N_ou)*100; 
        Y=[Y(:,1:N_in) Yo];
    end
    
    W=ones(3,size(X,2));                       
    [R0,t0]=IRLS_Welsch(X,Y,W,100);             % Welsch+IRLS method
    residual = R0*X+repmat(t0,1,N_in+N_ou)-Y;
    rms0 = norm(sqrt(sum(residual(1:N_in).^2,2)./N_in));  %计算正确匹配点的RMSE
        

    [R1,t1,i] = SA_Cauchy_IRLS_rigid_model(X,Y,1.3);
    residual = R1*X+repmat(t1,1,N_in+N_ou)-Y;
    rms1 = norm(sqrt(sum(residual(1:N_in).^2,2)./N_in));  %计算正确匹配点的RMSE
    if rms1<3*noise
        success = success+1;
    end

end
success = success/nTest









