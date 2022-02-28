
%%%% ����ͬ�������ά�㼯��׼ʾ�� %%%%
clc;clear;close all;warning off

addpath(genpath('RigidEstimation'));
N_in = 100;                    % �㼯����ȷƥ������
noise=0.1;                     % �㼯����������0.1m
nTest=100;                     % 100�ζ���ʵ��
O_rate = 0.9;                  % �ֲ���ʣ�����90%

success=0;                     %�ɹ���
for i=1:nTest
    i
    N_ou=round(N_in/(1-O_rate))-N_in;         % �ֲ�۲�����
    X=randn(3,N_in+N_ou)*100;                 % ������ɸ�˹�ֲ��㼯X ~ N(0,100^2)
    r=xrand(3,1,[-pi/2 pi/2]);                % ���������ת����
    R= rodrigues(r);                          % ������ת����
    t=xrand(3,1,[-100 100]);                  % �������ƽ������
    
    Ygt=R*X+repmat(t,1,N_in+N_ou);         % ������ֵƥ���Ygt
    randomvals=randn(3,N_in+N_ou)*noise;   % ���ɸ�˹�ֲ����� randomvals ~ N(0,0.1^2)
    Y=Ygt+randomvals;                         % ���ɴ�����ƥ���
    
    maxYn = max(Y,[],2);
    minYn = min(Y,[],2);
    
    if N_ou~=0                                % ������ɴֲ�㣬X��Yo֮��û���κι�ϵ���������κ�ģ��
        Yo=randn(3,N_ou)*100; 
        Y=[Y(:,1:N_in) Yo];
    end
    
    W=ones(3,size(X,2));                       
    [R0,t0]=IRLS_Welsch(X,Y,W,100);             % Welsch+IRLS method
    residual = R0*X+repmat(t0,1,N_in+N_ou)-Y;
    rms0 = norm(sqrt(sum(residual(1:N_in).^2,2)./N_in));  %������ȷƥ����RMSE
        

    [R1,t1,i] = SA_Cauchy_IRLS_rigid_model(X,Y,1.3);
    residual = R1*X+repmat(t1,1,N_in+N_ou)-Y;
    rms1 = norm(sqrt(sum(residual(1:N_in).^2,2)./N_in));  %������ȷƥ����RMSE
    if rms1<3*noise
        success = success+1;
    end

end
success = success/nTest









