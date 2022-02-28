function [R,t,i] = SA_Cauchy_IRLS_rigid_model(src,dst,tau)

prev_cost=10^15;                            %��ʼ��������
maxIter =100;                               %����������
n=size(src,2);                              %�۲�ֵ����
weights=ones(1,n);                          %Ȩֵ����
for i=1:maxIter
    [R,t,fit]=rigidTrans(src,dst,weights);  %����ƥ���src��dst������㼯����ĸ���任R,t;
    residuals = sqrt(sum((fit-dst).^2));    %����в�
    if i==1
        alpha = max(abs(residuals));        %alpha��ֵ��Ϊ���в�ֵ
    end
    cost = sum(weights.*residuals.^2);      %���㵱ǰ����������w*r^2
    inl_lable = residuals<3*alpha;            %�޳����ֲִ��
    src=src(:,inl_lable);                   %����ƥ���
    dst=dst(:,inl_lable);
   
    E=residuals(inl_lable);                 %����ƥ���в�
    weights=alpha^2./(alpha^2+E.^2);        %����SA-CauchyȨֵ
    cost_diff = abs(cost - prev_cost);      %������
    alpha = alpha / tau;                    %����alpha��
    prev_cost = cost;                       %��������
    
    if cost_diff < 1e-2||alpha<1
        break;
    end
end

