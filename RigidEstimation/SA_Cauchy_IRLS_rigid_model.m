function [R,t,i] = SA_Cauchy_IRLS_rigid_model(src,dst,tau)

prev_cost=10^15;                            %初始能量代价
maxIter =100;                               %最大迭代次数
n=size(src,2);                              %观测值数量
weights=ones(1,n);                          %权值向量
for i=1:maxIter
    [R,t,fit]=rigidTrans(src,dst,weights);  %根据匹配点src和dst，计算点集满足的刚体变换R,t;
    residuals = sqrt(sum((fit-dst).^2));    %计算残差
    if i==1
        alpha = max(abs(residuals));        %alpha初值设为最大残差值
    end
    cost = sum(weights.*residuals.^2);      %计算当前总能量代价w*r^2
    inl_lable = residuals<3*alpha;            %剔除部分粗差点
    src=src(:,inl_lable);                   %更新匹配点
    dst=dst(:,inl_lable);
   
    E=residuals(inl_lable);                 %更新匹配点残差
    weights=alpha^2./(alpha^2+E.^2);        %更新SA-Cauchy权值
    cost_diff = abs(cost - prev_cost);      %能量差
    alpha = alpha / tau;                    %更新alpha；
    prev_cost = cost;                       %更新能量
    
    if cost_diff < 1e-2||alpha<1
        break;
    end
end

