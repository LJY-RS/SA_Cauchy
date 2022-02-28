function [R,T]=IRLS_Welsch(x,y,W,iter)
w=W;
for i=1:iter
    [x_,R,T,W]=Welsch(x,y,2.985,W(1,:));
    if max(max(abs(w-W)))<1e-4
        break;
    end
    w=W;
end