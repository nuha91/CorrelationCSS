function [arr] = findminimum(inarr,n)
for i=1:n
    for j=1:n
        if(i<=j)
            inarr(i,j)=inf;
        end
    end
end;
arr=zeros(1,40);
for k=1:2:39
   [r,c]=find(inarr==min(min(inarr)));
   inarr(r,c)=inf;
   arr(k)=r(1);
   arr(k+1)=c(1);
end

