function z = mutualInformation(x, y)
    %x=abs(x);
    %y=abs(y);
    assert(numel(x) == numel(y));
    n = numel(x);
    x = reshape(x,1,n);
    y = reshape(y,1,n);
    
    l = min(min(x),min(y));
    x = x-l+1;
    y = y-l+1;
    k = max(max(x),max(y));
    idx = 1:n;
    Mx = sparse(idx,round(x),1,round(n),round(k),round(n));
    My = sparse(idx,round(y),1,round(n),round(k),round(n));
    Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
    Hxy = -dot(Pxy,log2(Pxy+eps));

    Px = mean(Mx,1);
    Py = mean(My,1);

    % entropy of Py and Px
    Hx = -dot(Px,log2(Px+eps));
    Hy = -dot(Py,log2(Py+eps));

    % mutual information
    z = Hx + Hy - Hxy;
end