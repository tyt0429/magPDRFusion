function [loocs3] = Zeros_finding_Rowdata_final(y,x,f,tol)
n = length(y);
ind = 1:(n-1);

% list of intervals with a zero crossing
k = find(((y(ind)<=0) & (y(ind+1)>0)) | ...
    ((y(ind)>=0) & (y(ind+1)<0)));

% list of zero crossings
xc = [];

% intervals where y is zero at both ends of the
% interval are exactly zero are indeterminate.
% the solution may be anywhere on that interval.
% I'll choose to return both endpoints of the
% interval.
L = (y(k)==0) & (y(k+1)==0);
if any(L)
    xc = x(k(L));
    k(L)=[];
end

% interpolate to find x. I've already removed
% the constant intervals at zero
if size(x) ~= size(y)
    x = x';
end
if ~isempty(k)
    s = (y(k+1)-y(k))./(x(k+1)-x(k));
    xc = [xc,x(k) - y(k)./s];
end

% patch for last element exactly zero
if y(end)==0
    xc(end+1) = x(end);
end
if 0 == 1
    if xc(1) == 0
        xc = xc(2:end);
    end
end
i = 1;
j = 2;
kk = 2;
xc_f = zeros();
xc_f_indx = zeros();
loocs3 = zeros();
xc_f(1) = xc(1);
xc_f_indx(1) = 1;
loocs3(1) = k(xc_f_indx(1));
while i < length(xc)
    delta = abs(xc(i)-xc(j));
    if delta < (1/((f+tol)*2))
        while delta < (1/((f+tol)*2))
            j = j + 1;
            delta = abs(xc(i)- xc(j));
        end
        xc_f(kk) = xc(j);
        xc_f_indx(kk) = j;
        loocs3(kk) = k(xc_f_indx(kk));
        i = j;
        kk = kk + 1;
    else
        xc_f(kk) = xc(j);
        xc_f_indx(kk) = j;
        loocs3(kk) = k(xc_f_indx(kk));
        i = i +1;
        j = i;
        kk = kk + 1;
    end
end
end

