function b=indicateursARPHMM(x,minormax)

b=zeros(size(x,1),1);
for i=1:size(x,1)
    if ~all(isnan(x(i,:)))
        f=find(~isnan(x(i,:)));
        r=x(i,f);
        if strcmp(minormax,'min')
            [~,e]=min(r);
            b(i) = f(e);
        elseif strcmp(minormax,'max')
            [~,e]=max(r);
            b(i) = f(e);
        else
            error('??')
        end
    end
end