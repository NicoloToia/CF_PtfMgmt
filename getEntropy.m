function e = getEntropy(x)
e = 0;
for i=1:length(x)
    if(x(i)<1e-6)
        
    else
        e = e + x(i)*log(x(i));
    end
end
e = -e;
end