
function s = logsumexp(x)
    %numerically stable implementation of LSE
    xstar = max(x);
    s=xstar+log(sum(exp(x-xstar)));
end
