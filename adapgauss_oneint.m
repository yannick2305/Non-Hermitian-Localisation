%{
    ------------------------------------------------------------------------------
    Author(s):    [Vladimir ROKHLIN , Hanwen ZHANG]
    Date:         [April 2024]
    Description:  [Adaped Gaussian integrator]
    ------------------------------------------------------------------------------
%}
function value = adapgauss_oneint(a,b,fun,x0,w0)
    u = (b + a)/2;
    v = (b - a)/2;

    x = v*x0 + u;
    w = w0*v;

    value = w*fun(x);
end
