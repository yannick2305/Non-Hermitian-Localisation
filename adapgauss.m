%{
    ------------------------------------------------------------------------------
    Author(s):    [Vladimir ROKHLIN , Hanwen ZHANG]
    Date:         [May 2024]
    Description:  [Gaussian integration]
    ------------------------------------------------------------------------------
%}



function [rint, ier, maxrec, numint] = adapgauss(fun,a,b,torr,N)

%                        input parameters:
% 
%   a,b - the ends of the interval on which the integral is to
%        be evaluated
%   fun - the user-supplied function to be integrated. the calling
%        sequence of fun must be
%  
%         fun(x,par1,par2).                            (1)
%  
%         in (1), x is a point on the interval [a,b] where
%         the function is to be evaluated, and par1, par2
%         are two parameters to be used by fun; they can be
%         variables or arrays, real or integer, as desired.
%         fun is assumed to be real *8.
%   par1, par2 - partameters to be used by the user-supplied
%        function fun (see above)
%   m - the order of the quadrature to me used on each subinterval
%   eps - the accuracy (absolute) to which the integral will be
%        evaluated
%  
%                        output parameters:
%  
%   ier - error return code.
%           ier=0 means normal conclusion
%           ier=8 means that at some point, one subinterval in the
%                 subdivision was smaller than (b-a)/2**200. this
%                 is a fatal error.
%           ier=16 means that the total number of subintervals in the
%                 adaptive subdivision of [a,b] turned out to be greater
%                 than 100000.  this is a fatal error.
%  
%   rint - the integral as evaluated
%   maxrec - the maximum depth to which the recursion went at its
%          deepest point. can not be greater than 200, since at that
%          point ier is set to 8 and the execution of the subroutine
%          terminated.
%   numint - the totla number of intervals in the subdivision. can not
%          be greater than 100000,  since at that
%          point ier is set to 16 and the execution of the subroutine
%          terminated.



    if nargin == 4
        N = 6; % default gauss order
    end
    
    nnmax = 100*100000; % max number of interval division
    maxdepth = 200; %maxdepth of recursion allowed
    maxrec = 0; % maxdepth of recursion reached
    stack = zeros(2,200); vals = zeros(200,1);
    
    j = 1; % stack depth parameter
    rint = 0.0; % integration value initialization


    % Calculate gauss points and weights on [-1,1]
    [x0,w0] = adapgauss_legpts(N);

    stack(1,j) = a; stack(2,j) = b;
    vals(1) = adapgauss_oneint(a,b,fun,x0,w0);
    
    % start recursion
    for i = 1:nnmax

        numint = i;
        if maxrec < j
            maxrec = j;
        end

        % splitting [a,b] into [a,c] [c,a]

        c = (stack(2,j) + stack(1,j))/2;


        value2 = adapgauss_oneint(stack(1,j),c,fun,x0,w0);

        value3 = adapgauss_oneint(c,stack(2,j),fun,x0,w0);

        dd = abs(value3 + value2 - vals(j));


        % if done, move up the stack
        ifdone = 0;
        if dd < torr
            ifdone = 1;
        end

        if ifdone == 1
            rint = rint + value2 + value3;
            j = j-1;
        elseif ifdone == 0
        % if not, move down the stack
        % [a,c] at j+1, [c,b] replaces [a,b] at j
            stack(1,j+1) = stack(1,j);
            stack(2,j+1) = (stack(1,j) + stack(2,j))/2;
            vals(j+1) = value2;

            stack(1,j) = (stack(1,j) + stack(2,j))/2;
            vals(j) = value3;
            j = j + 1;
        end

        if j == 0
            ier = 0;
            break;
        end

        if j >= maxdepth
            ier = 8;
            break;
        end

        ier = 16;
    end

end
