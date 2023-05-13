function d = angdiff2(th1, th2)


th1=degtorad(th1);
th2=degtorad(th2);
    if nargin < 2

        d = th1;
    else
        d = th1 - th2;
    end

    
    d = mod(d+pi, 2*pi) - pi;
    d = radtodeg(d);