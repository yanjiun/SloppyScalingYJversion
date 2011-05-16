import numpy

def Jac(x,func,args=(),eps = None,center_diff = False):
    if eps is None:
        eps = numpy.sqrt(2.2e-16)
    else:
        eps = numpy.sqrt(eps)
    dx = x*0.
    f0 = func(x,*args)
    J = numpy.empty((len(f0),len(x)))
    if not center_diff:
        for j in xrange(len(x)):
            h = eps*abs(x[j])
            if h == 0.:
                h = eps
            dx *= 0.
            dx[j] = h            
            J[:,j] = (func(x + dx,*args) - f0)/h
        return J
    else:
        for j in xrange(len(x)):
            h = eps*abs(x[j])
            if h == 0.:
                h = eps
            dx *= 0.
            dx[j] = h/2
            J[:,j] = (func(x + dx,*args) - func(x-dx,*args))/h
        return J
    
