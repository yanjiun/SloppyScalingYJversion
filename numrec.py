import numpy
from numpy import dot
from numpy import transpose as T
from numpy.linalg import inv
from numpy.linalg import norm
from fdjac import Jac as fdJac
from sys import stdout

def leastsq(func,Dfunc,x0,args=(),
            full_output = 0,callback = None,
            ## Convergence Criteria
            xtol=1.49012e-08,ftol=1.49012e-08,gtol = 1.49012e-08,Cgoal = -1,
            ## Give up criteria
            maxfev = 0, maxjev = 0,
            ## Other parameters
            lam=.001,lamup = 10.,lamdown = 10.,
            epsfcn=2.2e-16,
            min_dt = 0.1, threshold_dt = 0.9, alpha = 1.,
            d=None, flags = ['noa'],verbose = False):

    #temp = Dfunc(x0,*args)
    N = len(x0)
    info = {}  #full_output Dictionary
    dinfo = {} #diag_output Dictionary
    xs = []
    x = x0.copy()
    dx = numpy.ones(x.shape)

    class Counters:
        nfev = 0
        njev = 0

        @staticmethod
        def Continue():
            if maxfev:
                if Counters.nfev >= maxfev:
                    return False
            if maxjev:
                if Counters.njev >= maxjev:
                    return False
            if not maxjev and not maxfev:
                if Counters.njev >= 100*(N+1):
                    return False
            return True

    def F(x):
        Counters.nfev += 1
        return func(x,*args)

    def Jac(x):
        Counters.njev += 1
        if Dfunc is None:
            return fdJac(x,func,args,epsfcn,False)
        else:
            return Dfunc(x,*args)

    if d is None:
        d = numpy.eye(len(x))

    steps = 0
    smalls = 0 #Counts how many small steps in a row we have taken
    new_r = F(x)
    iterations = 0
    while Counters.Continue() and not(smalls > 1):
        r = new_r
        J = Jac(x)
        gradC = dot(T(J),r)
        cost = norm(r)
        cont = True
        #Find an acceptable value of lambda
        while cont:
            iterations += 1
            g = dot(T(J),J) + lam*d
            ginv = inv(g)
            if 'noa' not in flags:
                ##Add a 2nd order term
                v = -dot(ginv,gradC)
                h = epsfcn**.25/max(norm(v),1e-8)
                rfor = F(x + h*v)
                rback = F(x - h*v)
                D2r = (rfor - 2*r + rback)/h/h
                a =- dot(D2r,dot(J,ginv))
                ## dt = max( min(1.,alpha*norm(v)/norm(a)),min_dt)
                dt = 1.0
                dx = (v + .5*a*dt)*dt
            else:
                dx = -dot(ginv,gradC)
                v = dx
                a = v*0.
                dt = 1.
            new_r = F(x+dx)
            if verbose:
                print "Iterations: ", iterations, ' Lambda is: ', lam, "Grad is:", sum(abs(gradC)), 'Cost is: ', sum(r*r)/2, sum(new_r**2)/2
                print "|v| = ", norm(v)
                print "|a| = ", norm(a)
                print "a/v = ", norm(a)/norm(v)
                print "cos phi =", dot(a,v)/norm(a)/norm(v)
                stdout.flush()

            if (norm(new_r) > cost and norm(dx) > xtol) or norm(a)/norm(v) > alpha:
                lam *= lamup
            else:
                if dt < threshold_dt:
                    lam *= lamup
                else:
                    lam /= lamdown
                cont= False
        xs.append(x.copy())
        x += dx
        steps += 1
        if callback is not None:
            callback(x)

        ## Convergence criterion
        prered = norm(r) - norm(r + dot(J,dx))
        actred = norm(r) - norm(new_r)
        conv = 0
        if abs(actred) <= ftol and abs(prered) <= ftol and .5*abs(actred/prered) <= 1.:
            conv = 1
        if norm(dx) <= xtol:
            conv = 2
        if norm(dot(r,J)) <= gtol:
            conv = 3
        if sum(r*r)/2 < Cgoal:
            conv = 4
        if conv != 0:
            smalls += 1
        else:
            smalls = 0

    info['nfev'] = Counters.nfev
    info['njev'] = Counters.njev
    info['xs'] = xs
    info['steps'] = steps
    if full_output:
        return x, info
    else:
        return x

"""
from Method import Method

numrec = Method(leastsq,"F J x",(),{},'Numerical Recipes Levenberg-Marquardt','numrec')
numrec_accel = Method(leastsq,"F J x",(),{'flags':[]},'Numerical Recipes Levenberg-Marquardt with acceleration','numreca')
numrec_dg = Method(leastsq,"F J x",(),{'lamup':2.,'lamdown':10.},'Numerical Recipes Levenberg-Marquardt with Delayed Gratification','numrecdg')
numrec_dg_accel = Method(leastsq,"F J x",(),{'flags':[],'lamup':2.,'lamdown':10.},'Numerical Recipes Levenberg-Marquardt with acceleration and Delayed Gratification','numrecdga')
"""
