import scipy

class Rosenbrock:
    def __init__(self,R):
        self.R = R

    def Cost(self,params):
        x = params[0]
        y = params[1]
        cost = (1.-x)**2 + self.R*(y-x**2)**2
        return cost

    def AnalyJac(self,params):
        x = params[0]
        y = params[1]
        jac = scipy.array([0,0])
        jac[0] = 2*(1-x)-4*self.R*x*(y-x**2)
        jac[1] = 2*self.R*(y-x**2)
        return jac

    def hess(self,params):
        x = params[0]
        y = params[1]
        hess = scipy.zeros((2,2))
        hess[0,0] = 2-4*self.R*y +12*self.R*x**2
        hess[1,1] = 2*self.R
        hess[1,0] = -4*self.R*x
        hess[0,1] = hess[1,0]
        return hess
