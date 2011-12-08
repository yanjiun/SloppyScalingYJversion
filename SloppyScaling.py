import numpy
import scipy
import pylab
import copy
import time
from numpy import exp
from numpy import sin
from numpy import log

super_computer = False
if super_computer is False:
    import scipy.optimize
    import scipy.special
import WindowScalingInfo as WS
reload(WS)
import Utils
reload(Utils)
#import BoldAccel
#reload(BoldAccel)
#import Bold
#reload(Bold)
#import minpack
#reload(minpack)
#import numrec
#reload(numrec)
#import levmar_compare
#reload(levmar_compare)

f_counter=numpy.array([0])

class ScalingTheory:
    """
    A ScalingTheory's job is to provide a function Y(X) that predicts
    theory for an experiment, given a set of independent variables and
    parameters. The independent variables are those specifying the Data
    being described; the parameters describe universal critical
    exponents, universal scaling functions, and analytic
    and singular corrections to scaling.
    The theory is represented in a string consisting of a Python command.
    The variables are unpacked on the fly and the string is executed...
    For application convenience, you may use the natural variables for
    X and Y (say, 'S' and 'A') in the expressions, and set Xname and Yname
    appropriately.
    #
    Example of implementation:
    sizeHisto = ScalingTheory('S**(-tau)*numpy.exp((-(S*(R-Rc)**sigma)/XS)**nS)',
                'tau, sigma, XS, nS, Rc', (1.5,0.5,1.0,1.0,2.0),
                independentNames = 'R',
                scalingX = 'S*r**sigma', scalingY = 'D*S**tau',
                scalingXTeX = r'$S r^\sigma$',
                scalingYTeX = r'$D S^\tau$',
                title= 'Avalanche histo$'
                scalingTitle= 'Avalanche histo scaling plot'
                Xname = 'S', XscaledName='Ss', Yname = 'D', normalization = True)
    """
    def __init__(self, Ytheory, parameterNames, initialParameterValues, \
                 independentNames, \
                 scalingX = 'X', scalingY = 'Y', scalingW = None, \
                 scalingXTeX = r'${\cal{X}}$', \
                 scalingYTeX = r'${\cal{Y}}$', \
                 scalingWTeX = r'${\cal{X}}$',\
                 title = 'Fit', scalingTitle = 'Scaling Collapse',
                 Xname='X', XscaledName = 'Xs', Yname='Y', WscaledName = 'Ws',\
                 XnameTeX = r'${\cal{X}}$', YnameTeX = r'${\cal{Y}}$',\
                 fixParameter = False, fixedParameters = "", SetFixedParamsPass = False, 
                 normalization = None, priorList = None, priorCoeff=1.0, dYdp=None, Jac_dict=None):
        #YJC: added WscaledName = 'Ws' and scalingW =None to theory to make the scaling function easier to read, need to keep default none for scalingW, because some theories don't have this second variable
        #YJC: want to make the scaling variables a list?  So we can specify as many as we want
        #YJC: also added XnameTeX and YnameTeX to make axis labels prettier
        self.Ytheory = Ytheory
        self.parameterNames = parameterNames
        self.parameterNameList = parameterNames.split(",")
        self.initialParameterValues = initialParameterValues
        self.parameterNames0 = parameterNames 
        self.parameterNameList0 = self.parameterNameList
        self.initialParameterValues0 = initialParameterValues
        self.independentNames = independentNames
        self.Xname = Xname
        self.XnameTeX = XnameTeX
        self.XscaledName = XscaledName
        self.Yname = Yname
        self.YnameTeX = YnameTeX
        self.WscaledName = WscaledName
        self.scalingX = scalingX
        self.scalingY = scalingY
        self.scalingW = scalingW
        self.scalingXTeX = scalingXTeX
        self.scalingYTeX = scalingYTeX
        self.scalingWTeX = scalingWTeX
        self.title = title
        self.scalingTitle = scalingTitle
        self.normalization = normalization
        self.fixParameter = fixParameter
        self.SetFixedParamsPass = SetFixedParamsPass
        self.priorList = priorList
        self.priorCoeff = priorCoeff
        self.dYdp = dYdp
        self.Jac_dict=Jac_dict

    def Y(self, X, parameterValues, independentValues):
        """
        Predicts Y as a function of X
        """
        # Set values of parameters based on vector of current guess
        # Set values of independent variables based on which curve is being fit
        # Set up vector of independent variable from X
        # Warning: local variables in subroutine must be named
        # 'parameterValues', 'independentValues', and 'X'
        exec(self.parameterNames + " = parameterValues")
        exec(self.independentNames + " = independentValues")
        #if len(self.fixedParamNamesList) > 0:
        #    exec(self.fixedParamNames + "= self.fixedParamValues")
        if self.fixParameter:
            for par, val in self.fixedParameters:
                exec(par + " = " + str(val))
        exec(self.Xname + ' = X')
        if self.XscaledName:
            exec(self.XscaledName +'='+ self.scalingX)
        #YJC: added scalingW here
        if self.scalingW is not None:
            exec(self.WscaledName +"=" + self.scalingW)
        exec("Y = " + self.Ytheory)
        if self.normalization:
            fn = getattr(self, self.normalization)
            Y = fn(X, Y, parameterValues, independentValues)
        return Y

    def ScaleX(self, X, parameterValues, independentValues):
        """
        Rescales X according to scaling form
        """
        # Set values of parameters, independent variables, and X vector
        # Warning: local variables in subroutine must be named
        # 'parameterValues', 'independentValues', and 'X'
        exec(self.parameterNames0 + " = parameterValues")
        exec(self.independentNames + " = independentValues")
        #if len(self.fixedParamNamesList) > 0:
        #    exec(self.fixedParamNames + "= self.fixedParamValues")
        if self.fixParameter:
            for par, val in self.fixedParameters:
                exec(par + " = " + str(val))
        exec(self.Xname + " = X")
        #YJC: added scalingW here too
        if self.scalingW is not None:
            exec(self.WscaledName + '=' +self.scalingW)
        exec("XScale = " + self.scalingX)
        return XScale

    def ScaleW(self, X, parameterValues,independentValues):
        """
        rescales X acording to Wscaled form
        """
        if self.scalingW is not None:
            exec(self.parameterNames0 + " = parameterValues")
            exec(self.independentNames + " = independentValues")
            if self.fixParameter:
                for par, val in self.fixedParameters:
                    exec(par + " = " + str(val))
            exec(self.Xname + " = X")
            exec("Wscale = " + self.scalingW)
        else:
            print "WARNING: 2nd scaling variable does not exist, please define"
            exec("Wscale = X")

        return Wscale


    def ScaleY(self, X, Y, parameterValues, independentValues):
        """
        Rescales Y according to form
        """
        # Set values of parameters, independent variables, and X vector
        # Warning: local variables in subroutine must be named
        # 'parameterValues', 'independentValues', and 'X'
        exec(self.parameterNames + " = parameterValues")
        exec(self.independentNames + " = independentValues")
        #if len(self.fixedParamNamesList) > 0:
        #    exec(self.fixedParamNames + "= self.fixedParamValues")
        if self.fixParameter:
            for par, val in self.fixedParameters:
                exec(par + " = " + str(val))
        exec(self.Xname + " = X")
        #YJC: added scalingW here too
        if self.scalingW is not None:
            exec(self.WscaledName + '=' +self.scalingW)
        if self.XscaledName:
            exec(self.XscaledName + "="+self.scalingX)
        exec(self.Yname + " = Y")
        exec("YScale = " + self.scalingY)
        return YScale

    def SetFixedParams(self, fixedParameters=None):
        """
        Sets parameters to fixed values.
        fixedParameters is a list of tuple(s) of the type:
        [('param1', val1), ('param2', val2)]
        To "unfix" parameters, call with fixedParameters=None
        """
        # Set first the fixParameter to true and
        # pass the parameters to fix
        
        #if not self.fixParameter:
        #    self.fixParameter = True
        
        # Then adjust the Names and the Values
        
        if fixedParameters is not None:
            self.fixParameter = True
            pNames, pValues = \
                    Utils.reduceParameters(self.parameterNames0, \
                                           self.initialParameterValues0, \
                                           fixedParameters)
            self.parameterNames = pNames
            self.parameterNameList = pNames.split(",")
            self.initialParameterValues = pValues
            self.fixedParameters = fixedParameters
            self.fixParameter = True
            self.SetFixedParamsPass = True
        else:
            if self.SetFixedParamsPass:
                self.parameterNames=self.parameterNames0
                self.parameterNameList = self.parameterNameList0
                self.initialParameterValues = self.initialParameterValues0
                self.fixParameter = False
                self.fixedParameters = None

    def Prior(self, parameterValues, parameterNameList):
        """
        reads the priorlist to add to the residuals,
        if prior is None, return empty list []
        """
        exec(self.parameterNames +'= parameterValues')
        prior = []
        if self.priorList:
            priorList = self.priorList
            for p in priorList:
                # new version for more flexible priors
                exec('prior.append('+p+')')
                #
                #below is old version of code
                #
                #index = parameterNameList.index(p)
                #val = parameterValues[index]
                #val = self.priorCoeff*val
                #prior.append(val)
            prior=scipy.array(prior)*self.priorCoeff
        return prior

    #
    # Various options for normalization
    #
    def NormBasic(self, X, Y, parameterValues, independentValues):
        """
        Must guess at bin sizes for first and last bins
        """
        norm = Y[0]*(X[1]-X[0])
        norm += sum(Y[1:-1] * (X[2:]-X[:-2])/2.0)
        # GF: Why not this below?
        #norm += sum(Y[1:-2] * (X[2:-1]-X[:-3])/2.0)
        norm += Y[-1]*(X[-1]-X[-2])
        return Y/norm
    
    def NormIntegerSum(self, X, Y, parameterValues, independentValues, \
                xStart=1., xEnd=1024.):
        """
        Function summed over positive integers equals one; brute force
        up to xEnd
        """
        x = numpy.arange(xStart, xEnd)
        return Y/sum(self.Y(x, parameterValues, independentValues))
        
    def NormLog(self,X,Y,parameterValues, independentValues):
        """
        This kind of normalization is correct
        if the data are uniform in log scale,
        as prepared by our code toBinDistributions.py
        """
        lgX = numpy.log10(X)
        D = numpy.around(lgX[1] - lgX[0],2)
        bins = 10**(lgX+D/2.) - 10**(lgX-D/2.)
        return Y/sum(Y*bins)

    def DY(self,X,parameterValues,independentValues):
        """
        to put in dervatives for the Jacobian (need to combine for all data sets)
        """
        jac=numpy.array([])
        if self.dYdp is not None: 
            dYdp = self.dYdp
            exec(self.parameterNames + "=parameterValues")
            exec(self.independentNames + "= independentValues")
            exec(self.Xname+"=X")
            if self.WscaledName:
                exec(self.WscaledName + '=' +self.scalingW)
            if self.XscaledName:
                exec(self.XscaledName + "="+self.scalingX)
            if self.Jac_dict:#this is for the common factors to precalculate
                for name in self.Jac_dict.keys():
                    if name!='prior' and name!='priorvalues':
                        exec(name + "=" +  self.Jac_dict[name])      
            exec("jac=numpy.array("+dYdp+")")
        return jac


class Data:
    """
    A Data object contains a series of curves each for a set of independent
    control parameters. For example, the X values might be avalanche sizes
    (Xname = 'S'), the Y values might be percentage area covered by
    avalalanches of that size (Yname = 'A'),
    the sigmas the standard errors in the means, and an independent control
    parameters might be the demagnetizing field (independent = 'k'). If,
    as for A(S), the data plots are typically log-log set self.linlog = 'log';
    for things like V(t,T) set self.linlog = 'lin'.
    """ 
    def __init__(self, linlog = 'log'):
        self.experiments = []
        self.X = {}
        self.Y = {}
        self.linlog = linlog
        self.pointType = {}
        self.errorBar = {}
        self.fileNames = {}
        self.defaultFractionalError = {}
        self.initialSkip = {}
        self.finalSkip = {}

    def InstallCurve(self, independent, fileName, defaultFractionalError = 0.01,\
                     pointSymbol="o", pointColor="b", utilizeDefaultError=True,\
                     xCol=0, yCol=1, errorCol = 0, initialSkip = 0, finalSkip=0, factorError = 1.0):
        """
        Curves for independent control parameters given by "independent"
        loaded from "fileName". Plots use, for example, pointSymbol from 
            ['o','^','v','<','>','s', 'p', 'h','+','x']
        and pointColor from 
            ['b','g','r','c','m','burlywood','chartreuse']
        
        factorError is to artificially increase error bars for better fits
        """
        
        # check if independent is a tuple
        if not isinstance(independent, tuple):
            print "Warning: the independent variable is not a tuple"
            independent = tuple(independent)
        #
        self.experiments.append(independent)
        self.fileNames[independent] = fileName
        self.initialSkip[independent] = initialSkip
        self.finalSkip[independent]=finalSkip
        self.pointType[independent] = pointColor + pointSymbol
        self.defaultFractionalError[independent] = defaultFractionalError
        try:
            infile = open(fileName, 'r')
            lines = infile.readlines()
            infile.close()
            success = 1
            numbers = [line.split() for line in lines]
            self.X[independent] = numpy.array( \
                        [float(line[xCol]) for line in numbers])
            self.Y[independent] = numpy.array( \
                        [float(line[yCol]) for line in numbers])
            if errorCol:
                self.errorBar[independent] =  \
                        numpy.array([float(line[errorCol])*factorError for line in numbers])
                if utilizeDefaultError:
                    self.errorBar[independent]=numpy.array([max(e, emin) for (e,emin) in zip(self.errorBar[independent],self.Y[independent]*defaultFractionalError)])
            else:
                self.errorBar[independent] = \
                    self.Y[independent] * defaultFractionalError
        except IOError:
            print "File %s not found"%fileName
            success = 0
        except ValueError:
            print independent, fileName
            success = 0
        return success

class Model:
    """
    A Model object unites Theory with Data. It's primary task is to 
    calculate the residuals (the difference between theory and data)
    and the cost.
    """
    def __init__(self, theory, data, name, sorting):
        self.theory = theory
        self.data = data
        self.name = name
        self.sorting = sorting

    def AnalyJac(self, parameterValues, parameterNames=None):
        """
        uses the dYdp function in the theory class to calculate the Jacobian for these data sets
        (should weight them by errorbars)
        """
        model_num=1 #counter
        length =0 
        for independentValues in self.data.experiments:
            initialSkip = self.data.initialSkip[independentValues]
            finalSkip = self.data.finalSkip[independentValues]
            if finalSkip !=0:
                X = self.data.X[independentValues][initialSkip:int(-1*finalSkip)]
                errorBar = self.data.errorBar[independentValues][initialSkip:int(-1*finalSkip)]
            else:
                X = self.data.X[independentValues][initialSkip:]
                errorBar = self.data.errorBar[independentValues][initialSkip:]

            new_jac=self.theory.DY(X,parameterValues,independentValues)/errorBar
            new_jac=new_jac.transpose()

            if model_num==1:
                jac = new_jac
                model_num += 1
            else:
                jac = numpy.concatenate((jac,new_jac))

            #if there are fixed parameters, strike from jacobian
            #if self.theory.fixedParameters is not None:
            #    fixedParameters=self.theory.fixedParameters
            #    for fixedP in fixedParameters:
            #        index = self.theory.parameterNameList0.index(fixedP)
            #        jac.


            # last add priors (!!! these only now work for priors of the form n^2)
            # need to add flexible priors (take derivatives of priors and specify a new list?...)
        
        if self.theory.priorList:
            if parameterNames is None:
                parameterNames=self.theory.parameterNameList0
            exec(self.theory.parameterNames +"=parameterValues")
            priorList = self.theory.priorList[:]
            priorJac = numpy.zeros((len(priorList),len(parameterNames)))
            for pname in self.theory.Jac_dict['prior'].keys():
                i = priorList.index(self.theory.Jac_dict['priorvalues'][pname])
                p_index=parameterNames.index(pname)
                exec("temp_value="+self.theory.Jac_dict['prior'][pname])
                priorJac[i,p_index] = temp_value*self.theory.priorCoeff
            jac =numpy.concatenate((jac,priorJac))       
        return jac


    def Residual(self, parameterValues, dictResidual=False):
        """
        Calculate the weighted residuals,
        with the weights = 1 / errorbar
        """
        #YJC: want to add priors as extra rows to the residuals
        # put in simple prior on all parameters that start with U
        # we will want to define a more flexible function that 

        f_counter[0]+=1

        if dictResidual:
            residuals = {}
        else:
            residuals = numpy.array([])
            
        for independentValues in self.data.experiments:
            initialSkip = self.data.initialSkip[independentValues]
            finalSkip = self.data.finalSkip[independentValues]
            if finalSkip !=0:
                X = self.data.X[independentValues][initialSkip:int(-1*finalSkip)]
                Y = self.data.Y[independentValues][initialSkip:int(-1*finalSkip)]
                errorBar = self.data.errorBar[independentValues][initialSkip:int(-1*finalSkip)]
            else:
                X = self.data.X[independentValues][initialSkip:]
                Y = self.data.Y[independentValues][initialSkip:]
                errorBar = self.data.errorBar[independentValues][initialSkip:]
            Ytheory = self.theory.Y(X, parameterValues, independentValues)
            res = (Ytheory-Y)/errorBar
            if max(res) == numpy.inf:
                print "data set has infinite residual", independentValues
            if dictResidual:
                residuals[independentValues] = res
            else:
                residuals = numpy.concatenate((residuals,res))
                    
        #YJC: put this in joint Model too!
        priors = self.theory.Prior(parameterValues, self.theory.parameterNameList)

        #for p in pNames:
        #    if p.startswith('U'):
        #        index = pNames.index(p)
        #        val = pValues[index]
        #        priors.append(val)
        #priors = numpy.array(priors)
        if dictResidual is False:
            residuals = numpy.concatenate((residuals,priors))
        return residuals
        
    def Cost(self, parameterValues=None):
        """
        Sum of the squares of the residuals
        """
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        residuals = self.Residual(parameterValues)
        return sum(residuals*residuals)
    
    def SST(self, parameterValues=None):
        """
        SST is the sum of the squares about the mean
        """
        sst = 0.
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        for independentValues in self.data.experiments:
            initialSkip = self.data.initialSkip[independentValues]
            finalSkip = self.data.finalSkip[independentValues]
            if finalSkip !=0:
                Y = self.data.Y[independentValues][initialSkip:int(-1*finalSkip)]
                errorBar = self.data.errorBar[independentValues][initialSkip:int(-1*finalSkip)]
            else:
                Y = self.data.Y[independentValues][initialSkip:]
                errorBar = self.data.errorBar[independentValues][initialSkip:]
            sst_partial = (Y-numpy.mean(Y))/errorBar
            sst += sum(sst_partial*sst_partial)
        return sst
    
    def R_square(self,parameterValues=None):
        """
        Calculates the R-square = 1 - cost / SST
        where SST is the sum of the squares about the mean
        """
        sst = self.SST(parameterValues)
        cost = self.Cost(parameterValues)
        return 1.- cost/sst
        
    def getLabel(self, names, values, withRescale = False, pow10first=False, sigma = 0.45):
        """
        Get the Labels to be plotted. 
        """
        #lb_name = (names[-1] ==  ',') and names[:-1] or names[-1]
        lb = names + " = "
        lb += ",".join([str(i) for i in values])
        if len(values)==2:
            L, k = values
        if len(values)==3:
            L, k, W = values

        if withRescale:
            if len(values)==2:
                lb = names + "="
                lb += str(values[0])
                if pow10first:
                    lb += r", $10^{%d}$" %(int(round(numpy.log10(1.0*k/L))))
                else:
                    lb += ", %.3e" %(1.0*k/L)
            if len(values)==3:
                lb = r"$k, W_s =$"
                #lb += str(L)
                if pow10first:
                    lb += r"$10^{%d}$" %(int(round(numpy.log10(1.0*k/L))))
                else:
                    lb += "%.2e" %(1.0*k/L)
                lb += ",%.2f" %(W*(1.0*k/L)**sigma)
                #lb += str(W)
            #lb += ",%.3e" %((1.0*k/L)**(-sigma)/L)
            #for nm, val in zip(a,b):
            #    exec(nm + "= " + str(val))
            #if len(values) == 2:
            #    lb += str(1.0*k/L)**sigma
            #elif len(values) == 3:
            #    lb += str((1.0*k/L)**sigma*W)[0:5]
        return lb
    
    def getAxis(self,X,Y):
        """
        return the proper axis limits for the plots
        """
        out = []
        mM = [(min(X),max(X)),(min(Y),max(Y))]
        for i,j in mM:
            #YJC: checking if values are negative, if yes, return 0 and break
            if j <0 or i <0:
                return 0
            log_i = numpy.log10(i)
            d, I = numpy.modf(log_i)
            if log_i < 0:
                add = 0.5 *(numpy.absolute(d)<0.5)
            else:
                add = 0.5 *(numpy.absolute(d)>0.5)
            m = numpy.floor(log_i) + add
            out.append(10**m)
            log_j = numpy.log10(j)
            d, I = numpy.modf(log_j)
            if log_j < 0:
                add = - 0.5 *(numpy.absolute(d)>0.5)
            else:
                add = - 0.5 *(numpy.absolute(d)<0.5)
            m = numpy.ceil(log_j) + add
            out.append(10**m)
        return tuple(out)
        
    def PlotFunctions(self, parameterValues=None, plotCollapse = False, 
                fontSizeLabels = 24, fontSizeLegend=12, pylabLegendLoc=0, interactive=True, pow10first=False,plotPoints=100,theory=True, collapseAxis='Xscaled'):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        # XXX Having problems with pylab.ioff(), ipython command shell freezes
        #if not interactive:
        #    pylab.ioff()
        #    pylab.clf()
        ax0 = [1.e99,0,1.e99,0]
        if self.data.linlog == 'log':
            minY = 1.e99
            for independentValues in self.data.experiments:
                Y = self.data.Y[independentValues]
                if plotCollapse:
                    X = self.data.X[independentValues]
                    Y = self.theory.ScaleY(X,Y,parameterValues,\
                                           independentValues)
                minY = min(minY,min(Y))
        
        #pylab.plot([],label=r'$win (k/L)^{\sigma_k \zeta}$')
        
        if self.sorting:
            # preserve order of values as provided
            # by Utils.get_independent
            data_experiments = self.data.experiments
        else:
            # set sorted 
            data_experiments = sorted(self.data.experiments)
            
        for independentValues in data_experiments:
            X = self.data.X[independentValues]
            Y = self.data.Y[independentValues]
            if self.data.linlog == 'log':
                Xtheory = scipy.logspace(scipy.log10(min(X)),scipy.log10(max(X)),num=plotPoints)
            elif self.data.linlog == 'lin':
                Xtheory = scipy.linspace(scipy.log10(min(X)),scipy.log10(max(X)),num=plotPoints)
            else:
                Xtheory = X
            Ytheory = self.theory.Y(Xtheory, parameterValues, independentValues)
            pointType = self.data.pointType[independentValues]
            errorBar = self.data.errorBar[independentValues]
            
            if plotCollapse:
                # Scaled error bars and Y need not-rescaled X
                errorBar = self.theory.ScaleY(X, errorBar, parameterValues, \
                                              independentValues)
                Y = self.theory.ScaleY(X, Y, parameterValues, independentValues)
                Ytheory = self.theory.ScaleY(Xtheory, Ytheory, \
                                             parameterValues, independentValues)
                # Then rescale X too
                if collapseAxis=='Xscaled':
                    X = self.theory.ScaleX(X, parameterValues, independentValues)
                    Xtheory = self.theory.ScaleX(Xtheory,parameterValues,independentValues)
                elif collapseAxis=='Wscaled':
                    X = self.theory.ScaleW(X, parameterValues, independentValues)
                    Xtheory = self.theory.ScaleW(Xtheory,parameterValues,independentValues)
                else:
                    print "please specify which axis to collapse, we are plotting default x now"
            # Avoid error bars crossing zero on log-log plots
            if self.data.linlog == 'log':
                errorBarDown = errorBar * (errorBar < Y) + (Y -minY) * (errorBar > Y)
                y_error=[errorBarDown,errorBar]
            else:
                y_error=errorBar
                
            # Prepare the labels
            lb = self.getLabel(self.theory.independentNames, independentValues,pow10first=pow10first)
            #pylab.rcParams.update({'legend.fontsize':fontSizeLegend})
            pylab.rcParams.update({'xtick.labelsize':24,
                               'xtick.major.size':20,
                               'xtick.minor.size':10,
                               'ytick.labelsize':24,
                               'ytick.major.size':20,
                               'ytick.minor.size':10,
                               'lines.markersize':10,
                               'axes.labelsize':24,\
                               'legend.fontsize':20,
                               'legend.columnspacing':1.5,
                               'figure.figsize':[10.,8.],\
                               'text.usetex':False,
                               })
            # Note use figure.figsize:[10,10] for legend on bottom plots
            # [10,8] for regular plots
            font = {'family':'serif',
                    'serif':'Times New Roman'}
            pylab.rc('font',**font)
            pylab.axes([0.15,0.15,0.95-0.15,0.90-0.125])
            #pylab.axes([0.15,0.35,0.95-0.15,0.95-0.35]) # for legend on bottom
            #####################
            if self.data.linlog == 'log' or self.data.linlog == 'lin':
                if self.data.linlog == 'log':
                    plot_fn = getattr(pylab,'loglog')
                elif self.data.linlog == 'lin':
                    plot_fn = getattr(pylab,'plot')
                # Plot first data with their error
                plot_fn(X,Y,pointType[1])
                pylab.errorbar(X,Y, yerr=y_error, fmt=pointType,label=lb)
                axis_dep = self.getAxis(X,Y)
                # Get the current values of the axis
                # YJC: some values of binned data are negative, modified getAxis to check, and return 0 if negative values encountered 
                if axis_dep ==0:
                    print "this data set has negative values", independentValues
                    print "\n"
                else:
                    for i, Ax in enumerate(axis_dep):
                        ax0[i] = i%2 and max(ax0[i],Ax) or min(ax0[i],Ax)
                # Plot the theory function
                if theory:
                    plot_fn(Xtheory,Ytheory,pointType[0])
            else:
                print "Format " + self.data.linlog + \
                        " not supported yet in PlotFits"

        #print "range of axes are:", tuple(ax0)
        pylab.axis(tuple(ax0))
        pylab.legend(loc=pylabLegendLoc, ncol=2)
        #pylab.legend(loc=(-0.16,-0.52),ncol=3)
        if plotCollapse:
            if collapseAxis =='Xscaled':
                pylab.xlabel(self.theory.scalingXTeX, fontsize=fontSizeLabels)
            elif collapseAxis == "Wscaled":
                pylab.xlabel(self.theory.scalingWTeX, fontsize=fontSizeLabels)
            pylab.ylabel(self.theory.scalingYTeX, fontsize=fontSizeLabels)
            #pylab.title(self.theory.scalingTitle)
        else:
            pylab.xlabel(self.theory.XnameTeX, fontsize=fontSizeLabels)
            pylab.ylabel(self.theory.YnameTeX, fontsize=fontSizeLabels)
            #pylab.title(self.theory.title, fontsize=fontSizeLabels)
        # XXX Turn on if ioff used pylab.ion()
        #pylab.ion()
        #pylab.show()
        
    def PlotResiduals(self, parameterValues=None, \
                      fontSizeLabels = 18, pylabLegendLoc=0,pow10first=False):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        #pylab.ioff()
        #pylab.clf()
        residuals = self.Residual(parameterValues, dictResidual=True)
        x0 = 0
        pylab.rcParams.update({'figure.figsize':[10.,8.]})
        pylab.axes([0.15,0.15,0.95-0.15,0.90-0.125]) 
        for independentValues in self.data.experiments:
        #for independentValues in sorted(residuals):  
            res = residuals[independentValues]
            xStep = len(res)
            x = range(x0,x0+xStep)
            x0 += xStep
            pointType = self.data.pointType[independentValues]
            lb = self.getLabel(self.theory.independentNames, independentValues,pow10first=pow10first)
            pylab.plot(x,res,pointType, label=lb)
        pylab.ylabel("Weighted residuals")
        pylab.axhline(y=0,color='k')
        pylab.legend(loc=pylabLegendLoc, ncol=1)
        #pylab.ion()
        #pylab.show()

    def PlotEnsemblePredictions(self,ens, sampling_freq=10, fontSizeLabels=24, fontSizeLegend=20,rescale=1.,pow10first=False):
        """
        This function plots the fit and the range of predictions given by the parameter ensemble input.
        The sampling_freq option avoids sampling all of the ensemble points for faster plotting
        Rescale is an option for when the sampling is done at a lower temperature. 
        """
        pylab.figure()
	pylab.axes([0.15,0.35,0.95-0.15,0.95-0.35])
	data_experiments = self.data.experiments
        for independentValues in data_experiments:
            Xdata = self.data.X[independentValues]
            Ydata = self.data.Y[independentValues]
	    Xtheory = scipy.logspace(scipy.log10(min(Xdata)),scipy.log10(max(Xdata)),num=100)
	    pointType = self.data.pointType[independentValues]
	    errorBar = self.data.errorBar[independentValues]
	    mean_theory = scipy.zeros(len(Xtheory))
	    std_theory = scipy.zeros(len(Xtheory))
            for i in range(0, num_to_count):
		    ens_theory = self.theory.Y(Xtheory,ens[i*sampling_freq],independentValues)
		    mean_theory += ens_theory
		    std_theory += (ens_theory)**2
	    mean_theory = mean_theory/(1.0*num_to_count)
	    std_theory = scipy.sqrt((std_theory-num_to_count*mean_theory**2)/(num_to_count-1.))
	    pylab.loglog(Xdata,Ydata,pointType[1])
	    lb = self.getLabel(self.theory.independentNames,independentValues,pow10first=pow10first)
	    pylab.errorbar(Xdata,Ydata, yerr=errorBar, fmt=pointType,label=lb)
	    pylab.loglog(Xtheory,mean_theory,pointType[0])
	    axis_dep=self.getAxis(Xdata,Ydata)
	    #upper_bound = mean_theory+std_theory
	    #lower_bound = mean_theory-std_theory
	    upper_bound = scipy.exp(scipy.log(mean_theory) + scipy.log(1.+std_theory/mean_theory)*rescale)
	    lower_bound = scipy.exp(scipy.log(mean_theory)+scipy.log(1.-std_theory/mean_theory)*rescale)
	    for i in range(0, len(lower_bound)):
		    if lower_bound[i]<=0:
			    lower_bound[i]=10.**(-16)
	    pylab.fill_between(Xtheory,lower_bound,y2=upper_bound,color=pointType[0],alpha=0.2)

	    for i, Ax in enumerate(axis_dep):
		    ax0[i] =i%2 and max(ax0[i],Ax) or min(ax0[i],Ax)
	pylab.axis(tuple(ax0))
	pylab.legend(loc=(-0.15,-0.52),ncol=3)
        
    def BestFit(self,initialParameterValues = None, method = None, fixedParams=None):
        """
        uses various algorithms to fit scaling functions to data
        inputs:
        method: specifies different algorithms
                'lm' or None: minpack LM
                'lm_accel': LM + acceleration
                'bold': Cyrus' algorithm
                'boldAccel': Cyrus' algorithm +acceleration
        outputs:
        full parameter list, fitting algorithm output
        """

        if fixedParams:
            if not isinstance(fixedParams, list):
                fixedParams=[fixedParams]
            #Check now if the name is correct
            l_index=[]
            for index, par in enumerate(fixedParams):
                pName, pValue = par
                if pName not in self.theory.parameterNameList0:
                    print "%s is not a valid name. Ignored" %pName
                    l_index.append(index)
            if l_index:
                for i in l_index:
                    fixedParams.pop(i)

        self.theory.SetFixedParams(fixedParams)

        if initialParameterValues is None:
            initialParameterValues = self.theory.initialParameterValues
        #d = numpy.ones(len(initialParamaeterValues))
        start_time = time.time()
        if method is None or method == 'lm':
            out = scipy.optimize.minpack.leastsq(self.Residual,initialParameterValues,full_output=1, ftol=1.e-16)
        elif method == 'boldAccel':
            initialParameterValues=numpy.array(initialParameterValues)
            out = BoldAccel.leastsq(self.Residual,None,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,ibold=0,verbose=True)
        elif method == 'bold':
            initialParameterValues = numpy.array(initialParameterValues)
            out = Bold.leastsq(self.Residual,None,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,ibold=0,verbose=True)
        #out = minpack.leastsq(self.Residual,self.AnalyJac,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,Cgoal=4.e04)
        elif method == 'lm_accel':
            initialParameterValues=numpy.array(initialParameterValues)
            out = numrec.leastsq(self.Residual,self.AnalyJac,initialParameterValues,full_output=1,verbose=True, flags=[],maxfev=500)
        else:
            print "fitting method is not included"
            out = None
        end_time = time.time()
        print "fitting took (mins)", (end_time-start_time)/60.
        print "number of function evals:", f_counter
        
        if fixedParams:
            outputParameterValues = self.MergeFixedAndVariableParams(fixedParams,out[0])
            self.theory.SetFixedParams()
        else:
            outputParameterValues = out[0]


        return outputParameterValues, out
    
    def PlotBestFit(self, initialParameterValues = None, \
                    figFit = 1, figCollapse=2, figResiduals=3, fontSizeLabels=12, fixedParams = None, pow10first=False):
        
        #YJC: added abilitiy to set fixedParams for this, and also modified output to match what is done in composite theory

        if fixedParams:
            if not isinstance(fixedParams, list):
                fixedParams=[fixedParams]
            #Check now if the name is correct
            l_index=[]
            for index, par in enumerate(fixedParams):
                pName, pValue = par
                if pName not in self.theory.parameterNameList0:
                    print "%s is not a valid name. Ignored" %pName
                    l_index.append(index)
            if l_index:
                for i in l_index:
                    fixedParams.pop(i)

        # Call setfixedParams even if fixedParams=None to 
        # check if original Names and values have to be used

        self.theory.SetFixedParams(fixedParams)
   
        if initialParameterValues is None:
            initialParameterValues = self.theory.initialParameterValues
        
        print 'initial cost = ', self.Cost(initialParameterValues)
        optimizedParameterValues = self.BestFit(initialParameterValues)[1][0]
        covar = self.BestFit(initialParameterValues)[1][1]

        #YJC: if covariance matrix is singular still proceed
        if covar is None:
            print "covariance matrix is singular in some direction"
            errors = None
        else:
            errors = [covar[i,i]**0.5 for i in range(len(covar))]

        print 'optimized cost = ', self.Cost(optimizedParameterValues)
        #print 'R-value = ', self.R_square(optimizedParameterValues)
        
        if fixedParams:
            print "====== Fixed Parameters ======"
            for pName, pValue in fixedParams:
                print "%s = %2.2f" %(pName, pValue)

        if errors:
            print "====== Fitted Parameters (with one sigma error) ======"
            for name, val, error in \
                zip(self.theory.parameterNameList,optimizedParameterValues, errors):
                print name + "= %2.6f +/- %2.6f" %(val, error)
            print "======================================================"
        else:
            print "====== Fitted Parameters (error unavailable) ========="
            for name, val in \
                zip(self.theory.parameterNameList,optimizedParameterValues):
                print name + "= %2.6f" %val
            print "======================================================"

        pylab.figure(figFit)
        self.PlotFunctions(optimizedParameterValues, interactive=False, pow10first=pow10first)
        pylab.figure(figCollapse)
        self.PlotFunctions(optimizedParameterValues, plotCollapse = True, interactive=False,pow10first=pow10first)
        pylab.figure(figResiduals)
        self.PlotResiduals(optimizedParameterValues, pow10first=pow10first)

        #YJC: need to call these twice to show figures properly
        pylab.figure(figFit)
        pylab.figure(figCollapse)
        pylab.figure(figResiduals)

        #YJC: want to include the fixed Parameters in the output
        if fixedParams:
            outputParameterValues = self.MergeFixedAndVariableParams(fixedParams,optimizedParameterValues)
            self.SetFixedParams()
            #allParameterNames=self.theory.parameterNameList0
            #fittedParameterNames=self.theory.parameterNameList
            #outputParameterValues=numpy.zeros(len(allParameterNames))
            #for i in range(0,len(fittedParameterNames)):
            #    p_index=allParameterNames.index(fittedParameterNames[i])
            #    outputParameterValues[p_index]=optimizedParameterValues[i]
            #for pName, pValue in fixedParams:
            #    p_index =allParameterNames.index(pName)
            #    outputParameterValues[p_index]=pValue
        else:
            outputParameterValues=optimizedParameterValues

        return outputParameterValues

    def MergeFixedAndVariableParams(self,fixedParams,optimizedParameterValues):
        allParameterNames=self.theory.parameterNameList0
        fittedParameterNames=self.theory.parameterNameList
        outputParameterValues=numpy.zeros(len(allParameterNames))
        for i in range(0,len(fittedParameterNames)):
            p_index=allParameterNames.index(fittedParameterNames[i])
            outputParameterValues[p_index]=optimizedParameterValues[i]
        for pName, pValue in fixedParams:
            p_index =allParameterNames.index(pName)
            outputParameterValues[p_index]=pValue

        return outputParameterValues

class CompositeModel:
    """
    Class combining several Models into one
    The main job of CompositeModel is to combine the parameter lists and
    initial values into a single structure, and then to impose that structure
    on the individual theories.
    Also, plots and stuff should be delegated to the individual theories.
    """
    class CompositeTheory:
        def __init__(self):
            self.parameterNames = ""
            self.initialParameterValues = []
            self.parameterNameList = []
            
    def __init__(self, name):
        self.Models = {}
        self.theory = self.CompositeTheory()
        self.name = name
        self.SetFixedParamsPass = False
        self.OriginalParams={} #to store the original parameters associated with each model
        
    def InstallModel(self,modelName, model):
        self.Models[modelName] = model
        self.OriginalParams[model]=model.theory.parameterNameList[:]
        th = self.theory
        for param, init in zip(model.theory.parameterNameList, \
                                model.theory.initialParameterValues):
            if param not in th.parameterNameList:
                th.parameterNameList.append(param)
                th.initialParameterValues.append(init)
            else:
                # Check if shared param has consistent initial value
                # between models
                paramCurrentIndex = th.parameterNameList.index(param)
                paramCurrentInitialValue = \
                        th.initialParameterValues[paramCurrentIndex]
                if paramCurrentInitialValue != init:
                    print "Initial value %f"%(init,) \
                     + " for parameter " + param + " in model " + modelName \
                     + " \n disagrees with value %f"%(paramCurrentInitialValue)\
                     + " already stored for previous theory in " \
                     + " CompositeTheory. Ignoring new value."
                    
        th.parameterNames = ",".join(th.parameterNameList)
        #th.initialParameterValues =tuple(th.initialParameterValues)
        #
        # Update list of parameter names and values for all attached models
        #
        for currentModel in self.Models.values():
            currentModel.theory.parameterNames=th.parameterNames
            currentModel.theory.parameterNames0=th.parameterNames
            currentModel.theory.parameterNameList=th.parameterNameList
            currentModel.theory.parameterNameList0=th.parameterNameList
            currentModel.theory.initialParameterValues=tuple(th.initialParameterValues)
            currentModel.theory.initialParameterValues0=tuple(th.initialParameterValues)
        #
        # Remember original Names and values
        th.initialParameterValues0 = copy.copy(th.initialParameterValues)
        th.parameterNames0 = copy.copy(th.parameterNames)
        th.parameterNameList0 = copy.copy(th.parameterNameList)
        
    def SetFixedParams(self, fixedParameters):
        """
        Sets parameters in fixedParamNames to their initial values,
        and updates the parameter values, names of the composite model
        fixedParameters is a list of tuple(s) of the type: [('par1', val1)]
        """
        th = self.theory
        if fixedParameters:
            pNames, pValues = Utils.reduceParameters(th.parameterNames0,\
                                                 th.initialParameterValues0,\
                                                 fixedParameters)
            th.parameterNames = pNames
            th.parameterNameList = pNames.split(",")
            th.initialParameterValues = pValues
            for currentModel in self.Models.values():
                currentModel.theory.fixParameter = True
                currentModel.theory.SetFixedParams(fixedParameters)
            self.SetFixedParamsPass = True
        else:
            if self.SetFixedParamsPass:
                th.parameterNames = th.parameterNames0
                th.parameterNameList = th.parameterNameList0
                th.initialParameterValues = th.initialParameterValues0
                for currentModel in self.Models.values():
                    currentModel.theory.parameterNames=th.parameterNames0
                    currentModel.theory.parameterNameList=th.parameterNameList0
                    currentModel.theory.initialParameterValues=th.initialParameterValues0
                    currentModel.theory.fixParameter = False
                    currentModel.theory.fixedParameters = None
            
    def Residual(self, parameterValues):
        residuals = numpy.array([])
        for model in self.Models.values():
            modelResidual = model.Residual(parameterValues)
            residuals = numpy.concatenate((residuals,modelResidual))
        return residuals

    def AnalyJac(self,parameterValues):
        modelnum=1 
        #need to rewrite this to account for Shared Params
        # make a dictionary of original model parameters 
        for model in self.Models.values():
            modelJac = model.AnalyJac(parameterValues,parameterNames=self.OriginalParams[model][:])
            if modelnum==1:
                (ps1,ps2)=numpy.shape(modelJac) #shape with data x parameters
                total_param_num = len(self.theory.parameterNameList0)  
                temp_jac = numpy.zeros((ps1,total_param_num))
                current_param_list = self.OriginalParams[model]
                for c_index, param in enumerate(current_param_list):
                    p_index=self.theory.parameterNameList0.index(param)
                    temp_jac[:,p_index]=modelJac[:,c_index]
                modelnum+=1
                jacs = temp_jac
            else:
                (ps1,ps2)=numpy.shape(modelJac) #shape with data x parameters
                total_param_num = len(self.theory.parameterNameList0)
                temp_jac = numpy.zeros((ps1,total_param_num))
                current_param_list = self.OriginalParams[model]
                for c_index, param in enumerate(current_param_list):
                    p_index=self.theory.parameterNameList0.index(param)
                    temp_jac[:,p_index]=modelJac[:,c_index]
                jacs = numpy.concatenate((jacs, temp_jac),axis=0)
                modelnum+=1
        return jacs
    
    def Cost(self, parameterValues=None):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        residuals = self.Residual(parameterValues)
        return sum(residuals*residuals)
        #return sum(numpy.absolute(residuals))
    
    def SST(self, parameterValues=None):
        sst = 0.
        for model in self.Models.values():
            sst += model.SST(parameterValues)        
        return sst
        
    def R_square(self,parameterValues):
        """
        Calculates the R-square = 1 - cost / SST
        where SST is the sum of the squares about the mean
        """
        sst = self.SST(parameterValues)
        cost = self.Cost(parameterValues)
        return 1.- cost/sst
        
    def PlotFits(self, parameterValues=None, \
                 fontSizeLabels = 24, pylabLegendLoc=0, figNumStart=1,pow10first=False,interactive=False):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        figNum = figNumStart-1
        for model in self.Models.values():
            figNum+=1
            pylab.figure(figNum)
            model.PlotFunctions(parameterValues, fontSizeLabels=fontSizeLabels, pylabLegendLoc=pylabLegendLoc, interactive=interactive,pow10first=pow10first)
            # Weird bug: repeating figure needed to get to show
            pylab.figure(figNum)
    
    def PlotResiduals(self, parameterValues=None,fontSizeLabels = 18, pylabLegendLoc =0, figNumStart=1, pow10first=False):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        figNum = figNumStart-1
        for model in self.Models.values():
            figNum += 1
            pylab.figure(figNum)
            model.PlotResiduals(parameterValues, fontSizeLabels=fontSizeLabels, pylabLegendLoc=pylabLegendLoc, pow10first=pow10first)
            pylab.figure(figNum)

    def PlotCollapse(self, parameterValues=None, \
                 fontSizeLabels = 24, pylabLegendLoc=0, figNumStart=1, interactive=True,pow10first=False):
        if parameterValues is None:
            parameterValues = self.theory.initialParameterValues
        figNum = figNumStart-1
        for model in self.Models.values():
            figNum += 1
            pylab.figure(figNum)
            model.PlotFunctions(parameterValues, fontSizeLabels=fontSizeLabels, \
                                pylabLegendLoc=pylabLegendLoc, plotCollapse = True,interactive=interactive,pow10first=pow10first)
            pylab.figure(figNum)

    def PlotEnsemblePredictions(self,ens,figNumStart=1,fontSizeLabels=24,fontSizeLegend=18,sampling_freq=10.,rescale=1., pow10first=False):
        figNum=figNumStart-1
        for model in self.Models.values():
            figNum+=1
            pylab.figure(figNum)
            model.PlotEnsemblePredictions(self,ens,sampling_freq=sampling_freq, fontSizeLabels=fontSizeLabels,fontSizeLegend=fontSizeLegend,rescale=rescale,pow10first=False)
            pylab.figure(figNum)
            
    def BestFit(self,initialParameterValues=None, method=None , fixedParams=None):
        """
        uses various algorithms to fit function
        inputs:
        method: specifies different algorithms
                'lm' or None: minpack LM
                'lm_accel': LM + acceleration
                'bold': Cyrus' algorithm
                'boldAccel': Cyrus' algorithm +acceleration
        outputs:
            full parameter list, fitting algorithm outputs
        """

        if fixedParams:
            if not isinstance(fixedParams, list):
                fixedParams=[fixedParams]
            #Check now if the name is correct
            l_index=[]
            for index, par in enumerate(fixedParams):
                pName, pValue = par
                if pName not in self.theory.parameterNameList0:
                    print "%s is not a valid name. Ignored" %pName
                    l_index.append(index)
            if l_index:
                for i in l_index:
                    fixedParams.pop(i)
            self.SetFixedParams(fixedParams)

        if initialParameterValues is None:
            initialParameterValues = self.theory.initialParameterValues
        #d = numpy.ones(len(initialParameterValues))
        start_time = time.time()
        if method is None or method == 'lm':
            out = minpack.leastsq(self.Residual,initialParameterValues,full_output=1, ftol=1.e-16)
        elif method == 'boldAccel':
            initialParameterValues=numpy.array(initialParameterValues)
            out = BoldAccel.leastsq(self.Residual,None,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,ibold=0,verbose=True)
        elif method == 'bold':
            initialParameterValues = numpy.array(initialParameterValues)
            out = Bold.leastsq(self.Residual,None,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,ibold=0,verbose=True)
        #out = minpack.leastsq(self.Residual,self.AnalyJac,initialParameterValues,gtol=1e-8,xtol=1.49e-8,ftol=1e-16,full_output=1,Cgoal=4.e04)
        elif method == 'lm_accel':
            initialParameterValues=numpy.array(initialParameterValues)
            out = numrec.leastsq(self.Residual,self.AnalyJac,initialParameterValues,full_output=1,verbose=True, flags=[],maxfev=500)
        else:
            print "fitting method is not included"
        end_time = time.time()
        print "fitting took time (mins): ", (end_time-start_time)/60.
        print "number of function_calls:", f_counter
        
        if fixedParams:
            outputParameterValues = self.MergeFixedAndVariableParams(fixedParams,out[0])
        else:
            outputParameterValues = out[0]

        return outputParameterValues, out
        
    def PlotBestFit(self, initialParameterValues=None, \
                    figNumStart = 1, fixedParams = None,pow10first=False):
        
        if fixedParams:
            if not isinstance(fixedParams, list):
                fixedParams = [fixedParams]
            # Check now if the name is correct
            l_index = []
            for index, par in enumerate(fixedParams):
                pName, pValue = par
                if pName not in self.theory.parameterNameList0:
                    print "%s is not a valid name. Ignored" % pName
                    l_index.append(index)
            if l_index:
                for i in l_index:
                    fixedParams.pop(i)
                    
        # Call setfixedParams even if fixedParams = None to check
        # if original Names and values have to be used
        self.SetFixedParams(fixedParams)
        if initialParameterValues is None:
            initialParameterValues = self.theory.initialParameterValues

        print 'initial cost = ', self.Cost(initialParameterValues)
        out = self.BestFit(initialParameterValues)
        optimizedParameterValues = out[1][0]
        covar = out[1][1]
        #YJC: added check in case covariance is singular?
        if covar is not None:
            errors = [covar[i,i]**0.5 for i in range(len(covar))]
            #inv_t_student = scipy.special.stdtrit(len(errors),0.90)
        else:
            errors = None
        print 'optimized cost = ', self.Cost(optimizedParameterValues)
        #YJC: what are these? SST and R-value seem to be wrong...
        #print 'optimized SST = ', self.SST(optimizedParameterValues)
        #print 'R-value = ', self.R_square(optimizedParameterValues)
        print
        if fixedParams:
            if not isinstance(fixedParams, list):
                fixedParams = [fixedParams]
            print "=== Fixed parameters ================"
            for pName,pValue in fixedParams:
                print "%s = %2.2f" % (pName, pValue)
        # Print parameter values
        # YJC: changed printing here to print one sigma error instead of 95% confidence level
        if errors is not None:
            print "=== Fitting parameters (with one sigma error)=============="
            for name, val, error in \
                    zip(self.theory.parameterNameList,optimizedParameterValues, errors):
                print "%7s = %2.6f +/- %2.6f" %(name, val, error)
            print "=========================================================="
        else:
             print "=== Fitting parameters (no error available)=============="
             for name, val in \
                    zip(self.theory.parameterNameList,optimizedParameterValues):
                print "%7s = %2.6f " %(name, val)
             print "========================================================="

        #
        # Print plots
        #
        figNum = figNumStart-1
        for model in self.Models.values():
            for FT in [False,True]:
                figNum+=1
                pylab.figure(figNum)
                model.PlotFunctions(optimizedParameterValues, plotCollapse = FT, interactive=False,pow10first=pow10first)
                # Weird bug: repeating figure needed to get to show
                pylab.figure(figNum)
            figNum+=1
            pylab.figure(figNum)
            model.PlotResiduals(optimizedParameterValues,pow10first=pow10first)
            pylab.figure(figNum)

        #YJC: want to include the fixed Parameters in the output
        if fixedParams:
            outputParameterValues = self.MergeFixedAndVariableParams(fixedParams,optimizedParameterValues)
            #allParameterNames=self.theory.parameterNameList0
            #fittedParameterNames=self.theory.parameterNameList
            #outputParameterValues=numpy.zeros(len(allParameterNames))
            #for i in range(0,len(fittedParameterNames)):
            #    p_index=allParameterNames.index(fittedParameterNames[i])
            #    outputParameterValues[p_index]=optimizedParameterValues[i]
            #for pName, pValue in fixedParams:
            #    p_index =allParameterNames.index(pName)
            #    outputParameterValues[p_index]=pValue
        else:
            outputParameterValues=optimizedParameterValues

        return outputParameterValues


    def MergeFixedAndVariableParams(self,fixedParams,optimizedParameterValues):
        allParameterNames=self.theory.parameterNameList0
        fittedParameterNames=self.theory.parameterNameList
        outputParameterValues=numpy.zeros(len(allParameterNames))
        for i in range(0,len(fittedParameterNames)):
            p_index=allParameterNames.index(fittedParameterNames[i])
            outputParameterValues[p_index]=optimizedParameterValues[i]
        for pName, pValue in fixedParams:
            p_index =allParameterNames.index(pName)
            outputParameterValues[p_index]=pValue

        return outputParameterValues
