from numpy import exp
import os
import SloppyScaling
reload(SloppyScaling)
import WindowScalingInfo as WS
reload(WS)
import Utils
reload(Utils)

name = 'A00' # This is the name used in the files

Xname = 's'
XnameTeX = r'$s$'
XscaledName = 'Ss'
Xscaled = "(s*(1.0*k/L)**(sigma_k*(1.+zeta)))"
XscaledTeX = r'$S_s = s/L_k^{(1+\zeta)}$'

WscaledName = 'Ws'
Wscaled = '(W*(1.0*k/L)**(sigma_k))'

Yname = 'A00' # This must be the name of the module !!!!!!!!!
YnameTeX =r'$A_{00}$'
# asymptotes to A(s) distribution for large Ws
Ytheory = "(s*(k/L)**(sigma_k*(1.+zeta)))**(2.-tau) /s * exp(-((Ss**ns/Ws**nw)*Ixsw+Z0+Z1*Ss**0.5+Zn*Ss**n00))"

Yscaled = \
        "(1.0*k/L)**((sigma_k)*(tau-2.)*(1.+zeta)) * s**(tau-1.) * A00"
YscaledTeX = \
   r'${\cal{A}}_{00} = L_k^{(2-\tau) (1+\zeta)} s^{\tau-1} A_{00}$'

title = r'$A_{00}(s,k,W)$'+': Area covered by avalanches of size s in window of width W'
scalingTitle = r'${\cal A}_{00}(S_s,W_s)$'+' scaling function'
scalingTitle = None

#
# Include corrections
#

# for large Ws, this should just be the A(s,k) below 
#Ytheory = "(1.0*k/L)**(sigma_k*(2.-tau)*(1.+zeta)) * S**(2.-tau) / S *exp(-((S**(1./(1.+zeta))*(1.0*k/L)**sigma_k)*Ixs)**ns+U0)"
#Ytheory_corrections= "exp(-U1/S)"
Ytheory_corrections = "exp(A1/s+A2/s**2)"

parameterNames = "tau,sigma_k,zeta,ns,nw,n00,Ixsw,Z0,Z1,Zn"
#parameterNames_corrections = "U1"
parameterNames_corrections = "A1,A2"
initialParameterValues = (1.27,0.46,0.63,1.41,2.33,1.0,5.32,2.3,-0.4,0.57)
#initialParameterValues_corrections = (0.1,)
initialParameterValues_corrections = (-0.7,-0.9)

# priors
priorList=['ns','nw','n00']
#priorList=None

if WS.corrections_to_scaling:
    Ytheory = Ytheory + "*" + Ytheory_corrections
    parameterNames = parameterNames + "," + parameterNames_corrections
    initialParameterValues = initialParameterValues + initialParameterValues_corrections

Jac_dict={}

Jac_dict['Sh']="numpy.sqrt(Ss)"

Jac_dict['expfactor'] = "numpy.e**(A2/s**2 + A1/s - (Ixsw*Ss**ns)/Ws**nw - Z0 - Ss**0.5*Z1 - Zn*Ss**n00)"

Jac_dict['priorvalues']={'ns':'ns','nw':'nw','n00':'n00'}
Jac_dict['prior']={'ns':'1','nw':'1','n00':'1'}


dYdp = "[-((expfactor*Sh**(4 - 2*tau)*numpy.log(Sh**2))/s), -(expfactor*Sh**(4 - 2*tau)*(2*Ixsw*Sh**(2*ns)*(ns - nw + ns*zeta) + Ws**nw*(1 + zeta)*(-4 + 2*tau + Sh*Z1 + 2*n00*Sh**(2*n00)*Zn))*numpy.log(k/L))/(2.*s*Ws**nw), -(expfactor*Sh**(4 - 2*tau)*sigma_k*(2*Ixsw*ns*Sh**(2*ns) + Ws**nw*(-4 + 2*tau + Sh*Z1 + 2*n00*Sh**(2*n00)*Zn))*numpy.log(k/L))/(2.*s*Ws**nw), -((expfactor*Ixsw*Sh**(4 + 2*ns - 2*tau)*numpy.log(Sh**2))/(s*Ws**nw)), (expfactor*Ixsw*Sh**(4 + 2*ns - 2*tau)*numpy.log(Ws))/(s*Ws**nw), -((expfactor*Sh**(4 + 2*n00 - 2*tau)*Zn*numpy.log(Sh**2))/s), -((expfactor*Sh**(4 + 2*ns - 2*tau))/(s*Ws**nw)), -((expfactor*Sh**(4 - 2*tau))/s), -((expfactor*Sh**(5 - 2*tau))/s), -((expfactor*Sh**(4 + 2*n00 - 2*tau))/s), (expfactor*Sh**(4 - 2*tau))/s**2, (expfactor*Sh**(4 - 2*tau))/s**3]"


# If single independent parameter, must have comma after it -- makes it a tuple
theory = SloppyScaling.ScalingTheory(Ytheory, parameterNames, \
                initialParameterValues, WS.independentNames, \
                scalingX = Xscaled, scalingY = Yscaled, scalingW=Wscaled,\
                scalingXTeX = XscaledTeX, \
                scalingYTeX = YscaledTeX, \
                title = title, \
                scalingTitle = scalingTitle, \
                Xname=Xname, XnameTeX=XnameTeX, XscaledName=XscaledName,\
                Yname=Yname, YnameTeX=YnameTeX, WscaledName=WscaledName,\
                normalization = WS.normalization,priorList=priorList, dYdp=dYdp,Jac_dict=Jac_dict)

data = SloppyScaling.Data()

loaded = 0
for independent in WS.independentValuesA00:
    L, k, W = independent
    ext =  "_" + WS.simulType + ".bnd"
    k_string = "_k"
    fileName = "".join([WS.dataDirectory,name,"_W=",str(W),\
                        k_string,str(k), "_System_Size=",str(2*L), "x", str(L), ext])
    success = data.InstallCurve(independent, fileName, \
        pointSymbol=WS.SymbolA00[independent], \
        pointColor=WS.ColorA00[independent], \
        errorCol=2,initialSkip = WS.rows_to_skip)
    loaded += success

nFiles = len(WS.independentValuesA00)
if loaded ==  nFiles:
    print "Loaded %2d/%2d files (%s)" % (loaded, nFiles, name)
else:
    print "====================="
    print "Attention! %2d/%2d files are missing (%s)" %  (nFiles-loaded, nFiles, name)
    print "====================="

f = __file__
f = f.split("/")[-1]
thisModule = f.split('Module.py')[0]
exec(thisModule + "= SloppyScaling.Model(theory, data, name, WS.sortedValues)")
