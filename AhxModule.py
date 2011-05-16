from numpy import exp
import os
import SloppyScaling
reload(SloppyScaling)
import WindowScalingInfo as WS
reload(WS)
import Utils
reload(Utils)


name = 'A11' # This is the name used in the files

Xname = 'hx'
XnameTeX = r'$h_x$'
XscaledName = 'hs'
##Xscaled = '(s*(1.0*k/L)**(sigma_k*zeta)/W)'
#XscaledTeX = r'$s \kappa^{\sigma_k \zeta} / W$'
Xscaled = '(hx*(1.0*k/L)**(sigma_k*(zeta)))'
XscaledTeX = r'$h_s=h_x k^{\sigma_k (\zeta)}$'
WscaledName = 'Ws'
Wscaled = '(1*(1.0*k/L)**(sigma_k))'

Yname = 'Ahx' # This must be the name of the module !!!
YnameTeX = r'$A_{hx}$'
expansion = '(M1*(hs*Ws)**0.5+M2*(hs*Ws)**nd**2+U0+U1*(hs)**(m1**2)+U2*((hs*Ws)/Ws**p1)**(-m2**2))'
#expansion = '(U0+U1*(hs)**(m1**2))'
Ytheory = '(hs)**((2.-tau)*(1+zeta)/zeta)*(1./hx)*exp(-'+expansion+')'
Yscaled = '(hs)**((tau-2.)*(1.+zeta)/zeta)*hx*A11'       
#Yscaled = 'Ss**((tau-2.)*(1.+zeta)/zeta)*s*A11'
YscaledTeX = \
   r'${\cal{A}}_{hx}=(h/L_k^{\zeta})^{(\tau-2) (1+\zeta)/\zeta} hx A_{hx}$'

title = r'$A_{hx}(h_x,L_k,W)$'+': Area covered by local height $h_x$'

scalingTitle = None

#
# Include corrections
#

#
Ytheory_corrections = "exp(Ah1/hx+Ah2/hx**2)"

parameterNames = "tau,sigma_k,zeta,M1,M2,U0,U1,U2,m1,m2,nd,p1"
#parameterNames = "tau,sigma_k,zeta,U0,U1,m1"
parameterNames_corrections = "Ah1,Ah2"

# parametervalues to match original working function
initialParameterValues = (1.2636,0.4630,0.63,-0.5,0.21,0.42,0.54,0.86,1.21,1.27,1.05,1.653)
#initialParameterValues = (1.2636,0.4360,0.63,0.42,0.54,1.21)
initialParameterValues_corrections = (0.2,0.6)


# Correct if spaces are included in the parameters names
parameterNames = parameterNames.replace(" ","")
#parameterNames_corrections = parameterNames_corrections.replace(" ","")

# priors to include (parameters to be squared)
# 1st try just putting priors on powers in coefficients
priorList=['m2']
#priorList = None

corrections_to_scaling=True

if corrections_to_scaling:
    Ytheory = Ytheory + "*" + Ytheory_corrections
    parameterNames = parameterNames + "," + parameterNames_corrections
    initialParameterValues = initialParameterValues + initialParameterValues_corrections

theory = SloppyScaling.ScalingTheory(Ytheory, parameterNames, \
                initialParameterValues, 'L,k', \
                scalingX = Xscaled, scalingY = Yscaled,scalingW = Wscaled,\
                scalingXTeX = XscaledTeX, \
                scalingYTeX = YscaledTeX, \
                title = title, \
                scalingTitle = scalingTitle, \
                Xname=Xname, XnameTeX=XnameTeX, XscaledName=XscaledName, \
                Yname=Yname, YnameTeX=YnameTeX, \
                normalization = None, priorList=priorList, priorCoeff=WS.priorCoeff)

data = SloppyScaling.Data()

Values= [[4096,[0.4096,4.096,40.96]],\
             [8192,[0.08192]],\
             [16384,[0.016384,0.0016384]]] 

independentNames, independentValues = Utils.get_independent(Values,sorting=True)

Symbol, Color = Utils.MakeSymbolsAndColors(independentValues)

loaded = 0
for independent in independentValues:
    L, k = independent
    W = 1
    ext =  "_" + WS.simulType + ".bnd"
    k_string = "_k"
    fileName = "".join([WS.dataDirectory,name,"_W=1",\
                        k_string,str(k), "_System_Size=",str(2*L), "x", str(L), ext])
    success = data.InstallCurve(independent, fileName, \
        defaultFractionalError=0.1,\
        pointSymbol=Symbol[independent], \
        pointColor=Color[independent], \
        errorCol=2,initialSkip = WS.rows_to_skip)
    loaded += success

# YJC: declared nFiles
nFiles = len(independentValues) 

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
