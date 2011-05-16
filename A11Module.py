from numpy import exp
import os
import SloppyScaling
reload(SloppyScaling)
import WindowScalingInfo as WS
reload(WS)
import Utils
reload(Utils)


name = 'A11' # This is the name used in the files

Xname = 's'
XnameTeX = r'$s$'
XscaledName = 'Ss'
##Xscaled = '(s*(1.0*k/L)**(sigma_k*zeta)/W)'
#XscaledTeX = r'$s \kappa^{\sigma_k \zeta} / W$'
Xscaled = '(s*(1.0*k/L)**(sigma_k*(1.+zeta)))'
XscaledTeX = r'$S_s=s/L_k^{(1+\zeta)}$'
WscaledName = 'Ws'
Wscaled = '(W*(1.0*k/L)**(sigma_k))'

Yname = 'A11' # This must be the name of the module !!!
YnameTeX = r'$A_{11}$'
#expansion = '(M1*Ss**0.5+M2*Ss+T0+T1*Ws+U0*Ss**m0+U1*Ws**(p1**2)*Ss**(-m1**2)+U2*Ws**(p2**2)*Ss**(-m2**2)+U3*Ws**(-p3**2)*Ss**(m3**2)+U4*Ws**(-p4**2)*Ss**(m4**2)+U5*Ws**p5*Ss**m5)'
#expansion = '(M0+M1*Ss**0.5+U0*Ss**(1.+m0**2)+U1*(Ss/Ws)**(m1**2)+U2*(Ss/Ws**(1.+zeta))**(-m2**2))'
#expansion = '(Ws**(p1)*(U1*(Ss/Ws)**(m1**2)+U2*(Ss/Ws**(1.+zeta))**(-m2**2))+U0)'
expansion = '(M1*Ss**0.5+M2*Ss**nd**2+U0+U1*(Ss/Ws)**(m1**2)+U2*(Ss/Ws**p1)**(-m2**2))'

Ytheory = '(Ss/Ws)**((2.-tau)*(1+zeta)/zeta)*(1./s)*exp(-'+expansion+')'
Yscaled = '(Ss/Ws)**((tau-2.)*(1.+zeta)/zeta)*s*A11'       
#Yscaled = 'Ss**((tau-2.)*(1.+zeta)/zeta)*s*A11'
YscaledTeX = \
   r'${\cal{A}}_{11}=(s/W L_k^{\zeta})^{(\tau-2) (1+\zeta)/\zeta} s A_{11}$'

title = r'$A_{11}(s,k,W)$'+': Area covered by avalanches of size s in window of width W'
#scalingTitle = r'${\cal A}_{11}(S_s,W_s)$' + ' scaling function'
scalingTitle = None

#
# Include corrections
#

#YJC: note that these below aren't corrections, but the scaling functions
#Ytheory_corrections = "exp(U0-U1*"+Ws+"**n2/Ss**n3)"
#Ytheory_corrections = "exp(U0-U1*"+Ws+"**n2/Ss**n3-Ux*"+Ws+"**n4/Ss**n5)"
#Ytheory_corrections = "exp(U0-U1*"+Ws+"/(1.0*Ss)+U2*"+Ws+"/(1.0*Ss**2.))"
#Ytheory_corrections = "exp(Ah1/s+Ah2/s**2)*exp(Uh1*(s*(1.0*k/L)**sigma_k/win**(1.+zeta)) + Uh2/(s*(k/L)**sigma_k/win**(1.+zeta)))"
#
Ytheory_corrections = "exp(Ah1/s+Ah2/s**2)"

#parameterNames = "tau,sigma_k,zeta,M1,M2,T0,T1,U0,U1,U2,U3,U4,U5,p1,p2,p3,p4,p5,m0,m1,m2,m3,m4,m5"
parameterNames = "tau,sigma_k,zeta,M1,M2,U0,U1,U2,m1,m2,nd,p1"
#parameterNames = "tau,sigma_k,zeta,U0,U1,U2,m1,m2,p1"
parameterNames_corrections = "Ah1,Ah2"

# parametervalues to match original working function
# new best fit parametervalues
#initialParameterValues= (1.252,0.457,0.643,1.66,5.69,-1.39,-0.52,-6.3,1.08,0.0,0.19,-10.16,11.14,1.64,0.0,-1.40,-0.32,0.0085,0.97,1.3,-1.18,1.39,0.15,0.061)
#initialParameterValues=(1.252,0.457,0.643,1.66,5.69,-1.39,-0.52,-6.3,1.08,0.0,0.19,-10.16,11.14,1.0,0.0,-1.0,-0.32,  0.0085,0.97,0.8,-1.18,1.0,0.15,0.061) 
initialParameterValues_corrections = (0.0,0.0)

# new theory parameter values
initialParameterValues = (1.252,0.457,0.643,-0.073,0.14,0.48,0.34,0.40,1.4,1.4,1.19,1.58)
#initialParameterValues = (1.252,0.457,0.643,0.29,0.50,1.88,1.4,1.4,1.0)


# Correct if spaces are included in the parameters names
parameterNames = parameterNames.replace(" ","")
#parameterNames_corrections = parameterNames_corrections.replace(" ","")

# priors to include (parameters to be squared)
# 1st try just putting priors on powers in coefficients
priorList=['m1','m2','log(nd**2)','p1']
#priorList = None

if WS.corrections_to_scaling:
    Ytheory = Ytheory + "*" + Ytheory_corrections
    parameterNames = parameterNames + "," + parameterNames_corrections
    initialParameterValues = initialParameterValues + initialParameterValues_corrections

Jac_dict={}

Jac_dict['expfactor']= " exp(Ah2/s**2 + Ah1/s - M1*Ss**0.5 - M2*Ss**nd**2 - U0 - U1*(Ss/Ws)**m1**2 - U2/(Ss/Ws**p1)**m2**2)"

Jac_dict['prefactor']= "(Ss/Ws)**(((2 - tau)*(1 + zeta))/zeta)"

Jac_dict['prior'] = {'m1':'1','m2':'1','nd':'2/nd','p1':'1'}
Jac_dict['priorvalues']={'m1':'m1','m2':'m2','nd':'log(nd**2)','p1':'p1'}

dYdp = '[-((expfactor*prefactor*(1 + zeta)*scipy.log(((k/L)**(sigma_k*zeta)*s)/W))/(s*zeta)), (expfactor*prefactor*(-(m1**2*U1*(((k/L)**(sigma_k*zeta)*s)/W)**m1**2*zeta) - 0.5*M1*Ss**0.5*(1 + zeta) - M2*nd**2*Ss**nd**2*(1 + zeta) - (-2 + tau)*(1 + zeta) + (m2**2*U2*(1 - p1 + zeta))/(Ss/Ws**p1)**m2**2)*scipy.log(k/L))/s, (expfactor*prefactor*(sigma_k*(2 - 0.5*M1*Ss**0.5 - M2*nd**2*Ss**nd**2 - tau - m1**2*U1*(((k/L)**(sigma_k*zeta)*s)/W)**m1**2 + (m2**2*U2)/(Ss/Ws**p1)**m2**2 + 2/zeta - tau/zeta)*scipy.log(k/L) + ((-2 + tau)*scipy.log(((k/L)**(sigma_k*zeta)*s)/W))/zeta**2))/s, -((expfactor*prefactor*Ss**0.5)/s), -((expfactor*prefactor*Ss**nd**2)/s), -((expfactor*prefactor)/s), -((expfactor*(((k/L)**(sigma_k*zeta)*s)/W)**(m1**2 - ((-2 + tau)*(1 + zeta))/zeta))/s), -((expfactor*prefactor)/(s*(Ss/Ws**p1)**m2**2)), (-2*expfactor*m1*U1*(((k/L)**(sigma_k*zeta)*s)/W)**(m1**2 - ((-2 + tau)*(1 + zeta))/zeta)*scipy.log(((k/L)**(sigma_k*zeta)*s)/W))/s, (2*expfactor*m2*prefactor*U2*scipy.log(Ss/Ws**p1))/(s*(Ss/Ws**p1)**m2**2), (-2.*expfactor*M2*nd*prefactor*Ss**nd**2.*scipy.log(Ss))/s, -((expfactor*m2**2*prefactor*U2*scipy.log(Ws))/(s*(Ss/Ws**p1)**m2**2)), (expfactor*prefactor)/s**2, (expfactor*prefactor)/s**3]'


theory = SloppyScaling.ScalingTheory(Ytheory, parameterNames, \
                initialParameterValues, WS.independentNames, \
                scalingX = Xscaled, scalingY = Yscaled, scalingW=Wscaled,\
                scalingXTeX = XscaledTeX, \
                scalingYTeX = YscaledTeX, \
                title = title, \
                scalingTitle = scalingTitle, \
                Xname=Xname, XnameTeX=XnameTeX, XscaledName=XscaledName, \
                Yname=Yname, YnameTeX=YnameTeX, WscaledName=WscaledName, \
                normalization = WS.normalization, priorList=priorList, priorCoeff=WS.priorCoeff,dYdp=dYdp,Jac_dict=Jac_dict)

data = SloppyScaling.Data()

loaded = 0
for independent in WS.independentValues:
    L, k, W = independent
    ext =  "_" + WS.simulType + ".bnd"
    k_string = "_k"
    fileName = "".join([WS.dataDirectory,name,"_W=",str(W),\
                        k_string,str(k), "_System_Size=",str(2*L), "x", str(L), ext])
    success = data.InstallCurve(independent, fileName, \
        defaultFractionalError=0.1,\
        pointSymbol=WS.Symbol[independent], \
        pointColor=WS.Color[independent], \
        errorCol=2,initialSkip = WS.rows_to_skip)
    loaded += success

# YJC: declared nFiles
nFiles = len(WS.independentValues) 

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
