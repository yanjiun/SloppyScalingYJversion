from numpy import exp
import os
import SloppyScaling
reload(SloppyScaling)
import WindowScalingInfo as WS
reload(WS)
import Utils
reload(Utils)


name = 'A10' # This is the name used in the files

Xname = 's'
XnameTeX = r'$s$'
XscaledName = 'Ss'
Xscaled = "(s*(1.0*k/L)**(sigma_k*(1.+zeta)))"
XscaledTeX = r'$S_s = s/L_k^{(1+\zeta)}$'
WscaledName = 'Ws'
Wscaled = '(W*(1.0*k/L)**(sigma_k))'

Yname = 'A10' # This must be the name of the module !!!
YnameTeX = r'$A_{10}$'
Ytheory = "((1.0*k/L)**(sigma_k))**((2.-tau)*(1.+zeta))*s**(1.-tau+1./(1.+zeta))/W*exp(-(Iy0+Iy1*Ss**0.5+Iyn*Ss**n10+Iy11*(Ss**n2_10/Ws**n3_10)))"
Yscaled = \
        "((1.0*k/L)**(sigma_k))**(-(2.-tau)*(1.+zeta)) * s**(-1.+tau-1./(1.+zeta)) * W * A10"
YscaledTeX = \
   r'${\cal{A}}_{10} = L_k^{(2-\tau) (1+\zeta)} s^{-1+\tau-1/(1+\zeta)} W A_{10}$'

title = r'$A_{10}(s,k,W)$'+': Area covered by avalanches of size s in window of width W'
#scalingTitle = r'${\cal A}_{10}(S_s,W_s)$'+' scaling function'
scalingTitle = None

#
# Include corrections
#
Ytheory_corrections = "exp(Aw1/s+Aw2/s**2)"

parameterNames = "tau,sigma_k,zeta,Iy0,Iy1,Iyn,Iy11,n10,n2_10,n3_10"
parameterNames_corrections = "Aw1,Aw2"
initialParameterValues = (1.27,0.46,0.63,1.34,0.31,0.02,1.0,1.0,1.0,1.0)
initialParameterValues_corrections = (0.1,0.1)

#priorList
priorList = ['n10','n2_10','n3_10']
#priorList = None

# Correct if spaces are included in the parameters names
parameterNames = parameterNames.replace(" ","")
parameterNames_corrections = parameterNames_corrections.replace(" ","")

if WS.corrections_to_scaling:
    Ytheory = Ytheory + "*" + Ytheory_corrections
    parameterNames = parameterNames + "," + parameterNames_corrections
    initialParameterValues = initialParameterValues + initialParameterValues_corrections

Jac_dict={}
 
Jac_dict['expfactor']='numpy.e**(-Iy0 + Aw2/s**2 + Aw1/s - Iy1*Ss**0.5 - Iyn*Ss**n10 - (Iy11*Ss**n2_10)/Ws**n3_10)'

Jac_dict['Sh']='numpy.sqrt(Ss)'

Jac_dict['prefactor']='s**(1 - tau + 1/(1 + zeta))/((k/L)**sigma_k)**((-2 + tau)*(1 + zeta))'

Jac_dict['priorvalues']={'n10':'n10','n2_10':'n2_10','n3_10':'n3_10'}
Jac_dict['prior']={'n10':'1','n2_10':'1','n3_10':'1'}

dYdp="[-((expfactor*prefactor*((1 + zeta)*numpy.log((k/L)**sigma_k) + numpy.log(s)))/W), -(expfactor*prefactor*((Iy1*Sh + 2*(-2 + Iyn*n10*Sh**(2*n10) + tau))*Ws**n3_10*(1 + zeta) + 2*Iy11*Sh**(2*n2_10)*(n2_10 - n3_10 + n2_10*zeta))*numpy.log(k/L))/(2.*W*Ws**n3_10), (expfactor*prefactor*(-((-2 + tau)*numpy.log((k/L)**sigma_k)) + (-(Iy1*Sh*sigma_k)/2. - Iyn*n10*Sh**(2*n10)*sigma_k - (Iy11*n2_10*Sh**(2*n2_10)*sigma_k)/Ws**n3_10)*numpy.log(k/L) - numpy.log(s)/(1 + zeta)**2))/W, -((expfactor*prefactor)/W), -((expfactor*prefactor*Sh)/W), -((expfactor*prefactor*Sh**(2*n10))/W), -((expfactor*prefactor*Sh**(2*n2_10))/(W*Ws**n3_10)), -((expfactor*Iyn*prefactor*Sh**(2*n10)*numpy.log(Sh**2))/W), -((expfactor*Iy11*prefactor*Sh**(2*n2_10)*numpy.log(Sh**2))/(W*Ws**n3_10)), (expfactor*Iy11*prefactor*Sh**(2*n2_10)*numpy.log(Ws))/(W*Ws**n3_10), (expfactor*s**(-tau + 1/(1 + zeta)))/(((k/L)**sigma_k)**((-2 + tau)*(1 + zeta))*W), (expfactor*s**(-1 - tau + 1/(1 + zeta)))/(((k/L)**sigma_k)**((-2 + tau)*(1 + zeta))*W)]"


# If single independent parameter, must have comma after it -- makes it a tuple
theory = SloppyScaling.ScalingTheory(Ytheory, parameterNames, \
                initialParameterValues, WS.independentNames, \
                scalingX = Xscaled, scalingY = Yscaled, scalingW=Wscaled,\
                scalingXTeX = XscaledTeX, \
                scalingYTeX = YscaledTeX, \
                title = title, \
                scalingTitle = scalingTitle, \
                Xname=Xname, XscaledName=XscaledName, WscaledName=WscaledName,\
                Yname=Yname, XnameTeX=XnameTeX, YnameTeX=YnameTeX,\
                normalization = WS.normalization,priorList=priorList, priorCoeff=WS.priorCoeff,dYdp=dYdp, Jac_dict=Jac_dict)

data = SloppyScaling.Data()

loaded = 0
for independent in WS.independentValuesA10:
    L, k, W = independent
    ext =  "_" + WS.simulType + ".bnd"
    k_string = "_k"
    fileName = "".join([WS.dataDirectory,name,"_W=",str(W),\
                        k_string,str(k), "_System_Size=",str(2*L), "x", str(L), ext])
    success = data.InstallCurve(independent, fileName, \
        pointSymbol=WS.SymbolA10[independent], \
        pointColor=WS.ColorA10[independent], \
        errorCol=2,initialSkip = WS.rows_to_skip)
    loaded += success

nFiles = len(WS.independentValuesA10)
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
