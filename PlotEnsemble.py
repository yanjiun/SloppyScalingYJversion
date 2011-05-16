# PlotEnsemble.py
#
# Bryan Daniels
# 2007-05-21 - 2007-05-24
# 2007-05-30
# 2007-06-14 added plotContours2D for the actual cost
#
# 2007-12-07 debugging
# 2007-12-09 - 2007-12-12
#
# 2008.04.14 added support for arbitrary projections
#
# Yan-Jiun Chen 
# 2009.12.03 code added function to plot predictions, and histograms
# 2009.03.08 make plots nicer...


import pylab, scipy
import cPickle
import copy
from SloppyCell.ReactionNetworks import *


def scatterColors(xdata,ydata,colordata,size=50,alpha=0.75):
	"""
	Makes a scatter plot with colors specified by colordata.
	"""
	pylab.scatter(xdata,ydata,size*pylab.ones(len(xdata)),colordata,\
			alpha=alpha,faceted=False)

def plotEnsemble2D(ens,v1,v2,colordata=None,hess=None,\
		   size=50,labelBest=True,ensembleAlpha=0.75,contourAlpha=1.0):
	"""
	Plots a 2-dimensional projection of a given parameter
	ensemble, along given directions:
	     -- If v1 and v2 are scalars, project onto plane given by
		those two bare parameter directions.
	     -- If v1 and v2 are vectors, project onto those two vectors.
	
	When given colordata (either a single color, or an array
	of different colors the length of ensemble size), each point
	will be assigned a color based on the colordata.
	
	With labelBest set, the first point in the ensemble is
	plotted larger (to show the 'best fit' point for a usual 
	parameter ensemble).
	
	If a Hessian is given, cost contours will be plotted
	using plotContours2D.
	"""
	if pylab.shape(v1) is ():
		xdata = pylab.transpose(ens)[v1]
		ydata = pylab.transpose(ens)[v2]
		
		# label axes
		param1name, param2name = '',''
		try:
		    paramLabels = ens[0].keys()
		except:
		    paramLabels = None
		if paramLabels is not None:
		    param1name = ' ('+paramLabels[param1]+')'
		    param2name = ' ('+paramLabels[param2]+')'
		pylab.xlabel('Parameter '+str(v1)+param1name)
		pylab.ylabel('Parameter '+str(v2)+param2name)
	else:
		xdata = pylab.dot(ens,v1)
		ydata = pylab.dot(ens,v2)

	if colordata==None:
		colordata = pylab.ones(len(xdata))
		
	if labelBest: # plot first as larger circle
		if pylab.shape(colordata) is (): # single color
		    colordata0 = colordata
		    colordataRest = colordata
		else: # specified colors
		    colordata0 = [colordata[0]]
		    colordataRest = colordata[1:]
		scatterColors(xdata[1:],ydata[1:],colordataRest,		\
				size,alpha=ensembleAlpha)
		scatterColors([xdata[0]],[ydata[0]],colordata0,			\
				size*4,alpha=ensembleAlpha)
	else:
		scatterColors(xdata,ydata,colordata,size,alpha=ensembleAlpha)
		
	if hess is not None:
		plotApproxContours2D(hess,param1,param2,pylab.array(ens[0]),	\
			alpha=contourAlpha)


def plotApproxContours2D(hess,param1,param2,bestfit,plotPoints=100,		\
			xrange=None,yrange=None,numContours=20,alpha=1.0):
	"""
	Given the Hessian, plots the corresponding local 
	approximation of contours of constant cost, in the
	plane of the best-fit parameter point.
	
	param1,param2 : the two bare parameter directions of
			the two-dimensional plot
	bestfit	      : the vector location of the best fit
	"""
	#eigvals, eigvecs = Utility.eig(hess)
	
	costXY = lambda x,y: 0.5*hess[param1,param1]*(x-bestfit[param1])**2 	\
			       + hess[param1,param2]*(x-bestfit[param1])*	\
						     (y-bestfit[param2]) 	\
			   + 0.5*hess[param2,param2]*(y-bestfit[param2])**2
	
	# for debugging
	#return costXY
	
	contourFromFunction(costXY,plotPoints=plotPoints,xrange=xrange,		\
			    yrange=yrange,numContours=numContours,alpha=alpha)


def plotContours2D(cost,vec1,vec2,bestfit,plotPoints=100,			\
		xrange=None,yrange=None,bareParameters=False,			\
		numContours=20,alpha=1.0):
	"""
	Given a cost function that takes vectors
	of parameters, plots contours of constant cost, in the
	plane of the best-fit parameter point
	along vectors vec1 and vec2, with the bestfit
	parameters at the origin.
	
	bareParameters (True) : if True, instead plots along bare parameter
				directions given by integers vec1 and vec2,
				with origin at zero instead of at bestfit 
				parameters
	"""
	bestfit = copy.copy(pylab.array(bestfit))
	
	if bareParameters:
	   def costXY(x,y):
	   	currentLoc = bestfit
	   	currentLoc.itemset(vec1,x)
	   	currentLoc.itemset(vec2,y)
	   	return cost(currentLoc)
	else:	
	   costXY = lambda x,y: cost(bestfit + x*vec1 + y*vec2)
	   
	# for debugging
	#return costXY 
	   
	contourFromFunction(costXY,plotPoints=plotPoints,xrange=xrange,\
			    yrange=yrange,numContours=numContours,alpha=alpha)

# 6.14.07
# 12.10.07	
def contourFromFunction(XYfunction,plotPoints=100,\
			xrange=None,yrange=None,numContours=20,alpha=1.0, contourLines=None):
	"""
	Given a 2D function, plots constant contours over the given
	range.  If the range is not given, the current plotting
	window range is used.
	"""
	
	# set up x and y ranges
	currentAxis = pylab.axis()
	if xrange is not None:
		xvalues = pylab.linspace(xrange[0],xrange[1],plotPoints)
	else:
		xvalues = pylab.linspace(currentAxis[0],currentAxis[1],plotPoints)
	if yrange is not None:
		yvalues = pylab.linspace(yrange[0],yrange[1],plotPoints)
	else:
		yvalues = pylab.linspace(currentAxis[2],currentAxis[3],plotPoints)
	
	#coordArray = _coordinateArray2D(xvalues,yvalues)
	# add extra dimension to this to make iterable?
	# bug here!  need to fix for contour plots
	z = map( lambda y: map(lambda x: XYfunction(x,y), xvalues), yvalues)
	if contourLines:
		pylab.contour(xvalues,yvalues,z,contourLines,alpha=alpha)
	else:
		pylab.contour(xvalues,yvalues,z,numContours,alpha=alpha)
	
	
def _coordinateArray2D(xvalues,yvalues):
	"""
	Used to create a 2D array of coordinates for use in making
	contour plots.
	"""
	#xvalues2D = pylab.repeat([xvalues],len(yvalues))
	#yvalues2D = pylab.repeat([yvalues],len(xvalues)).transpose()
	
	#coordArray = []
	
	for x in xvalues:
		for y in yvalues:
			coordArray.append([x,y])
	
	return coordArray
	
def _map2D(func,array2D):
	return map( lambda row: map( func,row ) , array2D )
	
def getSloppyParameters(model,ens):
	"""
	Returns a list of the 
	sloppiest parameter for each member of the
	given parameter ensemble.  (here, the sloppiest parameter is
	defined as the one having the largest element in the
	sloppiest eigenvector)
	"""
	sloppy_param_list = []	
	for params in ens:
		J, JtJ = model.GetJandJtJInLogParameters(pylab.log(params))
		u, v = Utility.eig(JtJ)
		last_eigvect = v[:,len(v)-1]
		sloppy_param_list.append(abs(last_eigvect).argmax())
	return sloppy_param_list
	
	
def getStiffParameters(model,ens):
	"""
	Returns a list of the 
	stiffest parameter for each member of the
	given parameter ensemble.  (here, the stiffest parameter is
	defined as the one having the largest element in the
	stiffest eigenvector)
	"""
	stiff_param_list = []	
	for params in ens:
		J, JtJ = model.GetJandJtJInLogParameters(pylab.log(params))
		u, v = Utility.eig(JtJ)
		first_eigvect = v[:,0]
		stiff_param_list.append(abs(first_eigvect).argmax())
	return stiff_param_list

# the functions below work with SloppyScaling instead of Sloppy Cell...
# function to make histograms of params
def params_hist(m, ens,sampling_freq=1000):
    #ens_file=open(filename,'rb')
    #ens, ens_Fs, ratio=cPickle.load(ens_file)
    #ens_file.close()
    param_Names = m.theory.parameterNameList
    for i in range(0, len(param_Names)):
        exec(param_Names[i]+'s=[]')
    N_sample = len(ens)/sampling_freq
    for i in range(0, N_sample):
        for j in range(0, len(param_Names)):
            exec(param_Names[j]+'s.append(ens[i][j])')
    for i in range(0, len(param_Names)):
        exec('pylab.hist('+param_Names[i]+'s)')
        pylab.savefig('temp/'+param_Names[i]+'.png')
        pylab.clf()
    param_series=[]
    for i in range(0,len(param_Names)):
        exec('param_series.append('+param_Names[i]+'s)')
    return param_series

def params_error(ens):
	mean = scipy.zeros(scipy.shape(ens[0]))
	sum_of_squares = scipy.zeros(scipy.shape(ens[0]))
	N = len(ens)
	for i in range(0, N):
		mean += ens[i]
		sum_of_squares +=ens[i]**2
	mean = mean/N 
	std = scipy.sqrt((sum_of_squares-N*mean**2)/(N-1.))

	return mean, std

#make function to generate series in eigendirections...
def get_vector_series(ens, eigvecs):
	vector_series=[]
	for i in range(0, len(eigvecs)):
		temp_series=[]
		for j in range(0, len(ens)):
			temp_series.append(scipy.sum(ens[j]*eigvecs[i]))
		vector_series.append(temp_series)
	return vector_series
	
def plot_ensemble_traj(m, ens, sampling_freq=10, fontSizeLabels=24, fontSizeLegend=20,rescale=1.,paper_fig=False):
    """
    this function plots the range of ensemble predictions for a SloppyScaling 
    Model class, to speed up the calculation, one can specify the sampling         frequency of the ensemble
     if the ensemble was sampled at a lower temperature than
    desired, one can set a rescale factor "rescale = T_H/T_L"
    """
    
    if paper_fig:
	    fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
	    inches_per_pt = 1.0/72.27               # Convert pt to inch
	    golden_mean = (scipy.sqrt(5)-1.0)/2.0         # Aesthetic ratio
	    fig_width = fig_width_pt*inches_per_pt  # width in inches
	    fig_height = fig_width*golden_mean      # height in inches
	    fig_size =  [fig_width,fig_height]
	    pylab.rcdefaults()
	    params = {'axes.labelsize': 10,\
		      'text.fontsize': 10,\
		      'legend.fontsize': 8,\
		      'xtick.labelsize': 10,\
		      'ytick.labelsize': 10,\
		      'lines.markersize':3,
		      'text.usetex': True,\
		      'figure.figsize': fig_size}
	    pylab.rcParams.update(params)
	    font = {'family':'serif',
		    'serif':'Times New Roman'}
	    pylab.rc('font',**font)
    else:
	    pylab.rcdefaults()
	    pylab.rcParams.update({'backend':'ps',
			       'xtick.labelsize':24,
                               'xtick.major.size':20,
                               'xtick.minor.size':10,
                               'ytick.labelsize':24,
                               'ytick.major.size':20,
                               'ytick.minor.size':10,
                               'lines.markersize':10,
                               'axes.labelsize':24,\
                               'legend.fontsize':20,
                               'legend.columnspacing':1.5,
                               'figure.figsize':[10.,10.],\
                               'text.usetex':False,
                               })
            font = {'family':'serif',
                    'serif':'Times New Roman'}
            pylab.rc('font',**font)
           

    num_to_count = int(len(ens)/sampling_freq)

    # plot trajectories
    # want to calculate the mean of the ensemble and also the range of the ensemble... (use pylab.fill to plot)
    # but you can't do this by calculating mean of params, you have to do this by calculating mean of points along trajectory!!!!!!!
    # for using jointModules
    pylab.ioff()
    ax0 = [1.e99,0,1.e99,0]
    for model in m.Models.values():
        pylab.figure()
	#pylab.axes([0.2,0.2,0.95-0.2,0.95-0.20])
	#pylab.axes([0.125,0.10,0.95-0.125,0.95-0.10])
	pylab.axes([0.15,0.35,0.95-0.15,0.95-0.35])
	data_experiments = model.data.experiments
        #data_experiments = sorted(model.data.experiments)
        for independentValues in data_experiments:
            Xdata = model.data.X[independentValues]
            Ydata = model.data.Y[independentValues]
	    Xtheory = scipy.logspace(scipy.log10(min(Xdata)),scipy.log10(max(Xdata)),num=100)
	    pointType = model.data.pointType[independentValues]
	    errorBar = model.data.errorBar[independentValues]
	    mean_theory = scipy.zeros(len(Xtheory))
	    std_theory = scipy.zeros(len(Xtheory))
	    #max_theory = scipy.zeros(len(Xdata))
	    #min_theory = scipy.zeros(len(Xdata))
            for i in range(0, num_to_count):
		    ens_theory = model.theory.Y(Xtheory,ens[i*sampling_freq],independentValues)
		    mean_theory += ens_theory
		    std_theory += (ens_theory)**2
	    mean_theory = mean_theory/(1.0*num_to_count)
	    std_theory = scipy.sqrt((std_theory-num_to_count*mean_theory**2)/(num_to_count-1.))
	    pylab.loglog(Xdata,Ydata,pointType[1])
	    lb = model.getLabel(model.theory.independentNames,independentValues,pow10first=True)
	    pylab.errorbar(Xdata,Ydata, yerr=errorBar, fmt=pointType,label=lb)
	    pylab.loglog(Xtheory,mean_theory,pointType[0])
	    axis_dep=model.getAxis(Xdata,Ydata)
	    #upper_bound = mean_theory+std_theory
	    #lower_bound = mean_theory-std_theory
	    upper_bound = scipy.exp(scipy.log(mean_theory) + scipy.log(1.+std_theory/mean_theory)*rescale)
	    lower_bound = scipy.exp(scipy.log(mean_theory)+scipy.log(1.-std_theory/mean_theory)*rescale)
	    for i in range(0, len(lower_bound)):
		    if lower_bound[i]<=0:
			    lower_bound[i]=10.**(-16)
	    pylab.fill_between(Xtheory,lower_bound,y2=upper_bound,color=pointType[0],alpha=0.2)
	    #if model==m.Models.values()[1]:
	    #print Xdata, mean_theory-std_theory, mean_theory+std_theory, lower_bound
            #        break
	    for i, Ax in enumerate(axis_dep):
		    ax0[i] =i%2 and max(ax0[i],Ax) or min(ax0[i],Ax)
	pylab.axis(tuple(ax0))
	pylab.legend(loc=(-0.15,-0.52),ncol=3)

	if paper_fig:
		pylab.xlabel(model.theory.XnameTeX)
		pylab.ylabel(model.theory.YnameTeX)
	else:
		pylab.xlabel(model.theory.XnameTeX)
		pylab.ylabel(model.theory.YnameTeX)

    pylab.ion()
    pylab.show()
    
    return Xtheory,mean_theory, std_theory

def plot_error_traj(m, mean, std, fontSizeLabels=24, fontSizeLegend=12):
    # plot trajectories
    # use mean and std to calculate range of predictions
    pylab.ioff()
    ax0 = [1.e99,0,1.e99,0]
    for model in m.Models.values():
        pylab.figure()
        data_experiments = sorted(model.data.experiments)
        for independentValues in data_experiments:
            Xdata = model.data.X[independentValues]
            Ydata = model.data.Y[independentValues]
	    pointType = model.data.pointType[independentValues]
	    errorBar = model.data.errorBar[independentValues]
	    mean_theory = model.theory.Y(Xdata,mean,independentValues)
	    lower_bound = model.theory.Y(Xdata,mean-std,independentValues)
	    upper_bound = model.theory.Y(Xdata,mean+std,independentValues)
	    pylab.loglog(Xdata,Ydata,pointType[1])
	    lb = model.getLabel(model.theory.independentNames,independentValues)
	    pylab.errorbar(Xdata,Ydata, yerr=errorBar, fmt=pointType,label=lb)
	    pylab.loglog(Xdata,mean_theory,pointType[0])
	    axis_dep=model.getAxis(Xdata,Ydata)
	    for i in range(0, len(lower_bound)):
		    if lower_bound[i]<=0:
			    lower_bound[i]=10.**(-16)
		    if lower_bound[i] > upper_bound[i]:
			    upper = lower_bound[i]
			    lower = upper_bound[i]
			    upper_bound[i]=upper
			    lower_bound[i]=lower
	    pylab.fill_between(Xdata,lower_bound,upper_bound,color=pointType[0],alpha=0.2)
	    for i, Ax in enumerate(axis_dep):
		    ax0[i] =i%2 and max(ax0[i],Ax) or min(ax0[i],Ax)
	pylab.axis(tuple(ax0))
	pylab.rcParams.update({'legend.fontsize':fontSizeLegend})
	pylab.legend(loc=0)
	pylab.xlabel(model.theory.XnameTeX,fontsize=fontSizeLabels)
	pylab.ylabel(model.theory.YnameTeX,fontsize=fontSizeLabels)
	pylab.title(model.theory.title,fontsize=fontSizeLabels)
    pylab.ion()
    pylab.show()
