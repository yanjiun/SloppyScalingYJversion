import sys
import getopt
import cPickle

import numpy
import numpy.linalg

super_computer = False
if super_computer is False:
    import scipy.fftpack
    import scipy.stats

import copy, time, cPickle, os
import logging
logger=logging.getLogger('Ensembles')
hdlr = logging.FileHandler('ensemble_run.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
import shutil

from fdjac import Jac

import SloppyScaling
reload(SloppyScaling)
import Utils
import WindowScalingInfo as WS
reload(WS)

# make ensembles
# inputs: SloppyScaling model class m, params, steps, max_run_hours, temperature = 0.5 Cbf...

# YJC: Feb 11, 2010 Modified way to make sampling matrix, putting in option where you can add diagonal term to JtJ.

def autocorrelation(series):
    """
    Return the normalized autocorrelation of a series using the FFT.
    """
    # We need to de-mean the series. Also, we want to pad with zeros to avoid
    #  assuming our series is periodic.
    f = numpy.fftpack.rfft(numpy.asarray(series)-numpy.mean(series), 
                           n = 2*len(series))
    # The inverse fft of |f|**2 is the autocorrelation
    ac = numpy.fftpack.irfft(abs(f)**2)
    # But we padded with zeros, so it's too long
    ac = ac[:len(series)]
    # And we need to adjust to remove the effect of the zeros we added
    ac *= len(series)/(len(series) - 1.0*numpy.arange(len(series)))
    # Return the normalized ac
    return ac/ac[0]

def ensemble(m, bf_params, start_params=None,steps=numpy.inf, max_run_hours=numpy.inf, step_scale=1.,sing_val_cutoff=None, diag_num = None,recalc_hess_alg=False,save_hours=numpy.inf,save_to=None,temperature=None,hess=None):

    """
    Generate a Bayesian ensemble of parameter sets consistent with the data in the model.  This is Ryan's importance sampling code modified to work with SloppyScaling

    Inputs:
     m -- SloppyScaling Model to generate the ensemble for
     params -- Initial parameter (numpy) array to start from 
     steps -- Maximum number of Monte Carlo steps to attempt
     max_run_hours -- Maximum number of hours to run
     temperature -- Temperature of the ensemble
     step_scale -- Additional scale applied to each step taken. step_scale < 1
                   results in steps shorter than those dictated by the quadratic
                   approximation and may be useful if acceptance is low.
     sing_val_cutoff -- Truncate the quadratic approximation at eigenvalues
                        smaller than this fraction of the largest.
     diag_num -- Number for diagonal matrix to add to hessian.
     save_hours --- If save_to is not None, the ensemble will be saved to
                       that file every 'save_hours' hours.
     save_to --- Filename to save ensemble to.
     hess --- how to calculate hessian

    Outputs:
     ens, ens_fes, ratio
     ens -- List of parameter sets in the ensemble
     ens_fes -- List of free energies for each parameter set
     ratio -- Fraction of attempted moves that were accepted

    The sampling is done by Markov Chain Monte Carlo, with a Metropolis-Hasting
    update scheme. The canidate-generating density is a gaussian centered on the
    current point, with axes determined by the hessian. For a useful 
    introduction see:
     Chib and Greenberg. "Understanding the Metropolis-Hastings Algorithm" 
     _The_American_Statistician_ 49(4), 327-335
    """

    # assuming parameters we start from are best fit params
    # T = 0.5*(cost at best fit)

    #bf_params = m.BestFit()[0]

    if temperature is None:
        temperature = m.Cost(bf_params)*2.0/len(bf_params) # acceptance T # divide by numbers of parameters!
    temp_sampling = 1. # sampling T

    # specify which parameter set to start stepping from
    if start_params is None:
        curr_params = bf_params
    else:
        curr_params = start_params

    #1. Calculate Hessian (approximate as JtJ)

    if hess:
        hessian = m.hess(curr_params)
    else:
        curr_jac = m.AnalyJac(curr_params)
        #curr_jac = Jac(curr_params,m.Residual)
        hessian = numpy.dot(curr_jac.transpose(),curr_jac)

    #2 do singular valued decomposition of Hessian to get sampling matrix to generate candidate moves (define subroutine for this)
    # cuttoff for very very small singular values
    if diag_num:
        diag = diag_num*numpy.eye(len(hessian))
    else:
        diag=None

    samp_mat = _sampling_matrix(hessian, cutoff=sing_val_cutoff, diag=diag, temperature=temp_sampling, step_scale=step_scale)

    # start stepping in directions from sampling matrix
    #curr_params = copy.deepcopy(params)
    curr_cost = m.Cost(curr_params)

    ens, ens_Fs = [curr_params],[curr_cost]

    log_params = False # put this here for now
    
    start_time = time.time()
    accepted_moves=0
    moves=0
    while moves < steps+1:
        if (time.time()-start_time) >= max_run_hours*3600:
            break
            
        #Generate the trial move from the quadratic approximation
        deltaParams = _trial_move(samp_mat)

        # got rid of this step, since we already scale the step...
        #Scale the trial move by the step_scale and the temperature
        #scaled_step = step_scale * deltaParams
        scaled_step = deltaParams
        
        if log_params:
            next_params = curr_params * numpy.exp(deltaParams)
        else:
            next_params = curr_params + deltaParams
        
        next_cost = m.Cost(next_params)

        if recalc_hess_alg and not numpy.isinf(next_cost):
            #next_jac = Jac(next_params, m.Residual)
            if hess:
                next_hess = m.hess(next_params)
            else:
                next_jac = m.AnalyJac(next_params)
                #next_jac = Jac(next_params,m.Residual)
                next_hess = numpy.dot(next_jac.transpose(),next_jac)
            next_samp_mat = _sampling_matrix(next_hess, sing_val_cutoff,diag, temp_sampling, step_scale)
            accepted = _accept_move_recalc_alg(curr_cost, hessian, next_cost,next_hess, deltaParams, temperature, diag=diag, cutoff=sing_val_cutoff)
        else:
            accepted = _accept_move(next_cost-curr_cost,temperature)
            print accepted

        if accepted:
            accepted_moves +=1.
            curr_params = next_params
            curr_cost = next_cost
            if recalc_hess_alg:
                hessian = next_hess
                samp_mat = next_samp_mat
        ens.append(curr_params)
        ens_Fs.append(curr_cost)
        moves+=1
    
        if (time.time() - start_time) >= save_hours*3600:
            ratio = 1.0*accepted_moves/moves
            _save_ens(ens,ens_Fs,ratio,save_to)
            save_hours += save_hours #save at next save_hours 
            # clear lists after saving
            ens, ens_Fs =[curr_params],[curr_cost]

    if len(ens)>1:
        ratio = 1.0*accepted_moves/moves
    else:
        ratio = 0

    # put in option for saving, just in case 
    if save_to is not None:
        _save_ens(ens, ens_Fs, ratio, save_to)

    print "took time (hours):" ,(time.time()-start_time)/3600.

    return ens, ens_Fs, ratio

# using Ryan's definition of creating the sampling matrix and of accepting moves with metropolis and of various subroutines

def _save_ens(ens, ens_Fs, ratio, save_to):
    #temp_name = save_to + '_temporary'
    f = file(save_to, 'ab')
    cPickle.dump((ens, ens_Fs, ratio), f, 2)
    #shutil.move(temp_name, save_to)
    f.close()
    logger.debug('Ensemble of length %i saved to %s.' % (len(ens), save_to))
    logger.debug('Acceptance ratio so far is %f.' % ratio)


def _sampling_matrix(hessian, cutoff=None, diag = None, temperature=1, step_scale=1):
    # basically need SVD of hessian - singular values and eigenvectors
    # hessian = u * diag(singVals) * vh
    #u, sing_vals, vh = numpy.linalg.svd(hessian)

    # scroll through the singular values and find the ones whose inverses will
    # be huge and set them to zero also, load up the array of singular values 
    # that we store
    # cutoff = (1.0/_.singVals[0])*1.0e03
    # double cutoff = _.singVals[0]*1.0e-02
    # when cutoff is set to zero it means that all values are included
    # cutoff*(sloppiest eigenvalue)
    
    if cutoff:
        u, sing_vals, vh = numpy.linalg.svd(hessian)
        cutoff_sing_val = cutoff * max(sing_vals)
        #when cutoff is set to zero it means that all values are included
        D = 1.0/numpy.maximum(sing_vals, cutoff_sing_val)
        samp_mat = numpy.transpose(vh)*numpy.sqrt(D)
    # instead of cutoff use another method, adding diagonal term to hessian
    elif diag is not None:
        u, sing_vals, vh = numpy.linalg.svd(hessian+diag)
        D = 1.0/sing_vals
        samp_mat = numpy.transpose(vh)*numpy.sqrt(D)
        cutoff_sing_val = diag[0,0]
    else: 
        u, sing_vals, vh = numpy.linalg.svd(hessian)
        D = 1.0/sing_vals
        samp_mat = numpy.transpose(vh)*numpy.sqrt(D)
        cutoff_sing_val = 0

    ## now fill in the sampling matrix ("square root" of the Hessian)
    ## note that sqrt(D[i]) is taken here whereas Kevin took sqrt(D[j])
    ## this is because vh is the transpose of his PT -JJW
    #samp_mat = numpy.transpose(vh) * numpy.sqrt(D)

    # Divide the sampling matrix by an additional factor such
    # that the expected quadratic increase in cost will be about 1.
    cutoff_vals = numpy.compress(sing_vals < cutoff_sing_val, sing_vals)
    if len(cutoff_vals):
        scale = numpy.sqrt(len(sing_vals) - len(cutoff_vals)
                           + sum(cutoff_vals)/cutoff_sing_val)
    else:
        scale = numpy.sqrt(len(sing_vals))

    samp_mat /= scale
    samp_mat *= step_scale
    samp_mat *= numpy.sqrt(temperature)

    return samp_mat

def _accept_move(delta_F, temperature):
    """
    Basic Metropolis accept/reject step.
    """
    p = numpy.random.rand()
    return (p < numpy.exp(-delta_F/temperature))
    
def _accept_move_recalc_alg(curr_F, curr_hess, next_F, next_hess, 
                            step, T, diag=None, cutoff=None):
    """
    Accept/reject when each the sampling matrix is recalculated each step.
    """
    pi_x = numpy.exp(-curr_F/T)
    # This is the current location's covariance sampling matrix
    # try making the sigma_curr_inv = diagonal matrix of singular values
    
    # using this method is right now is problematic, let's make use of the SVD of the original thing... 
    #sigma_curr = numpy.dot(curr_samp_mat, numpy.transpose(curr_samp_mat))
    #sigma_curr_inv = numpy.linalg.inv(sigma_curr)

    if diag is not None:
        curr_jtj = curr_hess + diag
        next_jtj = next_hess + diag
    elif cutoff:
        u1, sv1, vh1 = numpy.linalg.svd(curr_hess)
        sv1 = max(sv1, sv1*cutoff)
        curr_jtj = numpy.dot(u1*sv1,vh1)
        u2, sv2, vh2 = numpy.linalg.svd(next_hess)
        sv2 = max(sv2, sv2*cutoff)
        next_jtj = numpy.dot(u2*sv2,vh2)
    else:
        curr_jtj = curr_hess
        next_jtj = next_hess

    # This is the transition probability from the current point to the next.
    #q_x_to_y = numpy.exp(-_quadratic_cost(step, sigma_curr_inv))\
    #        / numpy.sqrt(numpy.linalg.det(sigma_curr))
 
    # sampling matrix has problem of making the det(sigma_curr) really small, let's try using the hessian for the problem right now...
    q_x_to_y = numpy.exp(-_quadratic_cost(step,curr_jtj))*numpy.sqrt(numpy.linalg.det(curr_jtj))
    
    pi_y = numpy.exp(-next_F/T)
    
    # old stuff that causes problems
    #sigma_next = numpy.dot(next_samp_mat, numpy.transpose(next_samp_mat))
    #sigma_next_inv = numpy.linalg.inv(sigma_next)
    #q_y_to_x = numpy.exp(-_quadratic_cost(-step, sigma_next_inv))\
    #        / numpy.sqrt(numpy.linalg.det(sigma_next))

    q_y_to_x = numpy.exp(-_quadratic_cost(-step,next_jtj))*numpy.sqrt(numpy.linalg.det(next_jtj))

    p = numpy.random.rand()
    #accepted = (pi_y*q_y_to_x)/(pi_x*q_x_to_y)
    accepted = (numpy.exp(-(next_F-curr_F)/T))*(q_y_to_x/q_x_to_y)
    if accepted < 0:
        print "error, negative value in acceptance probability"
    did_accepted = p < accepted
    print numpy.exp(-(next_F-curr_F)/T), q_y_to_x/q_x_to_y, accepted, did_accepted
    import sys
    sys.stdout.flush()

    return p < abs(accepted)

def _trial_move(sampling_mat):
    randVec = numpy.random.randn(len(sampling_mat))
    trialMove = numpy.dot(sampling_mat, randVec)

    return trialMove

def _quadratic_cost(trialMove, hessian):
    """
    The cost from the quadratic approximation of a trialMove, given the hessian.

    (Note: the hessian here is assumed to be the second derivative matrix of the
     cost, without an additional factor of 1/2.)
    """
    quadratic = 0.5*numpy.dot(numpy.transpose(trialMove), 
                              numpy.dot(hessian, trialMove))
    return quadratic


def main(argv):
    random_seed = 0
    # get random number from command line input, and seed random number generator
    options, remainder = getopt.getopt(argv,'s:')
    for opt, arg in options:
        if opt == "-s":
            random_seed = int(arg)

    numpy.random.seed(seed=random_seed)
    
    # load best fit parameters to start from
    bf_file = open('bf_params_fixed_w_zeta.dat','rb')
    bf_params=cPickle.load(bf_file)
    bf_file.close()
    
    #start_file = open('bf_params_new.dat','rb')
    #start_params = cPickle.load(start_file)
    #start_file.close()
    start_params = bf_params

    # make jointModule
    for moduleName in WS.moduleNames:
        mod_name = moduleName + "Module"
        exec("import "+ mod_name+" as "+moduleName)
        exec("reload ("+moduleName+")")

    jointModuleName = "".join(WS.moduleNames)

    jointModule = SloppyScaling.CompositeModel(jointModuleName)
    for moduleName in WS.moduleNames:
        exec("obj_mod = "+moduleName)
        m = getattr(obj_mod,moduleName)
        jointModule.InstallModel(moduleName, m)

    #jointModule.SetFixedParams(fixedParameters=[('zeta',0.63)])

    size = 10**6
    fileName = "Ensemble_small_T_fixed_all_"+str(random_seed)
    ens, Fs, ratio=ensemble(jointModule,bf_params,start_params=bf_params,steps=size,step_scale=1.,save_hours = 1,recalc_hess_alg=True,save_to = fileName,temperature=1.)
    print "acceptance ratio:", ratio

if __name__ == "__main__":
    main(sys.argv[1:])
