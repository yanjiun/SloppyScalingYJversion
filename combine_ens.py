import glob
import cPickle

filenames = glob.glob('pruned_ens_jeff*')
i = 0

# to combine ensembles
total_ens=[]
total_Fs=[]

for name in filenames:
	i += 1
	file = open(name,'r')
	print "file opened:", str(i)
	eof=False
	while eof is not True:
		try:			
			ens,ensFs=cPickle.load(file)
			total_ens += ens
			total_Fs += ensFs
			print "Ensemble_loaded"
                except:
			eof = True
			print "file not long enough"
	
	file.close()
	

file = open('combined_ensemble_jeff','wb')
cPickle.dump((total_ens,total_Fs),file)
file.close()
		
