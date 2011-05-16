Gimport glob
import cPickle

filenames = glob.glob('ensembles_ept/Ensemble_0506_*')
i = 0

# to load and save new parameters

for name in filenames:
	i += 1
	file = open(name,'r')
	eof=False
	j = 1
	while eof not True:
		try:			
			exec("ens"+str(j)+", ens_Fs"+str(j)", ratio"+str(j)+" = cPickle.load(file)")
			j +=1     
                except:
			eof = True
			print "file not long enough"
	
	file.close()
	file = open('start_params'+str(i)+'.dat','wb')
	cPickle.dump(ens2[-1],file)
	file.close()

	print "file "+str(i)+" loaded and new start params saved"
	for l in range(1,j):
		print "acceptance ratio:"
		exec("print ratio"+str(j))
		
