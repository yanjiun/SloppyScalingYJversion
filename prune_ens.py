import glob
import cPickle

filenames = glob.glob('Ensemble_low_T')
i = 0
corr_length = 1000

# to prune

for name in filenames:
	i += 1
	file = open(name,'rb')
	print "file opened"
	eof=False
	j = 1
	total_ens=[]
	total_Fs=[]
	while eof is not True:
		try:			
			ens1,ens_Fs1,ratio1 = cPickle.load(file)
			total_ens+=ens1
			total_Fs +=ens_Fs1
			j +=1
			print "ensemble loaded", j
                except:
			eof = True
			print "file not long enough"

	file.close()
	
	pruned_ens = []
	pruned_Fs = []
	corr_num = 1

	while corr_length*corr_num < len(total_ens)-1:
		pruned_ens.append(total_ens[corr_length*corr_num])
		pruned_Fs.append(total_Fs[corr_length*corr_num])
		corr_num += 1

   	file = open('pruned_lowT_'+str(i),'wb')
	cPickle.dump((pruned_ens,pruned_Fs),file)
	file.close()
	
	print "pruned ensemble saved: ", str(i)
