import os

particle="e-"	
energy=[1,2,3,5,10,15,20,25,30]
seed=0

for i in range(0,len(energy)):
	print particle +" "+ str(energy[i]) +" "+ str(seed)
	os.system("/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/script/runShowerProcessor.sh "+particle+" "+str(energy[i])+" "+str(seed))
