import os

particle="pi-"	
energy=[ i*10 for i in range(1,4)]
seed=0

for i in range(0,len(energy)):
	print particle +" "+ str(energy[i]) +" "+ str(seed)
	os.system("/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/script/runShowerProcessor.sh "+particle+" "+str(energy[i])+" "+str(seed))
