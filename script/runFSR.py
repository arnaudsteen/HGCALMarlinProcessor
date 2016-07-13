import os

energy=[1,2,3]
prefix=["multi_fsr","single"]

for i in energy:
    for j in prefix:
        os.system("./script/runFSRProcessor.sh "+str(i)+" "+j)
