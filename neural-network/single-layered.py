#Neural Network
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt

names = ['scodeno','clump_th','un_cell_sz','un_cell_sh','marg_adhe','s_epi_cell_sz','bare_uc','bland_chrm','nor_nucl','mitoses','class']
df = pd.read_csv("C:\\Users\\Omkar\\Desktop\\bcw.data", names= names)


del df['scodeno']

cls=[]
for cl in df['class']:
    if cl == 2:
        cls.append(0)
    else:
        cls.append(1)
    
for miss in df['bare_uc']:
    if miss=='?':
        i=df[df.bare_uc == miss].index
        df['bare_uc'][i]=0
        


del df['class']

matr=[]
names = ['clump_th','un_cell_sz','un_cell_sh','marg_adhe','s_epi_cell_sz','bare_uc','bland_chrm','nor_nucl','mitoses']

for n in names:
    matr.append(df[n].tolist())

    
matrx =np.array(list(np.float_(matr)))
inp = matrx.T

clss =[cls]
clss = np.array(clss).T



synp0 = 2*np.random.random((9,1)) - 1



def sig(x,der=False):
    op = 1/(1+np.exp(-x*0.001))
    if der==False:
        return op
    else:
        return x*(1-x)

for i in range(100000):
    l0=inp
    l1 = sig(np.dot(l0,synp0))
    err = (clss - l1)
    l1_del = err * sig(l1,True)
    synp0 += np.dot(inp.T,l1_del)
    #print(synp0)

rng =l1.max()-l1.min()
result =[]  

for i in l1:
    if i > (0.6*rng + l1.min()):
        result.append(1)
    else:
        result.append(0)
cnt=0

for i,j in zip(result,cls):
    if i==j:
        cnt=cnt+1

print('Percentage correctness:',cnt*100.0/699)

x_a = sig(np.dot(inp,synp0))
y_a = np.array(result)

plt.scatter(x_a,y_a)        

