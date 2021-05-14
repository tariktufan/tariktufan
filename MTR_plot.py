import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch as psd
#import time
def imag2real(ms):
    msM=np.sqrt(ms*np.conj(ms))
    mst=np.zeros(np.shape(msM))
    rm=np.real(ms)       
    msTeta=np.arccos(rm/msM)
    boyut=np.size(np.shape(ms));
    if boyut==2:
        s1,s2=np.shape(ms)
        for k in range(s1):
            for i in range(s2):
                if msTeta[k][i]>np.pi/2:
                    mst[k][i]=-msM[k][i]
                elif msTeta[k][i]<np.pi/2:
                    mst[k][i]=msM[k][i]
                else:
                    mst[k][i]=msM[k][i] #not true        
    elif boyut==1:
        s1=np.size(ms,0)
        for k in range(s1):
            if msTeta[k]>np.pi/2:
                mst[k]=-msM[k]
            elif msTeta[k]<np.pi/2:
                mst[k]=msM[k]
            else:
                mst[k]=msM[k] #not true         
               
    return mst, msTeta

def MTR_Size(MTR):
    if np.size(MTR,0)>1:
        kl=np.real(MTR[1,0:16]-MTR[0,0:16])
        i=1
        while kl[i] == 0:
            sz=i+1
            i+=1
    else:
        sz=4
    prc=(np.size(MTR,1)-(sz+2))/sz
    return int(sz),int(prc)

def MTR_WSort(MTR,DfNm,Karmasik=0,Fsnr=[0,100],MSsnr=[-10,10]):
    if Karmasik==0:
        MTR=np.real(MTR)
    sz,prc=MTR_Size(MTR)
    f1c=[]
    ms1c=[]
    for i in range(np.size(MTR,0)):        
        for k in range(sz):
            if MTR[i,k]==DfNm:
                f1c.extend(np.real(MTR[i,sz+2:sz+2+prc-1]/2/np.pi))
                ms1c.extend(MTR[i,sz+2+prc+k:sz+2+prc+(prc-1)*(sz-1)+k:sz-1])            
    #breakpoint()  
    #start=time.time()
    frq=f1c.copy()
    MS1=ms1c.copy()
    r=np.argsort(f1c)
    kn=-1
    for i in range(np.size(r,0)):
        ff=f1c[r[i]]
        msms=ms1c[r[i]];
        if ff>=Fsnr[0] and ff<=Fsnr[1] and msms>=MSsnr[0] and msms<=MSsnr[1]:
            kn+=1
            frq[kn]=ff
            MS1[kn]=msms
    frq=frq[0:kn]
    MS1=MS1[0:kn]
    #end=time.time()
    #print(end-start)
    if Karmasik==1:
        if min(np.isreal(MS1))==False:
            MSC=MS1
            MS1,*MSA=imag2real(MSC) #MSA angle btw Comp and Real
        #else:
        #    MSC=MS1
        #    MSA=np.zeros(np.shape(MS1));    
    return frq , MS1

def Cizdir(i,x,y,f1,MS,xt='x',yt='y',titl='0',fxx=0,Pxx=0):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.scatter(x,y,s=2)
    #ax1.xlabel(xt)
    #ax1.ylabel(yt)
    #ax1.title(titl)   
    sz = 16    
    plt.rc('font', size=sz)         # controls default text sizes
    plt.rc('axes', titlesize=sz)     # fontsize of the axes title
    plt.rc('axes', labelsize=sz)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=sz)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=sz)    # fontsize of the tick labels
    plt.rc('legend', fontsize=sz)    # legend fontsize
    plt.rc('figure', titlesize=sz)  # fontsize of the figure title
    ax1.grid()
    ax1.set_xlim([-2, 2])
    ax1.set_ylim([0, 13])
    xsnr=6 #len(f1)
    for j in range(xsnr):
        ax1.plot(MS[i,j],f1[j],color='green', marker='o',markersize=6)
    #plt.figure(i+25)
    ax2.semilogy(fxx[0],Pxx)
    for j in range(xsnr):
        ax2.plot(f1[j]*[1,1],[min(Pxx),max(Pxx)],color='red', linestyle='dashed',linewidth=3)
    ax2.set_xlim([0, 13])
    figManager = plt.get_current_fig_manager()
    figManager.full_screen_toggle()
    plt.show()


## LOAD CALCULATED SYSTEM PARAMETERS
mat = sio.loadmat('D:\MatlabR2016a\Codes\YSG\Bildiri\OperationalEnergy\System(1.14Hz)\System\SystemProp.mat')
MS=mat['MS']
f1=mat['f1']
## LOAD ESTIMATED FREQUENCY AND MODE SHAPE PARAMETERS OF SYSTEM 
resp = sio.loadmat('D:\CozumlemeSonuc\MTR_Y_Inp_01-101-115.mat')#('D:\MatlabR2016a\Codes\YSG\Bildiri\OperationalEnergy\System(1.14Hz)\MTR\MTR_Y_Inp_01-101-115.mat')
#resp = sio.loadmat('D:\CozumlemeSonuc\MTR_Y_Inp_01-101-115.mat')
MTR=resp['MTR']
yp=resp['yp']
fs=resp['fs']
sz,prc=MTR_Size(MTR)
k=np.real(np.unique(MTR[:,0:sz-1]))
k=k.astype(int)
for i in range(np.size(k,0)):
    DFM=k[i]
    Fsnr=[0,13]
    MSsnr=[-2,2];
    F1,MS1=MTR_WSort(MTR,DFM,0,Fsnr,MSsnr);
    fxx, Pxx = psd(yp[i], fs, nperseg=1024)
    TotalPower=sum(Pxx)*(fxx[0][1]-fxx[0][0])
    print(TotalPower)
    Cizdir(i,MS1,F1,f1,MS,'Mode Shape','Frequency(Hz)',str(DFM),fxx,Pxx)
    


