import os
import numpy as np
from scipy.io import loadmat,savemat
from scipy.fft import fft, ifft
from scipy import signal
from itertools import chain
#import matplotlib.pyplot as plt
def MeanZero(yp,fs=200,ek=1024):
    adet,lm_nz=np.shape(yp)
    lm=ek+fs*2
    if lm>lm_nz:
        lm=lm_nz
    for j in range(adet):
        for i in range(0,lm_nz-lm,lm):
            vct=list(range(i,i+lm))
            nd=vct[-1]
            k=yp[j,vct]
            kmean=np.average(k)
            yp[j,vct]=k-kmean
        if lm<lm_nz:
            vct=list(range(nd+1,lm_nz))
            k=yp[j,vct]
            kmean=np.average(k)
            yp[j,vct]=k-kmean
    return yp
Klasor = 'D:\\MatlabR2016a\\Codes\\YSG\\Makale_2BlindTestCantileverExample\\OrjinalVeri\\Veri1\\'#'D:\\MatlabR2016a\\Codes\\YSG\\Bildiri\\OperationalEnergy\\System(1.14Hz)\\Analyzed\\'
KayitKlasor='D:\CozumlemeSonuc\\'
# Verilerin inceleneceði klasör (tüm mat dosyalarý)
Forcetype=1
Datatype=0
xDof=6
xInpDof=0
dmp=0.05 # kayýt yapýlýrken isim oluþturulurken kullanýlýr.
l2=[]
l1=[]
HnklMtSz=1000
ndgr=1700
dgr03=4096
ResponseType=3
hata=0
HATA=[]
for f in os.listdir(Klasor):
    if len(f)>5:
        if f[-4:-1]+f[-1]=='.mat':
            print(f)
            DosyaAdi=[Klasor+f]
            C=f.split('-')
            a1=int(C[1])
            a2=int(C[2][0:-4])
            SrDrc=list(range(a1,a2+1))
            mat = loadmat(DosyaAdi[0])
            y=mat['y']
            yp=y.T
            fs=200
            dt=1/fs
            if np.size(yp,0)>np.size(yp,1):
                yp=yp.T
            ysz1,ysz2=np.shape(yp)
            kn=-1
            MTR=[]
            ##Hankel Deneme
            xDOF=np.size(yp,0);
            xOutDof=list(range(0,xDOF))
            #ShowPlot=0
            RemoveMean=1
            #DeTrend=0
            if RemoveMean==1:
                x=list(range(np.size(yp[1,:],0)))
                #plt.scatter(x,yp[1,:])
                yp=MeanZero(yp,fs)
                #yp=signal.detrend(yp)
                #plt.scatter(x,yp[1,:])
                #plt.show()
            dn=500
            if ysz2<5000:
                dn=50
            elif ysz2<20000:
                dn=200
            nmdn=int(np.ceil((ysz2-4*HnklMtSz)/dn))
            mktr=ysz1-1
            nmdnMAX=nmdn
            if nmdn>500:
                nmdn=500
            mtn=list(range(0,ysz1))
            lnx=np.size(mtn,0)-1;
            z=len(SrDrc)
            MTR=np.zeros([nmdn,z*(dgr03+1)+2])
            for j in range(nmdn):
                if j==nmdnMAX:
                    vctn=list(range(j*dn,ysz2))
                else:
                    vctn=list(range(j*dn,j*dn+4*HnklMtSz))
                print(str(j)+':'+str(min(vctn))+'-'+str(max(vctn)))
                ypa=yp[:,vctn]
                NFFT=dgr03*2-2
                fs1=(np.array(range(dgr03))*fs/NFFT).T
                w1=2*np.pi*fs1
                ms1=np.ones([mktr,int(NFFT/2)+1]);
                for i1 in range(mktr):
                    vct=fft(ypa[i1,:],NFFT)/fft(ypa[-1,:],NFFT);
                    a=int(NFFT/2)+1
                    ms1[i1,range(a)]=vct[range(a)];
                kn=kn+1;        
                wdgr=np.size(w1,0);
                if wdgr<dgr03:
                    hata=hata+1;
                    HATA.extend([j])
                    np.append(w1,-np.ones([dgr03-wdgr,1]))
                    np.append(ms1,np.zeros([len(mtn),dgr03-wdgr]))
                MSL=np.reshape(ms1.T,[1,(lnx)*dgr03])
                v=[];
                v.extend(SrDrc)
                v.extend([min(vctn),max(vctn)])
                v.extend(w1[0:dgr03].T)
                v.extend(MSL[0])
                MTR[kn,:]=v                
                if j==nmdn-1:
                    fname= KayitKlasor+'MTR_'+f
                    data = {'MTR': MTR, 'Forcetype': Forcetype, 'Datatype': Datatype, 'fs': fs, 'xDof': xDof, 'xInpDof': xInpDof, 'dmp': dmp, 'yp': yp, 'j': j}
                    savemat(fname, data)

                    

