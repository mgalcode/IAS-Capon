import matplotlib.pyplot as plt
from subroutines import *
from obspy.core import read
import scipy as sp




# ==== USER INPUT PARAMETER ===
nsamp = 8000
smin = -50.0
smax = 50.0
sinc = 0.5
cap_find = 120
cap_fave = 10
dl = 1

st = read('WRA.2010.001.10.00.BHZ.mseed',format='MSEED')
dic_meta = get_metadata('WRA.metadata')
# =============================


nr = st.count()              
dt = st[0].stats.delta
rx,ry = metric_mseed(st,dic_meta,nr)         
         
            
fk,sxopt,syopt,vel,rp,baz = FK(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz = IAS_FK(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz,maa,pwe = Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,dl,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz,maa,pwe = IAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,dl,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz = CAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)




#print arrival stats
print_stats(fk,threshold=0.15)



#generating figure
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
im = ax.imshow(fk.T,extent=[smin,smax, smax, smin],cmap='gist_stern_r',interpolation='none')
plt.title('Slowness Spectrum at %.03f +- %.03f[Hz]' %(cap_find/(nsamp*dt),cap_fave/(nsamp*dt)))
ax.set_xlim([smin,smax])
ax.set_ylim([smin,smax])
ax.set_xlabel('East/West Slowness [s/deg]')
ax.set_ylabel('North/South Slowness [s/deg]')
circle=plt.Circle((0,0),sp.sqrt((0.3*111.19)**2),color='w',fill=False,alpha=0.4)
plt.gcf().gca().add_artist(circle)
circle=plt.Circle((0,0),sp.sqrt((0.24*111.19)**2),color='w',fill=False,alpha=0.4)
plt.gcf().gca().add_artist(circle)
cbar = fig.colorbar(im)
cbar.set_label('relative power (dB)',rotation=270)

            
plt.show()


