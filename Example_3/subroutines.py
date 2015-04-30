# -*- coding: utf-8 -*-
def Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,find,fave,delta,dl,overlap,taper):
    import util as ut
    import numpy as np
    import scipy as sp
    from scipy.weave import converters





    '''
    function to calculate the capon algorithm.
    variables used:
    ---------------
    nsamp ... amount of snapshots of the data sequence(in a window!)
    nr    ... amount of array stations/sensors
    rx    ... vector of r_x component for all stations. rx is in [deg]. reference station of the array is always first sensor in the sequence.
    ry    ... vector of r_y component ..
    nwin  ... amount of time windows that need to be extracted from the time series
    st    ... matrix of data. first element is sensor station and second element is time series. dimensions are st[nr][nsamp]
    smin  ... minimum slowness that is used to perform the slowness search
    smax  ... max slowness. values -50 and 50 are standard for these smin and smax
    sinc  ... the increment of the slowness. the lower the more accurate the finer the grid-search is. default value should be 0.5
    find  ... frequency at which the DOA is performed. warning: it is not [Hz]. It connected to the frequency by the equation freq = find/(nsamp*delta)
              where freq is in [Hz]. 
    fave  ... frequency to average over. DOA is performed at find +- fave respectively at freq +- df
    delta ... time steps. dependent on array.
    freq  ... frequency at which the DOA is performed in [Hz]
    df    ... averaging frequency
    kmin  ... minimum wavenumber. it is used later for the calculation as all angles need to be considered of the arriving signal.
              the steering vector is basically exp(ikr) where r is the interstation distance in [deg].
    kmax  ... max wavenumber
    kinc  ... wavenumber increment
    nk    ... amount of slowness steps/ respectively wavenumber steps from min to ma
    '''
    # Variable initialisation. If more than 1 time window is calculated nsamp is now the number of snapshots for each timewindow
    freq = find/(nsamp*delta)
    df = fave/(nsamp*delta)
    print 'Capon DOA estimation is performed at:','freq',freq,'+-',df
    kmin = 2*sp.pi*smin*freq
    kmax = 2*sp.pi*smax*freq
    kinc = 2*sp.pi*sinc*freq
    nk = int(((kmax-kmin)/kinc+0.5)+1) 
    
    if (overlap == True):  
        nwin = int(np.array(st[0].stats.npts/nsamp))*2-1
        xt = np.zeros((nr,nwin,nsamp))
        # splitting time sequence in nwin parts and storing result in xt[][][]
        # removing mean from time windows and applying tapper algorithm to surpress impact of windowing
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp/2:(j+2)*nsamp/2]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)
    else:
        nwin = int(np.array(st[0].stats.npts/nsamp))
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp:(j+1)*nsamp]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)


   
    
    # generating FFT for time series and calculating spectral matrix where smr is the real part and smi is the imag part
    # FFT has to be performed for real valued data!!! rfft from numpy will perform the job (rfft will return only positive freq.)
    # the loop over n averages the frequency components around the central frequency find
    # old code
    '''  
    rft = np.zeros((nr,nsamp/2+1))
    ift = np.zeros((nr,nsamp/2+1))
    for i in range(nwin):
         for j in range(nr):
             tp = np.fft.rfft(xt[j][i],nsamp)
             rft[j]=(tp.real)
             ift[j]=(tp.imag)

         for k in range(nr):
             for m in range(nr):
                 for n in range(find-fave,find+fave+1):
                     smr[k][m] += rft[k][n]*rft[m][n]+ift[k][n]*ift[m][n]
                     smi[k][m] += rft[m][n]*ift[k][n]-rft[k][n]*ift[m][n]
    '''
    
    # new code
    rft = np.zeros((nwin,nr,nsamp/2+1))
    ift = np.zeros((nwin,nr,nsamp/2+1))
    smr = np.zeros((nr,nr))
    smi = np.zeros((nr,nr))
    for i in range(nwin):
        for j in range(nr):
            tp = np.fft.rfft(xt[j][i],nsamp)
            rft[i][j]=(tp.real)
            ift[i][j]=(tp.imag)

    code = """
           int i,j,l,n;
               for(n=0;n<nwin;n++){
     	       for(i=0;i<nr;i++){
	           for(j=0;j<nr;j++){
		       for(l=find-fave;l<=find+fave;l++){
                           smr(i,j)+=rft(n,i,l)*rft(n,j,l)+ift(n,i,l)*ift(n,j,l);
                           smi(i,j)+=rft(n,j,l)*ift(n,i,l)-rft(n,i,l)*ift(n,j,l);
                }}}}
                """
    sp.weave.inline(code,['smr','smi','rft','ift','nr','find','fave','nwin'],type_converters=converters.blitz,compiler = 'gcc')

    smr /= nwin
    smi /= nwin
    
    
    pwe = 0.0
    for i in range(nr):
        pwe += np.abs(smr[i][i]+1j*smi[i][i])**2
    pwe /= float(nr)
    
    print 'Diagonal Loading On!'
    mi = np.identity(nr)
    smr += mi*smr.trace()/(nsamp)*1



    # applying weigths (check capon1969 for more info on the equation)
    wmean = 0.0
    w = np.zeros(nr)
    for i in range(nr):
        w[i] = (smr[i][i]*smr[i][i]+smi[i][i]*smi[i][i])**(-0.25)
        wmean += 1.0/(w[i]**2)
    wmean /= float(nr)*(2*fave+1)
    for i in range(nr):
        for j in range(nr):
            smr[i][j] *=w[i]*w[j]
            smi[i][j] *=w[i]*w[j]

    # calcilate invers of cross power spectral density / spectral covariance matrix
    itx = np.linalg.inv(smr+1j*smi)
    ismr = itx.real
    ismi = itx.imag
    

    
    # calculating the spectrum here that is stored in fk. grid search is performed over k_x and k_y
    # vectoring is applied to improve performance. rx for instance is a vector that is computed in the exp function. (thx scipy)
    ''' old code
    steer = np.zeros(nr,dtype=complex)
    fk = np.zeros((nk,nk),dtype=complex)
    for i in range(nk):
        kx=-(kmin+float(i*kinc))
        for j in range(nk):
            ky=-(kmin+float(j*kinc))
            steer=np.exp(1j*(kx*(rx[0]-rx)+ky*(ry[0]-ry)))
            fk[i][j]=1. / steer.T.conj().dot(itx).dot(steer)
    '''
        
    # new code
    # kx and ky are reverted so we can plot it North being the x axis
    # this revert will be later taken into account so we dont get false vaules
    fk = np.zeros((nk,nk))
    code1 = """
            int i,j,m,n,fcnt;
            float kx,ky,arg;
            fcnt=nr;
	    for(i=0;i<nk;i++){
	        kx=-(kmin+(double)i*kinc);
	        for(j=0;j<nk;j++){
	            ky=-(kmin+(double)j*kinc);
 	            fk(i,j)=0.0;
	    	    for(m=0;m<fcnt;m++){
		        fk(i,j)+=ismr(m,m);
		        for(n=m+1;n<fcnt;n++){
		            arg=kx*(rx(m)-rx(n))+ky*(ry(m)-ry(n));
			        fk(i,j)+=2.0*(ismr(m,n)*cos(arg)-ismi(m,n)*sin(arg));
		        }
		    }
 	        fk(i,j)=1.0/fk(i,j);
	        }
	    }
	    """
    sp.weave.inline(code1,['fk','ismr','ismi','nk','rx','ry','kinc','kmin','nr'],type_converters=converters.blitz,compiler = 'gcc',headers=['<math.h>'])

    # search for the maximum in order to normalize later. this is programmed very performance poor (should be revised later).
    # further recalculations are done to generate output

    max = 0.0

    for i in range(nk):
        for j in range(nk):
            if(fk[i][j].real>max):
                max=fk[i][j].real
                # x and y are set back again to get correct values
                sxopt=smin+i*sinc
                syopt=smin+j*sinc

    rp = sp.sqrt(sxopt**2+syopt**2)
    #6371*2*pi/360 -> 111.19 [km/deg]
    vel = 111.19/rp
    baz=np.math.atan2(sxopt,syopt)*180.0/3.1415926
    if(baz<0.0):
       baz+=360.0

    #return fk.real/max,sxopt,syopt,vel,rp,baz 
    return 10*sp.log10(fk.real/max),sxopt,syopt,vel,rp,baz,sp.log10(wmean),sp.log10(pwe)
    
    
    
def IAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,find,fave,delta,dl,overlap,taper):
    import util as ut
    import numpy as np
    import scipy as sp
    from scipy.weave import converters

    freq = find/(nsamp*delta)
    df = fave/(nsamp*delta)
    print 'IAS Capon DOA estimation is performed at:','freq',freq,'+-',df
    kmin = 2*sp.pi*smin*freq
    kmax = 2*sp.pi*smax*freq
    kinc = 2*sp.pi*sinc*freq
    nk = int(((kmax-kmin)/kinc+0.5)+1) 
    
    if (overlap == True):  
        nwin = int(np.array(st[0].stats.npts/nsamp))*2-1
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp/2:(j+2)*nsamp/2]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)
    else:
        nwin = int(np.array(st[0].stats.npts/nsamp))
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp:(j+1)*nsamp]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)

    smr = np.zeros((2*fave+1,nr,nr))
    smi = np.zeros((2*fave+1,nr,nr))

    rft = np.zeros((nwin,nr,nsamp/2+1))
    ift = np.zeros((nwin,nr,nsamp/2+1))
    for i in range(nwin):
        for j in range(nr):
            tp = np.fft.rfft(xt[j][i],nsamp)
            rft[i][j]=(tp.real)
            ift[i][j]=(tp.imag)
    
    code = """
           int i,j,l,n;
               for(n=0;n<nwin;n++){
     	       for(i=0;i<nr;i++){
	           for(j=0;j<nr;j++){
		       for(l=find-fave;l<=find+fave;l++){
                           smr(l-find+fave,i,j)+=rft(n,i,l)*rft(n,j,l)+ift(n,i,l)*ift(n,j,l);
                           smi(l-find+fave,i,j)+=rft(n,j,l)*ift(n,i,l)-rft(n,i,l)*ift(n,j,l);
                }}}}
                """
    sp.weave.inline(code,['smr','smi','rft','ift','nr','find','fave','nwin'],type_converters=converters.blitz,compiler = 'gcc')

    smr /= nwin
    smi /= nwin
    
    fw = 0.0
    fe = 0.0
    for m in range(2*fave+1):
        wmean = 0.0
        w = np.zeros(nr)
        for i in range(nr):
            w[i] = (smr[m][i][i]*smr[m][i][i]+smi[m][i][i]*smi[m][i][i])**(-0.25)
            wmean += 1.0/(w[i]**2)
            fw    += 1.0/(w[i]**2)
            fe    += np.abs(smr[m][i][i]+1j*smi[m][i][i])**2
        for i in range(nr):
            for j in range(nr):
                smr[m][i][j] *=w[i]*w[j]
                smi[m][i][j] *=w[i]*w[j]
    fw /= nr*(2*fave+1)
    fe /= nr
                
                
    print 'Diagonal Loading On!'
    mi = np.identity(nr)
    for i in range(2*fave+1):
        smr[i] += mi*smr[i].trace()/(nsamp)*dl
        
        

    tx = smr +1j*smi
    itx = np.zeros((2*fave+1,nr,nr),dtype=complex)
    for m in range(2*fave+1):
        itx[m] = np.linalg.inv(tx[m])
    ismr = itx.real
    ismi = itx.imag


    tfk = np.zeros((nk,nk))
    fk = np.zeros((nk,nk))

    code1 = """
            int i,j,m,n,fcnt,g;
            float kx,ky,arg,freq;
            fcnt=nr;
            for(g=find-fave;g<=find+fave;g++){
            freq=(double)g/((double)nsamp*delta);
            kmin=2*3.1415926*freq*smin;
            kinc=2*3.1415926*freq*sinc;
	    for(i=0;i<nk;i++){
	        kx=-(kmin+(double)i*kinc);
	        for(j=0;j<nk;j++){
	           ky=-(kmin+(double)j*kinc);
 	            fk(i,j)=0.0;
	    	    for(m=0;m<fcnt;m++){
		        fk(i,j)+=ismr(g-find+fave,m,m);
		        for(n=m+1;n<fcnt;n++){
		            arg=kx*(rx(m)-rx(n))+ky*(ry(m)-ry(n));
			        fk(i,j)+=2.0*(ismr(g-find+fave,m,n)*cos(arg)-ismi(g-find+fave,m,n)*sin(arg));
		        }
		    }
 	        fk(i,j)=1.0/fk(i,j);
	        }
	    }
	    tfk+=fk;
	    }
	    """
    sp.weave.inline(code1,['fk','tfk','ismr','ismi','nk','rx','ry','kinc','kmin','nr','nsamp','delta','smin','sinc','fave','find'],type_converters=converters.blitz,compiler = 'gcc',headers=['<math.h>'])
    fk=tfk

    max = 0.0
    for i in range(nk):
        for j in range(nk):
            if(fk[i][j].real>max):
                max=fk[i][j].real
                sxopt=smin+i*sinc
                syopt=smin+j*sinc
            if(fk[i][j].real<0):
                fk[i][j]=0
    rp = sp.sqrt(sxopt**2+syopt**2)
    vel = 111.19/rp
    baz=np.math.atan2(sxopt,syopt)*180.0/3.1415926
    if(baz<0.0):
       baz+=360.0
    
    res = 10*sp.log10(fk.real/max)
    return res,sxopt,syopt,vel,rp,baz,sp.log10(fw),sp.log10(fe)


    
    
def FK(nsamp,nr,rx,ry,st,smin,smax,sinc,find,fave,delta,overlap,taper):
    import util as ut
    import numpy as np
    import scipy as sp
    from scipy.weave import converters


    freq = find/(nsamp*delta)
    df = fave/(nsamp*delta)
    print 'FK DOA estimation is performed at:','freq',freq,'+-',df
    kmin = 2*sp.pi*smin*freq
    kmax = 2*sp.pi*smax*freq
    kinc = 2*sp.pi*sinc*freq
    nk = int(((kmax-kmin)/kinc+0.5)+1) 
    

    if (overlap == True):  
        nwin = int(np.array(st[0].stats.npts/nsamp))*2-1
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp/2:(j+2)*nsamp/2]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)
    else:
        nwin = int(np.array(st[0].stats.npts/nsamp))
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp:(j+1)*nsamp]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)


    smr = np.zeros((nr,nr))
    smi = np.zeros((nr,nr))


    rft = np.zeros((nwin,nr,nsamp/2+1))
    ift = np.zeros((nwin,nr,nsamp/2+1))
    for i in range(nwin):
        for j in range(nr):
            tp = np.fft.rfft(xt[j][i],nsamp)
            rft[i][j]=(tp.real)
            ift[i][j]=(tp.imag)
    
    code = """
           int i,j,l,n;
               for(n=0;n<nwin;n++){
     	       for(i=0;i<nr;i++){
	           for(j=0;j<nr;j++){
		       for(l=find-fave;l<=find+fave;l++){
                           smr(i,j)+=rft(n,i,l)*rft(n,j,l)+ift(n,i,l)*ift(n,j,l);
                           smi(i,j)+=rft(n,j,l)*ift(n,i,l)-rft(n,i,l)*ift(n,j,l);
                }}}}
                """
    sp.weave.inline(code,['smr','smi','rft','ift','nr','find','fave','nwin'],type_converters=converters.blitz,compiler = 'gcc')



    denom = nwin*(2*fave+1.0)
    smr /=denom
    smi /=denom


    wmean = 0.0
    w = np.zeros(nr)
    for i in range(nr):
        w[i] = (smr[i][i]*smr[i][i]+smi[i][i]*smi[i][i])**(-0.25)
        wmean += 1.0/(w[i]**2)
    wmean /= float(nr)
    for i in range(nr):
        for j in range(nr):
            smr[i][j] *=w[i]*w[j]
            smi[i][j] *=w[i]*w[j]


    fk = np.zeros((nk,nk))
    code1 = """
            int i,j,m,n,fcnt;
            float kx,ky,arg;
            fcnt=nr;
	    for(i=0;i<nk;i++){
	        kx=-(kmin+(double)i*kinc);
	        for(j=0;j<nk;j++){
	            ky=-(kmin+(double)j*kinc);
 	            fk(i,j)=0.0;
	    	    for(m=0;m<fcnt;m++){
		        fk(i,j)+=smr(m,m);
		        for(n=m+1;n<fcnt;n++){
		            arg=kx*(rx(m)-rx(n))+ky*(ry(m)-ry(n));
			        fk(i,j)+=2.0*(smr(m,n)*cos(arg)-smi(m,n)*sin(arg));
		        }
		    }
	        }
	    }
	    """
    sp.weave.inline(code1,['fk','smr','smi','nk','rx','ry','kinc','kmin','nr'],type_converters=converters.blitz,compiler = 'gcc',headers=['<math.h>'])
    

    max = 0.0
    for i in range(nk):
        for j in range(nk):
            if(fk[i][j].real>max):
                max=fk[i][j].real
                sxopt=smin+i*sinc
                syopt=smin+j*sinc
            if(fk[i][j].real<0):
                fk[i][j]=0
    rp = sp.sqrt(sxopt**2+syopt**2)
    vel = 111.19/rp
    baz=np.math.atan2(sxopt,syopt)*180.0/3.1415926
    if(baz<0.0):
       baz+=360.0
    return 10*sp.log10(fk.real/max),sxopt,syopt,vel,rp,baz  
    
    
    
def IAS_FK(nsamp,nr,rx,ry,st,smin,smax,sinc,find,fave,delta,overlap,taper):
    import util as ut
    import numpy as np
    import scipy as sp
    from scipy.weave import converters

    
    freq = find/(nsamp*delta)
    df = fave/(nsamp*delta)
    print 'IAS FK DOA estimation is performed at:','freq',freq,'+-',df
    kmin = 2*sp.pi*smin*freq
    kmax = 2*sp.pi*smax*freq
    kinc = 2*sp.pi*sinc*freq
    nk = int(((kmax-kmin)/kinc+0.5)+1) 


    if (overlap == True):  
        nwin = int(np.array(st[0].stats.npts/nsamp))*2-1
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp/2:(j+2)*nsamp/2]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)
    else:
        nwin = int(np.array(st[0].stats.npts/nsamp))
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp:(j+1)*nsamp]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)


    smr = np.zeros((2*fave+1,nr,nr))
    smi = np.zeros((2*fave+1,nr,nr))


    

    rft = np.zeros((nwin,nr,nsamp/2+1))
    ift = np.zeros((nwin,nr,nsamp/2+1))
    for i in range(nwin):
        for j in range(nr):
            tp = np.fft.rfft(xt[j][i],nsamp)
            rft[i][j]=(tp.real)
            ift[i][j]=(tp.imag)
    
    code = """
           int i,j,l,n;
               for(n=0;n<nwin;n++){
     	       for(i=0;i<nr;i++){
	           for(j=0;j<nr;j++){
		       for(l=find-fave;l<=find+fave;l++){
                           smr(l-find+fave,i,j)+=rft(n,i,l)*rft(n,j,l)+ift(n,i,l)*ift(n,j,l);
                           smi(l-find+fave,i,j)+=rft(n,j,l)*ift(n,i,l)-rft(n,i,l)*ift(n,j,l);
                }}}}
                """
    sp.weave.inline(code,['smr','smi','rft','ift','nr','find','fave','nwin'],type_converters=converters.blitz,compiler = 'gcc')


    smr /= nwin
    smi /= nwin


    for m in range(2*fave+1):
        wmean = 0.0
        w = np.zeros(nr)
        for i in range(nr):
            w[i] = (smr[m][i][i]*smr[m][i][i]+smi[m][i][i]*smi[m][i][i])**(-0.25)
            wmean += 1.0/(w[i]**2)
        wmean /= float(nr)
        for i in range(nr):
            for j in range(nr):
                smr[m][i][j] *=w[i]*w[j]
                smi[m][i][j] *=w[i]*w[j]


    tfk = np.zeros((nk,nk))
    fk = np.zeros((nk,nk))

    code1 = """
            int i,j,m,n,fcnt,g;
            float kx,ky,arg,freq;
            fcnt=nr;
            for(g=find-fave;g<=find+fave;g++){
            freq=(double)g/((double)nsamp*delta);
            kmin=2*3.1415926*freq*smin;
            kinc=2*3.1415926*freq*sinc;
	    for(i=0;i<nk;i++){
	        kx=-(kmin+(double)i*kinc);
	        for(j=0;j<nk;j++){
	           ky=-(kmin+(double)j*kinc);
 	            fk(i,j)=0.0;
	    	    for(m=0;m<fcnt;m++){
		        fk(i,j)+=smr(g-find+fave,m,m);
		        for(n=m+1;n<fcnt;n++){
		            arg=kx*(rx(m)-rx(n))+ky*(ry(m)-ry(n));
			        fk(i,j)+=2.0*(smr(g-find+fave,m,n)*cos(arg)-smi(g-find+fave,m,n)*sin(arg));
		        }
		    }
	        }
	    }
	    tfk+=fk;
	    }
	    """
    sp.weave.inline(code1,['fk','tfk','smr','smi','nk','rx','ry','kinc','kmin','nr','nsamp','delta','smin','sinc','fave','find'],type_converters=converters.blitz,compiler = 'gcc',headers=['<math.h>'])
    fk=tfk


    max = 0.0
    for i in range(nk):
        for j in range(nk):
            if(fk[i][j].real>max):
                max=fk[i][j].real
                sxopt=smin+i*sinc
                syopt=smin+j*sinc
            if(fk[i][j].real<0):
                fk[i][j]=0
    rp = sp.sqrt(sxopt**2+syopt**2)
    vel = 111.19/rp
    baz=np.math.atan2(sxopt,syopt)*180.0/3.1415926
    if(baz<0.0):
       baz+=360.0
    res = 10*sp.log10(fk.real/max)
    return res,sxopt,syopt,vel,rp,baz  


def CAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,find,fave,delta,overlap,taper):
    import util as ut
    import numpy as np
    import scipy as sp
    from scipy.weave import converters


    freq = find/(nsamp*delta)
    df = fave/(nsamp*delta)
    print 'CAS Capon DOA estimation is performed at:','freq',freq,'+-',df
    kmin = 2*sp.pi*smin*freq
    kmax = 2*sp.pi*smax*freq
    kinc = 2*sp.pi*sinc*freq
    nk = int(((kmax-kmin)/kinc+0.5)+1) 
    dt = delta
    
    if (overlap == True):  
        nwin = int(np.array(st[0].stats.npts/nsamp))*2-1
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp/2:(j+2)*nsamp/2]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)
    else:
        nwin = int(np.array(st[0].stats.npts/nsamp))
        xt = np.zeros((nr,nwin,nsamp))
        for i in range(nr):
            for j in range(nwin):
                xt[i][j] = st[i][j*nsamp:(j+1)*nsamp]
                xt[i][j] -= np.mean(xt[i][j])
                if(taper==True): xt[i][j] *= np.hanning(nsamp)

    
    tf = float()
    tr = np.empty(nr)
    ti = np.empty(nr)
    iismr = np.zeros((nr,nr))
    iismi = np.zeros((nr,nr))
    xsmr = np.zeros((2*fave+1,nr,nr))
    xsmi = np.zeros((2*fave+1,nr,nr))
    rft = np.zeros((nr,nsamp/2+1))
    ift = np.zeros((nr,nsamp/2+1))
    w = np.zeros(nr)
    fk = np.zeros((nk,nk))

    for i in range(nwin):
         for j in range(nr):
             tp = np.fft.rfft(xt[j][i],nsamp)
             rft[j]=tp.real
             ift[j]=tp.imag
         for k in range(nr):
             for m in range(nr):
                 for n in range(find-fave,find+fave+1):
                     xsmr[n-find+fave][k][m] += rft[k][n]*rft[m][n]+ift[k][n]*ift[m][n]
                     xsmi[n-find+fave][k][m] += rft[m][n]*ift[k][n]-rft[k][n]*ift[m][n]
    

    vari = ['xsmr','xsmi','fk','rft','ift','nr','find','fave','nwin',
            'nsamp','iismr','iismi','rx','ry','tr','ti','freq','tf',
            'dt','w','nk','smin','sinc']

    weave_info = { 'headers'       : ['<math.h>',
                                     '<gsl/gsl_blas.h>',
                                     '<gsl/gsl_complex.h>',
                                     '<gsl/gsl_complex_math.h>',
                                     '<gsl/gsl_matrix_complex_double.h>',
                                     '<gsl/gsl_permute_complex_double.h>',
                                     '<gsl/gsl_linalg.h>'],
                  'include_dirs'  : ['/sw/include/'],
                  'library_dirs'  : ['/sw/lib/'],
                  'libraries'     : ['gsl','gslcblas']}
    
            
    code = """
            int i,j,cn,l,ay,az,ix,iy,signum,m,n;
            double denom,wmean,sx,sy,arg,kx,ky;
            gsl_complex z;
            gsl_matrix_complex *mp,*minv;
            gsl_permutation *p;
            
            for(ix=0;ix<nk;ix++){
               sx = smin + ix*sinc;
               kx=-6.2832*(smin+(double)ix*sinc)*freq;
               for(iy=0;iy<nk;iy++){
                  sy = smin + iy*sinc;
                  ky=-6.2832*(smin+(double)iy*sinc)*freq;
                  for(i=0;i<nr;i++){
	             for(j=0;j<nr;j++){
	                iismr(i,j)=0.0;
	                iismi(i,j)=0.0;
	             }
	          }
		  for(l=find-fave;l<=find+fave;l++){
		     tf = double(l)/double(nsamp*dt);
		     for(az=0;az<nr;az++){
		        tr(az) = cos(2*3.14159265*(freq-tf)*(sx*(rx(0)-rx(az))+sy*(ry(0)-ry(az))));
		        ti(az) = -sin(2*3.14159265*(freq-tf)*(sx*(rx(0)-rx(az))+sy*(ry(0)-ry(az))));
		     }
     	             for(i=0;i<nr;i++){
	                for(j=0;j<nr;j++){
                           iismr(i,j) += -ti(i)*xsmi(l-find+fave,i,j)*tr(j) + ti(i)*xsmr(l-find+fave,i,j)*ti(j) + tr(i)*xsmi(l-find+fave,i,j)*ti(j) + tr(i)*xsmr(l-find+fave,i,j)*tr(j);
                           iismi(i,j) += +ti(i)*xsmi(l-find+fave,i,j)*ti(j) + ti(i)*xsmr(l-find+fave,i,j)*tr(j) + tr(i)*xsmi(l-find+fave,i,j)*tr(j) - tr(i)*xsmr(l-find+fave,i,j)*ti(j);
                        }
                     }
                  }
                  denom=(double)(nwin*(2*fave+1));
	          for(i=0;i<nr;i++){
	             for(j=0;j<nr;j++){
	                iismr(i,j)/=denom;
	                iismi(i,j)/=denom;
	             }
	          }
	          wmean=0.0;
                  for(i=0;i<nr;i++){
                      w(i)=pow(iismr(i,i)*iismr(i,i)+iismi(i,i)*iismi(i,i),-0.25);
          	    wmean+=1.0/(w(i)*w(i));
                  }
          	wmean/=(double)nr;
                  for(i=0;i<nr;i++){
                      for(j=0;j<nr;j++){
                          iismr(i,j)*=(w(i)*w(j));
                          iismi(i,j)*=(w(i)*w(j));
                      }   
                  }
                  
                  

       
                  mp=gsl_matrix_complex_alloc(nr,nr);
                  minv=gsl_matrix_complex_alloc(nr,nr);
                  p=gsl_permutation_alloc(nr);

                  for(i=0;i<nr;i++){
                      for(j=0;j<nr;j++){
                          z=gsl_complex_rect(iismr(i,j),iismi(i,j));
                          gsl_matrix_complex_set(mp,i,j,z);
                      }
                  }
                  
                  gsl_linalg_complex_LU_decomp(mp,p,&signum);
                  gsl_complex_abs(gsl_linalg_complex_LU_det(mp,signum));
                  gsl_linalg_complex_LU_invert(mp,p,minv);
        
                  for(i=0;i<nr;i++){
                      for(j=0;j<nr;j++){
                          z=gsl_matrix_complex_get(minv,i,j);
                          iismr(i,j)=GSL_REAL(z);
                          iismi(i,j)=GSL_IMAG(z);
                      }
                  }

                  gsl_matrix_complex_free(mp);
                  gsl_matrix_complex_free(minv);
                  gsl_permutation_free(p);
                  
                  fk(ix,iy)=0.0;
                  for(m=0;m<nr;m++){
		        fk(ix,iy)+=iismr(m,m);
		        for(n=m+1;n<nr;n++){
		            arg=kx*(rx(m)-rx(n))+ky*(ry(m)-ry(n));
			    fk(ix,iy)+=2.0*(iismr(m,n)*cos(arg)-iismi(m,n)*sin(arg));
		        }
		    }
 	        fk(ix,iy)=1.0/fk(ix,iy);
                  
                  
               }
            }      
            """
    sp.weave.inline(code,vari,type_converters=converters.blitz,compiler = 'gcc',**weave_info)            
   
    max = 0.0
    for i in range(nk):
        for j in range(nk):
            if(fk[i][j].real>max):
                max=fk[i][j].real
                sxopt=smin+i*sinc
                syopt=smin+j*sinc
            if(fk[i][j].real<0):
                fk[i][j]=0
    rp = sp.sqrt(sxopt**2+syopt**2)
    vel = 111.19/rp
    baz=np.math.atan2(sxopt,syopt)*180.0/3.1415926
    if(baz<0.0):
       baz+=360.0

    #return fk.real/max,sxopt,syopt,vel,rp,baz 
    res = 10*sp.log10(fk.real/max)
    #res = 10*sp.log10(fk.real)
    return res,sxopt,syopt,vel,rp,baz  


def metric(st):
    import util as ut
    import numpy as np
    import scipy as sp
    '''
    function takes data matrix and returns interstation distances rx and ry (vectors) in [deg].
    '''
    nr = len(st)
    rx = np.zeros(nr)
    ry = np.zeros(nr)
    for i in range(nr):
        decl,dist,az,baz = ut.grt(st[0].stats.sac.stla,st[0].stats.sac.stlo,st[i].stats.sac.stla,st[i].stats.sac.stlo)
        rx[i] = decl*sp.cos(0.017453*(90.0-az))
        ry[i] = decl*sp.sin(0.017453*(90.0-az))
    return rx,ry
    


def get_metadata(meta_f):
    d = dict()
    with open(meta_f) as f:
        for line in f:
            x = line.split('|')
            d[x[1]] = x[4],x[5]
    return d

def metric_mseed(st,d,nr):
    import util as ut
    import numpy as np
    import scipy as sp
    '''
    function takes data matrix and returns interstation distances rx and ry (vectors) in [deg].
    '''
    rx_0,ry_0 = d[st[0].stats.station]
    rx = np.zeros(nr)
    ry = np.zeros(nr)
    for i in range(nr):
        rx_i,ry_i = d[st[i].stats.station]
        decl,dist,az,baz = ut.grt(float(rx_0),float(ry_0),float(rx_i),float(ry_i))
        rx[i] = decl*sp.cos(0.017453*(90.0-az))
        ry[i] = decl*sp.sin(0.017453*(90.0-az))
    return rx,ry
    

    
def testfir(st,cb,ct,n):
    from scipy import signal
    import numpy as np
    
    sr = st[0].stats.sampling_rate/2.
    xx = np.empty([st.count(),n],)
    a = signal.firwin(n, cutoff = cb/sr, window = 'hamming')
    b = - signal.firwin(n, cutoff = ct/sr, window = 'hamming'); b[n/2] = b[n/2] + 1
    d = - (a+b); d[n/2] = d[n/2] + 1
    fft1 = np.abs(np.fft.fft(d))
    for i in range(st.count()):
        fft = np.fft.fft(st[i][:n])*fft1
        xx[i] = np.fft.ifft(fft)
    return xx
    
def print_stats(fk,threshold):
    import numpy as np
    import scipy as sp
    import scipy.ndimage.filters as filters

    tmp = []
    print '---------------------------'
    print '--- Arrival Information ---'
    print '---------------------------'
    print 
    print 'normalized power (dB)   ', 'velocity (km/s)   ', 'backazimuth (deg)'
    maxxi = (np.where(fk==filters.maximum_filter(fk, 5)))
    this=np.empty([2,len(maxxi[0])])
    lth = np.amin(fk)*threshold
    for i in range(len(maxxi[0])):
        this[0][i]=(maxxi[0][i]-100)*0.5
        this[1][i]=(maxxi[1][i]-100)*0.5
        if (fk[maxxi[0][i],maxxi[1][i]] > lth):
            baz=np.math.atan2(this[0][i],this[1][i])*180.0/3.1415926
            if(baz<0.0):
                baz+=360.0
            xvel = 111.19/sp.sqrt(this[0][i]**2+this[1][i]**2)
            xamp = fk[maxxi[0][i],maxxi[1][i]]
            tmp.append([xamp,xvel,baz])
    tmp.sort(reverse=True)
    for i in range(len(tmp)):
        print '%12.02f %19.02f %19.02f'%(tmp[i][0],tmp[i][1],tmp[i][2])


    
    
