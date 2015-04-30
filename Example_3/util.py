def grt(r1,r2,s1,s2):
    import numpy as np
    import math as m
    slat=s1*np.pi/180.
    slon=s2*np.pi/180.
    elat=r1*np.pi/180.
    elon=r2*np.pi/180.
    
    slat=m.atan(.996647*m.tan(slat))
    elat=m.atan(.996647*m.tan(elat))


    slat=np.pi/2.0-slat
    elat=np.pi/2.0-elat


    if(slon<0.0):
        slon+=2.0*np.pi
    if(elon<0.0):
        elon+=2.0*np.pi



    a=m.sin(elat)*m.cos(elon)
    b=m.sin(elat)*m.sin(elon)
    c=m.cos(elat)
    a1=m.sin(slat)*m.cos(slon)
    b1=m.sin(slat)*m.sin(slon)
    c1=m.cos(slat)

    cd=a*a1+b*b1+c*c1

    if(cd>1.0):
        cd=1.0
    if(cd<-1.0):
        cd=-1.0
    decl=m.acos(cd)*180.0/m.pi
    dist=decl*np.pi*6371.0/180.0

    tmp1=m.cos(elon)*m.cos(slon)+m.sin(elon)*m.sin(slon)
    tmp2a=1.0-cd*cd

    if tmp2a<=0.0:
	    tmp2=0.0
	    tmp3=1.0
    else:
        tmp2=m.sqrt(tmp2a)
        tmp3=(m.sin(elat)*m.cos(slat)-m.cos(elat)*m.sin(slat)*tmp1)/tmp2

    if(tmp3>1.0):
        tmp3=1.0
    if(tmp3<-1.0):
        tmp3=-1.0
    z=m.acos(tmp3)

    if((m.sin(slon)*m.cos(elon)-m.cos(slon)*m.sin(elon))<0.0):
        z=2.0*m.pi-z

    az=180.0*z/m.pi

    tmp1=m.cos(slon)*m.cos(elon)+m.sin(slon)*m.sin(elon)
    tmp2a=1.0-cd*cd
    if(tmp2a<=0.0):
        tmp2=0.0
        tmp3=1.0
    else: 
        tmp2=m.sqrt(tmp2a)
        tmp3=(m.sin(slat)*m.cos(elat)-m.cos(slat)*m.sin(elat)*tmp1)/tmp2

	

    if(tmp3>1.0):
        tmp3=1.0
    if(tmp3<-1.0):
        tmp3=-1.0
        
    bz=m.acos(tmp3)

    if((m.sin(elon)*m.cos(slon)-m.cos(elon)*m.sin(slon))<0.0):
        bz=2.0*m.pi-bz
        

    baz=180.0*bz/m.pi
    
    return decl,dist,az,baz
    

   
