IAS-Capon
=========

A collection of Direction of arrival (DOA) algorithms in Python for processing of passive seismic noise.
by Martin Gal, University of Tasmania
email: martin.gal@utas.edu.au
Date: 16.12.2014

----------------------------------------------------------------------------------------------------
Gal, M., Reading, a. M., Ellingsen, S. P., Koper, K. D., Gibbons, S. J., & N{"a}sholm, S. P. (2014). 
Improved implementation of the fk and Capon methods for array analysis of seismic noise. 
Geophysical Journal International, 198(2), 1045–1054. doi:10.1093/gji/ggu183
----------------------------------------------------------------------------------------------------

(I would like to thank Keith Koper who's conventional Capon code was very helpful and 
is partially present in the IAS implementation.)

Changelog
=========
30th April 2015
- Fixed bug with obspy future dependency

What you will need to run the code:
=========
Python 2.7.x
Numpy
Scipy
Matplotlib
Obspy
C-compiler 
(The code is mainly written in python but bottleneck calculations are scripted in 
scipy.weave which allows the inline implementation of C code in python. 
This boosts the performance considerably when processing multiple hours of data.)

For the use of CAS_Capon (which is an experimental code included in the script) you 
will need the GSL library as it makes heavy use of C to decrease the computational time.
This also requieres to edit library paths in the subroutines.py . In case you need help 
feel free to contact me.


Installation:
========
Compilation is not needed, only the above packages/modules are needed.
If you are new to python a simple solution is to download the anaconda package 
which comes with a lot of modules. Another option is to use pip, which allows 
easy installation of python modules.




Documentation:
========
General:
The script is divided into 2 parts. Git_Capon_mseed/sac.py call specific function 
to obtain the direction of arrival estimates and subroutine.py/util.py supply the 
algorithms/functions.
All control parameters i.e. user input is to be specified in Git_Capon_mseed/sac.py.

User input parameter in Git_Capon_mseed/sac.py:
example
# ==== USER INPUT PARAMETER ===
nsamp = 4000
smin = -50.0
smax = 50.0
sinc = 0.5
cap_find = 120
cap_fave = 15
dl = 1


st = read('location of mseed file',format='MSEED')
dic_meta = get_metadata('path to metadata from IRIS server (or check examples)')
# =============================


nsamp ... 	amount of snapshots of the data sequence(in a temporal subwindow!). 
			The script will take this number and extract as many temporal subwindows from a
			supplied time series as possible, i.e. if the time series has 40000 data points,
			the algorithm will divide the sequence in 10 temporal subwindows if nsamp=4000.

smin ... 	minimum boundary value of the projected slowness grid

smax ... 	minimum boundary value of the projected slowness grid

sinc ... 	grid spacing, i.e. step size on slowness grid

cap_find ... frequency at which the DOA is performed. warning: it is not [Hz]. 
			It is connected to the frequency by the equation freq = cap_find/(nsamp*delta) 
			where freq is in [Hz] and delta is the time difference between to data points.
			It is essentially the frequency at which the projection on the slowness grid is 
			carried out.

cap_fave ... if you desire a broadband representation of the spectrum cap_fave 
			is 1/2 the width of this frequency range. DOA estimation is performed at 
			cap_find +- cap_fave.

dl ... 		diagonal loading parameter. It is only implemented for the conventional 
			and IAS capon algorithm. Usually a value of 1 is enough to protect against bias 
			that arises from singular or close to singular cross spectral density matrices 
			(spectral covariance matrix). I you feel something is fishy with the resulting 
			spectrum turn it up to 10 or 100 and see if it make a difference, otherwise have 
			a look at the Gal et al. 2014 for more info.

read ... 	enter path to you mseed file (see examples)

get_metadata ... path to metadata file from iris (see examples)


For people that are interested in modifying the code to suite their needs, they will find 
some comments in the subroutines.py file. I have commented only the Capon function which 
should give the user enough information to modify any function as they are quite similar. 
Although the techniques are quite similar in their approach, I have made a separate function 
for each of them to make modifications more simple. 



Use:
=========
I would recommend to start with the examples and reproduce the included results. Try to
shuffle around the user supplied parameters to get a feeling on what the script is doing.

At the moment, you have the choice of 5 DOA methods:
conventional fk and Capon,
IAS fk and IAS Capon,
CAS_Capon which is an experimental code that performs coherent stacking of frequencies 
(which can be useful for short time series, low frequencies, to mitigate stability issues 
and others).

To select a certain algorithm, simply comment and uncomment the appropriate lines, i.e.
#fk,sxopt,syopt,vel,rp,baz = FK(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz = IAS_FK(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz,maa,pwe = Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,dl,overlap=True,taper=True)
fk,sxopt,syopt,vel,rp,baz,maa,pwe = IAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,dl,overlap=True,taper=True)
#fk,sxopt,syopt,vel,rp,baz = CAS_Capon(nsamp,nr,rx,ry,st,smin,smax,sinc,cap_find,cap_fave,dt,overlap=True,taper=True)
(here IAS Capon is selected.)

Further, you have the choice to use a 50% overlap for temporal subwindows (overlap=True).
The Hann-window function can be applied by taper=True to reduce spectral leakage.





Useful experiences with the code:
=======
Which algorithm to chose:
The conventional algorithms are the fastest. Fk is used by many people and offers a robust
estimate of the DOA analysis. The conventional Capon algorithm will give you a higher 
resolution, but can be prone to bias for very narrowband calculations. If these algorithms 
are used for a broad frequency range (i.e. >10% of the projection frequency) the resulting 
spectrum can show bias in slowness (more on this in Gal et al. 2014).
The IAS versions will give you a much cleaner slowness spectrum. If you have the 
computational time, you will be better of using it than the conventional option. 
As IAS generate small frequency bins and sums over the resulting spectra it will induce a 
bias into the Capon method. The bias can be reduced by applying diagonal loading. Here again 
IAS Capon shows a better resolution than IAS fk, but if you do not feel sure that your 
results are biased, compare your results with IAS fk. It is always wise to have more temporal 
subwindows than stations to ensure the cross power spectral density is not singular 
(this is also why windowing is beneficial as it increases the amount of temporal subwindows).
CAS Capon is very slow, but shows clear advantage for low frequencies, and low SNR. 
It will do a great job unless your sources are very closely spaced, then IAS will be the better choice.

Window length:
If too small, your frequency spectrum will be quite inaccurate and frequency leakage will 
influence your DOA. This means that sources from outside the frequency range of interest 
will enter your result. If the window is too long you will have less temporal windows to 
average over which will decrease the performance of fk and Capon.

Frequency width of analysis:
Depends on what you want to look at, the bigger it is the greater the computational cost.

Diagonal loading parameter:
It will be dependent on your array and other factors and the influence of it should 
be always studied. I had both cases where it was not necessary or a value of 1 was too 
small (Gal et al. 2014).

ARF:
Always be aware of the arrays ARF in the frequency band of interest! You can simply use 
Obspy to get the ARF.



Contact:
==========
For any questions or bug reports, feel free to contact me on the above email.
