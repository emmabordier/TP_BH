import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
f=fits.open('NVSSJ201943-364542.fits')  #on importe le fichier fits
f.info()
hdu=f[0]
data=hdu.data
header=hdu.header

#print (data)


data[260] #ou le 1230

plt.ion()
plt.figure(1)
xy=data[260]
plt.imshow(xy,interpolation='none', origin='lower',cmap='viridis',vmin=-0.5e-18) #toutes les valeurs plus basses que vmin sont affichées de la couleur la plus foncée
cb=plt.colorbar()
cb.set_label('Intensité')


#on construit le tableau des wcs
j=np.arange(0,2047,1)
wvl=header['CRVAL3']+header['CDELT3']*(j-header['CRPIX3'])


#on construit le spectre moyen de la galaxie
galaxy=data[:,44:53,44:53]   #on sélectionne les x et y où on a notre AGN
galaxy_x=np.mean(galaxy,axis=2) # on fait une moyenne sur chaque axe
spectrum=np.mean(galaxy_x,axis=1)


masque=np.where(~((wvl>17835) & (wvl<19950)))   #interbande des raies H de l'atmosphère
masque=masque[0]
spectre=spectrum[masque]
wvl1=wvl[masque]


plt.ion()
plt.figure(2)
x=wvl1
y=spectre
plt.plot(x,y,'r')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity')
plt.title('Spectre de la galaxie avec correction des raies H')

#on trace OHLINES pour savoir quelles valeurs enlever
oh_lines=np.loadtxt('OHlines.dat')
oh_wvl=oh_lines[:,0]
oh_i=oh_lines[:,1]

plt.ion()
plt.figure(3)
plt.plot(oh_wvl,oh_i,'g.-')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity')
plt.title('Intensités OH')

#print (oh_i[3])


bool_mask=(~((wvl>17835) & (wvl<19950)))            #on réalise le masque OHlines et on combine avec les raies H
for oh_wvl, oh_i in oh_lines:
    if (oh_i>70):
        line_bool_mask=(wvl>(oh_wvl+2.5))|(wvl<(oh_wvl-2.5))
        bool_mask=bool_mask & line_bool_mask

masque2=np.where(bool_mask)
masque2=masque2[0]
spectre2=spectrum[masque2]
wvl2=wvl[masque2]


plt.figure(4)                                         # spectre final de la galaxie (avec corrections)
plt.plot(wvl2,spectre2,'b')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity')
plt.title('Spectre de la galaxie avec les corrections des raies OH et H')



#fit gaussienne
bool_mask2=(wvl2<20439.6)|(wvl2>20486.0)     # on enlève le pic en Halpha
masque3=np.where(bool_mask2)
masque3=masque3[0]
spectre3=spectre2[masque3]
wvl3=wvl2[masque3]

plt.figure(5)
plt.plot(wvl3,spectre3,'g')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity')
plt.title('Spectre de la galaxie avec les corrections de l''atmosphère et le pic de Halpha en moins')




bool_mask3=(wvl3<15900.5)&(wvl3>15334.5)   #création tableau x et y pour raies O3
masque4=np.where(bool_mask3)
masque4=masque4[0]
spectreO3=spectre3[masque4]
wvlO3=wvl3[masque4]


bool_mask4=(wvl3<21300.9)&(wvl3>19810.7)   #création tableau x et y pour raies Halpha
masque5=np.where(bool_mask4)
masque5=masque5[0]
spectreHalpha=spectre3[masque5]
wvlHalpha=wvl3[masque5]


def gaussienne(lambda1,deltalambda,p,b,lambda0):
    return p*np.exp(-(lambda1-lambda0)**2/(2*deltalambda**2))+b  #c:pic

from scipy.optimize import curve_fit

popt, pcov = curve_fit(gaussienne, wvlO3, spectreO3, p0 = [200,28.57e-20,6e-20,15618.5])


plt.figure(6)                                        
plt.plot(wvlO3,spectreO3)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity (cgs)')
plt.title('Fit Spectre raies O3 ')
plt.plot(wvlO3, gaussienne(wvlO3, *popt), label='fit')

popt2, pcov2 = curve_fit(gaussienne, wvlHalpha, spectreHalpha, p0 = [200,10.57e-20,6e-20,20462.8])

plt.figure(7)                                        
plt.plot(wvlHalpha,spectreHalpha)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Intensity (cgs)')
plt.title('Fit Spectre raies Halpha')
plt.plot(wvlHalpha, gaussienne(wvlHalpha, *popt2), label='fit')  
