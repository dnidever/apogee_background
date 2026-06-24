import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dlnpyutils import utils as dln,plotting as pl,robust
from scipy.stats import binned_statistic
from astropy.table import Table,vstack


def plots():


    lcodr17 = Table.read('lco25m_medsky_cframe_results_dr17.062226.fits')
    apodr17 = Table.read('apo25m_medsky_cframe_results_dr17.062226.fits')

    apodr17['medsky'] *= 1.9  # convert to electrons
    lcodr17['medsky'] *= 3.0


    from astropy.coordinates import SkyCoord
    acoo=SkyCoord(apodr17['ra'],apodr17['dec'],unit='degree',frame='icrs')
    apodr17['glon']=acoo.galactic.l.degree
    apodr17['glat']=acoo.galactic.b.degree
    lcoo=SkyCoord(lcodr17['ra'],lcodr17['dec'],unit='degree',frame='icrs')
    lcodr17['glon']=lcoo.galactic.l.degree
    lcodr17['glat']=lcoo.galactic.b.degree

    lcodr17['medsky47'] = lcodr17['medsky']/lcodr17['exptime']*500  # to 47 reads
    apodr17['medsky47'] = apodr17['medsky']/apodr17['exptime']*500


    
    # --- APO DR17 ---

    hi, = np.where(np.abs(apodr17['glat'])>20)
    mid, = np.where(np.abs(apodr17['glat'])<5)

    
    # Sky vs. Moon fraction (high latitude)
    o=pl.hist2d(apodr17['moonfrac'][hi],apodr17['medsky47'][hi],log=True,xr=[0,1],yr=[-2,600],xtitle='Moon Fraction',nx=100,ny=100,
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO DR17 Sky Fiber Continuum Level (|b|>20)')
    res,xedge,binnumber = binned_statistic(apodr17['moonfrac'][hi],apodr17['medsky47'][hi],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(apodr17['medsky47'][hi])))
    plt.legend(loc='upper left')
    plt.savefig('apo25m_dr17_medsky47_electrons_moonfrac_highglat.png',bbox_inches='tight')

    # Sky vs. Moon fraction (midplane)
    o=pl.hist2d(apodr17['moonfrac'][mid],apodr17['medsky47'][mid],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO DR17 Sky Fiber Continuum Level (|b|<5)')
    res,xedge,binnumber = binned_statistic(apodr17['moonfrac'][mid],apodr17['medsky47'][mid],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(apodr17['medsky47'][mid])))
    plt.legend(loc='upper left')
    plt.savefig('apo25m_dr17_medsky47_electrons_moonfrac_midplane.png',bbox_inches='tight')

    # Sky vs. galactic latitude
    o=pl.hist2d(apodr17['glat'],apodr17['medsky47'],log=True,xr=[-90,90],yr=[-2,600],nx=100,ny=100,xtitle='Galactic Latitude',
                ytitle='APO DR17 Sky Fiber Continuum (electrons in 500s)',title='APO DR17 Sky Fiber Continuum Level')
    plt.axhline(40,c='red',linewidth=1)
    plt.axhline(100,c='red',linewidth=1)
    plt.axhline(200,c='red',linewidth=1)
    plt.axhline(300,c='red',linewidth=1)
    plt.savefig('apo25m_dr17_medsky47_electrons_glat.png',bbox_inches='tight')

    # Sky vs. airmass
    o=pl.hist2d(apodr17['airmass'],apodr17['medsky47'],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO DR17 Sky Fiber Continuum level')
    plt.savefig('apo25m_dr17_medsky47_electrons_airmass.png',bbox_inches='tight')

    # Sky vs. airmass (high latitude)
    o=pl.hist2d(apodr17['airmass'][hi],apodr17['medsky47'][hi],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO DR17 Sky Fiber Continuum level (|b|>20)')
    plt.savefig('apo25m_dr17_medsky47_electrons_airmass_highglat.png',bbox_inches='tight')

    

    # DR17 doesn't show the same increase in sky background that we see in the SDSS-V data.  WHY???
    # is the airmass distributiond different or moonphase?

    # I should really look at the sky fibers, not the sky model
    # maybe even in the ap1D files

    # it's possible that the DR17 plate sky values are just higher everywhere (bumped up) and the FPS
    # high-latitude values are lower

    # the Magellanic Clouds!!
    # we also don't have as many high-lat fields in DR17, maybe it's just harder to make out the details


    
    # --- LCO DR17 ---

    hi, = np.where(np.abs(lcodr17['glat'])>20)
    mid, = np.where(np.abs(lcodr17['glat'])<5)

    # Sky vs. Moon fraction (high latitude)
    o=pl.hist2d(lcodr17['moonfrac'][hi],lcodr17['medsky47'][hi],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO DR17 Sky Fiber Continuum Level (|b|>20)')
    res,xedge,binnumber = binned_statistic(lcodr17['moonfrac'][hi],lcodr17['medsky47'][hi],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(lcodr17['medsky47'][hi])))
    plt.legend(loc='upper left')
    plt.savefig('lco25m_dr17_medsky47_electrons_moonfrac_highglat.png',bbox_inches='tight')

    # Sky vs. Moon fraction (midplane)
    o=pl.hist2d(lcodr17['moonfrac'][mid],lcodr17['medsky47'][mid],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO DR17 Sky Fiber Continuum Level (|b|<5)')
    res,xedge,binnumber = binned_statistic(lcodr17['moonfrac'][mid],lcodr17['medsky47'][mid],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(lcodr17['medsky47'][mid])))
    plt.legend(loc='upper left')
    plt.savefig('lco25m_dr17_medsky47_electrons_moonfrac_midplane.png',bbox_inches='tight')

    # Sky vs. galactic latitude
    o=pl.hist2d(lcodr17['glat'],lcodr17['medsky47'],log=True,xr=[-90,90],yr=[-2,600],nx=100,ny=100,xtitle='Galactic Latitude',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO DR17 Sky Fiber Continuum Level')
    plt.axhline(40,c='red',linewidth=1)
    plt.axhline(100,c='red',linewidth=1)
    plt.axhline(200,c='red',linewidth=1)
    plt.axhline(300,c='red',linewidth=1)
    plt.savefig('lco25m_dr17_medsky47_electrons_glat.png',bbox_inches='tight')

    # Sky vs. airmass
    o=pl.hist2d(lcodr17['airmass'],lcodr17['medsky47'],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO DR17 Sky Fiber Continuum level')
    plt.savefig('lco25m_dr17_medsky47_electrons_airmass.png',bbox_inches='tight')

    # Sky vs. airmass (high latitude)
    o=pl.hist2d(lcodr17['airmass'][hi],lcodr17['medsky47'][hi],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO DR17 Sky Fiber Continuum level (|b|>20)')
    plt.savefig('lco25m_dr17_medsky47_electrons_airmass_highglat.png',bbox_inches='tight')
    
    

    
    # ======================================================================

    # SDSS-V


    lco = Table.read('lco25m_medsky_cframe_results_daily.062226.fits')
    apo = Table.read('apo25m_medsky_cframe_results_daily.062226.fits')

    apo['medsky'] *= 1.9  # convert to electrons
    lco['medsky'] *= 3.0


    from astropy.coordinates import SkyCoord
    acoo=SkyCoord(apo['ra'],apo['dec'],unit='degree',frame='icrs')
    apo['glon']=acoo.galactic.l.degree
    apo['glat']=acoo.galactic.b.degree
    lcoo=SkyCoord(lco['ra'],lco['dec'],unit='degree',frame='icrs')
    lco['glon']=lcoo.galactic.l.degree
    lco['glat']=lcoo.galactic.b.degree

    lco['medsky47'] = lco['medsky']/lco['exptime']*500  # to 47 reads
    apo['medsky47'] = apo['medsky']/apo['exptime']*500


    
    # --- APO SDSS-V ---

    hi, = np.where(np.abs(apo['glat'])>20)
    mid, = np.where(np.abs(apo['glat'])<5)

    # Sky vs. Moon fraction (high latitude)
    o=pl.hist2d(apo['moonfrac'][hi],apo['medsky47'][hi],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO SDSS-V Sky Fiber Continuum Level (|b|>20)')
    res,xedge,binnumber = binned_statistic(apo['moonfrac'][hi],apo['medsky47'][hi],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(apo['medsky47'][hi])))
    plt.legend(loc='upper left')
    plt.savefig('apo25m_daily_medsky47_electrons_moonfrac_highglat.png',bbox_inches='tight')

    # Sky vs. Moon fraction (midplane)
    o=pl.hist2d(apo['moonfrac'][mid],apo['medsky47'][mid],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO SDSS-V Sky Fiber Continuum Level (|b|<5)')
    res,xedge,binnumber = binned_statistic(apo['moonfrac'][mid],apo['medsky47'][mid],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(apo['medsky47'][mid])))
    plt.legend(loc='upper left')
    plt.savefig('apo25m_daily_medsky47_electrons_moonfrac_midplane.png',bbox_inches='tight')

    # Sky vs. Moon fraction (midplane)
    o=pl.hist2d(apo['glat'],apo['medsky47'],log=True,xr=[-90,90],yr=[-2,600],nx=100,ny=100,xtitle='Galactic Latitude',
                ytitle='APO DR17 Sky Fiber Continuum (electrons in 500s)',title='APO SDSS-V Sky Fiber Continuum Level')
    plt.axhline(40,c='red',linewidth=1)
    plt.axhline(100,c='red',linewidth=1)
    plt.axhline(200,c='red',linewidth=1)
    plt.axhline(300,c='red',linewidth=1)
    plt.savefig('apo25m_daily_medsky47_electrons_glat.png',bbox_inches='tight')

    # Sky vs. airmass
    o=pl.hist2d(apo['airmass'],apo['medsky47'],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO SDSS-V Sky Fiber Continuum level')
    plt.savefig('apo25m_daily_medsky47_electrons_airmass.png',bbox_inches='tight')

    # Sky vs. airmass (high latitude)
    o=pl.hist2d(apo['airmass'][hi],apo['medsky47'][hi],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='APO SDSS-V Sky Fiber Continuum level (|b|>20)')
    plt.savefig('apo25m_daily_medsky47_electrons_airmass_highglat.png',bbox_inches='tight')

    

    
    
    # --- LCO SDSS-V ---

    hi, = np.where(np.abs(lco['glat'])>20)
    mid, = np.where(np.abs(lco['glat'])<5)


    # Sky vs. Moon fraction (high latitude)
    o=pl.hist2d(lco['moonfrac'][hi],lco['medsky47'][hi],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500)',title='LCO SDSS-V Sky Fiber Continuum Level (|b|>20)')
    res,xedge,binnumber = binned_statistic(lco['moonfrac'][hi],lco['medsky47'][hi],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(lco['medsky47'][hi])))
    plt.legend(loc='upper left')
    plt.savefig('lco25m_daily_medsky47_electrons_moonfrac_highglat.png',bbox_inches='tight')

    # Sky vs. Moon fraction (midplane)
    o=pl.hist2d(lco['moonfrac'][mid],lco['medsky47'][mid],log=True,xr=[0,1],yr=[-2,600],nx=100,ny=100,xtitle='Moon Fraction',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO SDSS-V Sky Fiber Continuum Level (|b|<5)')
    res,xedge,binnumber = binned_statistic(lco['moonfrac'][mid],lco['medsky47'][mid],bins=np.linspace(0,1.1,12),statistic='median')
    xres = xedge[:-1]+0.5*(xedge[1]-xedge[0])
    plt.scatter(xres,res,c='r')
    gg,=np.where(np.isfinite(res))
    coef = np.polyfit(xres[gg],res[gg],2)
    xx = np.linspace(0,1.0,100)
    plt.plot(xx,np.polyval(coef,xx),c='r',label='Median={:.1f} electrons'.format(np.nanmedian(lco['medsky47'][mid])))
    plt.legend(loc='upper left')
    plt.savefig('lco25m_daily_medsky47_electrons_moonfrac_midplane.png',bbox_inches='tight')


    # Sky vs. galactic latitude
    o=pl.hist2d(lco['glat'],lco['medsky47'],log=True,xr=[-90,90],yr=[-2,600],nx=100,ny=100,xtitle='Galactic Latitude',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO SDSS-V Sky Fiber Continuum Level')
    plt.axhline(40,c='red',linewidth=1)
    plt.axhline(100,c='red',linewidth=1)
    plt.axhline(200,c='red',linewidth=1)
    plt.axhline(300,c='red',linewidth=1)
    plt.savefig('lco25m_daily_medsky47_electrons_glat.png',bbox_inches='tight')


    # Sky vs. airmass
    o=pl.hist2d(lco['airmass'],lco['medsky47'],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO SDSS-V Sky Fiber Continuum level')
    plt.savefig('lco25m_daily_medsky47_electrons_airmass.png',bbox_inches='tight')

    # Sky vs. airmass (high latitude)
    o=pl.hist2d(lco['airmass'][hi],lco['medsky47'][hi],log=True,nx=100,ny=100,xr=[1.0,2.25],yr=[-2,600],xtitle='Airmass',
                ytitle='Sky Fiber Continuum (electrons in 500s)',title='LCO SDSS-V Sky Fiber Continuum level (|b|>20)')
    plt.savefig('lco25m_daily_medsky47_electrons_airmass_highglat.png',bbox_inches='tight')
    
    # DR17 h igh latitude median
    # APO - 100.7
    # LCO - 50.8

    # SDSS-V  high latitude median
    # APO - 203.9
    # LCO - 49.1

    # We expect the APO level to be 2.4x the LCO level because of the different fiber diameters
    # 50*2.4=120, that's about right  (1.31" vs. 2.0", see Wilson+2019, pg.61)

    # APO is roughly a factor of 2x different, which is also the gain level
    # is it possible the gain is the difference
    
    # at moonfrac=0, the different is smaller, only about 40%
    

    #In [200]: len(apo)
    #Out[200]: 38057

    #In [201]: len(apodr17)
    #Out[201]: 60202

    #In [202]: len(lco)
    #Out[202]: 51155

    #In [203]: len(lcodr17)
    #Out[203]: 15622

    #o=pl.scatter(apodr17['mjd'].astype(float),apodr17['medsky47'],size=3)
    #o=pl.hist2d(apodr17['mjd'].astype(float),apodr17['medsky47'],yr=[0,600],log=True)

    # with the new scraping, the median APO SDSS-V sky level is 155 and about 129 for moonfrac<0.3

    # Make HTML page with table for all of the figures

    # combine the tables
    apo['telescope'] = 'apo'
    apo['survey'] = 'sdss5'
    apodr17['telescope'] = 'apo'
    apodr17['survey'] = 'dr17'
    apo = vstack((apodr17,apo))
    
    lco['telescope'] = 'lco'
    lco['survey'] = 'sdss5'
    lcodr17['telescope'] = 'lco'
    lcodr17['surey'] = 'dr17'
    lco = vstack((lcodr17,lco))

    # add per arcsec^2 of fiber
    apo['medsky47perarcsec2'] = apo['medsky47'] / 3.14
    lco['medsky47perarcsec2'] = lco['medsky47'] / 1.347
    
    #apo.write('apo25m_medsky_cframe_results_combined.fits',overwrite=True)
    #lco.write('lco25m_medsky_cframe_results_combined.fits',overwrite=True)

    
    import pdb; pdb.set_trace()
