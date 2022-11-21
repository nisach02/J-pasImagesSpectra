## setting the plot parameters
import matplotlib.pyplot as plt
from matplotlib import image
import pandas as pd
import numpy as np
from scipy.stats import sem 
from astropy.io import fits
from astropy.table import Table
#import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter

import glob

# set a grey background (use sns.set_theme() if seaborn version 0.11.0 or above) 
# sns.set(style="white")


#setting figure size
#plt.rcParams['figure.figsize'] = [8,6]
plt.rcParams['figure.dpi'] = 200
#setting x and y ticks to true
plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = 1
plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = 1
plt.rcParams["axes.linewidth"] = 1
plt.rcParams['font.family']='serif'
#plt.rcParams['text.usetex'] = True
#plt.rcParams['axes.labelsize'] = 'small'


jpasfilters = pd.read_csv('JpasNarrowFilters.csv')
jpaslambda = jpasfilters['lambda eff']
colors = jpasfilters['color_representation']

jpas_observed = pd.read_csv('/Users/nischal/Desktop/PhD/Neo/cigale_runs/21Nov2022/XrayIntegrated/xray_integratedV4.csv')

gal_spectra = fits.open('/Users/nischal/Desktop/PhD/Neo/cigale_runs/21Nov2022/XrayIntegrated/out/results.fits', memmap=True)
gal_spectra = Table(gal_spectra[1].data).group_by('best.universe.redshift')
jpasnames = gal_spectra['id']


# jpasnames = jpasnames[jpasnames == '2470-10043']   #to test for a single object in case of any problems.


#/Users/nischal/Desktop/PhD/Neo/cigale_runs/7Nov2022/outV7----Best-Fit
j = 0

for i in jpasnames:
    j +=1
    try:
        gal_spectra_s = gal_spectra[gal_spectra['id'] == i]
        jpas_observed_s = jpas_observed[jpas_observed['id'] == i]
        ##remember that the observed and fitted columns should be sliced according to the results.fits table.
        ##different runs have them in different tables. so be careful with it
        observed = jpas_observed_s[jpas_observed_s.columns[2::2][:-4]]
        fitted = gal_spectra_s.columns[360:].values()
        chisq = gal_spectra_s['best.reduced_chi_square']
        
        bestmodel = '/Users/nischal/Desktop/PhD/Neo/cigale_runs/21Nov2022/XrayIntegrated/out/' + str(i) + '_best_model.fits'
        
        gal_model = fits.open(bestmodel, memmap=True)
        gal_model = Table(gal_model[1].data)
        wavelength = gal_model['wavelength']
        wavelength_spec = wavelength * 10 #convert wavelength from nm to Å
    
        # assigning variable values... the summation is based on the models of CIGALE, and cigale-plots sed file
        agn = gal_model['agn.SKIRTOR2016_torus'] + gal_model['agn.SKIRTOR2016_polar_dust'] + gal_model['agn.SKIRTOR2016_disk']
        fnu = gal_model['Fnu']
        model = gal_model['L_lambda_total']
        stellar_unattenuated = gal_model['stellar.old'] + gal_model['stellar.young']
        stellar_attenuated = gal_model['stellar.old'] + gal_model['stellar.young'] + \
                gal_model['attenuation.stellar.old'] + gal_model['attenuation.stellar.young'] + \
            gal_model['nebular.absorption_old'] + gal_model['nebular.absorption_young']
    
        nebular = abs(gal_model['attenuation.nebular.lines_old'] + gal_model['attenuation.nebular.lines_young'] + \
            gal_model['attenuation.nebular.continuum_young'] + gal_model['attenuation.nebular.continuum_old'])
        dust = gal_model['dust']
    
        ### determining the scale factor to multiply the above variables by.
        ### because almost all of them are in the units of W/nm .... we need to convert it to mJy to plot
        c= 3* (10**8)
        surf = 4.0 * np.pi * (gal_spectra_s["best.universe.luminosity_distance"] ** 2)
        fact = 1e29 * 1e-9* (wavelength) ** 2 / c / surf
        minagn = np.median(agn*fact)
        minnebular = np.median(nebular*fact)
    
        jpasimage = '/Users/nischal/Desktop/PhD/ThumbnailsV2/XraysV4/' + str(i) + '.png'
        galaxy = image.imread(jpasimage)
        
    
        # ymin = gal_spectra_s['best.J0390'] - gal_spectra_s['best.J0390'] * 0.2
        # ymax = gal_spectra_s['best.J1007'] + gal_spectra_s['best.J1007'] * 0.2
        
        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(4,8), gridspec_kw={'height_ratios': [1,0.65],'width_ratios': [1]})
        ax1.imshow(galaxy)
        ax1.axis('off')
        
        ax2.plot(wavelength_spec, fact*stellar_attenuated, 
                                label = 'Stellar Emission', c='gold', linewidth = 0.75)
        ax2.plot(wavelength_spec, fact*stellar_unattenuated, 
                                label = 'Stellar Unattenuated', c='xkcd:deep sky blue', linestyle='--', linewidth = 0.75)
        ax2.plot(wavelength_spec, fact*nebular, 
                                label = 'Nebular Emission', c='green', linewidth = 0.75)        
        ax2.plot(wavelength_spec, fact*agn, 
                                label = 'AGN Emission', c='orange', linewidth = 0.75)
        ax2.plot(wavelength_spec, fact*model, 
                                label = 'Model Spectrum', c='black', linewidth = 0.75)
        ax2.scatter(jpaslambda[1:], fitted, label = 'Model Fluxes', facecolors='none', edgecolors='xkcd:strawberry')
        ax2.scatter(jpaslambda, observed, c=colors, label = 'Observed Fluxes')
        #ax2.title.set_text(str(chisq)[:3])




        ax2.set_xlim(3600,10000)

        ymin = np.min(observed.values) * 1e-1
        ymax = np.max(observed.values) * 1e1

        ax2.set_xscale('log')
        ax2.set_yscale('log')
        #ax2.set_yticks([0.01,1,2,3,4,5,6,7])
        #if np.max(agn*fact) > 0:    
        #    ax2.set_ylim(minagn,ymax)
        #else:
        #    ax2.set_ylim(ymin,ymax)
        ax2.set_xticks([4000,6000,8000,10000])
        ax2.set_ylim(ymin, ymax)
        
        ax2.set_xlabel(r'Observed $\lambda$ [$\AA$]', fontsize=6, labelpad=0)
        ax2.set_ylabel(r'S$_\nu$ [mJy]', fontsize=6, labelpad=0)
        #ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),1)))).format(y)))
        # ax2.ticklabel_format(useOffset=False, style='plain')
        ax2.xaxis.set_tick_params(labelsize=5)
        ax2.yaxis.set_tick_params(labelsize=5)
        #ax2.yaxis.set_major_formatter(pad=0)
        ax2.tick_params(axis='both', which='major', pad=1)
        
        plt.tick_params(axis='both', which='both', labelsize=5)
        plt.subplots_adjust(hspace=-0.6)
        plt.tight_layout()
        plt.legend(fontsize=5)
        
        """        
        if ymax < 0.1:
            plt.gca().yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
            ax2.yaxis.set_minor_formatter(FormatStrFormatter('%.2f'))
            plt.subplots_adjust(hspace=-0.2)

        elif ymax < 0.2:
            plt.gca().yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.03))
            ax2.yaxis.set_minor_formatter(FormatStrFormatter('%.2f'))
            
        elif ymax < 1:
            plt.gca().yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
            ax2.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
        
        elif ymax < 2:
            plt.gca().yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.4))
            ax2.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
        
        elif ymax < 4:
            plt.gca().yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.6))
            ax2.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
        else:
            pass
        
        """        
        
        savename = 'FinalImagesV1/' + str(j) + '.png'
        fig.savefig(savename,
                    bbox_inches='tight',
                    pad_inches=0,
                    format='png',
                    dpi=400)
        plt.show()
        
        # fig = plt.figure(num=2,figsize=(4,8))
        # gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[1,1])
        # ax0 = plt.subplot(gs[0])
        # ax1 = plt.subplot(gs[1])
    
        # ax0.plot(wavelength_spec, fact*model, 
        #                         label = 'Model Spectrum', c='black', linewidth = 0.75)
        # ax0.scatter(jpaslambda, observed, c=colors, label = 'Observed Points')
        # ax1.imshow(galaxy)
        # ax0.set_xlim(3700,10000)
        # ax0.set_ylim(ymin,ymax)
        # ax0.set_yscale('log')
    
        # plt.show()
    except:
        print ('Item not available with JPAS ID: ', str(i))
    

