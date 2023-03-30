import ROOT
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use([hep.style.ROOT, hep.style.firamath])
import numpy as np

_files_dict={
    '13_1_0_pre1':'orig.root',
    '13_1_0_pre1+PR':'pr.root',
}


def getSiPMonTileDigis(url,itSample=2):
    
    """reads all the SiPM-on-tile digis"""
    
    t=ROOT.RDataFrame('ana/hits',url)
    digi_coll=t.AsNumpy( ['event','detid','adc','isTOT','isSat','gain','isSci','layer'] )
    digis=pd.DataFrame(digi_coll)
    n=digis.shape[0]
    ceh=(digi_coll['layer']>26)
    sci=(digi_coll['isSci']==True)
    digis['digiid'] = np.where(ceh & sci, np.ones(n)*2,
                               np.where(ceh, np.ones(n), np.zeros(n)))
    digis=digis.rename(columns={'isTOT':'mode','adc':'ADC'})                              
    print(digis.head().T)
    return digis 


def processSimulations(overWrite=False):
    
    """calls the analyzer for all the DIGI files found in the base directory"""
   
    for name,f in _files_dict.items():
        try:
            outname=os.path.join(os.path.dirname(f),
                                 os.path.basename(f).replace('.root','.h5'))
            if not overWrite and os.path.isfile(outname):
                continue
            data=getSiPMonTileDigis(f)
            data.to_hdf(outname,key='digis')
            print(f,'->',outname,'with',data.shape)
        except Exception as e:
            print('Could not process',f)
            print(e)


def compareADCSpectra(data,key_list,outname,doSignal=True):

    """compare separately ADC / TOT spectra"""

    fig,ax=plt.subplots(1,2,figsize=(16,8))
    kwargs={'bins':np.linspace(0,1024,1025),
            'density':False,
            'linewidth':2,
            'histtype':'step'}

    for key in key_list:

        df=data[key]

        adc_mask=(df['mode']==False)
        tot_mask=~adc_mask
        mctruth_mask=np.ones_like(adc_mask.values,dtype=np.bool)
        if 'simE' in df.columns:
            mctruth_mask=(df['simE']>0)
            if not doSignal:
                mctruth_mask=~mctruth_mask
        adc_mask=adc_mask & mctruth_mask
        tot_mask=tot_mask & mctruth_mask

        ax[0].hist(df[adc_mask]['ADC']+0.5,label=getTitle(key), **kwargs)
        ax[1].hist(df[tot_mask]['ADC']+0.5,**kwargs) 
  
    for i in range(2):
        #ax[i].set_yscale('log')
        ax[i].text(0.05,0.95,
                   'ADC mode' if i==0 else 'TDC mode', 
                   horizontalalignment='left', 
                   verticalalignment='center', 
                   transform=ax[i].transAxes,
                   fontsize=14)
        ax[i].set_xlabel('ADC counts')
        ax[i].grid()
        if i==0: ax[i].legend()
    plt.tight_layout()
    plt.savefig('{}.png'.format(outname))

    #zoom
    for i in range(2):
        ax[i].set_xlim(0,50)
    plt.savefig('{}_zoom.png'.format(outname))
    plt.close() 


def compareDigisPerDetId(data,x,y,digiid,outname):

    """compare separately ADC / TOT spectra"""

    fig,ax=plt.subplots(1,2,figsize=(16,8))

    #merge and filter for true hits only
    df=data[x].merge(data[y],on=['event','detid','digiid'],how='inner',suffixes=('_x', '_y'))
    df=df[df['digiid']==digiid]
    print(df.head().T)

    common_mode=(df['mode_x']==df['mode_y'])
    common_gain=(df['gain_x']==df['gain_y'])
    adc_mode=(df['mode_x']==0)
    ax[0].scatter(df[common_mode & common_gain & adc_mode]['ADC_x'],
                  df[common_mode & common_gain & adc_mode]['ADC_y'],
                  label='Same mode/Same gain')
    ax[0].scatter(df[~common_mode & common_gain & adc_mode]['ADC_x'],
                  df[~common_mode & common_gain & adc_mode]['ADC_y'],
                  label='Mode switch/Same gain')
    ax[0].scatter(df[common_mode & ~common_gain & adc_mode]['ADC_x'],
                  df[common_mode & ~common_gain & adc_mode]['ADC_y'],
                  label='Same mode/Gain switch')
    ax[0].scatter(df[~common_mode & ~common_gain & adc_mode]['ADC_x'],
                  df[~common_mode & ~common_gain & adc_mode]['ADC_y'],
                  label='Mode switch/Gain switch')

    ax[1].scatter(df[common_mode & common_gain & ~adc_mode]['ADC_x'],
                  df[common_mode & common_gain & ~adc_mode]['ADC_y'],
                  label='Same mode/Same gain')
    ax[1].scatter(df[~common_mode & common_gain & ~adc_mode]['ADC_x'],
                  df[~common_mode & common_gain & ~adc_mode]['ADC_y'],
                  label='Mode switch/Same gain')
    ax[1].scatter(df[common_mode & ~common_gain & ~adc_mode]['ADC_x'],
                  df[common_mode & ~common_gain & ~adc_mode]['ADC_y'],
                  label='Same mode/Gain switch')
    ax[1].scatter(df[~common_mode & ~common_gain & ~adc_mode]['ADC_x'],
                  df[~common_mode & ~common_gain & ~adc_mode]['ADC_y'],
                  label='Mode switch/Gain switch')

    for i in range(2):    
        ax[i].plot([0,1024],[0,1024],'--',color='gray')
        ax[i].set_xlabel(x)
        ax[i].set_ylabel(y)
        ax[i].grid()
        ax[i].text(0.0,1.02,
                   'ADC mode' if i==0 else 'TDC mode', 
                   horizontalalignment='left', 
                   verticalalignment='center', 
                   transform=ax[i].transAxes,
                   fontsize=14)        
    ax[0].legend()
    plt.tight_layout()
    plt.savefig('{}_{}.png'.format(outname,digiid))


def makeValidationPlots():

    data={}
    for name,f in _files_dict.items():
        data[name]=pd.read_hdf(f.replace('.root','.h5'),key='digis')
    

    for digiid in [0,1,2]:
        #compareADCSpectra(data,data['digiid']==digiid,'adc')
        compareDigisPerDetId(data,'13_1_0_pre1','13_1_0_pre1+PR',digiid,'adcperdetid')


processSimulations(True)
makeValidationPlots()


