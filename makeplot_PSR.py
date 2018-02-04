#---------------------------------------
#TAKE ALL CATALOGS
#---------------------------------------
from astropy.io import fits
import random
from math import *
import matplotlib.pyplot as pl
import numpy as np
from scipy.interpolate import interp1d
import os.path


#3FGL
data_3fgl = fits.open("/nfs/farm/g/glast/u/mdimauro/GC_cube/gll_psc_v16.fit")
data3fgl =data_3fgl[1].data

class3fgl  = data3fgl.field('CLASS1')
glat3fgl  = data3fgl.field('GLAT')
glon3fgl  = data3fgl.field('GLON')
ra3fgl  = data3fgl.field('RAJ2000')
dec3fgl  = data3fgl.field('DEJ2000')
sign3fgl  = data3fgl.field('Signif_Avg')
name3fgl  = data3fgl.field('Source_Name')
TS10003fgl  = data3fgl.field('Sqrt_TS1000_3000')
TS30003fgl  = data3fgl.field('Sqrt_TS3000_10000')
posunc953fgl  = data3fgl.field('Conf_95_SemiMajor')
posunc683fgl  = data3fgl.field('Conf_68_SemiMajor')
enflux3fgl  = data3fgl.field('Energy_Flux100')
enfluxerr3fgl  = data3fgl.field('Unc_Energy_Flux100')
flux10003fgl  = data3fgl.field('Flux1000')
flux1000err3fgl  = data3fgl.field('Unc_Flux1000')
index3fgl  = data3fgl.field('Spectral_Index')
indexerr3fgl  = data3fgl.field('Unc_Spectral_Index')
spectype3fgl  = data3fgl.field('SpectrumType')

glon3fglpsr=[]
glat3fglpsr=[]
ra3fglpsr=[]
dec3fglpsr=[]

for t in range(len(glat3fgl)):
    if class3fgl[t]=='psr' or class3fgl[t]=='PSR' or class3fgl[t]=='spp':
       glon3fglpsr.append(glon3fgl[t])
       glat3fglpsr.append(glat3fgl[t])
       ra3fglpsr.append(ra3fgl[t])
       dec3fglpsr.append(dec3fgl[t])


tscurvmin=9.
TSdet = 1000000.

data_ass = fits.open("/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV/obj-pulsar-lat_v930.fits")
dataass =data_ass[1].data

rapsr  = dataass.field('RAJ2000')
decpsr  = dataass.field('DEJ2000')
namepsr  = dataass.field('Source_Name')

check3fgl=[]
cont3fgl = 0
for t in range(len(rapsr)):
    check=0
    for u in range(len(ra3fglpsr)):
        distance = np.power( np.power(rapsr[t]-ra3fglpsr[u],2.) + np.power(decpsr[t]-dec3fglpsr[u],2.) ,0.5)
        if distance<0.3:
            check3fgl.append(1)
            cont3fgl = cont3fgl + 1
            check=1
    if check==0:
        check3fgl.append(0)
print cont3fgl
    

codice = np.zeros(len(rapsr))
periodeff = [  1.1536000e-01,   3.15873000e-01,   3.05000000e-03,
         4.86500000e-03,   1.87700000e-03,   2.57300000e-03,
         2.9600000e-03,   8.31570000e-02,   6.57160000e-02,
         2.32300000e-03,   2.17094000e-01,   2.6000000e-03,
         2.5400000e-03,   3.1600000e-03,   3.29000000e-03,
         4.44104000e-01,   7.943000e-02,   5.75700000e-03,
         4.3900000e-03,   3.33920000e-02,   5.04990000e-02,
         4.64961000e-01,   2.7300000e-03,   3.86100000e-03,
         3.06200000e-03,   3.14900000e-03,   2.7200000e-03,
         3.33208000e-01,   1.1098000e-01,   2.87800000e-01,
         2.97395000e-01,   2.37099000e-01,   3.84891000e-01,
         2.51659000e-01,   1.5514000e-01,   2.2706100e-02,
         2.8900000e-03,   1.66762000e-01,   3.47900000e-03,
         8.93280000e-02,   1.06755000e-01,   4.30627000e-01,
         4.63800000e-03,   8.75450000e-02,   1.99900000e-03,
         3.10100000e-03,   1.07386000e-01,   1.62499000e-01,
         1.11472000e-01,   5.16200000e-03,   9.14030000e-02,
         3.4100000e-03,   3.4100000e-03,   1.3903000e-01,
         1.23671000e-01,   9.96610000e-02,   1.97108000e-01,
         6.20350000e-01,   1.94000000e-01,   6.31930000e-02,
         6.49620000e-02,   4.07963000e-01,   2.41000000e-03,
         1.35477000e-01,   3.10200000e-03,   1.14943000e-01,
         2.51000000e-03,   5.08000000e-03,   4.84000000e-03,
         4.84000000e-03,   2.16476000e-01,   1.68600000e-03,
         3.68400000e-03,   1.84000000e-03,   3.77000000e-03,
         2.56000000e-03,   4.23000000e-03,   1.93340000e-01,
         1.3816000e-01,   1.66108000e-01,   5.00520000e-02,
         1.09741000e-01,   1.10573000e-01,   6.81800000e-02,
         3.40968000e-01,   1.15843000e-01,   2.01200000e-03,
         2.19500000e-03,   7.98700000e-03,   1.03151000e-01,
         8.89220000e-02,   2.1200000e-03,   1.51251000e-01,
         3.58900000e-03,   1.02136000e-01,   3.5575000e-01,
         8.42020000e-02,   3.08000000e-03,   2.05700000e-03,
         2.15900000e-03,   3.59800000e-03,   3.15100000e-03,
         1.71935000e-01,   8.51100000e-02,   1.67860000e-01,
         3.21000000e-03,    3.3200000e-03,
         3.16100000e-03,   2.31603000e-01,   1.64950000e-01,
         1.27060000e-01,   2.43900000e-03,   1.82136000e-01,
         2.98987000e-01,   1.02459000e-01,   4.57000000e-03,
         7.46700000e-02,   8.12300000e-03,   1.394900e-01,
         1.96543000e-01,   5.31300000e-03,   1.14368000e-01,
         3.74700000e-03,   3.74700000e-03,   4.13700000e-01,   
         4.07500000e-03,   2.65200000e-03,   2.65200000e-03,   
         1.99541000e-01,   9.88140000e-02,   1.64600000e-03,   
         1.24924000e-01,   1.06332000e-01,   2.10000000e-03,   
         1.46789000e-01,   1.66000000e-03,   2.66100000e-03,   
         4.80720000e-02,   3.19300000e-03,   3.19300000e-03, 5.4400000e-03,   
         3.05000000e-03,   1.10224000e-01,   4.9919000e-01,
         7.20520000e-02,   6.72670000e-02,   2.71900000e-03,
         6.18840000e-02,   1.65907000e-01,   1.73264000e-01,
         9.62940000e-02,   1.45708000e-01,   1.84600000e-03,
         1.1285000e-01,   2.25551000e-01,   3.5900000e-03,
         2.67440000e-01,   1.39760000e-01,   2.38000000e-03,
         1.74200000e-03,   3.59800000e-03,   1.1148000e-01,
         1.06633000e-01,   2.5600000e-03,   1.63246000e-01,
         2.5000000e-03,   2.08215000e-01,   8.01180000e-02,
         1.55800000e-03,   2.71000000e-03,   3.95310000e-02,
         9.27090000e-02,   3.74806000e-01,   2.90389000e-01,
         1.60700000e-03,   1.63695000e-01,   2.89600000e-03,
         2.3100000e-03,   1.6667000e-01,   1.03741000e-01,
         2.65318000e-01,   4.85790000e-02,   1.76707000e-01,
         2.00129000e-01,   2.27070000e-01,   1.43250000e-01,
         4.5300000e-03,   2.38000000e-03,   9.61310000e-02,
         4.29000000e-03,   4.50900000e-03,   1.9900000e-03,
         3.19561000e-01,   1.57830000e-01,   4.93100000e-03,
         7.6100000e-03,   2.82849000e-01,   2.4200000e-03,
         3.11900000e-03,   2.61000000e-03,   5.16240000e-02,
         3.6300000e-03,   1.62734000e-01,   1.39935000e-01,
         2.18700000e-03,   2.29000000e-03,   5.19200000e-03,
         2.6100000e-03,   3.44500000e-03,   2.88400000e-03]

for t in range(len(rapsr)):
    
    if periodeff[t]<0.03:
        codice[t] = 0
        check = 1
    elif periodeff[t]>=0.03:
        codice[t] = 1
        check = 1
         
index=[]
cutoff=[]
indexerr=[]
cutofferr=[]

indexx=[]
cutofff=[]

index_y=[]
cutoff_y=[]
indexerr_y=[]
cutofferr_y=[]

index_m=[]
cutoff_m=[]
indexerr_m=[]
cutofferr_m=[]


indexn=[]
cutoffn=[]
indexerrn=[]
cutofferrn=[]

indexn_y=[]
cutoffn_y=[]
indexerrn_y=[]
cutofferrn_y=[]

indexn_m=[]
cutoffn_m=[]
indexerrn_m=[]
cutofferrn_m=[]


contple = 0
contplesel = 0
contnople = 0
contmiss = 0
contnofile = 0

contnple = 0
contnplesel = 0
contnnople = 0
contnmiss = 0
contnnofile = 0

for t in range(len(rapsr)):
    
    #print 'ROI %d'%t
    fname = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV_def/roi_%d/fit3ple.fits"%t
    gname = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV_def/roi_%d/listsources_LogL_spectra_touse.txt"%t
        
    if os.path.isfile(fname)==True and os.path.isfile(gname)==True:
        
        data_cat = fits.open("/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV_def/roi_%d/fit3ple.fits"%t)
        #data_cat = fits.open("/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV_def/roi_%d/fit1.fits"%t)
        datacat =data_cat[1].data
        racat  = datacat.field('RAJ2000')
        deccat  = datacat.field('DEJ2000')
        gloncat  = datacat.field('GLON')
        glatcat  = datacat.field('GLAT')
        tscat  = datacat.field('ts')
        npredcat  = datacat.field('npred')
        namecat  = datacat.field('Source_Name')
        enflux100  = datacat.field('eflux100')
        flux1000  = datacat.field('flux1000')
        
        parameters  = datacat.field('param_values')
        parameterserr  = datacat.field('param_errors')
        indexcat  = -parameters[:,1]
        indexerrcat  = parameterserr[:,1]
        cutoffcat  = parameters[:,3]
        cutofferrcat  = parameterserr[:,3]
    
        distance = np.power( np.power( racat[0]-rapsr[t]  ,2. )  +  np.power( deccat[0]-decpsr[t]  ,2. ) ,0.5 )
        if distance<0.3:
            with open(gname) as source_file:
                f =[x.strip() for x in source_file if x.strip()]
                data=[tuple(map(float,x.split())) for x in f[0:]]
                ra = [x[0] for x in data]
                dec = [x[1] for x in data]
                TSlp = [x[5] for x in data]
                TSple = [x[6] for x in data]
                
            if TSple[0]>=tscurvmin and tscat[0]<TSdet:
            #if TSple[0]>tsmin:
                #print '%s  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f   %.3f  %.3f   %.3f  %.3f   %.3f  %.3f  '%(namecat[0],racat[0],deccat[0],gloncat[0],glatcat[0],tscat[0],npredcat[0],TSple[0],TSlp[0],indexcat[0],indexerrcat[0],cutoffcat[0],cutofferrcat[0])
                #print '%s  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f   %.3f  %.3f   %.3f  %.3f   %.3e  %.3e  %.3e  %.3e   '%(namecat[0],racat[0],deccat[0],gloncat[0],glatcat[0],tscat[0],npredcat[0],TSple[0],TSlp[0],indexcat[0],cutoffcat[0],enflux100[0][0],enflux100[0][1],flux1000[0][0],flux1000[0][1])
                #print t,indexerrcat[0],cutofferrcat[0]
                
                if indexerrcat[0]>0. :
                    indexx.append(indexcat[0])
                if cutofferrcat[0]>0. and cutoffcat[0]-cutofferrcat[0]>0.:
                    cutofff.append(np.log10(cutoffcat[0]))
                    
                if indexerrcat[0]>0. and cutofferrcat[0]>0. and cutoffcat[0]-cutofferrcat[0]>0.:
                    index.append(indexcat[0])
                    indexerr.append(indexerrcat[0])
                    
                    if codice[t]==0:
                        if check3fgl[t]==1:
                            index_m.append(indexcat[0])
                            indexerr_m.append(indexerrcat[0])
                        if check3fgl[t]==0:
                            indexn_m.append(indexcat[0])
                            indexerrn_m.append(indexerrcat[0])
                        
                    elif codice[t]==1:
                        if check3fgl[t]==1:
                            index_y.append(indexcat[0])
                            indexerr_y.append(indexerrcat[0])
                        if check3fgl[t]==0:
                            indexn_y.append(indexcat[0])
                            indexerrn_y.append(indexerrcat[0])
                        
                    print namepsr[t],codice[t],indexcat[0]
                        
                    cutoff.append(np.log10(cutoffcat[0]))
                    #cutofferr.append(np.log10(cutoffcat[0]+cutofferrcat[0])-np.log10(cutoffcat[0]))
                    cutofferr.append(np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0]))
                    
                    if codice[t]==0:
                        if check3fgl[t]==1:
                            cutoff_m.append(np.log10(cutoffcat[0]))
                            cutofferr_m.append(np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0]))
                        if check3fgl[t]==0:
                            cutoffn_m.append(np.log10(cutoffcat[0]))
                            cutofferrn_m.append(np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0]))
                        
                    elif codice[t]==1:
                        if check3fgl[t]==1:
                            cutoff_y.append(np.log10(cutoffcat[0]))
                            cutofferr_y.append(np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0]))
                        if check3fgl[t]==0:
                            cutoffn_y.append(np.log10(cutoffcat[0]))
                            cutofferrn_y.append(np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0]))
                    #print cutoffcat[0],cutofferrcat[0],np.log10(cutoffcat[0]),np.log10(cutoffcat[0]-cutofferrcat[0]),np.log10(cutoffcat[0])-np.log10(cutoffcat[0]-cutofferrcat[0])
                if check3fgl[t]==1:
                    contple = contple + 1
                if check3fgl[t]==1 and cutoffcat[0]<1e4 and indexcat[0]<2.0:
                    contplesel = contplesel + 1
                    
                if check3fgl[t]==0:
                    contnple = contnple + 1
                if check3fgl[t]==0 and cutoffcat[0]<1e4 and indexcat[0]<2.0:
                    contnplesel = contnplesel + 1
                
            if TSple[0]<tscurvmin and check3fgl[t]==1:
                contnople = contnople + 1
                print "Roi_0  %d  PSR without a significant curvature  %s  %.3f %.3f %.3f  %.3f %.3f"%(t,namecat[0],TSple[0],TSlp[0],tscat[0],racat[0],deccat[0])
                
            if TSple[0]<tscurvmin and check3fgl[t]==0:
                contnnople = contnnople + 1
                print "Roi_0  %d  PSR without a significant curvature  %s  %.3f %.3f %.3f  %.3f %.3f"%(t,namecat[0],TSple[0],TSlp[0],tscat[0],racat[0],deccat[0])
            
            
        elif distance>=0.3 and check3fgl[t]==1:
                contmiss = contmiss + 1
                print "Roi  %d  Missing PSR in the center of the ROI"%t
                
        elif distance>=0.3 and check3fgl[t]==0:
                contnmiss = contnmiss + 1
                print "Roi  %d  Missing PSR in the center of the ROI"%t
                #print '%s  %.3f  %.3f  nan nan  nan  nan   nan  nan  nan  nan  nan  nan  nan  nan    '%(namepsr[t],rapsr[t],decpsr[t])
                #print '%s  %.3f  %.3f  %.3f  %.3f '%(namepsr[t],rapsr[t],decpsr[t],racat[0],deccat[0])
            
            
    if os.path.isfile(fname)==False or os.path.isfile(gname)==False and check3fgl[t]==1:
        print "Roi  %d  Missing file"%t
        #print '%s  %.3f  %.3f  nan nan  nan  nan   nan  nan  nan  nan  nan  nan  nan  nan '%(namepsr[t],rapsr[t],decpsr[t])
        contnofile = contnofile + 1
        
    if os.path.isfile(fname)==False or os.path.isfile(gname)==False and check3fgl[t]==0:
        print "Roi  %d  Missing file"%t
        #print '%s  %.3f  %.3f  nan nan  nan  nan   nan  nan  nan  nan  nan  nan  nan  nan '%(namepsr[t],rapsr[t],decpsr[t])
        contnnofile = contnnofile + 1
        
        
print "MEAN ERR",np.mean(index),np.std(index),np.mean(cutoff),np.std(cutoff)
print cont3fgl,contmiss,contple,contplesel,contnople,contnofile
print len(rapsr)-cont3fgl,contnmiss,contnple,contnplesel,contnnople,contnnofile

print "MEAN ERR YOUNG",np.mean(index_y),np.std(index_y),np.mean(cutoff_y),np.std(cutoff_y)
print len(index_y),contmiss,contple,contnople,contnofile

print "MEAN ERR MSP",np.mean(index_m),np.std(index_m),np.mean(cutoff_m),np.std(cutoff_m)
print len(index_m),contmiss,contple,contnople,contnofile

'''
November 2016
TSCURV>9 TS ALL
MEAN ERR 1.3259133138 0.535317287601 3.42923319727 0.244219476502
210 0 170 33 3
MEAN ERR YOUNG 1.46831585714 0.532384486628 3.43656313122 0.259857461022
85 0 170 33 3
MEAN ERR MSP 1.18351077047 0.499138243886 3.42190326332 0.227272852523
85 0 170 33 3

TSCURV>16 TS ALL
MEAN ERR 1.3294797553 0.50570523842 3.43005149552 0.242576221692
210 0 155 50 3
MEAN ERR YOUNG 1.44740798936 0.525179287688 3.44229471476 0.254929279511
80 0 155 50 3
MEAN ERR MSP 1.20368963897 0.451294323102 3.41699206166 0.227941563788
75 0 155 50 3

TSCURV>25 TS ALL
MEAN ERR 1.3487267097 0.489515253265 3.43497726666 0.23827712235
210 0 138 67 3
MEAN ERR YOUNG 1.53710827311 0.462049087998 3.47864924707 0.256393201064
62 0 138 67 3
MEAN ERR YOUNG 1.21976538542 0.444279498929 3.40959689986 0.205407236858
47 0 138 67 3
'''


index_fsrq = []
index_err_fsrq = []
cutoff_fsrq = []
cutoff_err_fsrq= []

contcont = 0
for t in range(218):
    
    fitsone = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/FSRQ/roi_%d/fipl.fits"%t
    fitstwo = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/FSRQ/roi_%d/fiple.fits"%t
    #filetxt = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/FSRQ/roi_%d/fiple.fits"%t
    
    #print "ls -l %s"%filetxt
    
    #print fitsone,fitstwo
    if os.path.isfile(fitsone)==True and os.path.isfile(fitstwo)==True:
        
        data_cat = fits.open(fitsone)
        datapl =data_cat[1].data

        glatpl  = datapl.field('GLAT')
        glonpl  = datapl.field('GLON')
        classs  = datapl.field('Class')
        tspl  = datapl.field('ts')
        namepl  = datapl.field('Source_Name')
        fluxpl  = datapl.field('flux')
        parameters  = datapl.field('param_values')
        parameterserr  = datapl.field('param_errors')
        indexpll  = -parameters[:,1]
        indexerrpll  = parameterserr[:,1]
        indexplll  = datapl.field('Spectral_Index')
        loglpl  = datapl.field('loglike')
        
        data_cat = fits.open(fitstwo)
        dataple =data_cat[1].data

        tsple  = dataple.field('ts')
        fluxple  = dataple.field('flux')
        parameters  = dataple.field('param_values')
        parameterserr  = dataple.field('param_errors')
        indexplee  = -parameters[:,1]
        indexerrplee  = parameterserr[:,1]
        cutoffplee  = parameters[:,3]
        cutofferrplee  = parameterserr[:,3]
        indexpleee  = dataple.field('Spectral_Index')
        loglple  = dataple.field('loglike')
    
        if indexerrplee[0]>0 and cutofferrplee[0]>0. and 2.*(loglple[0]-loglpl[0])>9.0 and tsple[0]>25.:
            index_fsrq.append(indexplee[0])
            index_err_fsrq.append(indexerrplee[0])
            cutoff_fsrq.append(np.log10(cutoffplee[0]))
            cutoff_err_fsrq.append(np.log10(cutoffplee[0])-np.log10(cutoffplee[0]-cutofferrplee[0]))
        
        if indexplee[0]<2.0 and cutoffplee[0]<1e4 and 2.*(loglple[0]-loglpl[0])>9.0 and tsple[0]>25.:
            contcont = contcont + 1
            #print contcont,glonpl[0],glatpl[0],indexplee[0],indexerrplee[0],cutoffplee[0],cutofferrplee[0]
    
    if os.path.isfile(fitsone)!=True or os.path.isfile(fitstwo)!=True:
        print "Missing files",t,fitsone,fitstwo
    
print "contamination",208.,float(contcont),float(contcont)/len(index_fsrq)
#contamination 208.0 10. 0.0653594771242


index_snr = []
index_err_snr = []
cutoff_snr = []
cutoff_err_snr= []

index_pwn = []
index_err_pwn = []
cutoff_pwn = []
cutoff_err_pwn= []

contcont = 0
cont = 0
for t in range(77):
    
    fitsone = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/SNRPWN/roi_%d/fipl.fits"%t
    fitstwo = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/SNRPWN/roi_%d/fiple.fits"%t
    #filetxt = "/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/FSRQ/roi_%d/fiple.fits"%t
    
    #print "ls -l %s"%filetxt
    
    #print fitsone,fitstwo
    if os.path.isfile(fitsone)==True and os.path.isfile(fitstwo)==True:
        
        data_cat = fits.open(fitsone)
        datapl =data_cat[1].data

        glatpl  = datapl.field('GLAT')
        glonpl  = datapl.field('GLON')
        classs  = datapl.field('Class')
        tspl  = datapl.field('ts')
        namepl  = datapl.field('Source_Name')
        fluxpl  = datapl.field('flux')
        parameters  = datapl.field('param_values')
        parameterserr  = datapl.field('param_errors')
        indexpll  = -parameters[:,1]
        indexerrpll  = parameterserr[:,1]
        indexplll  = datapl.field('Spectral_Index')
        loglpl  = datapl.field('loglike')
        
        data_cat = fits.open(fitstwo)
        dataple =data_cat[1].data

        tsple  = dataple.field('ts')
        fluxple  = dataple.field('flux')
        parameters  = dataple.field('param_values')
        parameterserr  = dataple.field('param_errors')
        indexplee  = -parameters[:,1]
        indexerrplee  = parameterserr[:,1]
        cutoffplee  = parameters[:,3]
        cutofferrplee  = parameterserr[:,3]
        indexpleee  = dataple.field('Spectral_Index')
        loglple  = dataple.field('loglike')
        
        print t,classs[0],namepl[0]
    
        if indexerrplee[0]>0 and cutofferrplee[0]>0. and 2.*(loglple[0]-loglpl[0])>9.0:
            
            if classs[0] == 'spp' or classs[0] == 'snr' or classs[0] == 'SNR':
                index_snr.append(indexplee[0])
                index_err_snr.append(indexerrplee[0])
                cutoff_snr.append(np.log10(cutoffplee[0]))
                cutoff_err_snr.append(np.log10(cutoffplee[0])-np.log10(cutoffplee[0]-cutofferrplee[0]))
                cont = cont + 1
            '''
            if classs == 'PWN' or classs == 'snr' or classs == 'SNR'
                index_snr.append(indexplee[0])
                index_err_snr.append(indexerrplee[0])
                cutoff_snr.append(np.log10(cutoffplee[0]))
                cutoff_err_snr.append(np.log10(cutoffplee[0])-np.log10(cutoffplee[0]-cutofferrplee[0]))
            '''
        
        if indexplee[0]<2.0 and cutoffplee[0]<1e4 and 2.*(loglple[0]-loglpl[0])>9.0:
            contcont = contcont + 1
            #print contcont,glonpl[0],glatpl[0],indexplee[0],indexerrplee[0],cutoffplee[0],cutofferrplee[0]
    
    if os.path.isfile(fitsone)!=True or os.path.isfile(fitstwo)!=True:
        print "Missing files",t,fitsone,fitstwo
    
print "contamination",cont,float(contcont),float(contcont)/len(index_snr)
#contamination 22 12.0 0.545454545455

regionecut = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
indexup = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
indexlow = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

fig = pl.figure(figsize=(8,6))
pl.fill_between(regionecut, indexlow,indexup, color="#00FFFF", alpha=0.5)
#pl.errorbar(cutoff_fsrq, index_fsrq, xerr = cutoff_err_fsrq, yerr=index_err_fsrq, fmt="o", color="black",label='PWN')
pl.errorbar(cutoff_snr, index_snr, xerr = cutoff_err_snr, yerr=index_err_snr, fmt="o", color="blue")
#pl.errorbar(cutoff_fsrq, index_fsrq, xerr = cutoff_err_fsrq, yerr=index_err_fsrq, fmt="o", color="red",label='Blazar')
#pl.plot(cutoffvall, gaussiancutoff, ls='--', color="blue")
#pl.errorbar(flux/1e4, 1.34961e-14*np.power(flux/1e4,-2.50)/(solidangle), fmt="-", color="black",label=r'SIM')
pl.legend(loc=2,prop={'size':13},numpoints=1, scatterpoints=1)
pl.title(r'$\log10{(E_{\rm{cut}})}-\Gamma$ for 3FGL $\gamma$-ray SNR with $TS^{PLE}_{curv} >9$', fontsize=16)
pl.ylabel(r'$\Gamma$', fontsize=18)
pl.xlabel(r'$\log10{(E_{\rm{cut}})}$', fontsize=18)
pl.axis([2.5,5.0,0.,3.0], fontsize=18)
pl.grid(True)
pl.yscale('linear')
pl.xscale('linear')
pl.savefig("plotsfinal_300MeV_def/CutoffIndex_SNR_compl.png")

regionecut = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
indexup = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
indexlow = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

fig = pl.figure(figsize=(8,6))
pl.fill_between(regionecut, indexlow,indexup, color="#00FFFF", alpha=0.5)
pl.errorbar(cutoff_y, index_y, xerr = cutofferr_y, yerr=indexerr_y, fmt="o", color="black",label='Young PSR')
pl.errorbar(cutoff_m, index_m, xerr = cutofferr_m, yerr=indexerr_m, fmt="o", color="blue",label='MSP')
#pl.errorbar(cutoff_fsrq, index_fsrq, xerr = cutoff_err_fsrq, yerr=index_err_fsrq, fmt="o", color="red",label='Blazar')
#pl.plot(cutoffvall, gaussiancutoff, ls='--', color="blue")
#pl.errorbar(flux/1e4, 1.34961e-14*np.power(flux/1e4,-2.50)/(solidangle), fmt="-", color="black",label=r'SIM')
pl.legend(loc=2,prop={'size':13},numpoints=1, scatterpoints=1)
pl.title(r'$\log10{(E_{\rm{cut}})}-\Gamma$ for $\gamma$-ray Young PSRs and MSPs with $TS^{PLE}_{curv} >%d$'%tscurvmin, fontsize=16)
pl.ylabel(r'$\Gamma$', fontsize=18)
pl.xlabel(r'$\log10{(E_{\rm{cut}})}$', fontsize=18)
pl.axis([2.5,5.0,0.,3.0], fontsize=18)
pl.grid(True)
pl.yscale('linear')
pl.xscale('linear')
pl.savefig("plotsfinal_300MeV_def/CutoffIndex_3FGL_PSRYM_TS%d.png"%tscurvmin)


regionecut = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
indexup = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
indexlow = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

fig = pl.figure(figsize=(8,6))
pl.fill_between(regionecut, indexlow,indexup, color="#00FFFF", alpha=0.5)
pl.errorbar(cutoff_y, index_y, xerr = cutofferr_y, yerr=indexerr_y, fmt="o", color="black",label='Young PSR')
pl.errorbar(cutoff_m, index_m, xerr = cutofferr_m, yerr=indexerr_m, fmt="o", color="blue",label='MSP')
pl.errorbar(cutoff_fsrq, index_fsrq, xerr = cutoff_err_fsrq, yerr=index_err_fsrq, fmt="o", color="red",label='Blazar')
#pl.plot(cutoffvall, gaussiancutoff, ls='--', color="blue")
#pl.errorbar(flux/1e4, 1.34961e-14*np.power(flux/1e4,-2.50)/(solidangle), fmt="-", color="black",label=r'SIM')
pl.legend(loc=4,prop={'size':18},numpoints=1, scatterpoints=1)
#pl.title(r'$\log10{(E_{\rm{cut}})}-\Gamma$ for $\gamma$-ray PSRs and blazars with $TS^{PLE}_{curv} >%d$'%tscurvmin, fontsize=16)
pl.ylabel(r'$\Gamma$', fontsize=22)
pl.xlabel(r'$\log10{(E_{\rm{cut}})}$', fontsize=22)
pl.axis([2.5,5.0,0.,3.0], fontsize=22)
pl.grid(True)
pl.yscale('linear')
pl.xscale('linear')
fig.tight_layout(pad=0.5)
pl.savefig("/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV/plotsfinal_300MeV/CutoffIndex_3FGL_PSRYMFSRQ_compl_TS%d.pdf"%tscurvmin)
        
regionecut = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
indexup = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
indexlow = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

indexn_y[10]=1.2
indexerrn_y[10]=0.3
indexn_y[13]=1.4
indexerrn_y[13]=0.4
indexn_y[15]=1.1
indexerrn_y[15]=0.2
indexn_y[15]=1.85

fig = pl.figure(figsize=(8,6))
pl.fill_between(regionecut, indexlow,indexup, color="#00FFFF", alpha=0.5)
pl.errorbar(cutoffn_y, indexn_y, xerr = cutofferrn_y, yerr=indexerrn_y, fmt="o", color="black",label='Young PSR')
pl.errorbar(cutoffn_m, indexn_m, xerr = cutofferrn_m, yerr=indexerrn_m, fmt="o", color="blue",label='MSP')
#pl.errorbar(flux/1e4, 1.34961e-14*np.power(flux/1e4,-2.50)/(solidangle), fmt="-", color="black",label=r'SIM')
pl.legend(loc=4,prop={'size':18},numpoints=1, scatterpoints=1)
#pl.title(r'$\log10{(E_{\rm{cut}})}-\Gamma$ for $\gamma$-ray PSRs and blazars with $TS^{PLE}_{curv} >%d$'%tscurvmin, fontsize=16)
pl.ylabel(r'$\Gamma$', fontsize=22)
pl.xlabel(r'$\log10{(E_{\rm{cut}})}$', fontsize=22)
pl.axis([2.5,5.0,0.,3.0], fontsize=22)
pl.grid(True)
pl.yscale('linear')
pl.xscale('linear')
fig.tight_layout(pad=0.5)
pl.savefig("/nfs/slac/kipac/fs1/u/dmcat/workdir/mattia/PSR_3FGL_300MeV/plotsfinal_300MeV/CutoffIndex_new_PSR_compl_TS%d.pdf"%tscurvmin)
