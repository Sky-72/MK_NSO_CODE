#!/usr/bin/env python
# coding: utf-8

# # Imports and Class

# In[ ]:


import numpy as np
import sunpy as sp
import cv2
import matplotlib.pyplot as plt
import matplotlib.colors
import scipy.io as sc
from scipy import ndimage
from astropy.io import fits 
#import matplotlib.animation as ani
import glob
from scipy.optimize import fmin
from scipy.ndimage.measurements import center_of_mass
from scipy.interpolate import UnivariateSpline
import matplotlib.patches as patches
get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


class SolarData:
    
    def getData(self, filename):
        '''Function to read a .SAV IDL file
        -----------------------------------
        Imput: .Sav file extension
        '''
        
        mData = sc.readsav(filename, verbose = True)
        return mData
    
    def NBplot(self, name, wavel, num = 0):
        '''Function for Plotting Narrowband Images
        -------------------------------------------
        Imputs: 
            name = Name of file \n
            wavel = wavelength of images \n
            num = image number (Set at 0 initially) \n
        -------------------------------------------- \n
        Outputs = Image of sun'''
        
        sData = name.nb_data_dstr
        plt.figure(figsize=(10,10))
        plt.title('Sun at {} angstroms'.format(wavel))
        plt.imshow(sData[num,:,:], interpolation='none', cmap= 'gist_gray', origin = 'lower')
        plt.colorbar()
        plt.show()
        
    def BBplot(self, name, wavel, num = 0):
        '''Function for Plotting BroadBand Images
        -------------------------------------------
        Imputs:
            name = Name of file \n
            wavel = wavelength of images \n
            num = image number (Set at 0 initially) \n
        -------------------------------------------- \n
        Outputs = Image of sun'''
        
        sData = name.bb_data_dstr
        plt.figure(figsize=(10,10))
        plt.title('Sun at {} angstroms'.format(wavel))
        plt.imshow(sData[num,:,:], interpolation='none', cmap= 'gist_gray', origin = 'lower')
        plt.colorbar()
        plt.show()
        
    def Createcontoursgraph(self, name, wavel, num = 0):
        '''Function for Plotting a Contour Image
        -------------------------------------------
        Imputs: \n
            name = Name of file \n
            wavel = wavelength of images \n
            num = image number (Set at 0 initially) \n
        -------------------------------------------- \n
        Outputs:  Image of sun contour lines'''
        
        mdata = name.nb_data_dstr[num,:,:]
        contours = ndimage.filters.gaussian_filter(mdata,5,mode='nearest')
        plt.figure(figsize=(10,10))
        plt.title('Sun at {} angstroms'.format(wavel))
        plt.imshow(mdata, interpolation='none', cmap= 'gist_gray', origin = 'lower')
        plt.colorbar()
        plt.contour(contours)
        plt.show()
        
    
    def rowplot(self, name, num = 0, y = 0):
        '''Function to create line plot from cut of image
        -------------------------------------------------
        Imputs: \n
            Name = name of file \n
            num = image number \n
            y = row in y \n
        Outputs: \n
            Line plot of light curve  \n
        '''
        gdata = name.nb_data_dstr
        plt.figure(figsize=(10,10))
        plt.title(f'Image {num}')
        plt.ylabel('Light Intesity')
        plt.xlabel('x position')
        plt.plot(gdata[num,y,:])
        plt.show()
                
    def findminx(self, name):
        
        xmdata = name.nb_data_dstr
        minx = []
        for i in range(len(xmdata[0,:,:])):
            minx.append([])
            for c in range(len(xmdata[0,:,:])):
                count = 4
                x = []
                y = []
                while count <= 8:
                    dat = xmdata[count,i,c]
                    x.append(count)
                    y.append(dat)
                    count += 1
                var = np.polyfit(x,y,2)
                a = var[0]
                b = var[1]
                d = var[2]
                eq = lambda x: a*x**2 + b*x + d 
                if a > 0:
                    mini = fmin(eq,6,disp=0)
                    minx[i].append(mini.tolist())
                else:
                    mini = fmin(lambda z: -eq(z), 6, disp = 0)
                    minx[i].append(mini.tolist())
                    continue
                
        xmin = minx
        return xmin
    
    def findminy(self, name):
        
        xmdata = name.nb_data_dstr
        miny = []
        for i in range(len(xmdata[0,:,:])):
            miny.append([])
            for c in range(len(xmdata[0,:,:])):
                count = 4
                x = []
                y = []
                while count <= 8:
                    dat = xmdata[count,i,c]
                    x.append(count)
                    y.append(dat)
                    count += 1
                var = np.polyfit(x,y,2)
                a = var[0]
                b = var[1]
                d = var[2]
                eq = lambda z: a*z**2 + b*z + d 
                if a > 0:
                    mini = fmin(eq,6,disp=0)
                    ysol = eq(mini)
                    miny[i].append(ysol.tolist())
                else:
                    mini = fmin(lambda z: -eq(z), 6, disp = 0)
                    ysol = eq(mini)
                    miny[i].append(ysol.tolist())
                    continue
                    
        ymin = miny
        
        return ymin
    
    def Mystats(self,name, x1 = 0, x2 = 999, y1 = 0, y2 = 999):
        sdata = name.bb_data_dstr
        mean = np.mean(sdata[:,x1:x2,y1:y2])
        rms = np.std(sdata[:,x1:x2,y1:y2])
        per = rms/mean
        
        return mean, rms, per
        
    def rms(self,name,x1 = 0, x2 = 999, y1 = 0, y2 = 999):
        rmsdata = name.bb_data_dstr
        rmsarray = []
        for i in range(len(rmsdata[:,0,0])):
            rms = np.std(rmsdata[i,y1:y2,x1:x2])
            rmsarray.append(rms)
        
        return rmsarray
    
    def brms(self,name,x1 = 0, x2 = 999, y1 = 0, y2 = 999):
        rmsdata = name.bb_data_dstr
        
        rms = np.std(rmsdata[:,y1:y2,x1:x2])

        return rms
    
    def linewidth(self,name,mp1 = 2, mp2 = 10):
        xmdata = name.nb_data_dstr
        width = []
        for i in range(len(xmdata[0,:,:])):
            width.append([])
            for c in range(len(xmdata[0,:,:])):
                count = 2
                x = []
                y = []
                while count <= 10:
                    dat = xmdata[count,i,c]
                    x.append(count)
                    y.append(dat)
                    count += 1
                var = UnivariateSpline(x,y)
                
                p1 = var(mp1)
                p2 = var(mp2)
                mini = fmin(var,6,disp=0)
                ysol = eq(mini)
                avemax = (p1+p2)/2
                if avemax > ysol:
                    midnum = (avemax - ysol) / 2
                midpnt = ysol + midnum
                newd = d - midpnt
                vari = [a,b,newd]
                points = np.roots(vari)
                lw = points[0] - points[1]
                width[i].append(lw)
            
               
        linewidth = np.array(width)
        
        return linewidth
    
    def cogfit(self, name, num = 6, relint = 2300):
        '''Function to Fit data with Center of Gravity
        -----------------------------------------------
        Inputs:
        - num = Image number'''
        cogdata = name.nb_data_dstr
        x = []
        y = []
        for i in range(len(cogdata[0,:,:])):
            com = center_of_mass(cogdata[num,i,:]>relint)
            y.append(i)
            x.append(com)
        
        return x, y 
    
    def cogfitmk2(self, data, num = 6, lim = 2500):
    
        cogdata = data.nb_data_crrt
        x = []
        y = []
        for i in range(len(cogdata[num,:,:])):
            com = center_of_mass(cogdata[num,i,:] * (cogdata[num,i,:] > lim))
            x.append(com)
            y.append(i)
            
        return x, y 
                
    def pre6173(self, name):
        
        predata = name.nb_data_dstr
        lnth = len(predata[:,500,500])
        precurve = []
        
        for i in range(len(predata[0,0,:])):
            
            anum = (predata[0,500,500]-predata[9,500,500])
            slope = anum/9
            eq3 = lambda x: -slope*x+data2.nb_data_dstr[0,500,500]
    
    
    def openall(sefl, wavel = 6563, band = 'nb'):
        dataname = glob.glob(f'/workspace/mkatilius/IDL/mData/20141025/results_000/{wavel}/' +band+ f'/{wavel}_'+ band+'nb*.sav')
        data = []
        for i in range(len(dataname)):
            datna = sc.readsav(dataname[i], verbose = True)
            data.append(datna)
        
        return data
    
 #   def prefiltcorrect(self, wavel = 6563)
        
            
    def Createmovie(self,wavel, num = 0):
        ''' WIP 
        Function to Create a movie of one image from whole data set
        -----------------------------------------------------------
        Input: \n 
            Wavel = wavelength \n
            num = image number \n
        Output: \n
            Movie of sun to see movement of gasses'
        -----------------------------------------------------------'''
      
        FFMpegWriter = ani.writers['ffmpeg']

        writer = FFMpegWriter(fps=10)

        fig = plt.figure(figsize=(10, 10))
        
        filename = glob.glob(F"/workspace/mkatilius/IDL/mData/20141025/results_000/{wavel}/nb/{wavel}_nb*.sav")
        
        with writer.saving(fig, "Test.mp4", 100):
            
            for i, f in enumerate(filename):
                file = sc.readsav(f)
                im = file.nb_data_dstr[num,:,:]
                
                plt.title(f'Sun at {wavel} angstroms')
                plt.imshow(im, interpolation='none', cmap= 'gist_gray', origin = 'lower')
                writer.grab_frame()


# # Solar Data Input

# In[ ]:


S = SolarData()


# In[ ]:


data1 = S.getData('/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6563_nbn/6563_nbn000.sav')
data2 = S.getData('/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6173_nbn/6173_nbn000.sav')
#data3 = S.getData('/workspace/mkatilius/IDL/mData/20141025/results_000/6563/nb/6563_nb173.sav')
#'/workspace/mkatilius/IDL/mData/20141025/results_000/6563/nb/
#'/workspace/mkatilius/IDL/mData/20141025/calibration/6563_prefilter.sav'
#'D:\\Coding\\Harris\\IDL\\DATA_Minda\\20141025\\results_000\\6563\\nb\\6563_nb000.sav'


# In[ ]:


xtickha = ['0',f'{200*.0954204:0.3f}',f'{400*.0954204:0.3f}',f'{600*.0954204:0.3f}',
         f'{800*.0954204:0.3f}',f'{1000*.0954204:0.3f}']
ytickha = ['0',f'{200*.0974302:0.3f}',f'{400*.0974302:0.3f}',f'{600*.0974302:0.3f}',
         f'{800*.0974302:0.3f}',f'{1000*.0974302:0.3f}']
xtickir = ['0',f'{200*.0953200:0.4f}',f'{400*.0953200:0.4f}',f'{600*.0953200:0.4f}',
         f'{800*.0953200:0.4f}',f'{1000*.0953200:0.4f}']
ytickir = ['0',f'{200*.0974150:0.3f}',f'{400*.0974150:0.3f}',f'{600*.0974150:0.3f}',
         f'{800*.0974150:0.3f}',f'{1000*.0974150:0.3f}']


# In[ ]:


f, ax = plt.subplots(nrows = 2, ncols = 3, figsize = (20,12), sharey = True, sharex = True)
f.subplots_adjust(wspace=0.0, hspace=0.2)
f.suptitle('H-Alpha Scans of the Sun', fontsize = 30)
waveha = data1.info_nb['WAVE'][0]

plt.setp(ax, xticks=[0, 200, 400, 600, 800, 999], xticklabels=ytickha ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickha)

ax[0,0].imshow(data1.nb_data_crrt[2,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax[1].imshow(data1.nb_data_crrt[3,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax[0,1].imshow(data1.nb_data_crrt[4,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax[3].imshow(data1.nb_data_crrt[5,:,:], origin = 'lowerleft', cmap = 'gray')
ax[0,2].imshow(data1.nb_data_crrt[6,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax[5].imshow(data1.nb_data_crrt[7,:,:], origin = 'lowerleft', cmap = 'gray')
ax[1,0].imshow(data1.nb_data_crrt[8,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax[7].imshow(data1.nb_data_crrt[9,:,:], origin = 'lowerleft', cmap = 'gray')
ax[1,1].imshow(data1.nb_data_crrt[10,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax[9].imshow(data1.nb_data_crrt[11,:,:], origin = 'lowerleft', cmap = 'gray')
ax[1,2].imshow(data1.nb_data_crrt[12,:,:].T, origin = 'lowerleft', cmap = 'gray')

ax[0,0].set_title(f'{waveha[2]/10:.3f} nm', fontsize = 14)
#ax[1].set_title(f'{waveha[3]/10:.2f} nm')
ax[0,1].set_title(f'{waveha[4]/10:.3f} nm', fontsize = 14)
#ax[3].set_title(f'{waveha[5]/10:.2f} nm')
ax[0,2].set_title(f'{waveha[6]/10:.3f} nm', fontsize = 14)
#ax[5].set_title(f'{waveha[7]/10:.2f} nm')
ax[1,0].set_title(f'{waveha[8]/10:.3f} nm', fontsize = 14)
#ax[7].set_title(f'{waveha[9]/10:.2f} nm')
ax[1,1].set_title(f'{waveha[10]/10:.3f} nm', fontsize = 14)
#ax[9].set_title(f'{waveha[11]/10:.2f} nm')
ax[1,2].set_title(f'{waveha[12]/10:.3f} nm', fontsize = 14)

#ax[1,0].set_xticklabels(ytickha)
#ax[1,1].set_xticklabels(ytickha)
#ax[1,2].set_xticklabels(ytickha)
#ax[0,0].set_yticklabels(xtickha)
#ax[1,0].set_yticklabels(xtickha)

ax[0,0].set_ylabel('y [arcsec]')
ax[1,0].set_ylabel('y [arcsec]')

ax[1,0].set_xlabel('x [arcsec]')
ax[1,1].set_xlabel('x [arcsec]')
ax[1,2].set_xlabel('x [arcsec]')
#ax[].set_xlabel('x [arcsec]')
#ax[4].set_xlabel('x [arcsec]')
#ax[5].set_xlabel('x [arcsec]')

f.savefig('Wavelength_points_HA.png')


# In[ ]:


f2, ax2 = plt.subplots(nrows = 2, ncols = 3, figsize = (20,12), sharey = True)
f2.subplots_adjust(wspace=0.0, hspace=0.2)
f2.suptitle('Iron I Scans of the Sun', fontsize = 25)
waveir = data2.info_nb['WAVE'][0]

plt.setp(ax2, xticks=[0, 200, 400, 600, 800, 999], xticklabels = ytickir ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickir)

ax2[0,0].imshow(data2.nb_data_crrt[0,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax2[1].imshow(data2.nb_data_crrt[1,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax2[0,1].imshow(data2.nb_data_crrt[2,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax2[3].imshow(data2.nb_data_crrt[3,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax2[0,2].imshow(data2.nb_data_crrt[4,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax2[5].imshow(data2.nb_data_crrt[5,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax2[1,0].imshow(data2.nb_data_crrt[6,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax2[7].imshow(data2.nb_data_crrt[7,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax2[1,1].imshow(data2.nb_data_crrt[8,:,:].T, origin = 'lowerleft', cmap = 'gray')
#ax2[9].imshow(data2.nb_data_crrt[9,:,:].T, origin = 'lowerleft', cmap = 'gray')

ax2[0,0].set_title(f'{waveir[0]/10:.3f} nm', fontsize = 16)
#ax2[1].set_title(f'Sun at Wavelength {waveir[1]/10:.2f}')
ax2[0,1].set_title(f'{waveir[2]/10:.3f} nm', fontsize = 16)
#ax2[3].set_title(f'Sun at Wavelength {waveir[3]/10:.2f}')
ax2[0,2].set_title(f'{waveir[4]/10:.3f} nm', fontsize = 16)
#ax2[5].set_title(f'Sun at Wavelength {waveir[5]/10:.2f}')
ax2[1,0].set_title(f'{waveir[6]/10:.3f} nm', fontsize = 16)
#ax2[7].set_title(f'Sun at Wavelength {waveir[7]/10:.2f}')
ax2[1,1].set_title(f'{waveir[8]/10:.3f} nm', fontsize = 16)
#ax2[9].set_title(f'Sun at Wavelength {waveir[9]/10:.2f}')

ax2[1,2].axis('off')
#ax2[1,0].set_xticklabels(ytickir)
#ax2[1,1].set_xticklabels(ytickir)
#ax2[0,0].set_xticklabels(ytickir)
#ax2[0,1].set_xticklabels(ytickir)
#ax2[0,2].set_xticklabels(ytickir)
#ax2[0,0].set_yticklabels(xtickir)
#ax2[0,1].set_yticklabels(xtickir)

ax2[0,0].set_ylabel('y [arcsec]')
ax2[1,0].set_ylabel('y [arcsec]')
ax2[0,0].set_xlabel('x [arcsec]')
ax2[0,1].set_xlabel('x [arcsec]')
ax2[0,2].set_xlabel('x [arcsec]')
ax2[1,0].set_xlabel('x [arcsec]')
ax2[1,1].set_xlabel('x [arcsec]')

f2.savefig('Wavelength_points_FEI.png')


# # Alignment Test

# In[ ]:


x = []
for i in range(len(data1.nb_data_dstr[0,0,:])):
    x.append([])
    for a in range(len(data1.nb_data_dstr[0,0,:])):
        if data1.nb_data_dstr[0,i,a] == data1.nb_data_dstr[0,0,0]:
            x[i].append(0)
            continue
        else:
            x[i].append(data1.nb_data_dstr[0,i,a])
y = []
for i in range(len(data2.nb_data_dstr[0,0,:])):
    y.append([])
    for a in range(len(data2.nb_data_dstr[0,0,:])):
        if data2.nb_data_dstr[0,i,a] == data2.nb_data_dstr[0,0,0]:
            y[i].append(0)
            continue
        else:
            y[i].append(data2.nb_data_dstr[0,i,a])
xmin, xmax = np.percentile(x,(5,95))
ymin, ymax = np.percentile(y,(5,95))


# In[ ]:


plt.figure(figsize = (10,10))
plt.imshow(data1.nb_data_dstr[6,:,:].T,origin = 'lowerleft', cmap = 'gray')
plt.axis('off')
plt.savefig('6563_nb000_aligntest.jpg',transparent=True)
#plt.figure(figsize = (10,10))
#plt.imshow(data2.nb_data_dstr[0,:,:],origin = 'lowerleft', cmap = 'gray')
#plt.axis('off')
#plt.savefig('6173_nb000_aligntest.jpg',transparent=True)


# # Masks and COG Fits

# In[ ]:


x,y = S.cogfitmk2(data1,num = 6, lim = 2300)
#x1,y1 = S.cogfit(data3, num = 6, relint = 2300)


# In[ ]:


plt.figure(figsize = (10,10))
plt.imshow(data1.nb_data_crrt[6,:,:].T * (data1.nb_data_crrt[6,:,:].T>2300),origin = 'lowerleft', cmap = 'gray')
plt.scatter(y,x, c = 'red', s = 0.5)
plt.show()
plt.figure(figsize = (10,10))
plt.imshow(data1.nb_data_crrt[6,:,:].T,origin = 'lowerleft', cmap = 'gray')
plt.scatter(y,x, c = 'red', s = 0.5)


# In[ ]:


xtickha = ['0',f'{200*.0954204:0.3f}',f'{400*.0954204:0.3f}',f'{600*.0954204:0.3f}',
         f'{800*.0954204:0.3f}',f'{1000*.0954204:0.3f}']
ytickha = ['0',f'{200*.0974302:0.3f}',f'{400*.0974302:0.3f}',f'{600*.0974302:0.3f}',
         f'{800*.0974302:0.3f}',f'{1000*.0974302:0.3f}']
plt.figure(figsize = (10,10))
plt.imshow(data1.nb_data_crrt[6,:,:].T * (data1.nb_data_crrt[6,:,:].T>2300),origin = 'lowerleft', cmap = 'gray')
#plt.scatter(x,y, c = 'red', s = 0.5)
plt.title('Flare Mask')
plt.xlabel('x [Arcsec]')
plt.ylabel('y [Arcsec]')
plt.xticks(np.arange(0,1001,200), xtickha)
plt.yticks(np.arange(0,1001,200), ytickha)
plt.savefig('Flare_mask.png')
plt.show()
plt.figure(figsize = (10,10))
plt.imshow(data1.nb_data_crrt[6,:,:].T,origin = 'lowerleft', cmap = 'gray')
plt.scatter(y,x, c = 'red', s = 0.5)
plt.title('Scatter Plot Over Flare Image')
plt.xlabel('x [Arcsec]')
plt.ylabel('y [Arcsec]')
plt.xticks(np.arange(0,1001,200), xtickha)
plt.yticks(np.arange(0,1001,200), ytickha)
plt.savefig('ScatterFlare.png')
plt.show()


# In[ ]:


plt.ylim((0,1000))
plt.xlim((0,1000))
plt.scatter(x,y , c = 'blue')
plt.scatter(x1,y1, c = 'red')


# 
# # Finding Minimums

# In[ ]:


minx = S.findminx(data1)


# In[ ]:


miny = S.findminy(data1)


# In[ ]:


count = 4
x = []
y = []
while count <= 8:
    dat = data1.nb_data_dstr[count,500,500]
    x.append(count)
    y.append(dat)
    count += 1
var = np.polyfit(x,y,2)
a = var[0]
b = var[1]
d = var[2]
eq = lambda x: a*x**2 + b*x + d
mini = fmin(eq,6,disp=0)
print(mini,eq(mini))


# # Finding Prefilter Curve

# In[ ]:


Ffdata = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6563_flat.sav')
#'/workspace/mkatilius/IDL/mData/20141025/results_000/6563/nb/'
#'D:\\Coding\\Harris\\IDL\\DATA_Minda\\20141025\\calibration\\6563_flat.sav'


# In[ ]:


ndata = []
for i in range(13):
    ndata.append([])
    for a in range(1):
        ndata[i].append([])
        for b in range(1):
            means = np.mean(Ffdata.flat[i,450:550,450:550])
            ndata[i][a].append(means)
numd = []
for i in range(13):
    numd.append(np.array(ndata)[i,0,0])
numd = np.array(numd)


# In[ ]:


hadata = S.getData('atlas_ha.sav')


# In[ ]:


#print(hadata.p[1200:1700])
#print(hadata.xl[1200:1700])

numd[2:11]


# In[ ]:


ndatha = [hadata.p[1223],hadata.p[1270],hadata.p[1320],hadata.p[1370],
           hadata.p[1420],hadata.p[1470],hadata.p[1520],hadata.p[1570],hadata.p[1610]]
ndathax = [hadata.xl[1223],hadata.xl[1270],hadata.xl[1320],hadata.xl[1370],
           hadata.xl[1420],hadata.xl[1470],hadata.xl[1520],hadata.xl[1570],hadata.xl[1610]]
print(ndatha)
print(ndathax)


# In[ ]:


plt.plot(ndathax,ndatha)
plt.plot(ndathax,numd[2:11])


# In[ ]:


prefilterc = ndatha/numd[2:11]
plt.plot(np.arange(0,9),prefilterc)


# In[ ]:


havar = np.polyfit(x,ndatha,2)
datavar = np.polyfit(x,numd[2:11], 2)
eq1 = lambda x: havar[0]*x**2 + havar[1]*x + havar[2]
eq2 = lambda z: datavar[0]*z**2 + datavar[1]*z + datavar[2]


# In[ ]:


plt.plot(data2.nb_data_dstr[:,500,500])
anum = (data2.nb_data_dstr[0,500,500]-data2.nb_data_dstr[9,500,500])
slope = anum/9

eq3 = lambda x: -slope*x+data2.nb_data_dstr[0,500,500]
x = np.arange(0,10)
plt.plot(x,eq3(x))


# In[ ]:


newline = data2.nb_data_dstr[:,500,500]/eq3(x)
plt.plot(newline)


# # RMS Calculations from BroadBand

# ## H-Alpha

# In[ ]:


dataname = glob.glob('/workspace/mkatilius/IDL/mData/20141025/results_000/6563/bb/6563_bb*.sav')
databbha = []

for i in range(len(dataname)):
    datax = S.getData(dataname[i])
    databbha.append(datax)


# In[ ]:


fitsname = sorted(glob.glob('/workspace/mkatilius/IDL/mData/20141025/whitelight/ScienceObservation/20141025_160750/s*.ScienceObservation.fits'))
times = []
for i in range(len(fitsname)):
    hdul = fits.open(fitsname[i])
    ti = hdul[1].header['DATE-OBS']
    timcut = ti[11:]
    times.append(timcut)


# In[ ]:


timeticks = []
for i in range(50):
    timeticks.append(f'{(i*20.31)/60:.03f}')


# In[ ]:


rmsvalha = []
for i in range(len(dataname)):
    rmsha = S.rms(databbha[i], x1 = 250, x2 = 500, y1 = 100, y2 = 250)
    rmsvalha.append(rmsha)
rmsvalha2 = []
for i in range(len(dataname)):
    rmsha2 = S.brms(databbha[i], x1 = 250, x2 = 500, y1 = 100, y2 = 250)
    rmsvalha2.append(rmsha2)


# In[ ]:


rmsdataha = []
for i in range(len(rmsvalha)):
    for c in range(len(rmsvalha[1])):
        rmsdataha.append(rmsvalha[i][c])


# In[ ]:


maxrmsha = np.amax(rmsdataha)
rmsdataha1 = rmsdataha/maxrmsha
maxrmsha2 = np.amax(rmsvalha2)
rmsdataha2 = rmsvalha2/maxrmsha2


# In[ ]:


plt.figure(figsize=(200,10))
plt.xlim(0,len(rmsdataha1))
plt.plot(np.arange(0,len(rmsdataha1)),rmsdataha1)
for i in range(len(rmsdataha1)):
    if i % 14 == 0:
        plt.axvline(i,marker = ',', color='r')
plt.title('H-Alpha RMS plot', fontsize = 25)
plt.ylabel('Normalized RMS Value')
plt.xlabel('Time [Minutes]')
plt.xticks(np.arange(0,len(rmsdataha1),56),(timeticks))
plt.savefig('rmsplotHalpha.png')


# In[ ]:


plt.figure(figsize=(20,10))
plt.xlim(0,len(rmsdataha2))
plt.plot(np.arange(0,len(rmsdataha2)),rmsdataha2)
plt.title('Basic H-Alpha RMS plot')
plt.ylabel('RMS Value percentage')
plt.xlabel('File Number')
plt.xticks(np.arange(0,len(rmsdataha2),5))
plt.savefig('rms2plotHalpha.png')


# In[ ]:


xtickha = ['0',f'{200*.0946605:0.4f}',f'{400*.0946605:0.4f}',f'{600*.0946605:0.4f}',
         f'{800*.0946605:0.4f}',f'{1000*.0946605:0.4f}']
ytickha = ['0',f'{200*.0967937:0.3f}',f'{400*.0967937:0.3f}',f'{600*.0967937:0.3f}',
         f'{800*.0967937:0.3f}',f'{1000*.0967937:0.3f}']
fig,ax = plt.subplots(1,figsize = (10,10))
plt.setp(ax, xticks=[0, 200, 400, 600, 800, 999], xticklabels=ytickha ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickha)
ax.set_title('Location of H-Alpha RMS value Calculations')
ax.set_ylabel('y [arcsec]')
ax.set_xlabel('x [arcsec]')
# Display the image
ax.imshow(databbha[0].bb_data_dstr[0,:,:].T,origin = 'lowerleft', cmap = 'gray',)

# Create a Rectangle patch
rect = patches.Rectangle((100,250),150,250,linewidth=1,edgecolor='r',facecolor='none')
    
# Add the patch to the Axes
ax.add_patch(rect)
plt.savefig('RMS_Value_Calc_Image_Halpha.png')


# In[ ]:


plt.imshow(databbha[0].bb_data_dstr[0,100:250,250:500],origin = 'lowerleft', cmap = 'gray',)


# ## Iron I

# In[ ]:


dataname2 = glob.glob('/workspace/mkatilius/IDL/mData/20141025/results_000/6173/bb/6173_bb*.sav')
databbir = []

for i in range(len(dataname2)):
    datax2 = S.getData(dataname2[i])
    databbir.append(datax2)


# In[ ]:


rmsvalir = []
for i in range(len(dataname2)):
    rmsir = S.rms(databbir[i], x1 = 250, x2 = 500, y1 = 100, y2 = 250)
    rmsvalir.append(rmsir)
rmsvalir2 = []
for i in range(len(dataname2)):
    rmsir2 = S.brms(databbir[i], x1 = 250, x2 = 500, y1 = 100, y2 = 250)
    rmsvalir2.append(rmsir2)


# In[ ]:


rmsdatair = []
for i in range(len(rmsvalir)):
    for c in range(len(rmsvalir[1])):
        rmsdatair.append(rmsvalir[i][c])


# In[ ]:


maxrmsir = np.amax(rmsdatair)
rmsdatair1 = rmsdatair/maxrmsir
maxrmsir2 = np.amax(rmsvalir2)
rmsdatair2 = rmsvalir2/maxrmsir2


# In[ ]:


plt.figure(figsize=(30,10))
plt.xlim(0,len(rmsdatair1))
plt.plot(np.arange(0,len(rmsdatair1)),rmsdatair1)
for i in range(len(rmsdatair1)):
    if i % 10 == 0:
        plt.axvline(i,marker = ',', color='r')
plt.title('Iron I RMS plot')
plt.ylabel('Normalized RMS Value')
plt.xlabel('Time [Minutes]')
plt.xticks(np.arange(0,len(rmsdatair1),40),timeticks)
plt.savefig('rmsplotIronI.png')


# In[ ]:


plt.figure(figsize=(20,10))
plt.xlim(0,len(rmsdatair2))
plt.plot(np.arange(0,len(rmsdatair2)),rmsdatair2)
plt.title('Basic Iron I RMS plot')
plt.ylabel('RMS Value percentage')
plt.xlabel('File Number')
plt.xticks(np.arange(0,len(rmsdatair2),5))
plt.savefig('rms2plotIronI.png')


# In[ ]:


xtickir = ['0',f'{200*.0946345:0.4f}',f'{400*.0946345:0.4f}',f'{600*.0946345:0.4f}',
         f'{800*.0946345:0.4f}',f'{1000*.0946345:0.4f}']
ytickir = ['0',f'{200*.0968169:0.3f}',f'{400*.0968169:0.3f}',f'{600*.0968169:0.3f}',
         f'{800*.0968169:0.3f}',f'{1000*.0968169:0.3f}']
fig,ax = plt.subplots(1,figsize = (10,10))
plt.setp(ax, xticks=[0, 200, 400, 600, 800, 999], xticklabels=ytickir ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickir)
ax.set_title('Location of Iron I RMS value Calculations')
ax.set_ylabel('y [arcsec]')
ax.set_xlabel('x [arcsec]')
# Display the image
ax.imshow(databbir[0].bb_data_dstr[0,:,:].T,origin = 'lowerleft', cmap = 'gray',)

# Create a Rectangle patch
rect = patches.Rectangle((100,250),150,250,linewidth=1,edgecolor='r',facecolor='none')
    
# Add the patch to the Axes
ax.add_patch(rect)
plt.savefig('RMS_Value_Calc_Image_IronI.png')


# # Testing Spline Coding

# In[ ]:


for i in range(len(data1.nb_data_dstr[0,500,:])):
    count = 2
    x = []
    y = []
    while count <= 10:
        dat = data1.nb_data_dstr[count,350,i]
        x.append(count)
        y.append(dat)
        count += 1
    var = UnivariateSpline(x,y,k=4)
    mini = fmin(var,6,disp = 0)
    mid = (var(2)+var(10))/2
    np.mean((mid,var(mini)))
    if i % 10:
        xnew = np.arange(2,10.1,.1)
        plt.plot(data1.nb_data_dstr[:,350,i])
        plt.plot(xnew,var(xnew))
        plt.show()


# In[ ]:


xnew = np.arange(2,11,.1)
plt.plot(data1.nb_data_dstr[:,350,350])
plt.plot(xnew,var(xnew))


# In[ ]:


count = 4
x = []
y = []
while count <= 8:
    dat = data1.nb_data_dstr[count,500,500]
    x.append(count)
    y.append(dat)
    count += 1
var = np.polyfit(x,y,2)
a = var[0]
b = var[1]
d = var[2]
eq = lambda x: a*x**2 + b*x + d 
mini = fmin(eq,6,disp=0)
eq(mini)


# In[ ]:


plt.plot(data1.nb_data_dstr[:,500,500])


# # Image Alignment and Other Things
# 

# In[ ]:


plt.figure(figsize = (15,15))
plt.imshow(data1.nb_data_dstr[6,:,:], cmap = 'gray', origin = 'lowerleft')
plt.show()


# In[ ]:


hdul = fits.open('/workspace/mkatilius/IDL/mData/20141025/scans/ScienceObservation/20141025_160750/s000.ScienceObservation.fits')


# In[ ]:


predata = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6563_prefilter.sav')


# In[ ]:


newdata = Ffdata.flat[:,200,200]/predata.prefilter_curve
plt.plot(data1.info_nb.WAVE[0],newdata)
plt.plot(hadata.xl,hadata.p)


# In[ ]:


x = []
for i in range(14):
    x.append(np.mean(data1.nb_data_dstr[i,450:550,450:550]))


# In[ ]:


plt.figure(figsize= (10,10))
ndata = data1.nb_data_dstr[:,350,350]/predata.prefilter_curve
plt.plot(data1.info_nb.WAVE[0],ndata*2.55)
plt.plot(hadata.xl,hadata.p)


# In[ ]:


data1.nb_data_dstr[0,10,:]


# In[ ]:


x = []
for i in range(len(data1.nb_data_dstr[0,0,:])):
    x.append([])
    for a in range(len(data1.nb_data_dstr[0,0,:])):
        if data1.nb_data_dstr[0,i,a] == data1.nb_data_dstr[0,0,0]:
            x[i].append(0)
            continue
        else:
            x[i].append(data1.nb_data_dstr[0,i,a])


# In[ ]:


plt.imshow(x, origin = 'lowerleft', cmap = 'gray')


# # Prefilter Graphs

# In[ ]:


hflats = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6563_flat.sav')
irflats = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6173_flat.sav')
hpref = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6563_prefilter.sav')
irpref = S.getData('/workspace/mkatilius/IDL/mData/20141025/calibration/6173_prefilter.sav')


# In[ ]:


avehflat = []
for i in range(14):
    means = np.mean(hflats.flat[i,450:550,450:550])
    avehflat.append(means)
avehflat = np.array(avehflat)
aveirflat = []
for i in range(10):
    means = np.mean(irflats.flat[i,450:550,450:550])
    aveirflat.append(means)
aveirflat = np.array(aveirflat)


# In[ ]:


maxh = max(avehflat)
iron1pref = aveirflat/irpref.prefilter_curve
halphapref = avehflat/hpref.prefilter_curve
iron1n = (iron1pref[0] + iron1pref[9])/2
halphan = (halphapref[2] + halphapref[10])/2


# In[ ]:


plt.plot(data1.info_nb['WAVE'][0]/10,avehflat/halphan, label = 'No Correction')
plt.plot(data1.info_nb['WAVE'][0]/10,halphapref/halphan, label = 'Corrected')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalized Light Intensity')
plt.title('H-Alpha Corrections')
plt.legend()
plt.savefig('Halpha Corrections')
plt.show()

plt.plot(data2.info_nb['WAVE'][0]/10,aveirflat/iron1n, label = 'No Correction')
plt.plot(data2.info_nb['WAVE'][0]/10,iron1pref/iron1n, label = 'Corrected')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalized Light Intensity')
plt.title('Iron I Corrections')
plt.legend()
plt.savefig('Iron I Corrections.jpg')
plt.show()


# # Pixel Walking Graphs

# In[ ]:


dataname = glob.glob('/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6563_nbn/6563_nbn*.sav')
datanbnha = []

for i in range(len(dataname)):
    datax = S.getData(dataname[i])
    datanbnha.append(datax)


# In[ ]:


plt.figure(figsize = (10,10))
plt.imshow(datanbnha[89].nb_data_crrt[6,:,:].T, origin = 'lowerleft', cmap = 'gray')
xtickha = ['0',f'{200*.0954204:0.3f}',f'{400*.0954204:0.3f}',f'{600*.0954204:0.3f}',
         f'{800*.0954204:0.3f}',f'{1000*.0954204:0.3f}']
ytickha = ['0',f'{200*.0974302:0.3f}',f'{400*.0974302:0.3f}',f'{600*.0974302:0.3f}',
         f'{800*.0974302:0.3f}',f'{1000*.0974302:0.3f}']


# In[ ]:


xtickha = ['0',f'{200*.0954204:0.3f}',f'{400*.0954204:0.3f}',f'{600*.0954204:0.3f}',
         f'{800*.0954204:0.3f}',f'{1000*.0954204:0.3f}']
ytickha = ['0',f'{200*.0974302:0.3f}',f'{400*.0974302:0.3f}',f'{600*.0974302:0.3f}',
         f'{800*.0974302:0.3f}',f'{1000*.0974302:0.3f}']
fig,ax = plt.subplots(1,figsize = (10,10))
ax.set_title('Location of Averaged Light Intensity')
plt.setp(ax, xticks=[0, 200, 400, 600, 800, 999], xticklabels=ytickha ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickha)
ax.set_ylabel('y [arcsec]')
ax.set_xlabel('x [arcsec]')
# Display the image
ax.imshow(datanbnha[89].nb_data_crrt[6,:,:].T, origin = 'lowerleft', cmap = 'gray')
ax.scatter(589,260, c = 'lightgreen', s = 75, marker = 'x')

# Create a Rectangle patch
rect1 = patches.Rectangle((675,400),50,50,linewidth=1,edgecolor='b',facecolor='none')
rect2 = patches.Rectangle((400,275),30,40,linewidth=1,edgecolor='orange',facecolor='none')

# Add the patch to the Axes
ax.add_patch(rect1)
ax.add_patch(rect2)
fig.savefig('FlarelightInt.png')


# In[ ]:


meanf = []
meanq = []
for i in range(14):
    avef = np.mean(datanbnha[89].nb_data_crrt[i,400:430,275:315])
    aveq = np.mean(datanbnha[89].nb_data_crrt[i,675:725,400:450])
    meanf.append(avef)
    meanq.append(aveq)


# In[ ]:


nmal = (meanq[2]+meanq[10])/2
wave = datanbnha[89].info_nb['WAVE'][0]


# In[ ]:


f, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, sharey = True, figsize = (20,8))
ax1.plot(wave/10, meanq/nmal)
ax2.plot(wave/10, meanf/nmal)


# In[ ]:


plt.figure(figsize = (10,8))
plt.plot(wave/10, meanq/nmal, label = 'Outside Flare')
plt.plot(wave/10, meanf/nmal, label = 'Inside Flare')
plt.title('Averaged Light Intesity Inside and Outside Flare')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Light Intensity')
plt.savefig('AveLightFlare')
plt.show()


# In[ ]:


nnal1 = (datanbnha[0].nb_data_crrt[2,589,260] + datanbnha[0].nb_data_crrt[10,589,260])/2
for i in range(len(datanbnha)):
    plt.plot(wave/10,datanbnha[i].nb_data_crrt[:,589,260]/nnal1)
    plt.ylim((0.6,1.3))
    plt.title(f'File Number {i}')
    plt.show()


# In[ ]:


hdul = fits.open('/workspace/mkatilius/IDL/mData/20141025/scans/ScienceObservation/20141025_160750/s060.ScienceObservation.fits')
hdul.info()
plt.imshow(hdul[1].data)


# In[ ]:


plt.figure(figsize = (10,8))
plt.plot(wave/10,datanbnha[1].nb_data_crrt[:,589,260]/nnal1, label = '16:08 UTC')
plt.plot(wave/10,datanbnha[60].nb_data_crrt[:,589,260]/nnal1, label = '16:13 UTC')
plt.plot(wave/10,datanbnha[135].nb_data_crrt[:,589,260]/nnal1, label = '16:19 UTC')
plt.plot(wave/10,datanbnha[160].nb_data_crrt[:,589,260]/nnal1, label = '16:21 UTC')
plt.ylabel('Normalize Light Intensity')
plt.xlabel('Wavelength [nm]')
plt.title('Change in Spectrum with Pixel by Flare')
plt.legend()
plt.savefig('Flarechangespec')
plt.show()


# # Bisector Graphs

# In[ ]:


bisec = S.getData('/workspace/mkatilius/IDL/mData/20141025/bisec/6563_bisec_089.sav')


# In[ ]:


av56 = (data1.info_nb['WAVE'][0][6] - data1.info_nb['WAVE'][0][5])/10
av67 = (data1.info_nb['WAVE'][0][7] - data1.info_nb['WAVE'][0][6])/10
nu59 = (data1.info_nb['WAVE'][0][6]-av56)
nu60 = data1.info_nb['WAVE'][0][6]
nu61 = (data1.info_nb['WAVE'][0][6] + av67)
nu62 = (data1.info_nb['WAVE'][0][6] + 2*av67)
nu63 = (data1.info_nb['WAVE'][0][6] + 3*av67)
nu64 = (data1.info_nb['WAVE'][0][6] + 4*av67)
nu65 = (data1.info_nb['WAVE'][0][6] + 5*av67)
cbar01w = (["{0:.3f}".format(nu59),"{0:.3f}".format(nu60), "{0:.3f}".format(nu61), "{0:.3f}".format(nu62),
           "{0:.3f}".format(nu63), "{0:.3f}".format(nu64), "{0:.3f}".format(nu65)])
cbar10w = (["{0:.3f}".format(nu61),"{0:.3f}".format(nu62),"{0:.3f}".format(nu63),"{0:.3f}".format(nu64),
           "{0:.3f}".format(nu65)])
print(cbar01w)
print(cbar10w)


# In[ ]:


intc1, intc2 = np.percentile(bisec.int_core,(5,95))
posc1, posc2 = np.percentile(bisec.pos_core,(5,95))
posf1, posf2 = np.percentile(bisec.pos_fft, (5,95))
fwhm1, fwhm2 = np.percentile(bisec.fwhm, (5,95))
fig, axes = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True, figsize = (10,10))

plt.setp(axes, xticks=[0, 200, 400, 600, 800, 999], xticklabels=ytickha ,
         yticks=[0, 200, 400, 600, 800, 999], yticklabels = xtickha)

im00 = axes[0,0].imshow(bisec.int_core.T, origin = 'lowerleft', cmap = 'gray', vmin = intc1)
im01 = axes[0,1].imshow(bisec.pos_core.T, origin = 'lowerleft', cmap = 'gray', vmin = posc1, vmax = posc2)
im10 = axes[1,0].imshow(bisec.pos_fft.T, origin = 'lowerleft', cmap = 'gray', vmin = posf1, vmax = posf2)
im11 = axes[1,1].imshow(bisec.fwhm.T, origin = 'lowerleft', cmap = 'gray', vmin = fwhm1, vmax = fwhm2)

cbar00 = fig.colorbar(im00, ax = axes[0,0])
cbar01 = fig.colorbar(im01, ax = axes[0,1])
cbar10 = fig.colorbar(im10, ax = axes[1,0])
cbar11 = fig.colorbar(im11, ax = axes[1,1])
cbar01.ax.set_yticklabels(cbar01w)
cbar10.ax.set_yticklabels(cbar10w)
axes[0,0].set_title('Light Intensity of Core')
axes[0,1].set_title('Position of Line Core')
axes[1,0].set_title('Interpolated Position of Line Core')
axes[1,1].set_title('Full Width at Half Max')

axes[1,0].set_xlabel('x [arcsec]')
axes[1,1].set_xlabel('x [arcsec]')
axes[0,0].set_ylabel('y [arcsec]')
axes[1,0].set_ylabel('y [arcsec]')

fig.savefig('BisectorStuff.png')


# # Flare Movement

# In[ ]:


dataname = glob.glob('/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6563_nbn/6563_nbn*.sav')
datanbnha = []

for i in range(len(dataname)):
    datax = S.getData(dataname[i])
    datanbnha.append(datax)


# In[ ]:


xs = []
ys = []

for i in range(len(datanbnha)):
    datax, datay = S.cogfitmk2(datanbnha[i], num = 6, lim = 2200)
    xs.append(datax)
    ys.append(datay)


# In[ ]:


plt.figure(figsize = (10,10))
for i in range(len(datanbnha)):
    plt.scatter(xs[i],ys[i], s = 0.5 )
    


# # Sun Images

# In[ ]:


flats = glob.glob('/home/mkatilius/Python/Solar_jpg/Sun_Images/20141025_*_4096_HMIIF.jpg')
magneto = glob.glob('/home/mkatilius/Python/Solar_jpg/Magnetograms/20141025_*_1024_HMIB.jpg')
aia131 = glob.glob('/home/mkatilius/Python/Solar_jpg/AIA_131/20141025_*_2048_0131.jpg')
aia1700 = glob.glob('/home/mkatilius/Python/Solar_jpg/AIA_1700/20141025_*_2048_1700.jpg')


# In[ ]:


newflat = flats[52:77:4]
imagef = []
for i in range(len(newflat)):
    imagef.append(plt.imread(newflat[i]))
    
newmag = magneto[52:77:4]
imagem = []
for i in range(len(newmag)):
    imagem.append(plt.imread(newmag[i]))
    
new131 = aia131[52:77:4]
image13 = []
for i in range(len(new131)):
    image13.append(plt.imread(new131[i]))
    
new1700 = aia1700[52:77:4]
image17 = []
for i in range(len(new1700)):
    image17.append(plt.imread(new1700[i]))


# In[ ]:


fig, axes = plt.subplots(nrows = 1, ncols = 6, figsize = (50,6))
fig.subplots_adjust(wspace=0.1, hspace=0)
fig.suptitle('Zoomed in Flatfields', fontsize = 25)
axes[0].imshow(imagef[0][2400:2900,2500:3200,:])
axes[1].imshow(imagef[1][2400:2900,2500:3200,:])
axes[2].imshow(imagef[2][2400:2900,2500:3200,:])
axes[3].imshow(imagef[3][2400:2900,2500:3200,:])
axes[4].imshow(imagef[4][2400:2900,2500:3200,:])
axes[5].imshow(imagef[5][2400:2900,2500:3200,:])
#axes[6].imshow(imagef[6])
axes[0].axis('off')
axes[1].axis('off')
axes[2].axis('off')
axes[3].axis('off')
axes[4].axis('off')
axes[5].axis('off')
#axes[6].axis('off')
axes[0].set_title('13:00 UTC', fontsize = 14)
axes[1].set_title('14:00 UTC', fontsize = 14)
axes[2].set_title('15:00 UTC', fontsize = 14)
axes[3].set_title('16:00 UTC', fontsize = 14)
axes[4].set_title('17:00 UTC', fontsize = 14)
axes[5].set_title('18:00 UTC', fontsize = 14)


# In[ ]:


fig2, axes2 = plt.subplots(nrows = 1, ncols = 6, figsize = (40,7))
fig2.subplots_adjust(wspace=0.1, hspace=0)
fig2.suptitle('Zoomed in Magnetograms', fontsize = 25)
axes2[0].imshow(imagem[0][550:775,600:825], cmap = 'gray')
axes2[1].imshow(imagem[1][550:775,600:825], cmap = 'gray')
axes2[2].imshow(imagem[2][550:775,600:825], cmap = 'gray')
axes2[3].imshow(imagem[3][550:775,600:825], cmap = 'gray')
axes2[4].imshow(imagem[4][550:775,600:825], cmap = 'gray')
axes2[5].imshow(imagem[5][550:775,600:825], cmap = 'gray')
#axes[6].imshow(imagem[6])
axes2[0].axis('off')
axes2[1].axis('off')
axes2[2].axis('off')
axes2[3].axis('off')
axes2[4].axis('off')
axes2[5].axis('off')
#axes[6].axis('off')
axes2[0].set_title('13:00 UTC', fontsize = 14)
axes2[1].set_title('14:00 UTC', fontsize = 14)
axes2[2].set_title('15:00 UTC', fontsize = 14)
axes2[3].set_title('16:00 UTC', fontsize = 14)
axes2[4].set_title('17:00 UTC', fontsize = 14)
axes2[5].set_title('18:00 UTC', fontsize = 14)


# In[ ]:


fig2, axes2 = plt.subplots(nrows = 1, ncols = 6, figsize = (40,7))
fig2.subplots_adjust(wspace=0.1, hspace=0)
fig2.suptitle('Zoomed in 131 Angstroms', fontsize = 25)
axes2[0].imshow(image13[0][1000:1550,1125:1600])
axes2[1].imshow(image13[1][1000:1550,1125:1600])
axes2[2].imshow(image13[2][1000:1550,1125:1600])
axes2[3].imshow(image13[3][1000:1550,1125:1600])
axes2[4].imshow(image13[4][1000:1550,1125:1600])
axes2[5].imshow(image13[5][1000:1550,1125:1600])
#axes[6].imshow(imagem[6])
axes2[0].axis('off')
axes2[1].axis('off')
axes2[2].axis('off')
axes2[3].axis('off')
axes2[4].axis('off')
axes2[5].axis('off')
#axes[6].axis('off')
axes2[0].set_title('13:12 UTC', fontsize = 14)
axes2[1].set_title('14:12 UTC', fontsize = 14)
axes2[2].set_title('15:12 UTC', fontsize = 14)
axes2[3].set_title('16:12 UTC', fontsize = 14)
axes2[4].set_title('17:12 UTC', fontsize = 14)
axes2[5].set_title('18:12 UTC', fontsize = 14)


# In[ ]:


fig3, axes3 = plt.subplots(nrows = 1, ncols = 6, figsize = (40,7))
fig3.subplots_adjust(wspace=0.1, hspace=0)
fig3.suptitle('Zoomed in 1700 Angstroms', fontsize = 25)
axes3[0].imshow(image17[0][1000:1550,1125:1600])
axes3[1].imshow(image17[1][1000:1550,1125:1600])
axes3[2].imshow(image17[2][1000:1550,1125:1600])
axes3[3].imshow(image17[3][1000:1550,1125:1600])
axes3[4].imshow(image17[4][1000:1550,1125:1600])
axes3[5].imshow(image17[5][1000:1550,1125:1600])
#axes[6].imshow(imagem[6])
axes3[0].axis('off')
axes3[1].axis('off')
axes3[2].axis('off')
axes3[3].axis('off')
axes3[4].axis('off')
axes3[5].axis('off')
#axes[6].axis('off')
axes3[0].set_title('13:03 UTC', fontsize = 14)
axes3[1].set_title('14:03 UTC', fontsize = 14)
axes3[2].set_title('15:03 UTC', fontsize = 14)
axes3[3].set_title('16:03 UTC', fontsize = 14)
axes3[4].set_title('17:03 UTC', fontsize = 14)
axes3[5].set_title('18:03 UTC', fontsize = 14)


# # Changing Flare

# In[ ]:


dataname = glob.glob('/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6563_nbn/6563_nbn*.sav')
datanbha = []

for i in range(len(dataname)):
    datax = S.getData(dataname[i])
    datanbha.append(datax)


# In[ ]:


plt.imshow(datanbha[173].nb_data_crrt[6,:,:].T,origin = 'lowerleft',cmap = 'gray')


# In[ ]:


waveha = data1.info_nb['WAVE'][0]
datarr = []
for a in range(len(datanbha)):
    datarr.append([])
    for b in range(len(datanbha[0].nb_data_crrt[6,:,:])):
        datarr[a].append(datanbha[a].nb_data_crrt[6,400,b])


# In[ ]:


yticklabel = [times[0][0:8],times[50][0:8],times[100][0:8],times[150][0:8]]
plt.figure(figsize = (10,10))
plt.imshow(np.array(datarr), origin = 'lowerleft', cmap = 'gray')
plt.title(f'Solar Flare Width in Time at {waveha[6]:.02f} nm')
plt.xlabel('x [arcsec]')
plt.xticks(np.arange(0,1001,200), xtickha)
plt.ylabel('Time in UTC')
plt.yticks(np.arange(0,200,50),yticklabel)
plt.savefig('Flaretime.png')


# In[ ]:


yticklabel


# In[ ]:


wave = datanbha[0].info_nb['WAVE'][0][6]/10


# # Contour Plots

# In[ ]:


bisec = S.getData('/workspace/mkatilius/IDL/mData/20141025/bisec/6563_bisec_089.sav')
xtickha = ['0',f'{200*.0954204:0.3f}',f'{400*.0954204:0.3f}',f'{600*.0954204:0.3f}',
         f'{800*.0954204:0.3f}',f'{1000*.0954204:0.3f}']
ytickha = ['0',f'{200*.0974302:0.3f}',f'{400*.0974302:0.3f}',f'{600*.0974302:0.3f}',
         f'{800*.0974302:0.3f}',f'{1000*.0974302:0.3f}']


# In[ ]:


levelsc = [0,400,800,1200,1600,2000,2100]
intc1, intc2 = np.percentile(bisec.int_core,(5,95))
plt.figure(figsize = (20,20))
plt.imshow(bisec.int_core.T, cmap = 'gray', origin = 'lower', vmin = intc1)
cp = plt.contour(bisec.int_core.T, levels = levelsc)

plt.title('Contour plot of the Filaments for the Intensity in the Core of H-Alpha')
plt.xticks(np.arange(0,1001,200), xtickha)
plt.yticks(np.arange(0,1001,200), ytickha)
plt.xlabel('x [arcsec]')
plt.ylabel('y [arcsec]')

norm= matplotlib.colors.Normalize(vmin = cp.cvalues.min(), vmax = cp.cvalues.max())
sm = plt.cm.ScalarMappable(norm=norm, cmap = cp.cmap)
sm.set_array([])
clb = plt.colorbar(sm, ticks=cp.levels)
clb.ax.set_title('Light Intensities')

plt.savefig('Filament_Contour.png')


# # Changing Flare
# 

# In[ ]:


dataname = glob.glob('D:\\Coding\\Harris\\IDL\\DATA_Minda\\20141025\\prefiltcalib\\6563_nbn\\6563_nbn*.sav')
datanbha = []

for i in range(len(dataname)):
    datax = S.getData(dataname[i])
    datanbha.append(datax)


# In[ ]:


plt.imshow(datanbha[173].nb_data_crrt[6,:,:],origin = 'lowerleft')


# In[ ]:


datarr = []
for a in range(len(datanbha)):
    datarr.append([])
    for b in range(len(datanbha[0].nb_data_crrt[6,:,:])):
        datarr[a].append(datanbha[a].nb_data_crrt[6,400,b])


# In[ ]:


plt.figure(figsize = (10,10))
plt.imshow(np.array(datarr), origin = 'lowerleft', cmap = 'gray')
plt.title('Solar Flare Width in Time')
plt.savefig('Flaretime')

