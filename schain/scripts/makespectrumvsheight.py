#READING FILE
import numpy as np
#import matplotlib.pyplot as plt
import os
from scipy import signal
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Almacenar espectros y tiempos
spectra_time = []
utc_time = []

#Write directory of .r files here
#dirr='/home/roberto/DATA/HAE_2009_12/data/'
#dirr='/media/roberto/DISCO_13/PUMA oblique/2009/2009_02/Hybrid_LP_1064NSA/'
#dirr='/home/roberto/puma/2009_12/d2009348/'
#dirr='/home/roberto/Minotaur/2024_01/MP_ISR/main_radar/rawdata/d2024027/'
dirr = '/media/cportilla/HDD/DAVID-OP/data-sevag/7/d2021301/'


#Reading

SAVE = 1
PROC = 1

if SAVE:
    outputfile = open("SpectrumDec2024027.dat","wb")
    outputfile_utc = open("TimeStampsDec2024027.dat","wb")



nprofiles= 128 #200
nsamp=800 #980
nreceivers=2
osamp=4
navg=5 #nb block averages
filepath=os.listdir(dirr) #Get the list of all files in directory
filepath=sorted(filepath) #Sort the files in ascending numerical order
nlag=16
nrange=176
hsize=244 #LP = 344

lag16=[]
x=0
filecount=0
iavg=0 #27-40

fileread=True
dirr_ = 'D2021301'

for file1 in filepath:
 filecount=filecount+1
 #if file1=='D2009345004.r': fileread=True
 if fileread==True and file1[-1]=='r' and file1[:8]==dirr_:
    print("file1")
    print(file1)
    inputfile = open(dirr+file1,"rb")
    print('inputfile',inputfile)
    fsize=os.stat(dirr+file1).st_size
    # 160321468
    # print('fsize',fsize)
    nblocks=int(((fsize)-hsize)/(nprofiles*nsamp*nreceivers*4))

    #from the manual nblocks per file
    nblocks=50

    #print('nblocks', nblocks)
    for k in range(nblocks):#nblocks): #only read one block, acording time step
        iavg=iavg+1
        if k==0:
            #hdr=hsize# 244 #492 #First header???
            nSize = np.fromfile(inputfile,np.int32,1,"") # int32 4 bytes
            #print(nSize)
            version = np.fromfile(inputfile,np.int16,1,"") # int16 2 bytes version
            dataBlock = np.fromfile(inputfile,np.int32,1,"")
            utc = np.fromfile(inputfile,np.int32,1,"")
            miliSecond = np.fromfile(inputfile,np.int16,1,"")
            timeZone = np.fromfile(inputfile,np.int16,1,"")
            dstFlag = np.fromfile(inputfile,np.int16,1,"")
            errorCount = np.fromfile(inputfile,np.int32,1,"")

            '''tempheader=np.fromfile(inputfile,np.int32,int(hdr/4)-6,"") #Skip header
            print("tempheader: ", tempheader)'''
            #print("utc: ",datetime.datetime.utcfromtimestamp(utc))

            #first skip of long header
            size1LongHeader = np.fromfile(inputfile,np.int32,1,"")

            nSamples = np.fromfile(inputfile,np.int32,1,"")[0]
            #print(nSamples)
            nProfiles = np.fromfile(inputfile,np.int32,1,"")[0]
            #print(nProfiles)
            nChannels = np.fromfile(inputfile,np.int32,1,"")[0]
            #print(nChannels)

            #print("size1LongHeader", size1LongHeader[0])
            tmpSkip1 = np.fromfile(inputfile,np.int32,int(size1LongHeader[0]/4)-3-1,"") #4 for 3 parameters

            #second skip of long header
            size2LongHeader = np.fromfile(inputfile,np.int32,1,"")
            #print("size2LongHeader", size2LongHeader[0])
            tmpSkip2 = np.fromfile(inputfile,np.int32,int(size2LongHeader[0]/4)-1,"")

            #n1 = np.fromfile(inputfile,np.int32,1,"")
            

            #third skip of long header
            size3LongHeader = np.fromfile(inputfile,np.int32,1,"")
            #print("size3LongHeader", size3LongHeader[0])
            n2 = np.fromfile(inputfile,np.int32,1,"")
            sizeBlock = np.fromfile(inputfile,np.int32,1,"")
            #print(n2)
            #print('nDataBlock', sizeBlock )

            tmpSkip3 = np.fromfile(inputfile,np.int32,int(size3LongHeader[0]/4)-2-1,"")
        
        

        
        else:
            #print("hola else hdr")
            #hdr=24  #Short header
            nSize = np.fromfile(inputfile,np.int32,1,"")
            version = np.fromfile(inputfile,np.int16,1,"")
            dataBlock = np.fromfile(inputfile,np.int32,1,"")
            utc = np.fromfile(inputfile,np.int32,1,"")
            miliSecond = np.fromfile(inputfile,np.int16,1,"")
            timeZone = np.fromfile(inputfile,np.int16,1,"")
            dstFlag = np.fromfile(inputfile,np.int16,1,"")
            errorCount = np.fromfile(inputfile,np.int32,1,"")

            print("utc: ",datetime.datetime.utcfromtimestamp(utc[0]))
        #exit(1)
         #print("utc: ",datetime.datetime.utcfromtimestamp(utc))
         #print("UTC: ", datetime.datetime.utcfromtimestamp(UTC))
        #tempheader=np.fromfile(inputfile,np.int32,int(hdr/4),"") #Skip header
        #print("tempheader: ", np.fromfile(inputfile,np.int32,int(hdr/4),""))
        #print("UTC: ", datetime.datetime.utcfromtimestamp(self.utc))
        '''
        print("nSize: ", nSize)
        print("version: ", version)
        print("dataBlock: ", dataBlock)
        print("utc: ", utc)
        print("miliSecond: ", miliSecond)
        print("timeZone: ", timeZone)
        print("dstFlag: ", dstFlag)
        print("errorCount: ", errorCount)
        '''

        
        

        

        #blockdata=np.fromfile(inputfile,np.int16,2*nsamp*nprofiles*nreceivers,"")
        blockdata=np.fromfile(inputfile,np.int16, 2*nSamples*nProfiles*nChannels,"") # *2 components I/Q and int16 for complete entire block
        #print(blockdata.shape)
        #blockdata = np.fromfile(inputfile,np.int16,int(sizeBlock[0]/2),"")
        #print('blockdata.shape', blockdata.shape)

        r=blockdata[0::2]#real
        #print('r.shape',r.shape)
        #print(r)
        #print('r.shape',r.shape)
        i=blockdata[1::2]#imag
        #print('i.shape',i.shape)
        complex_samp=r+1j*i #complex
        complex_samp1=complex_samp[0::nChannels]#nreceivers] #Channel 0
        complex_samp2=complex_samp[1::nChannels] #nreceivers] #Channel 1

        #print(complex_samp1.shape)
        #print(complex_samp2.shape)

        #tempsamp1=np.reshape(complex_samp1,(nprofiles,nsamp))
        #tempsamp2=np.reshape(complex_samp2,(nprofiles,nsamp))

        #ch0
        tempsamp1=np.reshape(complex_samp1,(nProfiles,nSamples))
        #ch1
        tempsamp2=np.reshape(complex_samp2,(nProfiles,nSamples))

        #print('tempsamp1',tempsamp1)
        #print('tempsamp1.shape',tempsamp1.shape)
        #print('tempsamp2.shape',tempsamp2.shape)

        #select profiles and heights, review of the manual
        #tempsamp1 = tempsamp1[0:128,0:800]
        #tempsamp2 = tempsamp2[0:128,0:800]

        print('tempsamp1[profile 1,last]',tempsamp1[1,:])

        print('tempsamp1.shape',tempsamp1.shape)
        print('tempsamp2.shape',tempsamp2.shape)


        #FFT
        n=40
        #added
        nFFTs = 200
        #freqs=np.fft.fftfreq(n,d=25e-6) #*2 is the result of zeropadding
        freqs=np.fft.fftfreq(nFFTs,d=25e-6)
        freqs=np.fft.fftshift(freqs)

        #Switch antenna
        #tempsamp1=tempsamp2
        if PROC:
        #'''
            for jj in range(500,3500,10):
            #for jj in range(0,3000,10):
              l=int(jj/3.75)
              r=int(jj/3.75)+n

              #Process Spectra (from antenna #1)
              #tempsamp1a=tempsamp1[:128,l:r]
              
              #tempsamp1a=tempsamp1[:128,l:r]
              tempsamp1a=tempsamp2[:128,l:r]
              
              # Zero padding: 40 muestras + 160 ceros
              padding = np.zeros((tempsamp1a.shape[0], nFFTs - n),dtype=tempsamp1a.dtype)
              tempsamp1a_padded = np.hstack((tempsamp1a, padding))

              # FFT de 200 puntos
              spectr1 = np.fft.fft(tempsamp1a_padded, axis=1)

            


              
              #spectr1=np.apply_along_axis(np.fft.fft,1,np.hstack((tempsamp1a,np.zeros((tempsamp1a.shape[0],tempsamp1a.shape[1])) ))).astype(np.complex)
              #spectr2=np.apply_along_axis(np.fft.fft,1,tempsamp1a).astype(np.complex)
              #spectr1=np.fft.fftshift(spectr1)[:,::2] #adding axis=1
              spectr1=np.fft.fftshift(spectr1,axes=1)
              spectr1=np.abs(spectr1)**2
              spectr1=np.sum(spectr1,axis=0)
              spectr1sliced=spectr1 #Only the relevant are selected
              if jj==500:
                spectrslicedallranges=spectr1sliced
              else:
                spectrslicedallranges=np.vstack((spectrslicedallranges,spectr1sliced))

                #print('spectr1sliced',spectr1sliced.shape)

            #Incoherent integration over blocks
            if iavg==1:
               spectrslicedtot= spectrslicedallranges
            else: #Add a third dimension (depth->dstack)
                spectrslicedtot=np.dstack((spectrslicedtot, spectrslicedallranges))
            print(iavg,file1)
            if iavg==navg: #We have reached the number of blocks required
                if navg>1:
                    spectrslicedtotnavg=np.sum(spectrslicedtot,axis=2)#Sum along depth which corresponds to blocks
                else:
                    spectrslicedtotnavg=spectrslicedtot
                if SAVE:
                    spectrslicedtotnavg.astype(np.float32).tofile(outputfile)
            
                    utc.astype(np.int32).tofile(outputfile_utc)
                print(spectrslicedtotnavg.shape)
                iavg=0
                #

    x=x+1
    inputfile.close()

if SAVE:
    outputfile.close()
    outputfile_utc.close()

