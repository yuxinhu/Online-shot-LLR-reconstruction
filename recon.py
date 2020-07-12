#!/usr/local/lib python3.5

import subprocess as sp
import tempfile as tmp
import numpy as np
import os
import cfl
from scipy.misc import imresize
import getopt,sys

# This file is saved at multiple places for different scanners:
#    alfa for LPCH: /lpchNAS/scratch/bin/yuxinh
#    pinin for others: /hospNAS/scratch/bin/yuxinh
#    como1 for 3t2 and 3t3: /usr/local/recon/yuxinh
# The scripts may be slightly different. This version is from como1, and I try to make the code organized in this version to make
# it eaiser for the readers.

# The main function is called llr. There you can find llr21, llr22, llr23, llr24. They are about the same way, except how the sensitivity
# map is calculated (three ways: using ESPIRiT/BART on b0 data, using low resolution b0 data, or use the one extracted from scanArchive), 
# and whether llr reconstruction is applied to the b0 data. And I wrote comments carefully for llr21. At the beginning of each function, 
# there are explanations about how the sensitivity map is calculated. I am not sure if my current modified Orchestra can extract correct
# sensitivity information (I was struggling with the shift and some transformation), so I will not suggest to use that (llr22), but feel
# free to explore. As for whether to do LLR on b0 data, I think do the LLR will help with the denoising (though the regularization parameter 
# probably needs to be tuned), but considering that will slow down the algorithm, I am not using it on the clinical scanners. 
# I would suggest to use low resolution b=0 image as the sensitivity map for breast reconstruction (llr23), since ESPIRiT/BART will crop 
# the sensitivity map based on the eigenvalues which may not look very "natural" (this strategy works well on brain though).
# For all these llrs, I am doing LLR reconstruction first without zero-filled the unacquired k-space due to partial Fourier (cropped instead),
# since this helps reduce the matrix size and accelerate the algorithm a bit, and won't change the performance much based on my observation. 
# After LLR, I would transform the image into kspace and zero-filled the kspace, then do homodyne reconstruction.

def bart(nargout, cmd, *args):
    # This functions comes directly from BART. 
    if type(nargout) != int or nargout < 0:
        print("Usage: bart(<nargout>, <command>, <arguements...>)");
        return None

    bart_path = os.environ['TOOLBOX_PATH'];

    if not bart_path:
        if os.path.isfile('/usr/local/bin/bart'):
            bart_path = '/usr/local/bin'
        elif os.path.isfile('/usr/bin/bart'):
            bart_path = '/usr/bin'
        else:
            raise Exception('Environment variable TOOLBOX_PATH is not set.')

    name = tmp.NamedTemporaryFile().name

    nargin = len(args);
    infiles = [name + 'in' + str(idx) for idx in range(nargin)]
    in_str = ' '.join(infiles)

    for idx in range(nargin):
        cfl.writecfl(infiles[idx], args[idx])

    outfiles = [name + 'out' + str(idx) for idx in range(nargout)]
    out_str = ' '.join(outfiles)

    #TODO: Windows option.
    ERR = os.system(bart_path + '/bart ' + cmd + ' ' + in_str + ' ' + out_str);

    for elm in infiles:
        if os.path.isfile(elm + '.cfl'):
            os.remove(elm + '.cfl')
        if os.path.isfile(elm + '.hdr'):
            os.remove(elm + '.hdr')

    output = []
    for idx in range(nargout):
        elm = outfiles[idx]
        if not ERR:
            output.append(cfl.readcfl(elm))
        if os.path.isfile(elm + '.cfl'):
            os.remove(elm + '.cfl')
        if os.path.isfile(elm + '.hdr'):
            os.remove(elm + '.hdr')

    if ERR:
       raise Exception("Command exited with an error.")

    if nargout == 1:
        output = output[0]

    return output


class Header:
    # Useful acquisition parameters
    nx = 0; 
    ny = 0; # Ny size of the acquisition matrix
    acqy = 0; # Actually acquired number of ky lines (consiering partial Fourier)
    nc = 0; # Number of coils
    nslice = 0; # Number of slices
    np = 0; # Number of phases. This is a term from the PSD. For example, for a normal DWI scan (with b=0 and b=1000), the phase would be 2 (1 + 1).
    nshot = 0; # Number of shots
    asset = 0; # Whether ASSET is used or not. I haven't tested if ASSET (further acceleration) is supported (maybe not).
    itera = 100; # Number of iterations
    lam = 0.002; # Regularization parameter for the locally low-rank term.
    flag = 0; # Whether or not partial Fourier is used.
    
def usage():
    print('main.py -p <inputfilepath> -d <dicompath> -o <oxpath> -b <bartpath> -i <iteration> -r <lambda>')

def loadHeader(filepath, itera, lam):
    # All useful information are saved in headerIn when reading the rawdata. This function loads the information and saves them into header.
    # Other information could be added if necessary. I would suggest you to read the modified Orchestra page about how we achieve it.
    header2 = cfl.readcfl(filepath + "/cfls/headerIn")
    h = Header()
    h.nx = int(abs(header2[0]))
    h.acqy = int(abs(header2[1]))
    if int(abs(header2[7])) == 0:
        h.ny = h.acqy;
        h.flag = 0;
    else:
        h.ny = h.acqy*2 - 2*int(abs(header2[7]))
        h.flag = 1; # flag for homodyne
    h.resX = int(abs(header2[2]))
    h.resY = int(abs(header2[3]))
    h.nc = int(abs(header2[4]))
    h.nslice = int(abs(header2[5]))
    h.np = int(abs(header2[6]))
    h.nshot = int(abs(header2[8]))
    h.asset = int(abs(header2[9]))
    h.itera = itera;
    h.lam = lam;
    return h  


def getScale(filepath, h):
    # Load the first nex of the middle slice to get the scaling factor for all slices.
    filename = filepath + "/cfls/P" + str(0) + "_S" + str(h.nslice//2) + "_Nex" + str(0)
    kspace = cfl.readcfl(filename)
    ktemp = np.zeros((h.nx,h.acqy,1,h.nc),dtype = np.complex64)
    ktemp[:,0:h.acqy:1,0,:] = kspace # zero-filled k-space
    image = bart(1,'rss 8',bart(1,'fft -i -u 7',ktemp)) # Fourier transform then root-sum-of-square along the coil dimension.
    scale = float(abs(np.amax(image,(0,1)))) 
    return scale

def llr21(Bartpath, filepath,itera,lam, saveRes = 0, use_L1 = 0):
# Given the BART path, this function do the LLR (with homodyne reconstruction) for all files in the given file folder.
# The data in the filepath should be extracted by our modified Orchestra. Please refer the that page for details. Basically, phase 0 should be 
# the non-diffusion-weighted data which will be used for sensitivity map calculation, shot-LLR reconstruction will be applied to all the following
# phases.

# Use sensitivity map from b0 (BART/ESPIRiT)
# No LLR recon on b0 data
# Zero-filled after LLR recon then homodyne (for acceleration)

    os.environ["TOOLBOX_PATH"] = Bartpath
    h = loadHeader(filepath, itera, lam)
    sens = np.zeros((h.nx,h.acqy,h.nc,h.nslice),dtype = np.complex64) # To store the sensitivity map
    scale = getScale(filepath, h)
    res = np.zeros((h.resX,h.resY,h.nslice,h.np),dtype = np.complex64) # To store the final result
    com_homo2 = 'homodyne -C 1 ' + str(1.0*h.acqy/h.ny) # Homodyne command using BART


    for p in range(h.np): # loop for phase
        for s in range(h.nslice): # loop for slice
            n = 0; # Count the nex. Unfortunately, I do not know if the header has the nex for each phase.
            # This script counts the nex by checking if the corresponding file exists.

            temp = np.zeros((h.nx,h.ny),dtype = np.complex64) # Stores the reconstructed image for this slice and phase.

            while True: # loop for nex
                filename = filepath + "/cfls/P" + str(p) + "_S" + str(s) + "_Nex" + str(n) # Check if the next nex exists.
                if (os.path.exists(filename + ".cfl") == False):
                    temp = temp / n # Divide by the nex, since all averages are simplied added together before.
                    break

                kspace = cfl.readcfl(filename)/scale # The k-space to be reconstruted.
                sens2 = np.zeros((h.nx,h.acqy,1,h.nc),dtype = np.complex64) # To store the corresponding sensitivity map.

                if p == 0:
                    image2 = bart(1,'fft -i -u 3',kspace)
                    kspace = np.swapaxes(np.expand_dims(kspace, axis=3), 2, 3)
                    if n == 0:
                        # If this is the first nex and first phase (b=0) for the current slice, calculate the sensitivity map using BART (ESPIRiT).
                        sens[:,:,:,s] = np.squeeze(bart(1,'ecalib -m 1 -t 0.005', np.roll(kspace, (h.acqy-h.ny)//2, axis = 1)));
                        # It is tricky here, because we are doing the recon on the non-zero-padded k-space for speed, which means the signal peak (k=0)
                        # is not in the center of the matrix, which requires a shift of the kspace (by np.roll).

                    image2 = np.sum(image2*np.conjugate(sens[:,:,:,s]),2) # Coil combination
                    image = np.zeros((h.nx,h.ny), dtype = np.complex64)
                    image[:,0:h.acqy] = bart(1,'fft -u 3',image2) # This is actually the zero-padded k-space after coil combination.
                    image = bart(1, com_homo2, image) # Homodyne reconsruction.

                else:
                    # Do shot-LLR then homodyne for phase > 1 data (non-diffusion-weighted images).

                    sens2[:,:,0,:] = sens[:,:,:,s] # Extract the previously calculated sensitivity map (based on b=0 data)

                    ktemp = np.zeros((h.nx,h.acqy,1,h.nc,1,h.nshot),dtype = np.complex64)
                    for kk in range(h.nshot):
                        # Allocate the kspace data for different shots (dimension: nx - ny (nacqy) - 1 - nc - 1 - nshot).
                        ktemp[:,kk:h.acqy:h.nshot,0,:,0,kk] = kspace[:,kk:h.acqy:h.nshot,:]

                    if use_L1 > 0:
                        com = 'pics -l1 -r ' + str(use_L1) + ' -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                        # This is too slow to add another L1 constraint (have to use ADMM to solve the problem)
                    else:       
                        com = 'pics -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)

                    print('reconstructing phase ' + str(p) + ' slice ' + str(s) + ' nex ' + str(n) + '\n')

                    image2 = bart(1,com,ktemp,sens2) # Reconstructed multi-shot images using shot-LLR
                    if h.flag == 0:
                        image = np.average(np.squeeze(image2),2)
                    else:
                        image3 = np.zeros((h.nx,h.ny,h.nshot),dtype = np.complex64)
                        image3[:,0:h.acqy,:] = bart(1,'fft -u 3',np.swapaxes(image2,2,5))
                        for kk in range(h.nshot):
                            image3[:,:,kk] = bart(1,com_homo2,image3[:,:,kk]) # Partial Fourier reconstruction.

                        image = np.average(abs(image3),2) # Using amplitude average for different shots (not a big difference after homodyne recon?)
                
                temp = temp + image * scale; # Adding images from different nex/averages.
                n = n + 1 # counts the nex/averages

            if saveRes == 1 and p > 0:
                # Save the temporary reconstruction results if required.
                cfl.writecfl(filepath + "/cfls/res21_P" + str(p) + "_S" + str(s), temp)   

            temp = imresize(abs(temp/10.0),(h.resX,h.resY),mode = 'F')
            res[:,:,s,p] = np.array(temp, dtype = np.complex64)

    if saveRes == 1:
        cfl.writecfl(filepath + "/cfls/sens_b0", sens)  

    cfl.writecfl(filepath + "/cfls/res",res) # Save the final reconstruction results for DICOM generation.
    return 0

def llr22(Bartpath, filepath,itera,lam, saveRes = 0, use_L1 = 0):
# Use sensitivity map from scanArchives
# Do LLR recon on b0 data (similar to product MUSE)
# Zero-filled after LLR recon then homodyne (for acceleration)

    os.environ["TOOLBOX_PATH"] = Bartpath
    h = loadHeader(filepath, itera, lam)

    #sens = np.zeros((h.nx,h.ny,h.nc,h.nslice),dtype = np.complex64)
    sens = cfl.readcfl(filepath + "/cfls/sens")
    # silly python
    res = np.zeros((h.resX,h.resY,h.nslice,h.np),dtype = np.complex64)
    com_homo = 'homodyne -C 1 ' + str(1.0*h.acqy/h.ny)
    com_homo2 = 'homodyne -C 1 ' + str(1.0*h.acqy/h.ny)
    #print(com_homo2)
    scale = getScale(filepath, h)
  
    for p in range(h.np): # loop for phase
        for s in range(h.nslice): # loop for slice
            n = 0;
            temp = np.zeros((h.nx,h.ny),dtype = np.complex64)
            while True: # loop for nex
                filename = filepath + "/cfls/P" + str(p) + "_S" + str(s) + "_Nex" + str(n)
                if (os.path.exists(filename + ".cfl") == False):
                    #res[:,:,s,p] = res[:,:,s,p] / n 
                    temp = temp / n
                    break
                kspace = cfl.readcfl(filename)
                #print(filename)

                if 1:
                    # do LLR then homodyne if phase > 1 (DWIs)
                    ktemp = np.zeros((h.nx,h.acqy,1,h.nc,h.nshot,1,1,1,1,1,1),dtype = np.complex64)
                    for kk in range(h.nshot):
                        ktemp[:,kk:h.acqy:h.nshot,0,:,kk,0,0,0,0,0,0] = kspace[:,kk:h.acqy:h.nshot,:]
                    ktemp = np.swapaxes(ktemp,4,10) # split different shots into different frames
                    if use_L1 > 0:
                        com = 'pics -l1 -r ' + str(use_L1) + ' -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    else:       
                        com = 'pics -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    sens2 = np.zeros((h.nx,h.acqy,1,h.nc),dtype = np.complex64)
                    sens2[:,:,0,:] = sens[:,:,:,s]
                    #print(com)

                    image2 = bart(1,com,ktemp/scale,sens2) # reconstructed multi-shot images
                    if h.nshot > 1:
                        image2 = np.swapaxes(image2,2,10)
                        image2 = image2[:,:,:,0,0,0,0,0,0,0,0];
                        image3 = np.zeros((h.nx,h.ny,h.nshot),dtype = np.complex64)
                        image3[:,0:h.acqy,:] = bart(1,'fft -u 3',image2)
                        if h.flag == 0:
                            image = np.average(abs(image2),2)
                        elif p == 0:
                            image = bart(1, com_homo2, np.average(image3,2))
                        else:
                            for kk in range(h.nshot):
                                image3[:,:,kk] = bart(1,com_homo2,image3[:,:,kk])
                                image = np.average(abs(image3),2)
                    else:
                        if h.flag == 0:
                            image = abs(image2)
                        else:
                            image = abs(bart(1,com_homo2,image3))
                    
                temp = temp + image * scale;
                #res[:,:,s,p] = res[:,:,s,p] + image*scale;
                n = n + 1
            temp = imresize(abs(temp/10.0),(h.resX,h.resY),mode = 'F')
            if saveRes == 1:
                cfl.writecfl(filepath + "/cfls/res_P" + str(p) + "_S" + str(s), temp)
            res[:,:,s,p] = np.array(temp,dtype = np.complex64)
    cfl.writecfl(filepath + "/cfls/res",res)
    return 0

def llr23(Bartpath, filepath,itera,lam, saveRes = 0, use_L1 = 0):
# Use sensitivity map from b0 (low res instead of BART)
# No LLR recon on b0 data
# Zero-filled after LLR recon then homodyne (for acceleration)

    os.environ["TOOLBOX_PATH"] = Bartpath
    h = loadHeader(filepath, itera, lam)
    sens = np.zeros((h.nx,h.acqy,h.nc,h.nslice),dtype = np.complex64)
    scale = getScale(filepath, h)
    # silly python
    res = np.zeros((h.resX,h.resY,h.nslice,h.np),dtype = np.complex64)
    com_homo2 = 'homodyne -C 1 ' + str(1.0*h.acqy/h.ny)
    #print(com_homo2)
    for p in range(h.np): # loop for phase
        for s in range(h.nslice): # loop for slice
            n = 0;
            temp = np.zeros((h.nx,h.ny),dtype = np.complex64)
            while True: # loop for nex
                filename = filepath + "/cfls/P" + str(p) + "_S" + str(s) + "_Nex" + str(n)
                if (os.path.exists(filename + ".cfl") == False):
                    #res[:,:,s,p] = res[:,:,s,p] / n 
                    temp = temp / n
                    break
                kspace = cfl.readcfl(filename)/scale
                #print(filename)
                sens2 = np.zeros((h.nx,h.acqy,1,h.nc),dtype = np.complex64)
                if p == 0:
                    image2 = bart(1,'fft -i -u 3',kspace)
                    #kspace = np.swapaxes(np.expand_dims(kspace,axis=3),2,3)
                    if n == 0:
                        #sens[:,:,:,s] = np.squeeze(bart(1,'ecalib -m 1 -t 0.005', np.roll(kspace, (h.acqy-h.ny)//2, axis = 1)));
                        kspace_temp = np.zeros((h.nx,h.acqy,h.nc),dtype = np.complex64)
                        kspace_temp[:,h.acqy//2 - 24: h.acqy//2 + 23,:] = kspace[:,h.ny//2 - 24: h.ny//2 + 23,:]
                        sens_temp = bart(1, 'fft -i -u 3', kspace_temp)
                        sens[:,:,:,s] = sens_temp / np.tile(np.expand_dims(np.sqrt(np.sum(abs(sens_temp)**2, axis = 2)), axis = 2),(1,1,h.nc))
                    #sens2[:,:,0,:] = sens[:,:,:,s]
                    #image2 = bart(1,'pics -r 0 -i 1 -w 1',kspace,sens2)
                    #image2 = np.sum(image2*sens[:,:,:,s],2)
                    image2 = np.sum(image2*np.conjugate(sens[:,:,:,s]),2)
                    image = np.zeros((h.nx,h.ny), dtype = np.complex64)
                    image[:,0:h.acqy] = bart(1,'fft -u 3',image2)
                    image = bart(1, com_homo2, image)
                else:
                    # do LLR then homodyne if phase > 1 (DWIs)
                    sens2[:,:,0,:] = sens[:,:,:,s]
                    ktemp = np.zeros((h.nx,h.acqy,1,h.nc,1,h.nshot),dtype = np.complex64)
                    for kk in range(h.nshot):
                        ktemp[:,kk:h.acqy:h.nshot,0,:,0,kk] = kspace[:,kk:h.acqy:h.nshot,:]
                    # split different shots into different frames
                    if use_L1 > 0:
                        com = 'pics -l1 -r ' + str(use_L1) + ' -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    else:       
                        com = 'pics -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    print('reconstructing phase ' + str(p) + ' slice ' + str(s) + ' nex ' + str(n) + '\n')
                    image2 = bart(1,com,ktemp,sens2) # reconstructed multi-shot images
                    if h.flag == 0:
                        image = np.average(np.squeeze(image2),2)
                    else:
                        image3 = np.zeros((h.nx,h.ny,h.nshot),dtype = np.complex64)
                        image3[:,0:h.acqy,:] = bart(1,'fft -u 3',np.swapaxes(image2,2,5))
                        for kk in range(h.nshot):
                            image3[:,:,kk] = bart(1,com_homo2,image3[:,:,kk])
                        image = np.average(abs(image3),2)
                temp = temp + image * scale;
                #res[:,:,s,p] = res[:,:,s,p] + image*scale;
                n = n + 1
            if saveRes == 1 and p > 0:
                cfl.writecfl(filepath + "/cfls/res22_lam_P" + str(p) + "_S" + str(s), temp)            
            temp = imresize(abs(temp/10.0),(h.resX,h.resY),mode = 'F')
            res[:,:,s,p] = np.array(temp, dtype = np.complex64)
    if saveRes == 1:
        cfl.writecfl(filepath + "/cfls/sens22_b0", sens)  
    cfl.writecfl(filepath + "/cfls/res",res)
    return 0



def llr24(Bartpath, filepath,itera,lam, saveRes = 0, use_L1 = 0):
# Use sensitivity map from b0 (low res instead of BART)
# Do LLR recon on b0 data (similar to product MUSE)
# Zero-filled after LLR recon then homodyne (for acceleration)

    os.environ["TOOLBOX_PATH"] = Bartpath
    h = loadHeader(filepath, itera, lam)
    sens = np.zeros((h.nx,h.acqy,h.nc,h.nslice),dtype = np.complex64)
    scale = getScale(filepath, h)
    # silly python
    res = np.zeros((h.resX,h.resY,h.nslice,h.np),dtype = np.complex64)
    com_homo2 = 'homodyne -C 1 ' + str(1.0*h.acqy/h.ny)
    #print(com_homo2)
    for p in range(h.np): # loop for phase
        for s in range(h.nslice): # loop for slice
            n = 0;
            temp = np.zeros((h.nx,h.ny),dtype = np.complex64)
            while True: # loop for nex
                filename = filepath + "/cfls/P" + str(p) + "_S" + str(s) + "_Nex" + str(n)
                if (os.path.exists(filename + ".cfl") == False):
                    #res[:,:,s,p] = res[:,:,s,p] / n 
                    temp = temp / n
                    break
                kspace = cfl.readcfl(filename)/scale
                #print(filename)
                sens2 = np.zeros((h.nx,h.acqy,1,h.nc),dtype = np.complex64)
                if p == 0:
                    image2 = bart(1,'fft -i -u 3',kspace)
                    #kspace = np.swapaxes(np.expand_dims(kspace,axis=3),2,3)
                    if n == 0:
                        #sens[:,:,:,s] = np.squeeze(bart(1,'ecalib -m 1 -t 0.005', np.roll(kspace, (h.acqy-h.ny)//2, axis = 1)));
                        kspace_temp = np.zeros((h.nx,h.acqy,h.nc),dtype = np.complex64)
                        kspace_temp[:,h.acqy//2 - 24: h.acqy//2 + 23,:] = kspace[:,h.ny//2 - 24: h.ny//2 + 23,:]
                        sens_temp = bart(1, 'fft -i -u 3', kspace_temp)
                        sens[:,:,:,s] = sens_temp / np.tile(np.expand_dims(np.sqrt(np.sum(abs(sens_temp)**2, axis = 2)), axis = 2),(1,1,h.nc))

                    #image2 = np.sum(image2*np.conjugate(sens[:,:,:,s]),2)
                    #image = np.zeros((h.nx,h.ny), dtype = np.complex64)
                    #image[:,0:h.acqy] = bart(1,'fft -u 3',image2)
                    #image = bart(1, com_homo2, image)
                if 1:
                    # do LLR then homodyne if phase > 1 (DWIs)
                    sens2[:,:,0,:] = sens[:,:,:,s]
                    ktemp = np.zeros((h.nx,h.acqy,1,h.nc,1,h.nshot),dtype = np.complex64)
                    for kk in range(h.nshot):
                        ktemp[:,kk:h.acqy:h.nshot,0,:,0,kk] = kspace[:,kk:h.acqy:h.nshot,:]
                    # split different shots into different frames
                    if use_L1 > 0:
                        com = 'pics -l1 -r ' + str(use_L1) + ' -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    else:       
                        com = 'pics -R L:7:7:' + str(h.lam) + ' -w 1 -i ' + str(h.itera)
                    print('reconstructing phase ' + str(p) + ' slice ' + str(s) + ' nex ' + str(n) + '\n')
                    image2 = bart(1,com,ktemp,sens2) # reconstructed multi-shot images
                    if h.flag == 0:
                        image = np.average(np.squeeze(image2),2)
                    else:
                        image3 = np.zeros((h.nx,h.ny,h.nshot),dtype = np.complex64)
                        image3[:,0:h.acqy,:] = bart(1,'fft -u 3',np.swapaxes(image2,2,5))
                        for kk in range(h.nshot):
                            image3[:,:,kk] = bart(1,com_homo2,image3[:,:,kk])
                        image = np.average(abs(image3),2)
                temp = temp + image * scale;
                #res[:,:,s,p] = res[:,:,s,p] + image*scale;
                n = n + 1
            if saveRes == 1 and p > 0:
                cfl.writecfl(filepath + "/cfls/res23_lam_P" + str(p) + "_S" + str(s), temp)            
            temp = imresize(abs(temp/10.0),(h.resX,h.resY),mode = 'F')
            res[:,:,s,p] = np.array(temp, dtype = np.complex64)
    if saveRes == 1:
        cfl.writecfl(filepath + "/cfls/sens23_b0", sens)  
    cfl.writecfl(filepath + "/cfls/res_23", res)
    return 0





