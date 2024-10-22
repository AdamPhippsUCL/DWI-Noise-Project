import numpy as np
from scipy.io import savemat
import SimpleITK as sitk
import glob
import os
import sys


def img2mat(folder, matfolder, extension = '.img'):
    
    '''
    Functions to convert all img files to mat files 
    
    '''
    
    
    # Read all img filenames in folder
    paths = glob.glob(f'{folder}/**/*{extension}')
    
    # Find relative paths (from folder) and file names
    relpaths = [path[len(folder):] for path in paths]   
    names = [(os.path.split(relpath)[-1])[0:-len(extension)] for relpath in relpaths]
    
    # Save each as .mat file
    for indx, path in enumerate(paths):
        
        
        # Attempt to read image
        try:
            img = sitk.GetArrayFromImage(sitk.ReadImage(path))
            print(f'Image read successfully: ...{relpaths[indx]}')
            print(f'Image shape: {np.shape(img)}')
        except:
            print(f'Unable to read image: ...{relpaths[indx]}')
            continue
            
        
        # Reshape image
        img = np.moveaxis(img , 0, -1)
           
        # Make folder for mat file
        newfolder = os.path.split(f'{matfolder}/{relpaths[indx][:-len(extension)]}.mat')[0]
        try:
            os.makedirs(
                newfolder
            )
        except:
            None
        
        # Save as mat
        savemat(
            f'{matfolder}/{relpaths[indx][:-len(extension)]}.mat',
            dict(zip(
                [names[indx]], [img]
            ))
        )
            
            
        
        
    
img2mat(
   imgfolder, 
   matfolder,
   extension)