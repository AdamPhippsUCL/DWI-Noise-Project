% Script to save .img files as .mat

pyfile = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Functions\Python\img2mat.py";

% Define img folder
imgfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20240319\ROIs";

% Define file extension
extension = '.img';

% Define .mat folder
matfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Imaging Data\ROIs";


pyrunfile(pyfile, imgfolder = imgfolder, matfolder = matfolder, extension=extension);

