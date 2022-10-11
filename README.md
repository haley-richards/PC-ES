# PC-ES

MATLAB App developed to localize hd-EEG electrodes from 3D Head Scans 
What you need:
• Point clouds of the EEG cap (.ply files)
• Reference data (optional)
Compatibility Requirements:
• Windows 10 OS
• Any hdEEG 256 Electrode cap 

******* For a detailed user guide with corresponding images, download the 'PC-ES' pdf document in the zip file**************

What it does and why:
Accurate localization of scalp electrodes is instrumental to mapping brain activity using high- density electroencephalography (hdEEG). Existing technologies for electrode localization can be costly and difficult to transport, which can restrict access to hdEEG monitoring. This software enables a novel and affordable electrode digitization technique by leveraging the 3D scanning feature in modern iPhones (and other 3D scanning devices). This interactive semi-automated MATLAB GUI intuitively allows the user to visualize, select, and label all electrodes in a 256 EEG cap using point clouds alone. To compare this method to other methods (photogrammetry, electromagnetic tracking, etc) there is an option to upload a set of reference data to compare to the results from PC-ES.

PC-ES was first designed to be compatible with 3D scanning features on modern iPhones using the Heges app. Since the Heges app can not take individual scans of the whole head, the app was originally designed to label and merge a variable number of partial scans of the head. However, PC-ES is compatible with any scanning device that produces a .ply file, as well as a single and complete scan of the head. While this app currently only supports an 256 electrode HD-EEG cap, future versions will be compatible with any EEG cap. 

 The processing approach is semi-automated, involving both user-driven and automatic processing steps.
After uploading the scan, the user crops the scan for improved processing speed if needed and deidentifies the facial region for HIPAA compliance. The user then selects and manually labels fiducial points (electrodes with known labels that are used as reference landmarks). These may be marked with selected colored stickers during data collection to ease in recognition. Lastly, the user selects the remaining electrodes visible in the scan without assigning labels, concluding the manual processing steps. Future versions of this app will allow for fully automated electrode selection using AI image classification algorithms. 

After electrodes and fiducial points are marked, the program automatically labels the selected points, and uses them to sequentially merge the point clouds. This results in full a point cloud of the whole head in a single reference frame, as well as a set of labeled coordinates. The program then identifies electrodes that the user has missed and interpolates them from the merged point cloud.
