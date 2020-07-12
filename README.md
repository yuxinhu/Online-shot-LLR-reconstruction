# Online shot-LLR reconstruction (including data reading and DICOM generation)

[Shot-LLR](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27488) is a multi-shot diffusion-weighted MRI reconstruction method. Here we provide one example about how to make it run online on GE MRI scanners at [our institution (Lucas center, Stanford University)](https://med.stanford.edu/rsl/about/lucas.html).

## Usage
"sllr" is an excutable python file, which you can put some command like the following in the Terminal:

        /sllr -p RAW_DATA_PATH -d OUTPUT_PATH
        
Sorry I would love to provide example input data here but they are all huge and there might be some privacy concerns. If you have access to [our group](https://med.stanford.edu/bmrgroup.html) server, you should be able to find some data under my folder (/bmrNAS/people/yuxinh).

## Pipeline
There are five steps from the raw data to the DICOM image on the scanner as below.

Step 1: Transfer the files (data and the acquisition information) from the scanner to the reconstruction server (since usually “fancy” reconstructions are computationally expensive and needs to run on some powerful servers). 

Step 2: Read the file and get the raw k-space data. We are using GE's Orchestra (can be found [here](https://collaborate.mr.gehealthcare.com/welcome
)) (SDK 1.9-1) in this example, which also provide Nyquist artifact correction and ramp sampling correction. Some modifications are needed to make Orchestra output the k-space data, and we provide one way to do this [here](). We are using cfl format from [BART](https://mrirecon.github.io/bart/) to save the data.

Step 3: Shot-LLR reconstruction. We are doing this slice-by-slice. For each slcie, we calculate the sensitivity map first, and then apply shot-LLR for each nex/average/repetition. This is written in Python and the reconstruction is mainly based on [BART](https://mrirecon.github.io/bart/) (version 0.4.03). Homodyne reconstruction is applied after shot-LLR reconstruction if Partial Fourier acquisition is used. Detailed comments in recon.py.

Spte 4: Generate DICOMs. Again, we use Orchestra in this step to make sure the header information is correct.

Step 5: Transfer the DICOMs back to the scanner.


Step 1 and Step 5 about transferring files are achieved with the help of [Dr. Marcus Alley](https://med.stanford.edu/profiles/marcus-alley), which may be implemented differently for different institutions and scanners. Steps 2 - 4 are provided in this example. The input is a folder containing the raw data from step 1 and the ouput is another folder containing the generated DICOMs for step 5. We provide the compiled BART and Orchestra (for Linux), while re-compiling may be needed. Also notice that we have modified Orchestra to save the k-space data and some header information. Before you dive into the scripts, I would suggest you to read [how we save the data]() in step 2.

### To be improved
(1) Parallel computing to accelerate the reconstruction.

(2) Orchestra and BART may be updated.
