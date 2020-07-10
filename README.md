# Online shot-LLR reconstruction (including data reading and DICOM generation)

[Shot-LLR](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27488) is a multi-shot diffusion-weighted MRI reconstruction method. Here we provide one example about how to make it run online on GE MRI scanners at [our institution (Lucas center, Stanford University)](https://med.stanford.edu/rsl/about/lucas.html).

## Pipeline
There are five steps from the raw data to the DICOM image on the scanner as below.

Step 1: Transfer the files (data and the acquisition information) from the scanner to the reconstruction server (since usually “fancy” reconstructions are computationally expensive and needs to run on some powerful servers). 

Step 2: Read the file and get the raw k-space data. We are using GE's Orchestra (can be found [here](https://collaborate.mr.gehealthcare.com/welcome
)) in this example, which also provide Nyquist artifact correction and ramp sampling correction. Some modifications are needed to make Orchestra output the k-space data, and we provide one way to do this [here](). We are using cfl format from [BART](https://mrirecon.github.io/bart/) to save the data.

Step 3: Shot-LLR reconstruction. We are doing this slice-by-slice. For each slcie, we calculate the sensitivity map first, and then apply shot-LLR for each nex/average/repetition. This is written in Python and the reconstruction is mainly based on [BART](https://mrirecon.github.io/bart/). Homodyne reconstruction is applied after shot-LLR reconstruction if Partial Fourier acquisition is used.

Spte 4: Generate DICOMs. Again, we use Orchestra in this step to make sure the header information is correct.

Step 5: Transfer the DICOMs back to the scanner.


Step 1 and Step 5 about transferring files are achieved with the help of [Dr. Marcus Alley](https://med.stanford.edu/profiles/marcus-alley), which may be implemented differently for different institutions and scanners. Steps 2 - 4 are provided in this example. The input is a folder containing the raw data from step 1 and the ouput is another folder containing the generated DICOMs for step 5.

### To be improved
(1) Parallel computing to accelerate the reconstruction.

(2) Orchestra and BART may be updated.
