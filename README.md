# shot-LLR-on-GE-s-scanners

Here we show one example of making your own reconstruction running online (right after the data collection) and putting images back to the scanner. The whole pipeline is as following, 
1) extract all required files (pfiles, scanArchives)to our own server,
2) use Orchestra to load the raw data, and apply some post-correction (ramp sampling correction and EPI odd-even correction in this example)
3) do your own reconstruction (shot-LLR for multi-shot DWI reconstruction in this example)
4) load the reconstructed images in Orchestra and call function "" to generate DICOMs 
5) put DICOMs back to the scanner.
