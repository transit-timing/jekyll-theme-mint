


# Transit timing




<<<<<<< HEAD
## Processing a single *TESS* light curve (either 2-minute or 30-minute data)

1. To process a single *TESS* light curve and extract mid-transit times, first go to the directory with the code:

```
cd 0_tt
```
2. In line 1185, specify the name of the target that you would like to study. For example, to run the code on WASP-012, we write
```
favorite = 'WASP-012'
```
*1_target_list.csv* contains a list of the targets studied in our work. You may only select a target name from this list (as some information from this file (such as transit depth, etc.) is used in the script). You should write the name of the target as it appears in the table (for example, WASP-012 instead of WASP-12).

3. You then need to decide which data you would like to analyse. You can check what kinds of *TESS* data are available for a given target using the *30-min TESS w/ lightkurve.ipynb* notebook. The notebook contains instructions on how to run it. After running the notebook, you should note the index of the data product that you would like to use; it will be shown in the SearchResult output of the *search_lightcurve* function from the *lightkurve* package. Specify this index in line 1186:
```
idx = 0
```
4. Now you are ready to run the script:
=======
## Processing *TESS* light curves 

1. To process *TESS* light curves stored in the *~/data/* folder and extract mid-transit times, first go to the directory with the code:

```
cd 0_tt
```
 Then run the script:
```
python3 main.py
```
### Output
The results are stored in the *~/output* folder. A typical folder that this script produces for a given target contains the following:

1. Figure of the  *TESS* light curve that was processed (transits that were processed to extract mid-transit times are shown in red).

![Original *TESS* light curve](/5_data/WASP-012_dir/WASP-012_a_TimeSeries.png)

2. Figure with the individual transits that were processed (a polynomial that was fitted to out-of-transit data during de-trending is shown in red).

![Individual transits](/5_data/WASP-012_dir/WASP-012_b_IndividualTransits.png)

3. Figure with the de-trended individual transits.

![De-trended individual transits](/5_data/WASP-012_dir/WASP-012_c_IndividualTransitsDetrended.png)

4. Folded light curve.

![Folded light curve at the first iteration of our procedure](/5_data/WASP-012_dir/WASP-012_d_FoldedLightCurve.png)

![Folded light curve at the second iteration of our procedure](/5_data/WASP-012_dir/WASP-012_e_FoldedLightCurve.png)

5. Individual events with the transit model fitted.

![Individual events with the transit model fitted](/5_data/WASP-012_dir/WASP-012_f_IndividualTransitsWithFit.png)

6. O-C diagram showing timing residuals. The timing residuals were calculated as follows: first, a linear model was fitted to the transit times extracted from *TESS* data. Then, the linear model was subtracted from the transit times at each epoch.

![Timing residuals](/5_data/WASP-012_dir/WASP-012_g_TimingResiduals.png)

7. *WASP-012_results.txt* contains extracted mid-transit times and their uncertainties.

8. *WASP-012_log.txt* contains intermediate outputs of the code, such as the found best-fit transit model parameters and different statistics of the produced fits.


**Note:** in order to process *TESS* light curves, you may want to first download the relevant *.fits* files. To download a given light curve, you may do the following:

```
cd 0_download_data
```
In the  *download_single_target.py* script, specify the name of the target and the light curve that you would like to download. Then run the script:
>>>>>>> 1bfeb1ffa1aa9c09458be4082e49ef8945366c30
```
python3 download_single_target.py
```
<<<<<<< HEAD
### Output
The results are be stored in *~/5_data* folder. A typical folder that this script produces for a given target contains the following:

1. Figure of the  *TESS* light curve that was processed (transits that were processed to extract mid-transit times are shown in red).

![Original *TESS* light curve](/5_data/WASP-012_dir/WASP-012_a_TimeSeries.png)

2. Figure with the individual transits that were processed (a polynomial that was fitted to out-of-transit data during de-trending is shown in red).

![Individual transits](/5_data/WASP-012_dir/WASP-012_b_IndividualTransits.png)

3. Figure with the de-trended individual transits.

![De-trended individual transits](/5_data/WASP-012_dir/WASP-012_c_IndividualTransitsDetrended.png)

4. Folded light curve at the first iteration of our procedure.

![Folded light curve at the first iteration of our procedure](/5_data/WASP-012_dir/WASP-012_d_FoldedLightCurve.png)

5. Folded light curve at the second iteration of our procedure.

![Folded light curve at the second iteration of our procedure](/5_data/WASP-012_dir/WASP-012_e_FoldedLightCurve.png)

6. Individual events with the transit model fitted.

![Individual events with the transit model fitted](/5_data/WASP-012_dir/WASP-012_f_IndividualTransitsWithFit.png)

7. O-C diagram showing timing residuals. The timing residuals were calculated as follows: first, a linear model was fitted to the transit times extracted from *TESS* data. Then, the linear model was subtracted from the transit times at each epoch.

![Timing residuals](/5_data/WASP-012_dir/WASP-012_g_TimingResiduals.png)

8. *WASP-012_results.txt* contains extracted mid-transit times and their uncertainties.

9. *WASP-012_log.txt* contains intermediate outputs of the code, such as the found best-fit transit model parameters and different statistics of the produced fits.


**Note:** in order to process *TESS* light curves in batches for multiple targets, you may first download the relevant *.fits* files from MAST and add a few lines in the *tt.py* script to read in the files in a loop.

=======
>>>>>>> 1bfeb1ffa1aa9c09458be4082e49ef8945366c30


# Figures

All of the O-C diagrams are available on Google Drive: 

# Tables

1. A single CSV file with all of the collected transit times for each target is available on Google Drive: 
2. A single CSV file with all of the target ephemerides is available on Google Drive: 

# Scraping literature
1. Download bibcodes of all papers for a given object on ADSABS website:
    - Go to https://ui.adsabs.harvard.edu/search
    - Do cone search using RaDec coordinates with 3" radius; put the .bib file into ~/article_database/{planet}/{planet}.bib. Alternatively, you can type *object:"planet_name"* in the search bar to retrieve articles for a planet with a given name *planet_name*
    
```
object:"20h30m54.13s 6d25m46.4s:3"
```
2. Run **bibcode2arxiv_id.py** to create a .txt file containing a list of arXiv IDs of papers for a given object with the name *planet_name*, using the arXiv IDs reported in the bibtex file downloaded in the previous step:
```  
python3 bibcode2arxiv_id.py -p planet_name
```
3. Run **run.py** to download *.tex* source files of papers with the arXiv IDs retrieved in the previous step. Untar the files and transfer all .tex files to a different directory.
    For a given planet, all tar.gz are stored in ~/arxiv_dir/{planet_name} directory.
    For a given planet, all untarred files are stored in ~/untar_dir/{planet_name} directory.
4. Run **read_tex.py** to find numbers > 2,000,000; if such numbers and the name of our planet are present in a given paper, a PDF of the paper is transferred to ~/{planet_name} directory for final manual review.
```
arxiv % python3 read_tex.py -s wasp -id 2
```
# Citation
```
DOI
```
