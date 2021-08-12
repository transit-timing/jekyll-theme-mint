


# Transit timing

## Processing *TESS* 2-minute cadence light curves in batches for multiple targets

### Step 1: Download .fits files for each target using X script
```
python3 X.py
```


### Step 2: Run tt.py script 
```
python3 tt.py
```
A typical folder that this script produces for a given target contains the following:

## Processing a single *TESS* light curve (either 2-minute or 30-minute data)

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
a one-line code block
```
