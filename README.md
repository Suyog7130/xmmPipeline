[![PDF](https://img.shields.io/badge/xmmPipeline-PDF-brightgreen?style=plastic)](https://github.com/Suyog7130/xmmPipeline/blob/main/docs/_build/latex/xmmpipeline.pdf)


Here are markdown README instructions for setting up your environment and installing HEASARC and XMMSAS:

---

# `epicPipeline` Setup Instructions

## 1. Clone the Repository

```bash
git clone https://github.com/yourusername/xmmPipeline.git
cd xmmPipeline
```

## 2. Create the Conda Environment

Create a new conda environment with all required dependencies:

```bash
conda env create -f environment.yml
conda activate xmm_epic_env
```

If you want to try the latest Python version and fallback automatically, use the provided setup script:

```bash
bash setup_env.sh
conda activate xmm_epic_env
```

## 3. Install HEASoft (HEASARC tools)

HEASoft is required for X-ray data analysis. Download and install it from the official site:

- [HEASoft Download & Installation Guide](https://heasarc.gsfc.nasa.gov/lheasoft/download.html)

**Quick steps:**

```bash
# Download the latest HEASoft source package
wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft.tar.gz
tar -xzvf lheasoft.tar.gz
cd heasoft-*
./configure
make
make install
```

- Follow the prompts and instructions on the HEASoft website for your OS.
- After installation, add the following to your `.bashrc` or `.zshrc`:

  ```bash
  source /path/to/heasoft/heasoft-init.sh
  ```

## 4. Install XMMSAS

XMMSAS is the XMM-Newton Science Analysis System. Download and install from:

- [XMMSAS Download & Installation Guide](https://www.cosmos.esa.int/web/xmm-newton/sas-download)

**Quick steps:**

- Download the appropriate installer for your OS.
- Follow the [official installation instructions](https://www.cosmos.esa.int/web/xmm-newton/sas-installation).
- Set up the environment variables in your `.bashrc` or `.zshrc`:

  ```bash
  export SAS_DIR=/path/to/xmmsas
  source $SAS_DIR/setsas.sh
  ```

## 5. Verify Installation

Check that the following commands work in your terminal:

```bash
heainit
sasversion
```

If both commands return version information, your setup is correct.

## 6. Run the Pipeline

You can now run the pipeline scripts, for example:

```bash
python epicPipeline.py --method reduceData --ra 192.0625 --dec 17.7739 --workdir /your/data/dir
```

## 7. Troubleshooting

- If you encounter issues with missing libraries, ensure your conda environment is activated and all dependencies are installed.
- For HEASoft and XMMSAS, always source their setup scripts before running the pipeline.

---

**Note:**  

- Replace `/path/to/heasoft` and `/path/to/xmmsas` with your actual installation paths.
- For more details, refer to the official documentation for [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) and [XMMSAS](https://www.cosmos.esa.int/web/xmm-newton/sas).

---

### Basic Description of pipeline modules

- **epicObj.py**, contains the main **epicObj** class which would be used by the rest of the programs. 

	- This would also contain two new functions, **writeCCDcoordsPickle** and **readCCDcoordsPickle**, for saving and reading the Pickle file contains the SrcBkg coordinates. 
	
	- Note that the presently available **readPickleFile** function in **spectra.py** actually reads the whole **xmmObj** from the saved Pickle of the whole class. This was probably done because I didn’t know back then how to do Class Inheritance and will be altered together with these modifications.*
	
- **epicPipeline.py**, contains separate methods for sequentially running different functions from the **epicObj** class.

	- **reduceData** method, to find obsIDs, download and reduce the EPIC data.
	
	- **combineAndFind** method, to find the Source CCDs, SrcBkg circle coordinates and to combine SrcBkg Event lists and instrumental and Flare GTIs.
	
		- This would run the **getBackgroundCircles** function which internally executes **findOverlap.py** to find the Background circle coordinates and the Source circle radius.
		
		- The found SrcBkg coordinate and circle parameters will be saved to **ccd\_coords\_info.pickle** file through **writeCCDcoordsPickle**.
		
		- After completion of this method, if Pile-up correction is required, it can be done using **epicPileup.py** which would read the pickle file using **readCCDcoordsPickle**.
		
	- **extractProds** method, to extract light curves from the Event Lists and create PNG images using matplotlib.
	
		- The method would start with reading the saved **ccd\_coords\_info.pickle** file contains a dictionary of the SrcBkg coordinates using **readCCDcoordsPickle**.
		
		- The **extractProds** method can be run in continuous manner directly after the **combineAndFind** method for cases when Pile-up correction is not wanted, or it can be run after correcting for Pile-up using **epicPileup.py** which would have rewritten **ccd\_coords\_info.pickle** after changing the Source coordinate parameters to an Annulus instead of a circle.
		
		- If any particular obsID is Piled-up then running **extractProds** would know from the content of **ccd\_coords\_info.pickle**
		
		- This would also run the **saveResults** function to copy the resulting files to a higher level **results** directory in **workdir**.
		
- **findOverlap.py**, internally called routine to find the Background circles.

- **epicPileup.py**, contains the **epicPileup** class which has two functions **checkPileUp** and **correctPileUp** for checking and correcting the Pile-up issue.

	- The **epicPileup** class inherits the **epicObj** from **epicObj.py** and reads the SrcBkg circle parameters from **ccd\_coords\_info.pickle** file already saved by **combineAndFind** method, using **readCCDcoordsPickle**.
	
	- If Pile-up correction is done, this routine would modify the **ccd\_coords\_info.pickle** file and rewrite it to the directory using **writeCCDcoordsPickle**. This rewritten pickle file would then be used by **spectra.py** and the **extractProds** method.
	
- **spectra.py**, contains the **spectra** class which inherits the **epicObj** from **epicObj.py** and contains additonal functions to extract Spectra, fit models to the Spectra and save plots and other results.

	- This can be renamed to **epicSpectra.py** since it’s only used for extracting EPIC spectra.
	
	- There are numerous modes, like Image mode, Timing mode, Burst mode etc. for which Spectra can be extracted. For the present work only the Image mode is important and thus, that’s the only mode for which **spectra** class has a function currently available. For Timing mode SrcBkg coordinates need to be in RAWX and RAWY anyway.

	- Contains the **xspec_fitSpectra** function which uses the **pyXspec.py** routine to fit models to the Spectra.

	- Since this will be executed after all the methods in **epicPipeline** and the Pile-up correction using **epicPileup.py** are finished, there shouldn’t be any difficulties arising from Pile-up correction related procedures.
	
- **pyXspec.py**, internally called routine for using **pyxspec** package and fitting models to the Spectra.

- **rgsPipeline.py**, for RGS data reduction and light curve and Spectra generation.
	- Nothing to be changed in this presently. Modifications may perhaps be required later on.
	
- **plotAnal.py**, contains utility functions for *matplotlib* plot beautification.

- **convert.py**, to convert format of RA and DEC.
