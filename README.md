# Fiberphotometry Analysis

Python tools for preprocessing, analyzing, and visualizing fiber photometry recordings, with integration of behavioral and physiological data.

## Overview

This repository contains scripts and reusable Python modules for processing fiber photometry recordings and combining them with behavioral and physiological measurements.

## Repository structure

```text
Fiberphotometry_analysis/
│
├── modules/
│   ├── behaviour/
│   └── common/
│
├── scripts/
│
├── tutorial/
│   └── Photometry data preprocessing.ipynb
│
├── requirements.txt
├── activate.bat
├── video_alignment_pipeline.bat
├── .gitignore
│
└── README.md
```

The analysis pipeline includes:

- Fiber photometry preprocessing :
    - Signal cleaning and filtering.
    - Delta F/F calculation.

- Behavioral data alignment and analysis (depending on task):
    - DeepLabCut tracking analysis and behavioral data extraction
    - Lick and airpuff detection 
    - Fear conditioning analysis.
    - Elevated plus maze analysis.
    - Whole-body plethysmography analysis.

- Final outputs :
    - Peri-event time histograms (PETH)
    - Mean dFF and AUC during behaviors or zones
    - Transient detection and quantification
    - Correlational analyses for dual color data

- Alignment of fiberphotometry data and behaviour on video

## Installation

### Requirements
The complete list of dependencies is available in `requirements.txt`.

### Install dependencies

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Windows:

```bash
.venv\Scripts\activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

## Data organization

The analysis scripts use an experiment directory containing the raw `Data` and an `Analysis` directory for generated results.

The main path is configured in:

```text
scripts/loader.py
```

For every analysis, you need to provide the list of mice you want to analyse and the name of the task. The code will retrieve the data corresponding to the task and ID of the mouse. Thus all subjects and raw data paths must be written in excels in the experiment directory (before running `scripts/loader.py`) :

- `experiment_path \ subjects.xlsx` :

| Subject | Batch | Sex | Group    |
| ------- | ----- | --- | -------- |
| 1009    | 1     | M   | Control  |
| 1008    | 1     | F   | Treated  |
| 1003    | 2     | M   | Treated  |
| 1001    | 2     | F   | Control  |
| 1000    | 2     | F   | Treated  |

Preprocessing and all individual analyses will be done on the mice listed in the `Included` sheet. For every task listed in  `experiment_path \ protocol.xlsx`, an exclusion sheet must be generated in `experiment_path \ subjects.xlsx`, with the name `Excluded_{task}`. All the mice listed in the exclusion sheet will not be taken into account in group analyses.

- `experiment_path \ protocol.xlsx` :

| Batch | Task          | Data_path           |
| ----- | ------------- | ------------------- |
| 1     | EPM           | 20260218_EPM        |
| 1     | FearConditioning| 20260219_FC       |
| 2     | EPM           | 20260520_EPM        |


## Analysis workflow

### 1. Preprocessing

Main script:

```text
scripts/script_preprocess.py
```

The preprocessing stage prepares raw fiber photometry data for downstream analysis.

It includes:

* Loading and organizing raw recordings.
* Deinterleaving and cleaning signals.
* Manual artifact flagging and correction.
* Calculation of delta F/F.
* Removing residual low-frequency components (photobleaching)
* Manual flagging of residual corrupted data (initial exponential decay, lost data due to patch cord detachment)
* Saving preprocessed data and diagnostic plots.
    * Deinterleaved data.
    * Cleaned data.
    * Processed dFF data.
    * Raw signal plots.
    * Cleaned signal plots.
    The exact filenames are defined in `scripts/script_preprocess.py`.

A tutorial Jupyter notebook details the preprocessing steps that are used in the pipeline :

```text
tutorial/Photometry data preprocessing.ipynb
```

### 2. Behavioral analysis

The main behavioral analysis scripts are:

| Script                     | Purpose                                           |
| -------------------------- | ------------------------------------------------- |
| `script_fiberboris.py`     | Analysis using BORIS manula behavioral scoring    |
| `script_fiberEPM.py`       | Elevated plus maze analysis using DeepLabCut data |
| `script_fiberFC.py`        | Fear conditioning analysis                        |
| `script_fiberlicks.py`     | Lick and airpuff analysis                         |
| `script_fiberplethysmo.py` | Fiber photometry and plethysmography analysis     |

The general purpose of those scripts is :

* Processing behavioral events (including extracting and processing DeepLabCut tracking data if relevant)
* Aligning behavior and photometry, and plotting.
* Extracting behavioral data and plotting occupation heatmaps (if relevant)
* Extracting mean dFF during behaviours / zones and plotting mean dFF on the arena heatmap (if relevant)

A few parameters have to be chosen before the analysis :

```python
# If data was recorded via Bonsai, with a camera_flashes csv file, set to True (else camera flashes were directly recorded in the Doric console)
bonsai_setup = True

# low-pass filter characteristics to remove artefactual jittering
ORDER = 4
CUT_FREQ = 20 #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 2
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0.5
```

The corresponding analyses will be stored at this location :

```python
# Create repository path where fiberbehav data will be stored
exp = 'EPM'
exp_path = experiment_path / 'Analysis' / exp
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
```

In the above example, the photometry data will be low-pass filetered with a 20Hz butterworth filter with order 4. The behavioural bouts will be fused if they are less than 2 seconds apart. If, in the resulting bouts, some are shorter than 1/2 second, do not included it in further analysis.

These parameters have to be adapted to your task.

You have to generate at least one analysis with THRESH_S = 0 and EVENT_TIME_THRESHOLD = 0 for some of the further analyses (like behavioural metrics and mean dFF extraction during behavior).

#### BORIS analysis (Deprecated)

```text
scripts/script_fiberboris.py
```

This script combines fiber photometry data with manual behavioral scoring from BORIS.

#### Elevated plus maze

```text
scripts/script_fiberEPM.py
```

This script combines fiber photometry with DeepLabCut tracking data.

It uses manually defined maze boundaries to classify the animal's position in the maze.
Then it aligns fiber photometry with behavioural data.

The script generates:

* Aligned fiber photometry and behavioral data and plots :
```text
repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
repo_path / f'{batch}_{mouse}_fiberbehav.csv'
repo_path / f'{batch}_{mouse}_fiberbehav.pdf'
repo_path / f'{batch}_{mouse}_fiberbehav.png'
```
* Behavioral data : 
    * Occupation heatmaps and pie charts:
    ```text
    behavioural_analysis_path / 'Grouped figures' / f'Group_{group}' / figure names
    ```
    * Quantification of time in zones, plotting individual pie charts and heatmaps:
    ```text
    behavioural_analysis_path / 'behav_summary.xlsx'
    save_dir / f'{batch}_{mouse}_epm_pie_chart.png / pdf'
    save_dir / f'{batch}_{mouse}_epm_heatmap_{n_bins}bins.png / pdf'
    '''
* Fiberphotometry data : 
    * Mean dFF and AUC during behaviours, with z-scored and non z-scored data. Only closed arm is included in the calculation of F0 and std0 for z-scoring.
    ```text
    repo_path / 'dFF_summary_raw.xlsx'
    repo_path / 'dFF_summary_zscored.xlsx'
    '''


#### Fear conditioning

```text
scripts/script_fiberFC.py
```

This script analyzes fear conditioning experiments.

It detects freezing events based on dlc tracking data, then aligns fiber photometry with behavioral data (CS+ events, CS− events, Shock events, Freezing).

It also supports behavioral bout processing and dFF analysis.

#### Lick and airpuff analysis

```text
scripts/script_fiberlicks.py
```

This script combines fiber photometry with lick and airpuff recordings.

It processes:

* Lick events.
* Airpuff events.
* Animal position.
* Port entry.
* Approach behavior.
* Head orientation.

#### Plethysmography

```text
scripts/script_fiberplethysmo.py
```

This script combines fiber photometry with whole-body plethysmography data.

It supports:

* Loading plethysmography data.
* Identifying sniff events.
* Aligning photometry and respiratory signals.
* Generating plots.
* PETH analysis around sniff events.

## 6. PETH analysis

Peri-event time histograms are used to analyze photometry signals around behavioral events.

Main modules:

```text
modules/common/behavplot.py
modules/common/plethyplot.py
modules/common/correlation.py
```

Main scripts:

```text
scripts/script_PETH.py
scripts/script_jointPETH.py
scripts/script_correlation.py
```

### PETH analysis

```text
scripts/script_PETH.py
```

This script computes PETHs for individual animals and behavioral events.

It loads combined fiber-behavior data:

```text
{batch}_{mouse}_fiberbehav.csv
```

The script then:

1. Selects a behavioral event.
2. Identifies event onsets or offsets.
3. Extracts photometry windows around the event.
4. Computes the PETH.
5. Saves the resulting data and plots.

The time windows are configured in the script.

For example:

```python
TIME_WINDOWS = [[5, 10], [5, 10]]
```

These values define the pre-event and post-event windows.

### Joint PETH analysis

```text
scripts/script_jointPETH.py
```

This script performs PETH analysis across multiple animals.

It loads individual animal fiber-behavior files, computes PETHs, and combines them for group analysis.

### Correlation analysis

```text
scripts/script_correlation.py
```

This script performs correlation analysis between photometry signals.

It includes cross-correlation analysis between the 465 nm and 560 nm channels.

The correlation functions are implemented in:

```text
modules/common/correlation.py
```

The module computes normalized cross-correlations and supports baseline cross-correlation analysis.

## 7. Transient detection

Main module:

```text
modules/common/transients.py
```

The transient analysis module detects and quantifies changes in the dFF signal.

It includes:

* Bandpass filtering.
* Peak detection.
* Transient analysis.
* Signal visualization.

The main function for filtering is:

```python
bandpass_filter(...)
```

It estimates the sampling rate from the DataFrame and applies a Butterworth bandpass filter.

The filtered signal is stored in a column called:

```text
Filtered dFF
```

## 8. Quantification

Main script:

```text
scripts/script_quantification.py
```

Main module:

```text
modules/common/quantification.py
```

The quantification pipeline extracts signal measurements during behavioral events or zones.

The main function is:

```python
extract_dff_summary(...)
```

It calculates metrics such as:

* Mean dFF.
* Area under the dFF curve (AUC).
* Signal measurements during specific behavioral zones.
* Signal measurements during specific behaviors.

The AUC is calculated using:

$$
AUC = \frac{\sum dFF}{fps}
$$

where `fps` is the sampling rate.

The function can also calculate z-scored dFF values using a baseline period.


## Data formats

The repository uses several data formats.

| Format       | Use                                                       |
| ------------ | --------------------------------------------------------- |
| `.csv`       | Processed fiber photometry, behavioral, and tracking data |
| `.txt`       | Raw capacitance and behavioral recordings                 |
| `.doric`     | Doric fiber photometry recordings                         |
| `.json`      | Maze, port, and video calibration coordinates             |
| `.xlsx`      | Quantification results                                    |
| `.png`       | Figures                                                   |
| `.pdf`       | Figures and reports                                       |

## Author

Alice Fermigier

Written within NutriNeuro Lab - UMR INRAE 1286 (2020-2024) and Neurocentre Magendie's INSERM U1215 SNaP Lab (2024-2026)

Repository:
https://github.com/AliceFermigier/Fiberphotometry_analysis
