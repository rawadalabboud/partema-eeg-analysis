# 🧠 partema-eeg-analysis

This repository contains the code and structure for preprocessing and analyzing **EEG**, **voice**, and **questionnaire** data collected as part of the PARTEMA study on **rTMS treatment response in treatment-resistant depression**.

---

## 📁 Project Structure

```plaintext
partema-eeg-analysis/
├── data/                        # Raw data (not tracked by Git)
│   ├── eeg/                    # Raw .bdf files
│   ├── voice/                  # Audio recordings
│   └── questionnaires/        # Raw CSV or Excel responses
│
├── preprocessing/
│   ├── python/                 # Scripts for BDF → BIDS conversion using MNE
│   └── matlab/                 # EEGLAB / TESA / FieldTrip preprocessing scripts
│
├── analysis/
│   ├── eeg/
│   │   ├── python/             # Python-based EEG analysis (e.g., MNE)
│   │   └── matlab/             # MATLAB analysis scripts (TESA, FieldTrip)
│   ├── voice/                  # Acoustic & linguistic feature extraction
│   └── questionnaires/        # Statistical analyses
│
├── notebooks/                  # Jupyter or Live Scripts for exploration
├── results/                    # Output figures, tables, stats
├── requirements.txt            # Python dependencies
├── .gitignore                  # Files/folders excluded from version control
├── .gitkeep                    # Placeholder files to track empty folders
└── README.md                   # This file
⚙️ Getting Started

🔬 Requirements
Python ≥ 3.8
mne
mne-bids
numpy, pandas, scipy
MATLAB ≥ R2021a
EEGLAB
TESA
FieldTrip
You can install Python dependencies with:
pip install -r requirements.txt
🔄 Setup (local dev)
git clone https://github.com/rawadalabboud/partema-eeg-analysis.git
cd partema-eeg-analysis
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
🧪 Pipeline Overview

EEG Data
Collected in .bdf format
Converted to BIDS with mne-bids
Preprocessed using EEGLAB + TESA (MATLAB)
Later analysis using TESA, FieldTrip, or MNE
Voice Data
Acoustic and linguistic features extracted with Python (openSMILE, pyAudioAnalysis, etc.)
Questionnaires
CSV/XLS data from patient self-reports
Analyzed in Python or R (e.g., Likert scoring, correlations)
🧠 Study Context

This project supports the PARTEMA study, which explores the use of EEG biomarkers, voice phenotyping, and digital assessments to predict response to rTMS in treatment-resistant depression.

🤝 Collaboration

This work is part of a CIFRE PhD project between:

🏥 NeuroStim & Clariane (clinical data collection)
🧪 Paris Brain Institute (analysis and methodology)
