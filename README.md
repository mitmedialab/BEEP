# BEEP

## Base Editing Evaluation Program

Assessing CRISPR-mediated base editing efficiency from Sanger sequencing ab1 files

### Dependencies

1. pandas
```
pip install pandas
```
2. Biopython
```
pip install Biopython
```

### What You Need

1. beep.py
2. Negative control of amplified region: control.ab1 file 
3. Sample ab1 files
4. .csv file containing ab1 file names, target sequences, base position in target, and desired base conversion (template provided)

### Installation
```
git clone https://github.com/mitmedialab/BEEP
```

### Running BEEP

To run multiple samples at once:
```
python beep.py ./example/template.csv
```
To run a single sample:
```
python beep.py folder_with_ab1s control.ab1 sample.ab1 target_sequence base_position base_conversion
```
Example:
```
python beep.py example control.ab1 sample.ab1 TCGGCCACCACAGGGAAGCT 6 CT
```
### Output
* Efficiences will be outputted to the provided .csv file under "Efficiency" header, or in the terminal for single sample usage. 
* .png files are created with edits visualized as chromatograms.

### Author

**Pranam Chatterjee** 