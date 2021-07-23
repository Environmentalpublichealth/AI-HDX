### public HDX datasets contain 60 min / 1 hour time point
- For all the exposure time in the files, if in second, keep 3600, if in min, keep 60, if in hour, keep 1.
- For all files, what should be output: "Sequence", "Maxuptake", "Exposure", "Update". There is an "Update SD" column, pay attension to it, we don't need uptake SD.
- We only care about stand alone protein for now, so here is what we want to keep in the "state":
    - `20181213_AriH2_StateData`: AriH2
    - `HDX_19884_kawasaki`: unbound
    - `HDX_19884_MIOO1`: unbound
    - `HDX_19884_Vietnam`: unbound
    - `MASSIVE-COMPLETE-9dffe2a9-display_quant_results-main`: default state
    - `state_19810`: apo
    - `State_data_Bcd1p_profile`: bcd
    - `State_data_M5_vs_M5-L18-RNA`: M5_alone
    - `State_data_RAR_profile_DR0_binding`: RAR_RXR
    - `State_data_RAR_profile_DR5_binding`: RAR_RXR
    - `State_data_Rtt106p(65-320)_profile`: Rtt
    - `State_data_RXR_profile_DR0_binding`: RAR_RXR
    - `State_data_RXR_profile_DR5_binding`: RAR_RXR
    - `TableS1_HDXMS_Data_Results_19199`: Unbound (in tab 'Uptake_summary') 

Need a python script to filter the tables.      

### Task 1
1. read csv file
2. select the desired state
3. select the desired time
4. output wanted columns into a new csv file

### Task 2
calculate MaxUptake for tables don't have it.
MaxUptake = (length of sequence) - 1 - (number of P)

### Task 3
Calculate tables not have 60 min or 1 hour or 3600 seconds.

### Task 4
We encoded the protein sequence from a different way. I used to use a fixed matrix to represent each AA acids in the sequence. Now, I changed to a methods used in alphafold, to align the entire protein against a database, which standardize the AA encoding. Then I will need the information of the full-length protein sequences for each of the data table.          
Please create a script that takes 'start' and 'end' positions and the short sequences, to assemble them into the full length. It is OK to have small gaps (ex. 1-7 amino acid + 10-20 amino acid) in the sequence. The alignment program can take care of it. The output could be in either these two formats:
```fasta
>20181213_AriH2_StateData
protein sequence MWEADAVSSVLSLLSILL
>HDX_19884_kawasaki
sequence
>filename
sequence
...
```
Or just put them as a table format like:
```txt
20181213_AriH2_StateData    sequence
HDX_19884_kawasaki  sequence
filename1    sequence
...
```
Some HDX data tables do not have start and end positions for the sequence, then we can't use them, so let's not worry about those.
