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

### Task 5
Build a simple model!            
I am working on some really complex deep learning models to learn the sequence and HDX rates. What you can do now, is to help me build a very very simple model as a prediction baseline!
* 1. Convert a peptide sequence into a vector. The machine learning model does not take text, so we need to encode protein sequence into a set of numbers. Here, let's use the simplest way: protein sequences is consist of 20 amino acides, we just count how many [A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y] in each peptide, and make it into a vector. For example, a peptide sequence 'GTAGSAEEPS' should encode as [2,0,0,2,0,2,0,0,0,0,0,0,1,0,0,2,1,0,0,0]. Do this encoding for all peptides in the HDX tables.
* 2. For the selected HDX result tables, if a D% is given, use that number, if not, calculate the D% with this equation `D% = uptake / (Maxuptake * 0.7) * 100`. To make things easy, we just use a constent back-exchange rate for all tables for now. But remember it is not always true for all experiments, back-exchange is ranged from 10% to 40%, we take the middle value for now.
* 3. Remove the overlap sequence, to avoid overfitting issue. I created a function for that, feel free to use it or you can creat one yourself.
```python
def removeOverlap(start_list, end_list):
  uniq_start = []
  uniq_end = []
  uniq_idx = []
  for idx1, val1 in enumerate(start_list[:-1]):
    overlap = False
    uniq_idx.append(idx1)
    uniq_start.append(val1)
    uniq_end.append(end_list[idx1])
    for val2 in start_list[idx1+1:]:
      if val1 >= val2 and end_list[idx1] <= end_list[idx1+1]:
        overlap = True
        break
    if overlap == True:
        uniq_start.pop()
        uniq_end.pop()
        uniq_idx.pop()
  return uniq_start, uniq_end, uniq_idx
```
`start_list` - a list of start positions for the peptides.         
`end_list` - a list of end positions for the peptides.       
Outputs: a list of unique start positions, a list of end positions and a list of the index, for filtering labels later on.      
* 4. Now start to build a model!
    - You need to know, the 20 numbers vector is your input, and the HDX rate is your label. Remember to convert the percentage numbers into (0, 1).
    - To do the training, we need to split the data into training set and validation set, follow the 7:3 rule. Randomly select 70% of the data as training and 30% of the data as testing. (Tip, use `train_test_split()` function in `sklearn.model_selection`)
    - Select a classifier and train it! A simple classifier can be `sklearn.svm.SVR`. 
    - Evaluate your prediction, by R squred and mean_squared_error. (Tips: `r2_score()` and `mean_squared_error()` functions in sklearn.metrics)
    - Ajust the classifier parameters to get a performance as good as possible, need to keep trying different parameters.
