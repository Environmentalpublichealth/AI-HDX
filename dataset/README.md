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
