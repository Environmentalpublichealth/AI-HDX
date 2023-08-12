### Deploy the codes as a web tool
Weekly meeting: 7 pm -7:30 pm CT Thursday June 8th to June 22th.

1. find a platform: wix.com


2. Try simple, example codes
3. Build our own!         

Type what you want to write. 

First step:
* user need to prepare their own inputs, and can upload their input files
  * Text boxes to paste input
  * A bottom to upload files 
* user run the code by submitting their inputs
  * a submit bottom
  * Redirect to a result page
  * a download bottom

Second step:
Users still need to prepare their own inputs.          
* input a protein structure file in PDB format
* The result page can draw the structure and match the dynamic data on the structure


### Getting shiny install and run
1. need to have python installed
2. intall shiny either way: mannual in here https://shiny.posit.co/py/docs/install.html
```bash
pip install shiny
```
or 
```bash
conda install conda-forge shiny
```
3. Test if shiny runs
```bash
shiny create .
```
This will creat a `app.py` in the current working directory. If you see this file in your folder, then do
```bash
shiny run --reload
```
This will output a bunch of test in your terminal, it will not automatic open a webpage in your browser.
```basth
INFO:     Will watch for changes in these directories: ['/Users/jiali/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/test']
INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [78432] using StatReload
INFO:     Started server process [78434]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
To make it in the browser, you need to copy the url 'http://127.0.0.1:8000' it generated and paste it into the web browser you use.

4. Use visual studio code to have your code and the webpage side by side (optional).
- Download and install VS code: https://code.visualstudio.com/. It is free! Don't search visual studio, it is a different app.
- install shiny extension in VS code: View -> extension, search 'shiny for python' and click install.
- run `shiny run --reload` in your terminal (not in the VS code terminal!)
- Get simple browser in VS code: View -> Command Pallette..., search 'simple browser' and choose it. Paste the URL shiny generated in the box and hit enter. A browser window will be open in the VS code.
- click the 'split editor right' button at the upper right corner, like '[|]' shape, and you can put the shiny app in simple browser on one side and open the app.py on the other side at the same time. The elements in the browser should be changing when editing the codes and save in app.py. 

5. Stop shiny       
Do control+C to end shiny in the terminal.

### Meeting schedule
* biweekly - Saturday, 2:30 pm EST

### Functions needed in app
1. input windows
   - :white_check_mark:one for the vector.txt (tab separated text file)
   - :white_check_mark:one for fragment.csv (comma separated text file)
   - :white_check_mark:user name text box
   - future feature: input protein sequences instead of the fragment, write a function to cut proteins and write fragments into csv

:white_check_mark:2. save the input files with username and time stamp - Anna's idea

3. migrate the codes from the AI-HDX model google colab to app.py
   - :white_check_mark:add a submit button. In the button action:
   - :white_check_mark:one function to preprocess the two inputs
   - :white_check_mark:one function to run machine learning models
   - :white_check_mark:Install tensorflow and need to make tensorflow compatible with shiny
   - A download button to download result
5. display the results
   - :white_check_mark:one function to simply print out the table
   - one function to visualize results on protein structure - need to find out
6. Plotting
   - :white_check_mark:one function to plot the result into a line/bar/scatter plot
   - integrate the plot to shiny
7. Web design
   - Font, size of the text
   - color of background
   - color of the button
   - color of the text box and input window
   - descriptions and instructions for the app usage
   - copyright, authorship claim
   - warnings about prediction accuracy
8. Auto-fill examples
   - a button to load example inputs
 
