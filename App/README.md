### Python library used
* python = 3.11
* shiny = 0.2.9
* tensorflow = 2.12.1
* numpy = 1.25.0
* pandas = 1.5.3
* jinja2 = 3.1.2

Create a virtual environment, use either pip or conda install the required libraries
```bash
conda install shiny
pip install tensorflow==2.12.* # we don't need cuda, because we are not using GPU
pip install pandas
pip install Jinja2 # a package required for table output
```
or
```bash
conda install -c conda-forge shiny
conda install -c conda-forge tensorflow
conda install pandas
conda install -c anaconda jinja2
```

### Dependency python programs
* `embedding.py` - a sequencing embedding function used in the model
* `AIHDX.py` - the function to run AI-HDX model
  * need to manually change the path to find pre-trained models now.
  * we should make it read the user file system and find path on its own. (to do)

### AI-HDX app
* `app.py` - the shiny app built to run AI-HDX model        
Make sure the dependencies are in the same folder of the `app.py`.

### Pre-trained models
The `models.zip` in this repository. Need to download and unzip it to use. 
