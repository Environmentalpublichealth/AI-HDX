### Install HHBlits via conda
```bash
module load Anaconda3/2020.07
conda create -n hhblits
source activate hhblits
conda install -c conda-forge -c bioconda hhsuite
```
Download database
```bash
mkdir databases
cd databases
nohup wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz &
tar -xvfz UniRef30_2020_06_hhsuite.tar.gz
```
The download takes ~6 hours, because it is connected to a German server.

### Run HHBlits
```bash
module load Anaconda3/2020.07

source activate hhblits

cd /scratch/user/jialiyu/BlitSearch
for file in ./*fasta
do
BASE=$(basename $file | sed 's/.fasta//g')
hhblits -i $file -ohhm $BASE.hhm -d databases/UniRef30_2020_06/UniRef30_2020_06
mv $BASE.hhm hhm_data
done
```
Each protein is stored in a single fasta file, so I used a for loop to run all files.

Transfer `.hhm` files to local.
```bash
rsync -a grace.tamu:/scratch/user/jialiyu/BlitSearch/hhm_data .
```
