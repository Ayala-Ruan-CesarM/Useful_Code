# Useful_Code
Code for specific task that save time

## How to mount and external hard drive in a WSL using ubuntu
### First, Insert drive in the Pc and check the letter assignd (E,G,F, etc) 
### Second, open your terminal and type 

``` 
cd /mnt/
sudo mkdir g (or what ever letter your disk is, in this example i am suing G)
sudo mount -t drvfs g: /mnt/g
``` 
And that's it!


## How to perform a task inside various folders with a bash script
### It requires to know the path to the folders 
``` 
folder1="PATH/TO/FOLDER1"
folder2="PATH/TO/FOLDER2"
folder3="PATH/TO/FOLDER3"
# Create the list with the folders as variables
folders_list=("$folder1" "$folder2" "$folder3")
# for cycle to iterate over each folder

for folder in "${folder_list[@}": do
  cd "$folder" || exit 1 # To check whether the folder exist
  # Do something. For example, create a file in each one of the three folders
  echo hello world > test.txt
  # this will create a file containing the words "hello world" in each folder 
done
``` 

## How to remove a set of sequences from a multi fasta file by its header
### It require to have the header IDs of the sequences without the ">" in a textfile
### This python script has three arguments and BioPython

* Sequence_1  
* Sequence_156  
* Sequence_11684 
* Sequence_20546
* Sequence_553218 

```
python Remove_Seqs_By_Header.py input.fasta headers.txt > output.fasta
``` 
## How to extract a set of sequences from a multi fasta file by its header
### It require to have the header IDs of the sequences without the ">" in a textfile

```
seqkit grep -f headers.txt -o Wanted_seqs.fasta input.fasta
``` 
