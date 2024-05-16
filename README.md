# Useful_Code
Code for specific task that save time

## How to mount and external hard drive in a WSL using ubuntu
### First, Insert drive in the Pc and check the letter assignd (E,G,F, etc) 
### Second, open your terminal and type 

``` 
cd /mnt
sudo mkdir g (or what ever letter your disk is, in this example i am suing G)
sudo mount -t drvfs g: /mnt/g
``` 
And that it!
