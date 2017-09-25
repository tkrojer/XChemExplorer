#!/bin/bash

echo "Please enter your fedID:"
echo " "
read fedid
echo "Now generating ssh key..."
echo " "
rm ~/.ssh/*
ssh-keygen -t rsa
ssh-copy-id $fedid@nx.diamond.ac.uk


echo "Please enter the full path to the directory where you want to mount the diamond file system:"
echo " "

read mountdir

string="sshfs -o Ciphers=arcfour -o Compression=no -o allow_other -o reconnect -o workaround=rename ${fedid}@nx.diamond.ac.uk:/dls/labxchem/data ${mountdir}"

echo $string > mount_diamond.sh

string2="cd ${mountdir}" 

echo $string2 >> mount_diamond.sh
