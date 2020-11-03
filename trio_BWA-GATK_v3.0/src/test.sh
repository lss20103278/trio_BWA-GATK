#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

#########################################################################
# File Name: src/test.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 22 Jan 2019 07:56:36 AM CST
# Usage:
#########################################################################

basedir=`cd $(dirname $0); pwd -P`
echo $basedir
