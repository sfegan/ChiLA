#!/bin/bash
#$ -e $HOME/logs
#$ -o $HOME/logs

exec /usr/bin/perl $HOME/ChiLA/VSScripts/vs_run_sge.pl $@
