#!/bin/bash
#$ -e $HOME/logs
#$ -o $HOME/logs

exec /usr/bin/perl $HOME/ChiLA/VSScripts/vse_run_sge.pl $@
