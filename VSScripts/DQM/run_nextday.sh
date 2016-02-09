#!/bin/sh

CHILADIR=/home/mdwood/ChiLA /home/mdwood/raid/nextday/run_nextday.pl /home/mdwood/raid/nextday/list.txt /home/mdwood/raid/nextday/cfg_2009/mlm432_pl35_2tel_s300.txt /home/mdwood/raid/nextday/cfg_2009/mlm432_pl25_2tel_s550.txt /home/mdwood/raid/nextday/cfg_2009/rbm432_pl35_2tel_s300.txt --lck  1>> /home/mdwood/raid/nextday/logs/log.`date +%Y%m%d`.txt 2>> /home/mdwood/raid/nextday/logs/log.`date +%Y%m%d`.txt 
