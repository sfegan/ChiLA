#!/bin/bash
/usr/bin/ssh -v -o TCPKeepAlive=true -L login4:3306:gamma5.astro.ucla.edu:3306 localhost -N &> ssh_log.`date +%Y%m%d%H%M%S`
echo ssh died with signal $? > ssh.`date +%Y%m%d%H%M%S`
