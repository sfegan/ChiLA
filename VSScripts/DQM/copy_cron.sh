#!/bin/bash

#
# Copy files from primary data directory to scratch. Only copy files
# that have not been modified since the last time we executed.
#
# Stephen Fegan -- 2006-09-20
#

# $Id: copy_cron.sh,v 1.1 2008/06/06 15:44:55 sfegan Exp $

# Variables
utc=`date -u +%Y%m%d`
srcdir=/net/gamma3/srv/raid7/data/d${utc}
dstdir=/home/sfegan/scratch/VBF/current
lastrun=${dstdir}/.last_iteration
utcfile=${dstdir}/.utc_${utc}
lockfile=${dstdir}/.lockfile

# Make destination directory if it does not exist
if test ! -d ${dstdir}
then
	mkdir ${dstdir}
fi

# Initialization
if test ! -e ${lastit}
then
	touch -t 197001010000 ${lastrun}
fi

# Test if the source directory is there
if test -d ${srcdir}
then
	# Lock the lockfile
	lockfile -l 3600 ${lockfile}

	# Delete old data
	if test ! -e ${utcfile}
	then
		rm -f ${dstdir}/.utc*
		rm -f ${dstdir}/*.cvbf
		touch ${utcfile}
	fi

	# Copy data
	for srcfile in ${srcdir}/*.cvbf
	do
		basefile=`basename ${srcfile}`
		dstfile=${dstdir}/${basefile}
		if test ${srcfile} -ot ${lastrun} -a ! -e ${dstfile}
		then
			cp ${srcfile} ${dstfile}
			touch ${lockfile}
		fi
	done
fi

# Update time guard
touch ${lastrun}
rm -f ${lockfile}
