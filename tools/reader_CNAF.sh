#!/bin/sh
# PATCH FOR BX_READER
# given by Alessandro Razeto (2012-11-30)

files=`echo $* | sed s.http://..g`
cat $files
