#!/bin/bash
# make a tarball of maphys
dir=$(dirname $0)
v=$( head -n1 ${dir}/../VERSION )

# Assumes
if git rev-parse 2>/dev/null ; then # is under git

    b=$( git status | head -n1 | cut -d ' ' -f 1,2,3  --complement)
    git archive --format=tar --prefix=maphys-${v}-${b}/ HEAD | gzip > maphys-${v}-${b}.tar.gz 

elif [ -d ${dir}/../.svn ] ; then # is under svn

    d=$PWD
    b=$( basename $(readlink -e ${dir}/..))
    [ -d /tmp/maphys-${v}-${b}/ ] && rm -rf /tmp/maphys-${v}-${b}/
    svn export ${dir}/.. /tmp/maphys-${v}-${b}/
    cd /tmp/ && tar zcvf ${d}/maphys-${v}-${b}.tar.gz maphys-${v}-${b}
    rm -rf /tmp/maphys-${v}-${b}/

else

    echo "version control not supported, cannot make tarball"
    exit 1

fi

