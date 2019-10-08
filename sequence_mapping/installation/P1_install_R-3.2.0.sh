#!/bin/bash -e

echo " Start compile R-3.2.0 (will take ~3 min)"

cd /home/jh7x3/CaTrace2Seq//tools/R-3.2.0

make clean

./configure --prefix=/home/jh7x3/CaTrace2Seq//tools/R-3.2.0  --with-readline=no --with-x=no

make

make install

echo "installed" > /home/jh7x3/CaTrace2Seq//tools/R-3.2.0/install.done

