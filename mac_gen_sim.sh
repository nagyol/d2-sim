#/bin/bash

#ID=`date -Isecond`.data

g++-13 -O3 ./main.cpp
./a.out > ./data/${@}
./plot.sh ./data/${@}
mv ./data/${@}.pdf ./pdfs/

