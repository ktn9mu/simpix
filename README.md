Work in progress: 
g++ -O3 -std=c++17 simpix.cpp $(root-config --cflags --libs) -lASImage -lGui -o simpix

code:
./simpix imageA.png imageB.png out_AB.png collage_AB.png 20000000
./simpix imageB.png imageA.png out_BA.png collage_BA.png 20000000

code:(slurm)
sbatch run_simpixAB.slurm
sbatch run_simpixBA.slurm

time A → B is 1756.689s
time B → A is 1857.061s





# simpix

C++ starter code
* simpix_start.cpp
use make to build this example

Usage: simapix_start image1 image2 <output=out.png>

Python starter code
* simpix_start.py

Usage: simapix_start image1 image2 <output=out.png>


how to run:

Build:
c++ -O3 -std=c++17 simpix.cpp \
  $(root-config --cflags) \
  -L$(root-config --libdir) \
  $(root-config --libs) \
  -lASImage \
  -o simpix

RUn:
./simpix A.png B.png out.png --batch

