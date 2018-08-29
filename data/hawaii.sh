#!/bin/bash
rm b18t1800.eps
#pscoast -R-161.76/-151.69/14.21/23.07 -JM6i -B2/2 -K -W0.5p >> b18t1800.eps
#pscontour b18T1800.dat  -R -JM -K -O -Cpac.cpt -I >> b18t1800.eps
psxy -R -JM -K -O -Sc0.10 -G0 b18T1800.dat >> b18t1800.eps
