#!/bin/bash
rm pacific.eps
pscoast -R-10.0/30.0/-30.0/24 -JM -B -K -W0.5p >> pacific.eps
#pscontour b18T1800.dat  -R -JM -K -O -Cpac.cpt -I >> b18t1800.eps
#psxy -R -JM -K -O -Sc0.10 -G0 b18T1800.dat >> b18t1800.eps
