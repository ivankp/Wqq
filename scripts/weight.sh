#!/bin/bash

./bin/weight_dist -i data/jul17.lhe -o weighted/opang/energy_cm/nocut2/ -v energy-cm -b 100:1000:25
#./bin/weight_dist -i data/wwnocuts_v2.lhe -o weighted/opang/energy_cm/nocut/  -v energy-cm  -b 100:1000:25

#./bin/weight_dist -i data/wwnocuts_v2.lhe -o weighted/opang_norm/energy/nocut/ -v energy -b 100:1000:25 --div_by_min -p
