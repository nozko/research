#!/bin/sh

interval="1"

python road_simulation.py $interval p
python road_simulation.py $interval s
python road_simulation.py $interval 1
python road_simulation.py $interval 2.0
python road_simulation.py $interval 2.1
python road_simulation.py $interval 2.2
python road_simulation.py $interval 3.0
python road_simulation.py $interval 3.1
python road_simulation.py $interval 4.0
python road_simulation.py $interval 4.1
python road_simulation.py $interval 5.0
