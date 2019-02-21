#launch many simulations of rising droplet

echo "Start submission"

#go to main direcotory
cd ~/Documents/C++/edge_state

#compile source code
make

#cd ~/work/scratch/edge_state

COUNTERshape=0
while [  $COUNTERshape -lt 10 ]; do

echo "Beginning of loop on shape $(($COUNTERshape + 1))"

let COUNTERshape=COUNTERshape+1
COUNTERround=0

while [ $COUNTERround -lt 100 ]; do

echo "Beginning of loop on round $(($COUNTERround + 1))"

let COUNTERround=COUNTERround+1
COUNTERdelta=0
while [  $COUNTERdelta -lt 3 ]; do

        let COUNTERdelta=COUNTERdelta+1
        nohup ./edge_state.o $COUNTERshape $COUNTERdelta $COUNTERround>edge$COUNTERdelta.out 2>&1 &

done

wait

done

#wait

done

#nohup ./rising_droplet.o $COUNTER>drop$COUNTER.out 2>&1 &

echo "End of submission"
