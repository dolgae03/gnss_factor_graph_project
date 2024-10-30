# cd ..
# mkdir build
cd build
# cmake ..
make



# ./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v1" --tau 0 --disable-tau --label "TEST"
./ceres_solver_gnss  --T 100 --dataset "constSig_v1" --tau 0 --disable-tau --label "TEST"



