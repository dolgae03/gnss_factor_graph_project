# cd ..
# mkdir build
cd build
# cmake ..
make



# ./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v1" --tau 0 --disable-tau --label "TEST"
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 0  --disable-tdcp --disable-tau --label "TEST"
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50  --disable-tdcp --disable-tau --label "TEST" --clock-const-weight 100.00
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50  --disable-tdcp  --clock-const-weight 100.00  --label "TEST"
./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50   --disable-tdcp --label "TEST" --clock-const-weight 10.00
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50  --disable-tdcp  --label "TEST"  --disable-clock-const 