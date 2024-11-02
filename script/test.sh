# cd ..
# mkdir build
cd build
# cmake ..
make



# ./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v1" --tau 0 --disable-tau --label "TEST"
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v1" --tau 0 --disable-tau --disable-clock-const 


# 
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 0  --disable-tau --disable-imu --disable-clock-const
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 0  --disable-tau --disable-tdcp 
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 0  --disable-tau --disable-clock-const
# ./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50  
./ceres_solver_gnss  --T 100 --dataset "constSig_v3" --tau 50  --disable-tdcp --disable-imu


# --disable-clock-const

