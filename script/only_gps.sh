# cd ..
# mkdir build
cd build
# cmake ..
make

#  rooftop4
# ./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 100 --start-epoch 600 --T 100


# ./ceres_solver_gnss  --tau 50 --df-pr-weight 0.5 --tdcp-weight 1250 --disable-clock-const --T 100 
# ./ceres_solver_gnss  --tau 50 --df-pr-weight 0.5 --disable-tau --disable-clock-const --T 100 


./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v5" --tau 50  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v5" --tau 50 --disable-tau
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v5" --tau 0  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v5" --tau 0 --disable-tau

./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v6" --tau 50  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v6" --tau 50 --disable-tau
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v6" --tau 0  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v6" --tau 0 --disable-tau

./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v7" --tau 50  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v7" --tau 50 --disable-tau
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v7" --tau 0  
./ceres_solver_gnss --disable-clock-const --T 100 --dataset "constSig_v7" --tau 0 --disable-tau



# --disable-df-pr --df-pr-weight 0.5
# --disable-tdcp --tdcp 0.5
# --disable-clock-const --clock-const-weight 0.5
# --disable-tau --tau-weight 25.51
