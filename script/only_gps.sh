# cd ..
# mkdir build
# cd build
# cmake ..
make

#  rooftop4
# ./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 100 --start-epoch 600 --T 100



#  monte-carlo simulation
# ./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 50 --start-epoch 1 --T 100


# ./ceres_solver_gnss --constellation gps --disable-tau --start-epoch 1 --T 100

# ./ceres_solver_gnss --constellation gps --tau-weight 100 --start-epoch 1 --T 100  
# ./ceres_solver_gnss --constellation gps  --start-epoch 1 --T 100  # time correlation factor 고려

# ./ceres_solver_gnss --disable-tau --tau 0
# ./ceres_solver_gnss --disable-tau --tau 100
./ceres_solver_gnss --tau 100