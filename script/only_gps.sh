# cd ..
# mkdir build
# cd build
# cmake ..
make

#  rooftop4
# ./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 100 --start-epoch 600 --T 100



#  monte-carlo simulation
# ./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 50 --start-epoch 1 --T 100


./ceres_solver_gnss --constellation gps --disable-clock-const --start-epoch 1 --T 100

