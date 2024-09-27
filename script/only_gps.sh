cd ..
mkdir build
cd build
cmake ..
make

./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 100 --start-epoch 646 --T 2
# ./ceres_solver_gnss --constellation gps --disable-clock-const --disable-tdcp --df-pr-weight 1 --start-epoch 630 --T 20

