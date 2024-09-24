cd ..
mkdir build
cd build
cmake ..
make

./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 1 --tdcp-weight 10 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps --disable-clock-const --disable-tdcp --df-pr-weight 1 --start-epoch 600 --T 100