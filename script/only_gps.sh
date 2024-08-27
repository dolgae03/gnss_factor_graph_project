cd ..
mkdir build
cd build
cmake ..
make

./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 50 --tdcp-weight 10 --start-epoch 600 --T 40
./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 50 --tdcp-weight 1 --start-epoch 600 --T 40
./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 50 --tdcp-weight 0.01 --start-epoch 600 --T 40
./ceres_solver_gnss --constellation gps --disable-clock-const --df-pr-weight 50 --tdcp-weight 0.001 --start-epoch 600 --T 40