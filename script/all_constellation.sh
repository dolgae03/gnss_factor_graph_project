cd ..
mkdir build
cd build
cmake ..
make


# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 100 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 10 --start-epoch 600 --T 100

# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 5 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 4 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 3 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 2 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 1 --start-epoch 600 --T 100

# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 0.9 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 0.8 --start-epoch 600 --T 100
# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --df-pr-weight 50 --tdcp-weight 0.7 --start-epoch 600 --T 100


# ./ceres_solver_gnss --constellation gps bds gal --disable-clock-const --disable-tdcp --df-pr-weight 50 --start-epoch 600 --T 100

./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.006 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.005 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.004 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.003 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.002 --start-epoch 600 --T 100
./ceres_solver_gnss --constellation gps bds gal --disable-tdcp --df-pr-weight 50 --clock-const-weight 0.001 --start-epoch 600 --T 100