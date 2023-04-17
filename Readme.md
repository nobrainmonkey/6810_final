# 2-D Ising Model Simulation 

Author: Xihe Han

This program can graph the spin-flipping process of a micro state and calculate the related macroscopic quantities (energy, specific heat, entropy, magnetization, magnetic susceptibility) with plotting capabilities following the metropolis algorithm for the 2-dimensional Ising model problem. 

## Dependencies:

This project requires:

1. Functional `C++` compiler that support `C++11` standard. 
2. `Gnu Make` for using the make file.
3. [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), a BLAS library for `C++`.
4. `OpenMP` for parallel computing.
5. `gnuplot` for plotting.

### Installing `Eigen`:

#### Arch Linux:

```
#	pacman -S eigen
```

#### Ubuntu:

```
#	apt-get install libeigen3-dev
```

#### MacOS:

```
$	brew install eigen
```

Note: `Eigen` is a header file only library.

 ## Building:

1. Clone the following GitHub repository and `cd` into the source files' location.

   ```
   $	git clone https://github.com/nobrainmonkey/6810_final && cd 6810_final/src
   ```

   

2. Use `Gnu Make` with the makefile `make_main` to build the project.

   ```
   $	make -f make_main
   ```

The project should now be built in the directory `6810_final/bin`.

### Potential Build Error:

#### Eigen Error:

If you have not downloaded `Eigen` from your package manager, it would not be placed in the system default `include` directory.

To fix this, either install `eigen` from your package manager, or do the following: 

1.  Locate the directory of `eigen` source file
2. Use the best editor- vim to open `/6810_final/src/make_main`
3. Edit line 6 , put the **absolute** directory of the `eigen` source file after `=`
4. Change line 9 to `-I (EIGEN_PATH) FLAGS= -Wall -O3 -fopenmp -std=c++11`
5. Make the project again.

## Usage:

To run the program, run the following command:

```
$	./calculate.x
```

A CLI interface should prompt you to enter simulation values. Here are some explanation of each option. The user will input the value they want to modify with their corresponding numbered label into the terminal emulator.

### Macroscopic Calculation:

This mode returns your a `.dat` file containing every macroscopic quantities' calculation. 

1. `thread` determines how many physical threads does the program utilize.
2. `sample` determines the number of Metropolis algorithm to perform. The final data will be an average over all samples. For good-quality graph, make `sample` to be 1000.  
3. `iterations` determines the number of spin flips for our system to perform in order to reach thermal equilibrium. **This value will be suggested automatically** to be the size of the micro state * 1000.
4. `J` is the spin interaction parameter of the Ising model Hamiltonian.
5. `h` is the field interaction parameter of the Ising model Hamiltonian.
6. `row` determines the number of rows of our spin system.
7. `col` determines the number of columns of our spin system. This does not need to equal to `row`.
8. `Tmin` determines the minimum temperature in which we start our evaluation. I would highly recommend not to decrease this value as for extremely low temperature, our micro state will tend to a local minimum more often.
9. `Tmax` determines the maximum temperature where we end our evaluation.
10. `dT` determines the temperature steps starting from `Tmin` to `Tmax`. 

After you are satisfied with the values, enter `0` to run the calculation. The data will be stored under `/6810_final/data`. A sample data is provided for a 20 x 20 system. under `/6810_final/good_data`.

After we each run, a window of `gnuplot` should automatically pop up with every graph of the calculated values with their corresponding labels.

Sample plot:

![sample_plot](/home/xihe/6810_final/sample_plot.svg)

System size = 10 x 10, sample size = 1000, Tmin = 0.5, Tmax= 6, dT=0.1

### Graphing Mode:

This mode gives you a terminal based visualization of how a spin system evolve using the metropolis algorithm.

#### Caution:

1. Currently,  this mode only work on Unix-like system. If you run a NT system, please refer to `/6810_final/src/microstate.cpp` and look for the `graph_evolve` method. I have provided instructions of how to make it work on NT systems.
2. To draw the spins correctly, make sure your terminal emulator have `UTF-8` encode support. 
3. Because different terminal emulators handle `clear` differently, your result may vary depending on  your terminal emulator. Note that there could potentially be **flashing image** for some terminal emulators. If you experience this, set `fps` to be around 10 might help.
   * Known working terminal emulators:
     * `Konsole`
     * `Alacritty` (flashing on high frame rate)
     * `Kitty`
   * Known none-working terminal emulators:
     * `gnome-terminal` 

1. The `row` and `col` variables are inherited from the calculation mode.
2. `graphing Temperature` determines the temperature where we bring our system to thermal equilibrium.
3. `graphing time` is the time in second where we graph the evolution of our spin system.
4. `graphing iteration per frame` determines the number of iteration passed each frame.
5. `graphing fps` determines how many frames we output on our terminal per second
6. `start graphing!` enters the graphing mode.

If for some reasons you cannot run this mode, here is what it looks like in real time for a 20 x 20 system:

 ![evolve](/home/xihe/6810_final/evolve.gif)

For performance and error analysis, please refer to `/6810_final/docs`. 

