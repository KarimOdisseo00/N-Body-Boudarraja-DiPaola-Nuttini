#ifndef CONSTANTS_H
#define CONSTANTS_H

//Dimensione del problema 1D 2D 3D
const unsigned int dim=3;

//Costante Gravitazionale Universale G


//Delta temporale espresso in secondi
const double dt=1;
        
//Tempo totale della simulazione espresso in secondi
const double totalTime = 10;

//Numeri di cicli
const int cycles = static_cast<int>(std::floor(totalTime/dt));

//Numero di particelle generate
const int numberOfParticles = 100;  

#endif  // CONSTANTS_H
