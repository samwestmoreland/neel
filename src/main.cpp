#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "./data.hpp"
#include "./io.hpp"
#include "./neel.hpp"

int main(int argc, char* argv[]) {

   read_inputfile();

   /* read in unit cell coordinates */
   read_coordinates();

   /* generate super cell for periodic boundaries */
   generate_supercell();

   populate_atom_eij_tensors();

   /* calculate neel energies using vector method */
   // calculate_k();

   /* calculate easy axis for each atom */
   calculate_easy_axes();

   return EXIT_SUCCESS;

}
