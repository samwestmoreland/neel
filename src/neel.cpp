#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "./data.hpp"
#include "./classes.hpp"
#include "./io.hpp"

int generate_supercell() {

   std::cout << "replicating unit cell..." << std::endl;;

   int global_id = 0;

   /* populate supercell using unitcell info */
   for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
         for (int k=0; k<3; ++k) {

            for (int atom=0; atom<unitcell.size(); ++atom) {

               atom_t tmp;

               tmp.type = unitcell[atom].type;

               tmp.gid = global_id;
               global_id ++;

               tmp.id = unitcell[atom].id;

               tmp.pos.x = unitcell[atom].pos.x + i * mat.ucd.x;
               tmp.pos.y = unitcell[atom].pos.y + j * mat.ucd.y;
               tmp.pos.z = unitcell[atom].pos.z + k * mat.ucd.z;

               supercell.push_back(tmp);

            }

         }
      }
   }

   return EXIT_SUCCESS;

}

int populate_atom_eij_tensors() {

   std::cout << "populating eij tensors..." << std::endl;

   /* determine beginning and end of central cell */
   int start = (supercell.size() - unitcell.size())/2;
   int end = (supercell.size() + unitcell.size())/2;

   double pair_count = 0;

   /* loop through central cell */
   for (int i=start; i<end; ++i) {

      /* loop through all atoms in supercell */
      for (int j=0; j<supercell.size(); ++j) {

         /* create a pair with these two atoms */
         pair_t temp_pair;

         temp_pair.atomi = supercell[i];
         temp_pair.atomj = supercell[j];

         if ((temp_pair.rij() <= sim.rcut) && (!temp_pair.are_same_atom())) {

            vec_t eij = temp_pair.eij();

            unitcell[i-start].ktensor[0] = eij.x * eij.x * temp_pair.lij();
            unitcell[i-start].ktensor[1] = eij.x * eij.y * temp_pair.lij();
            unitcell[i-start].ktensor[2] = eij.x * eij.z * temp_pair.lij();
            unitcell[i-start].ktensor[3] = eij.y * eij.y * temp_pair.lij();
            unitcell[i-start].ktensor[4] = eij.y * eij.z * temp_pair.lij();
            unitcell[i-start].ktensor[5] = eij.z * eij.z * temp_pair.lij();

         }

      }
   }

   return EXIT_SUCCESS;

}

int calculate_k() {

   std::cout << "calculating neel anisotropy using vector method...\n\n";

   for (int atom=0; atom<unitcell.size(); ++atom) {

//      if (unitcell[atom].is_re()) {

         std::ostringstream id_string;
         id_string << atom;

         std::string filename = unitcell[atom].type + id_string.str() + ".dat";

         std::ofstream fout (filename);

         /* polar angle */
         for (int phi=0; phi<=180; phi+=sim.angle_increment) {

            /* azimuthal angle */
            for (int theta=0; theta <=360; theta+=sim.angle_increment) {

               /* convert angle to radians */
               double phi_rad = phi / 180.0 * 3.14159;
               double theta_rad = theta / 180.0 * 3.14159;

               /* convert to cartesian */
               vec_t cart;
               cart.x = cos(theta_rad) * sin(phi_rad);
               cart.y = sin(theta_rad) * sin(phi_rad);
               cart.z = cos(phi_rad);

               vec_t temp_vector;

               /* col 1 */
               temp_vector.x = cart.x * unitcell[atom].ktensor[0]
                             + cart.y * unitcell[atom].ktensor[1]
                             + cart.z * unitcell[atom].ktensor[2];

               /* col 2 */
               temp_vector.y = cart.x * unitcell[atom].ktensor[1]
                             + cart.y * unitcell[atom].ktensor[3]
                             + cart.z * unitcell[atom].ktensor[4];

               /* col 3 */
               temp_vector.z = cart.x * unitcell[atom].ktensor[2]
                             + cart.y * unitcell[atom].ktensor[4]
                             + cart.z * unitcell[atom].ktensor[5];

               double k = temp_vector.x * cart.x +
                          temp_vector.y * cart.y +
                          temp_vector.z * cart.z;

               fout << phi << "\t"
                    << theta << "\t"
                    << k << "\n";

            }

            fout << "\n";

         }

         std::cout << "data output to " << filename << std::endl;

//      }
   }


   return EXIT_SUCCESS;

}

int calculate_easy_axes() {

   std::cout << "calculating neel anisotropy using vector method...\n\n";

   std::ofstream fout ("easy_axes.dat");

   for (int atom=0; atom<unitcell.size(); ++atom) {

      double min_energy = 10000;
      double min_phi;
      double min_theta;

         /* polar angle */
         for (int phi=0; phi<=180; phi+=sim.angle_increment) {

            /* azimuthal angle */
            for (int theta=0; theta <=360; theta+=sim.angle_increment) {

               /* convert angle to radians */
               double phi_rad = phi / 180.0 * 3.14159;
               double theta_rad = theta / 180.0 * 3.14159;

               /* convert to cartesian */
               vec_t cart;
               cart.x = cos(theta_rad) * sin(phi_rad);
               cart.y = sin(theta_rad) * sin(phi_rad);
               cart.z = cos(phi_rad);

               vec_t temp_vector;

               /* col 1 */
               temp_vector.x =
                  cart.x * unitcell[atom].ktensor[0] +
                  cart.y * unitcell[atom].ktensor[1] +
                  cart.z * unitcell[atom].ktensor[2];

               /* col 2 */
               temp_vector.y =
                  cart.x * unitcell[atom].ktensor[1] +
                  cart.y * unitcell[atom].ktensor[3] +
                  cart.z * unitcell[atom].ktensor[4];

               /* col 3 */
               temp_vector.z =
                  cart.x * unitcell[atom].ktensor[2] +
                  cart.y * unitcell[atom].ktensor[4] +
                  cart.z * unitcell[atom].ktensor[5];

               double k =
                  temp_vector.x * cart.x +
                  temp_vector.y * cart.y +
                  temp_vector.z * cart.z;

               if (k < min_energy) {
                  min_energy = k;
                  min_phi = phi;
                  min_theta = theta;
               }

            }

         }

         double phi_rad = min_phi / 180.0 * 3.14159;
         fout << atom << "\t" << unitcell[atom].pos.z << "\t" << cos(phi_rad) << "\n";
   }

   std::cout << "data output to 'easy_axes.dat'\n";


   return EXIT_SUCCESS;

}
