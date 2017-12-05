#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "../hdr/classes.hpp"

/* prototypes */
namespace cal
{
   double dot(vec_t i, vec_t j);
   double lij(int itype, int jtype, double r);
   void k_tensor(std::vector<atom_t> super, int ucsize, double rcut);
   void k_vec(std::vector<atom_t> uc, std::vector<atom_t> super, int ucsize, double rcut);
}

int main(int argc, char* argv[])
{
   /* set unit cell dimensions */
   vec_t ucd;
   ucd.x = 26.181;
   ucd.y = 26.181;

   /* set cut-off range */
   double rcut = 5.0;

   /* set name of coordinate file */
   std::string c_file = "coordinate_files/run00d1_final_frame.coords";

   std::ifstream coordfile;
   coordfile.open (c_file.c_str());

   if (!coordfile.is_open())
   {
      std::cout << "coordinate file \"" << c_file << "\" not found. exiting.\n";
      exit(EXIT_FAILURE);
   }

   std::string elstring;
   atom_t tmp;
   int aidcount = 0;
   std::vector<atom_t> uc;

   int nd_count = 0;
   int fe_count = 0;

   /* system dimensions */
   vec_t fe_min, fe_max;
   vec_t nd_min, nd_max;
   fe_min = 0;
   fe_max = 0;
   nd_min = 0;
   nd_max = 0;

   std::cout << "reading in coordinates...";

   while (coordfile >> elstring >> tmp.pos.x >> tmp.pos.y >> tmp.pos.z)
   {
      tmp.aid = aidcount;
      tmp.gid = aidcount;
      ++aidcount;
      tmp.k = 0;

      /* determine min and max coordinates for each species */
      if (elstring == "Fe")
      {
         tmp.type = 0;
         fe_count ++;

         if (tmp.pos.x < fe_min.x) fe_min.x = tmp.pos.x;
         else if (tmp.pos.x > fe_max.x) fe_max.x = tmp.pos.x;
         if (tmp.pos.y < fe_min.y) fe_min.y = tmp.pos.y;
         else if (tmp.pos.y > fe_max.y) fe_max.y = tmp.pos.y;
         if (tmp.pos.z < fe_min.z) fe_min.z = tmp.pos.z;
         else if (tmp.pos.z > fe_max.z) fe_max.z = tmp.pos.z;
      }
      else if (elstring == "Nd")
      {
         tmp.type = 1;
         nd_count ++;

         if (tmp.pos.x < nd_min.x) nd_min.x = tmp.pos.x;
         else if (tmp.pos.x > nd_max.x) nd_max.x = tmp.pos.x;
         if (tmp.pos.y < nd_min.y) nd_min.y = tmp.pos.y;
         else if (tmp.pos.y > nd_max.y) nd_max.y = tmp.pos.y;
         if (tmp.pos.z < nd_min.z) nd_min.z = tmp.pos.z;
         else if (tmp.pos.z > nd_max.z) nd_max.z = tmp.pos.z;
      }

      uc.push_back(tmp);
   }

   std::cout << "done" << "\n\n";

   std::cout << "system info" << std::endl;
   std::cout << "===============================" << std::endl;
   std::cout << std::endl;
   std::cout << "number of Fe atoms: " << fe_count << std::endl;
   std::cout << "number of Nd atoms: " << nd_count << std::endl;
   std::cout << std::endl;
   std::cout << "Fe atoms dimensions:" << std::endl;
   std::cout << "x: " << fe_min.x << " - " << fe_max.x << std::endl;
   std::cout << "y: " << fe_min.y << " - " << fe_max.y << std::endl;
   std::cout << "z: " << fe_min.z << " - " << fe_max.z << std::endl;
   std::cout << std::endl;
   std::cout << "Nd atoms dimensions:" << std::endl;
   std::cout << "x: " << nd_min.x << " - " << nd_max.x << std::endl;
   std::cout << "y: " << nd_min.y << " - " << nd_max.y << std::endl;
   std::cout << "z: " << nd_min.z << " - " << nd_max.z << std::endl;
   std::cout << std::endl;

   /* replicate uc */

   std::cout << "replicating unit cell...";
   std::vector<atom_t> super;
   int gidcount = 0;

   for (int i=0; i<3; ++i)
   for (int j=0; j<3; ++j)
   for (int site=0; site<uc.size(); ++site)
   {
      atom_t tmp;

      tmp.type  = uc[site].type;
      tmp.aid   = uc[site].aid;
      tmp.pos.x = uc[site].pos.x + i*ucd.x;
      tmp.pos.y = uc[site].pos.y + j*ucd.y;
      tmp.pos.z = uc[site].pos.z;

      tmp.gid = gidcount;
      ++gidcount;

      tmp.k = 0;

      super.push_back(tmp);
   }

   std::cout << "done\n";
   std::cout << "total atoms in super cell: " << super.size() << std::endl;
   std::cout << std::endl;
   std::cout << "outputting super cell coordinates to supercell.xyz\n";

   std::ofstream superxyz ("supercell.xyz");
   superxyz << super.size() << "\n\n\n";
   for (int i=0; i<super.size(); ++i)
   {
      if (super[i].type == 0) superxyz << "Fe\t";
      else if (super[i].type == 1) superxyz << "Nd\t";

      superxyz << super[i].pos.x << "\t"
               << super[i].pos.y << "\t"
               << super[i].pos.z << "\n";
   }

   std::cout << "\ncalculating using vector method...\n\n";

   std::clock_t start_vector;
   double vector_duration;

   start_vector = std::clock();

   cal::k_vec(uc, super, uc.size(), rcut);

   vector_duration = (std::clock() - start_vector) / (double) CLOCKS_PER_SEC;

   std::cout << "\ncalculating using tensor method...\n\n";

   std::clock_t start_tensor;
   double tensor_duration;

   start_tensor = std::clock();

   cal::k_tensor(super, uc.size(), rcut);

   tensor_duration = (std::clock() - start_tensor) / (double) CLOCKS_PER_SEC;

   std::cout << "vector calculation time: " << vector_duration << "s" << std::endl;
   std::cout << "tensor calculation time: " << tensor_duration << "s" << std::endl;
   std::cout << std::endl;

   return EXIT_SUCCESS;
}

namespace cal
{
   double lij(int itype, int jtype, double r)
   {
      double lij;
      double r0;
      double l0;

      r0 = 5.0;
      l0 = 1.0;

      /* fe */
      if (itype == 0) lij = l0*exp(-r/r0);

      /* nd */
      else if (itype == 1) lij = l0*exp(-r/r0);

      else lij = 0;

      return lij;
   }

   double dot(vec_t i, vec_t j)
   {
      return i.x*j.x + i.y*j.y + i.z*j.z;
   }

   void k_vec(std::vector<atom_t> uc, std::vector<atom_t> super, int ucsize, double rcut)
   {
      /* determine beginning and end of centre cell */
      int start = (super.size() - ucsize)/2;
      int end = (super.size() + ucsize)/2;

      double knd_hard = 0;
      double knd_easy = 0;

      double kfe_hard = 0;
      double kfe_easy = 0;

      double knd;
      double kfe;

      /* hard direction */
      vec_t hard;
      hard.x = 1;
      hard.y = 0;
      hard.z = 0;

      /* easy direction */
      vec_t easy;
      easy.x = 0;
      easy.y = 0;
      easy.z = 1;

      /* initialise arrays for atom resolved hard and easy energies */
      std::vector<double> k_hard_atom(ucsize, 0);
      std::vector<double> k_easy_atom(ucsize, 0);

      /* loop through centre cell */
      for (int i=start; i<end; ++i)

         /* loop through all atoms */
         for (int j=0; j<super.size(); ++j)
         {
            /* calculate neighbour vector */
            vec_t eij = super[j].pos - super[i].pos;

            /* calculate neighbour distance */
            double rij = eij.length();

            /* check if atom is within cut off range */
            if (rij <= rcut && rij > 1e-35)
            {
               /* add energy to hard energy array */
               k_hard_atom[i-start] += cal::lij(super[i].type, super[j].type, rij) * dot(hard, eij) * dot(hard, eij);

               /* add energy to easy energy array */
               k_easy_atom[i-start] += cal::lij(super[i].type, super[j].type, rij) * dot(easy, eij) * dot(easy, eij);

                  //  std::cout <<
                  //      cal::lij(super[i].type, super[j].type, rijk)
                  //      * dot(hard, eij) * dot(hard, eij) << std::endl;

               /* if atom is Nd */
               if (super[i].type == 1)
               {
                  /* add interaction energy to total hard energy */
                  knd_hard += cal::lij(super[i].type, super[j].type, rij) * dot(hard, eij) * dot(hard, eij);

                  /* add interaction energy to total easy energy */
                  knd_easy += cal::lij(super[i].type, super[j].type, rij) * dot(easy, eij) * dot(easy, eij);
               }

               /* if atom is Fe */
               else if (super[i].type == 0)
               {
                  kfe_hard += cal::lij(super[i].type, super[j].type, rij) * dot(hard, eij) * dot(hard, eij);

                  kfe_easy += cal::lij(super[i].type, super[j].type, rij) * dot(easy, eij) * dot(easy, eij);
               }
            }
         }

      /* divide by factor 2, as per neel expression */   
      knd_hard /= 2.0;
      knd_easy /= 2.0;
      kfe_hard /= 2.0;
      kfe_easy /= 2.0;

      /* k is then the difference between hard and easy energies */
      knd = (knd_hard - knd_easy);
      kfe = (kfe_hard - kfe_easy);

      double ktotal = knd + kfe;

      /* output stream for atom resolved energies */
      std::ofstream kout ("k.dat");

      for (int i=0; i<ucsize; ++i)
      {
         if (uc[i].type == 1)
         kout << i << "\t"
              << uc[i].pos.z << "\t"
              << k_hard_atom[i] << "\t"
              << k_easy_atom[i] << "\t"
              << k_hard_atom[i] - k_easy_atom[i]
              << std::endl;
      }

      std::cout << "knd_hard = " << knd_hard << "\n";
      std::cout << "knd_easy = " << knd_easy << "\n";
      std::cout << "knd = " << knd << "\n\n";

      std::cout << "kfe_hard = " << kfe_hard << "\n";
      std::cout << "kfe_easy = " << kfe_easy << "\n";
      std::cout << "kfe = " << kfe << "\n\n";

      std::cout << "k = " << ktotal << "\n";
   }

   void k_tensor(std::vector<atom_t> super, int ucsize, double rcut)
   {
      int start = (super.size() - ucsize)/2;
      int end = (super.size() + ucsize)/2;

      /* initialise system tensor */
      std::vector<double> tensor (6, 0);

      /* hard direction */
      vec_t hard;
      hard.x = 1;
      hard.y = 0;
      hard.z = 0;

      /* easy direction */
      vec_t easy;
      easy.x = 0;
      easy.y = 0;
      easy.z = 1;

      int pairs_within_range = 0;

      /* loop through centre cell */
      for (int i=start; i<end; ++i)

      /* loop through all atoms */
      for (int j=0; j<super.size(); ++j)
      {
         /* check atom is within cut off */
         vec_t eij = super[j].pos - super[i].pos;
         double rij = eij.length();
         if (rij < rcut && rij > 1e-35)
         {
            pairs_within_range ++;
            double lij = cal::lij(super[i].type, super[j].type, rij);

            tensor[0] += eij.x * eij.x * lij;
            tensor[1] += eij.x * eij.y * lij;
            tensor[2] += eij.x * eij.z * lij;
            tensor[3] += eij.y * eij.y * lij;
            tensor[4] += eij.y * eij.z * lij;
            tensor[5] += eij.z * eij.z * lij;
         }
      }

      /*
       * because we're just considering hard and easy directions,
       * lots of terms cancel and we end up only needing
       * two tensor elements
       */

      double k_hard = tensor[0] / 2.0;
      double k_easy = tensor[5] / 2.0;
      double k = (k_hard - k_easy);

      std::cout << "k_hard = " << k_hard << "\n";
      std::cout << "k_easy = " << k_easy << "\n";
      std::cout << "k = " << k << "\n\n";
      std::cout << "pairs in range = " << pairs_within_range << "\n";
   }

} /* end of namespace cal */
