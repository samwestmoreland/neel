#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"

int read_inputfile() {

   std::cout << "\nreading input file..." << std::endl;;
   std::ifstream inputfile ("neel.in");

   /* exit if no input file found */
   if (!inputfile.good()) {
      std::cout << "input file \"neel.in\" not found. exiting." << std::endl;
      exit(EXIT_FAILURE);
   }

   std::string line;

   /* read each line of file */
   while (std::getline(inputfile, line))
   {
      /* skip comments or blank lines */
      if ((line[0]=='#')||line.empty()) continue;

      /* remove spaces from line */
      line.erase(remove(line.begin(), line.end(), ' '), line.end());

      /* identify delimiters and initialise key string */
      const char* eq = "=";
      const char* cm = ",";

      std::string key = "";
      std::string val = "";

      /* point in line where key ends */
      int endvar;

      /* loop through characters in line to identify key */
      for (int i=0; i<line.length(); ++i) {

         if (line.at(i) != *eq) key.push_back(line.at(i));
         else {
            endvar = i+1;
            break;
         }
      }

      /* loop through characters after equals sign */
      for (int i=endvar; i<line.length(); ++i) val.push_back(line.at(i));

      /* material name */
      if (key == "material") mat.name = val;

      /* cut off radius */
      else if (key == "cutoffradius") sim.rcut = stof(val);

      /* angle increment */
      else if (key == "angleincrement") sim.angle_increment = stof(val);

      /* else unrecognised */
      else {
         std::cout << "input file error: i don't know what " << "\"" << key << "\" means. exiting.\n";
         exit(EXIT_FAILURE);
      }

   } /* end of while loop */

   get_unitcell_dimensions();

   return EXIT_SUCCESS;

}

int read_coordinates() {

   /* coordinate file name */
   std::string coord_file;
   if (mat.name == "interface") coord_file = "./coordinates/mirror_shifted.coords";
   else coord_file = "./coordinates/" + mat.name + ".coords";

   /* initialise filestream for coordinate file */
   std::ifstream coordstream (coord_file.c_str());

   int fe_count = 0;
   int re_count = 0;

   /* check that coordinate file opened correctly */
   if (coordstream.is_open()) {

      std::cout << "reading in coordinates from \"" << coord_file << "\"...";

      atom_t tmp;
      int atom_count = 0;
      unitcell_count = 0;

      /* loop through coordinate file */
      while (coordstream >> tmp.type >> tmp.pos.x >> tmp.pos.y >> tmp.pos.z) {

         /* assign atom ids */
         tmp.id = atom_count;
         atom_count ++;

         /* scale coordinates to unit cell dimensions */
         if (mat.name != "interface") tmp.pos = tmp.pos * mat.ucd;

         /* initialise tensors to 0 */
         for (int i=0; i<6; ++i) tmp.ktensor[i] = 0.0;

         unitcell.push_back(tmp);
         unitcell_count ++;

         /* update counters */
         if (tmp.is_fe()) fe_count ++;
         else if (tmp.is_re()) re_count ++;

      }

   }

   /* otherwise quit with error message */
   else {

      std::cout << "coordinate file \"" << coord_file << "\" not found. exiting.\n";
      exit(EXIT_FAILURE);

   }

   std::cout << "\nsystem information:" << std::endl;
   std::cout << std::endl;
   std::cout << "number of fe atoms: " << fe_count << std::endl;
   std::cout << "number of re atoms: " << re_count << std::endl;

   std::ofstream fout ("unitcell.xyz");

   fout << unitcell.size() << "\n\n";

   for (int i=0; i<unitcell.size(); ++i) {

      fout << unitcell[i].type << "\t"
           << unitcell[i].pos.x << "\t"
           << unitcell[i].pos.y << "\t"
           << unitcell[i].pos.z << "\n";
   }


   return EXIT_SUCCESS;

}

int get_unitcell_dimensions() {

//   if (mat.name == "ndfe12_0") {
//
//      mat.ucd.x = 17.025013;
//      mat.ucd.y = 17.025013;
//      mat.ucd.z = 4.842164;
//
//   }
//
//   else if (mat.name == "ndfe12_1") {
//
//      mat.ucd.x = 17.040534;
//      mat.ucd.y = 17.0406127;
//      mat.ucd.z = 4.8449;
//
//   }
//
//   else if (mat.name == "ndfe12_2") {
//
//      mat.ucd.x = 17.061906;
//      mat.ucd.y = 17.061906;
//      mat.ucd.z = 4.847393;
//
//   }
//
//   else if (mat.name == "ndfe12_3") {
//
//      mat.ucd.x = 17.07819;
//      mat.ucd.y = 17.082273;
//      mat.ucd.z = 4.849857;
//
//   }
//
//   else if (mat.name == "ndfe12_4") {
//
//      mat.ucd.x = 17.098785;
//      mat.ucd.y = 17.098785;
//      mat.ucd.z = 4.85208;
//
//   }

   if (mat.name == "ndfe12_0") {

      mat.ucd.x = 17.025013;
      mat.ucd.y = 17.025013;
      mat.ucd.z = 4.842164;

   }

   else if (mat.name == "ndfe12_1") {

      mat.ucd.x = 17.025013;
      mat.ucd.y = 17.025013;
      mat.ucd.z = 4.842164;

   }

   else if (mat.name == "ndfe12_2") {

      mat.ucd.x = 17.025013;
      mat.ucd.y = 17.025013;
      mat.ucd.z = 4.842164;

   }

   else if (mat.name == "ndfe12_3") {

      mat.ucd.x = 17.025013;
      mat.ucd.y = 17.025013;
      mat.ucd.z = 4.842164;

   }

   else if (mat.name == "ndfe12_4") {

      mat.ucd.x = 17.025013;
      mat.ucd.y = 17.025013;
      mat.ucd.z = 4.842164;

   }

   else if (mat.name == "bccfe") {

      mat.ucd.x = 2.856;
      mat.ucd.y = 2.856;
      mat.ucd.z = 2.856;

   }

   else if (mat.name == "interface") {

      mat.ucd.x = 26.403;
      mat.ucd.y = 26.403;
      mat.ucd.z = 255.93180;

   }

   else if (mat.name == "ndfeb") {

      mat.ucd.x = 8.8;
      mat.ucd.y = 8.8;
      mat.ucd.z = 12.2;

   }

   return EXIT_SUCCESS;

}
