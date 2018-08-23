#ifndef CLASSES_H_
#define CLASSES_H_

#include <string>
#include <cmath>
#include <iostream>

class vec_t {
public:
    double x;
    double y;
    double z;

    // overload '-' operator to subtract two vec_t objects
    vec_t operator-(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x - v.x;
        vec.y = this->y - v.y;
        vec.z = this->z - v.z;
        return vec;
    }

    vec_t operator*(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x * v.x;
        vec.y = this->y * v.y;
        vec.z = this->z * v.z;
        return vec;
    }

    vec_t operator+(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x + v.x;
        vec.y = this->y + v.y;
        vec.z = this->z + v.z;
        return vec;
    }
};

/* struct to hold material parameters */
struct material_t {
   std::string name;
   vec_t ucd;
};

/* simulation parameters */
struct simulation_parameters_t {
   double rcut;
   double angle_increment;
};

class atom_t {

   public:
      vec_t pos;
      std::string type;

      int id;     /* unique unit cell id */
      int gid;    /* global id for use in super cell */

      double ktensor[6];

      bool is_fe() {

         if (
               type == "Fe8i" ||
               type == "Fe8j" ||
               type == "Fe8f" ||
               type == "Fe" ||
               type == "Fe16k1" ||
               type == "Fe16k2" ||
               type == "Fe8j1" ||
               type == "Fe8j2" ||
               type == "Fe4e" ||
               type == "Fe4c")

            return true;

         else return false;

      }

      bool is_re() {

         if (
               type == "Nd" ||
               type == "Nd4f" ||
               type == "Nd4g" ||
               type == "Sm")

            return true;

         else return false;

      }

};

class pair_t {

   public:
      atom_t atomi;
      atom_t atomj;

      double rij() {

         double dx, dy, dz;

         dx = atomj.pos.x - atomi.pos.x;
         dy = atomj.pos.y - atomi.pos.y;
         dz = atomj.pos.z - atomi.pos.z;

         return sqrt(dx*dx + dy*dy + dz*dz);

      }

      /* calculate eij vector between atoms */
      vec_t eij() {

         vec_t eij = atomj.pos - atomi.pos;

         double dist = rij();

         /* normalise to unity */
         eij.x /= dist;
         eij.y /= dist;
         eij.z /= dist;

         return eij;
      }

      /* calculate lij for interaction */
      double lij() {
         if (atomi.type == "Ti" || atomj.type == "Ti")
            return 0.0;
         else
            return 1.0;
      }

      /* self interaction check */
      bool are_same_atom() {
         if (atomi.gid == atomj.gid)
            return true;
         else
            return false;
      }
};

#endif
