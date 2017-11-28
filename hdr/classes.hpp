/* vector class */
class vec_t
{
    public:
        double x, y, z;

        double length(void) {
            return sqrt(x*x+y*y+z*z);
        }

        vec_t operator - (const vec_t& v)
        {
            vec_t vec;
            vec.x = this->x - v.x;
            vec.y = this->y - v.y;
            vec.z = this->z - v.z;
            return vec;
        }

        void operator = (const double v)
        {
            x = v;
            y = v;
            z = v;
        }
};

struct atom_t
{
    int aid, gid;           /* atom id (uc / global) */
    int type;               /* species */
    vec_t pos;              /* position */
    double k;               /* anisotropy energy */
};
