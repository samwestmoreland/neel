#include<cmath>
#include<fstream>
#include<iostream>
#include<vector>

/* vector class */
class vec_t
{
    public:
        double x, y, z;

        double length(void)
        {
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
    vec_t pos;                /* position */
    double k;               /* anisotropy energy */
};

/* prototypes */
namespace cal
{
    double dot(vec_t i, vec_t j);
    vec_t eij(vec_t ipos, vec_t jpos);
    double lij(int itype, int jtype, double r);
    void k_tensor(std::vector<atom_t> super, int ucsize, double rcut);
    void k_vec(std::vector<atom_t> super, int ucsize, double rcut);
}

int main(int argc, char* argv[])
{
    vec_t ucd;
    if (argc == 3)
    {
         ucd.x = atof(argv[1]);
         ucd.y = atof(argv[1]);
         ucd.z = atof(argv[2]);
    }

    else
    {
        ucd.x = 8.8;
        ucd.y = 8.8;
        ucd.z = 12.2;
    }

    /* set cut-off range */
    double rcut = 10.0;

    std::string c_file = "coordinate_files/interface.coords";
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
        ++aidcount;
        tmp.k = 0;

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
    std::cout << "Fe atom dimensions:" << std::endl;
    std::cout << "x: " << fe_min.x << " - " << fe_max.x << std::endl;
    std::cout << "y: " << fe_min.y << " - " << fe_max.y << std::endl;
    std::cout << "x: " << fe_min.z << " - " << fe_max.z << std::endl;
    std::cout << std::endl;
    std::cout << "Nd atom dimensions:" << std::endl;
    std::cout << "x: " << nd_min.x << " - " << nd_max.x << std::endl;
    std::cout << "y: " << nd_min.y << " - " << nd_max.y << std::endl;
    std::cout << "x: " << nd_min.z << " - " << nd_max.z << std::endl;
    std::cout << std::endl;

    /* replicate uc */

    std::vector<atom_t> super;
    int gidcount = 0;

    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k)
                for (int site=0; site<uc.size(); ++site)
                {
                    atom_t tmp;

                    tmp.type = uc[site].type;
                    tmp.aid = uc[site].aid;
                    tmp.pos.x = uc[site].pos.x + i*ucd.x;
                    tmp.pos.y = uc[site].pos.y + j*ucd.y;
                    tmp.pos.z = uc[site].pos.z + k*ucd.z;

                    tmp.gid = gidcount;
                    ++gidcount;

                    tmp.k = 0;

                    super.push_back(tmp);
                }

    std::cout << "\ncalculating using vector method...\n\n";
    cal::k_vec(super, uc.size(), rcut);

    std::cout << "\ncalculating using tensor method...\n\n";
    cal::k_tensor(super, uc.size(), rcut);

    return EXIT_SUCCESS;
}

namespace cal
{
    vec_t eij(vec_t ipos, vec_t jpos)
    {
        return jpos - ipos;
    }

    double lij(int itype, int jtype, double r)
    {
        double lij;
        double r0 = 5.0;
        double l0 = -1.0;

        /* fe-fe */
        // if (itype == 0 && jtype == 0) lij = l0*exp(-r/r0);

        /* nd-* */
        if (itype == 1) lij = l0*exp(-r/r0);
        else lij = 0;

        return lij;
    }

    double dot(vec_t i, vec_t j)
    {
        return i.x*j.x
            +  i.y*j.y
            +  i.z*j.z;
    }

    void k_vec(std::vector<atom_t> super, int ucsize, double rcut)
    {
        int start = (super.size() - ucsize)/2;
        int end = (super.size() + ucsize)/2;

        double k_hard = 0;
        double k_easy = 0;
        double k;

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

        /* loop through centre cell */
        for (int i=start; i<end; ++i)

            /* loop through all atoms */
            for (int j=0; j<super.size(); ++j)
            {
                /* check atom is within cut off */
                vec_t eij = cal::eij(super[i].pos, super[j].pos);
                double rij = eij.length();
                if (rij < rcut && rij > 1e-35)
                {
                    k_hard += cal::lij(super[i].type, super[j].type, rij)
                        * dot(hard, eij) * dot(hard, eij);

                    k_easy += cal::lij(super[i].type, super[j].type, rij)
                        * dot(easy, eij) * dot(easy, eij);
                }
            }

        k_hard /= 2.0;
        k_easy /= 2.0;
        k = (k_hard - k_easy) / (double) ucsize;

        std::cout << "k_hard = " << k_hard << "\n";
        std::cout << "k_easy = " << k_easy << "\n";
        std::cout << "k = " << k << "\n";
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

        int total_pairs = 0;
        int pairs_within_range = 0;

        /* loop through centre cell */
        for (int i=start; i<end; ++i)

            /* loop through all atoms */
            for (int j=0; j<super.size(); ++j)
            {
                /* check atom is within cut off */
                vec_t eij = cal::eij(super[i].pos, super[j].pos);
                double rij = eij.length();
                total_pairs ++;
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

        /* because we're just considering hard and easy directions,
         * lots of terms cancel and we end up only needing
         * two tensor elements */

        double k_hard = tensor[0] / 2.0 / (double) ucsize;
        double k_easy = tensor[5] / 2.0 / (double) ucsize;
        double k = (k_hard - k_easy);

        std::cout << "k_hard = " << k_hard << "\n";
        std::cout << "k_easy = " << k_easy << "\n";
        std::cout << "k = " << k << "\n";
        std::cout << "pairs in range = " << pairs_within_range << "\n";
    }

} /* end of namespace cal */
