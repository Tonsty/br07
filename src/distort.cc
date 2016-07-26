/*
Copyright 2007 Benedict Brown
Princeton University

Compute histgram of distances of range scan points to final mesh

-----

This file is part of tps_alignment.

tps_alignment is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

tps_alignment is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <assert.h>
#include "XForm.h"
#include "TriMesh.h"
#include "KDtree.h"

using namespace std;

struct opts_t {

	float k1;
	float k2;
	float k3;
	float k4;
	float k5;
	float c0;

	char* tgt_name;

	opts_t(int *argc, char** argv) : k1(0.0001f), k2(0.0001f), k3(0.0001f), k4(0.0001f), k5(0.0f), c0(0.0001) {

		bool have_tgt = false;

		for (int i = 1; i < *argc; i++) {
			if (argv[i][0] != '-') {
				if(!have_tgt) {
					have_tgt = true;
					tgt_name = argv[i];
					fprintf(stderr, "tgt_name: %s\n", tgt_name);
				} else {
					fprintf(stderr, "Unexpected argument %d: %s\n", i, argv[i]);
					exit(-1);
				}
			} else if (!strcmp(argv[i], "-k1")) {
				assert(i < (*argc - 1));
				k1 = atof(argv[++i]);
			} else if (!strcmp(argv[i], "-k2")) {
				assert(i < (*argc - 1));
				k2 = atof(argv[++i]);
			} else if (!strcmp(argv[i], "-k3")) {
				assert(i < (*argc - 1));
				k3 = atof(argv[++i]);
			} else if (!strcmp(argv[i], "-k4")) {
				assert(i < (*argc - 1));
				k4 = atof(argv[++i]);
			} else if (!strcmp(argv[i], "-k5")) {
				assert(i < (*argc - 1));
				k5 = atof(argv[++i]);
			} else if (!strcmp(argv[i], "-c0")) {
				assert(i < (*argc - 1));
				c0 = atof(argv[++i]);
			} else {
				fprintf(stderr, "Unknown option %s\n", argv[i]);
				exit(-1);
			}
		}
	}
};

int main(int argc, char **argv) {

  opts_t opts(&argc, argv);

  // read in master mesh
  TriMesh *mesh = TriMesh::read(opts.tgt_name);    
  char tgt_tf_name[1024];
  strcpy(tgt_tf_name, opts.tgt_name);
  strcpy(tgt_tf_name + strlen(tgt_tf_name) - 3, "tf");
  xform xfin;
  xfin.read(tgt_tf_name);

  const float &k1 = opts.k1;
  const float &k2 = opts.k2;
  const float &k3 = opts.k3;
  const float &k4 = opts.k4;
  const float &k5 = opts.k4;
  const float &c0 = opts.c0;

  for (unsigned int i = 0; i < mesh->vertices.size(); i++) {
	  point &p = mesh->vertices[i];
	  float &xc = p[0];
	  float &yc = p[1];
	  float &zc = p[2];

	  float xn = xc/zc;
	  float yn = yc/zc;

	  float r2 = xn * xn + yn * yn;
	  float r4 = r2 * r2;
	  float r6 = r2 * r4;

	  float xg = 2 * k3 * xn * yn + k4 * (r2 + 2 * xn * xn);
	  float yg = k3 * (r2 + 2 * yn * yn) + 2 * k4 * xn * yn;

	  float scale = 1.0f + k1 * r2 + k2 * r4 + k5 * r6;
	  std::cout << "scale = " << scale << " shift = " << xg << ", " << yg << std::endl;
	  float xk = scale * xn + xg; // distorted x
	  float yk = scale * yn + yg; // distorted y

	  zc = zc / ( 1.0f - zc * c0);
	  xc = xk * zc;
	  yc = yk * zc;

	  p = xfin * p;
  }

  char output_name[1024];
  sprintf(output_name, "distorted/%s", opts.tgt_name);

  mesh->write(output_name);

  return 0;
}
