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


#if 0
#define MAX_KD_DIST (10 * 10)
#else
// for sea floor
#define MAX_KD_DIST (5 * 5)
#endif

// Find closest point to p on segment from v0 to v1
point closest_on_segment(const point &v0, const point &v1, const point &p) {
  vec v01 = v1 - v0;
  float d = (p - v0) DOT v01;
  d /= len2(v01);
  if (d < 0.0f)
    d = 0.0f;
  else if (d > 1.0f)
    d = 1.0f;
  return v0 + d * v01;
}


// Find closest point to p on face i of mesh
point closest_on_face(const TriMesh *mesh, int i, const point &p) {
  const TriMesh::Face &f = mesh->faces[i];
  const point &v0 = mesh->vertices[f[0]];
  const point &v1 = mesh->vertices[f[1]];
  const point &v2 = mesh->vertices[f[2]];
  vec a = v1 - v0, b = v2 - v0, p1 = p - v0, n = a CROSS b;

  float A[3][3] = { { a[0], b[0], n[0] },
                    { a[1], b[1], n[1] },
                    { a[2], b[2], n[2] } };
  float x[3] = { p1[0], p1[1], p1[2] };
  int indx[3];
  ludcmp<float, 3>(A, indx);
  lubksb<float, 3>(A, indx, x);

  if (x[0] >= 0.0f && x[1] >= 0.0f && x[0] + x[1] <= 1.0f)
    return v0 + x[0] * a + x[1] * b;

  point c01 = closest_on_segment(v0, v1, p);
  point c12 = closest_on_segment(v1, v2, p);
  point c20 = closest_on_segment(v2, v0, p);
  float d01 = dist2(c01, p);
  float d12 = dist2(c12, p);
  float d20 = dist2(c20, p);
  if (d01 < d12) {
    if (d01 < d20) return c01; else return c20;
  } else {
    if (d12 < d20) return c12; else return c20;
  }
}


// Find (good approximation to) closest point on mesh to p.
// Finds closest vertex, then checks all faces that touch it.
// If the mesh has no faces, returns the closest vertex
// Returns nearest vertex num, with nearest point in out
// Return < 0 for no match
int closest_pt(const TriMesh *mesh, const KDtree *kd, float mdist, const point &p, point &out) {
  // float maxdist2 = sqr(dist(p, mesh->bsphere.center) + mesh->bsphere.r);
  const float *match = kd->closest_to_pt(p, mdist);
  if (!match) return -1;

  int ind = (match - (const float *) &(mesh->vertices[0][0])) / 3;
  if (ind < 0 || ind >= (int) mesh->vertices.size()) {
    fprintf(stderr, "This can't happen - bogus index\n");
    return -1;
  }

  if (mesh->faces.size() == 0)  {
    out = mesh->vertices[ind];
    return ind;
  } else assert(mesh->adjacentfaces.size());

  const vector<int> &a = mesh->adjacentfaces[ind];
  if (a.empty()) {
    fprintf(stderr, "Found match to unconnected point\n");
    out =  mesh->vertices[ind];
    return ind;
  }

  float closest_dist2 = 3.3e33;
  for (unsigned int i = 0; i < a.size(); i++) {
    point c1 = closest_on_face(mesh, a[i], p);
    float this_dist2 = dist2(c1, p);
    if (this_dist2 < closest_dist2) {
      closest_dist2 = this_dist2;
      out = c1;
    }
  }
  return ind;
}

struct opts_t {

	float max_kd_dist;

	char* tgt_name;

	opts_t(int *argc, char** argv) : max_kd_dist(MAX_KD_DIST) {

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
			} else if (!strcmp(argv[i], "-max_kd_dist")) {
				assert(i < (*argc - 1));
				max_kd_dist = atof(argv[++i]);
				max_kd_dist = max_kd_dist * max_kd_dist;
			} else {
				fprintf(stderr, "Unknown option %s\n", argv[i]);
				exit(-1);
			}
		}
	}
};

struct DistanceInfo {
	int num_of_distances;
	float sum_of_distances;
	float mean_distance;
	DistanceInfo(int nod = 0, float sod = 0.0f) : num_of_distances(nod), sum_of_distances(sod) {
		if(num_of_distances != 0) mean_distance = sum_of_distances / num_of_distances;
		else mean_distance = 0.0f;
	}
};

int main(int argc, char **argv) {

  opts_t opts(&argc, argv);

  // read in master mesh
  TriMesh *mesh = TriMesh::read(opts.tgt_name);

  vector<DistanceInfo> distanceInfos(mesh->vertices.size());

  while (!feof(stdin)) {
    char fname[1024];
    scanf(" %s ", fname);

	if (!strcmp(opts.tgt_name, fname)) continue;
	
    TriMesh *m = TriMesh::read(fname);
	m->need_neighbors();
	m->need_adjacentfaces();

	KDtree *kd = new KDtree(m->vertices);
                                                                                                                                         
    for (unsigned int i = 0; i < mesh->vertices.size(); i++) {
#if 0
      const float *p = kd->closest_to_pt(mesh->vertices[i], MAX_KD_DIST);
      if (!p) continue; // skip outliers

      int j = (p - (const float *) &m->vertices[0]) / 3;
#else
      point p;
      int j = closest_pt(m, kd, opts.max_kd_dist, mesh->vertices[i], p);    
      if (j < 0) continue;
#endif

      // skip boundary points
      if (m->neighbors[j].size() != m->adjacentfaces[j].size()) continue;

      double d2 = dist2(mesh->vertices[i], p);
      double d = sqrt(d2);
	  
	  distanceInfos[i].num_of_distances++;
	  distanceInfos[i].sum_of_distances += d;
    }

    delete m;
  }

  char output_name[1024];
  sprintf(output_name, "distances/%s_DistanceInfo.txt", opts.tgt_name);
  FILE *h = fopen(output_name, "w");
  for (int i = 0; i < distanceInfos.size(); i++) {
	if(distanceInfos[i].num_of_distances != 0) distanceInfos[i].mean_distance = distanceInfos[i].sum_of_distances / distanceInfos[i].num_of_distances;
	else distanceInfos[i].mean_distance = 0.0f;
    fprintf(h, "%d %f %f\n", distanceInfos[i].num_of_distances, distanceInfos[i].sum_of_distances, distanceInfos[i].mean_distance);
  }
  fclose(h);

  return 0;
}
