/*
Copyright 2007 - 2011 Benedict Brown
Katholieke Universiteit Leuven

global_reg.cc
TPS Global Alignment

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

#include "global_reg.h"
#include "jama_lu.h"

#include "timestamp.h"
#include "KDtree.h"

#include "rand48.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>

#include <string.h>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>

#include <GMNR/MultiTPS.h>

using namespace std;

#define BREAK_FACTOR 5u

int pruned_points = 0, thinned_points = 0, conf_points = 0;

#undef COMPUTE_RMS_ERR

// #define CONF_THRESH 0.8f
#define CONF_THRESH -1
// #define MAX_SPRING_LEN 500.0f
#define MAX_SPRING_LEN 0.025f
// #define MIN_TARGET_DIST2 64.0f
// #define MIN_TARGET_DIST2 8.0f // for sea floor

#ifndef SQ
#define SQ(x) ((x) * (x))
#endif

FILE *outfile = stdout;

#define MAX_RELAX_ITERS 100
#define MIN_RELAX_ERR    0.0000000005f
// #define MIN_RELAX_ERR    0.000005f
#define MIN_RELAX_MOTION 0.005f

struct input_corr {
  int tgt;
  point p;
  float stability, max_stability, err;
  enum { STABLE, UNSTABLE, DISABLED } status;
};
typedef vector<input_corr> input_corr_vector;

// minimize the error for point tgt
// skip_unstable says whether unstable springs should be ignored
//
// we return the new point position
static point relax(int tgt, vector<float> &errors, vector<point> &targets,
                  const vector<spring_vector> &springs, const vector<bool> &use_points,
                  bool skip_unstable) {
  point old_pos = targets[tgt];
  point new_pos = targets[tgt];
  float new_err = errors[tgt];


  for (int iters = 0; iters < MAX_RELAX_ITERS && errors[tgt] > MIN_RELAX_ERR; iters++) {
    // calculate move downhill
    point dir;
    float total_weight = 0;

    // for each spring
    for (unsigned int i = 0; i < springs[tgt].size(); i++) {
      // skip springs that are invalid
      if (!use_points[springs[tgt][i].tgt]) continue;
      if (skip_unstable && !springs[tgt][i].stable) continue;

      // skip springs that are of 0 length (more or less), because
      // these have very high weight and are numerically unstable
      point p = targets[springs[tgt][i].tgt] - targets[tgt];
      float l = len(p);
      if (likely(l == 0)) continue;

      // calculate how far we'd need to move towards this point to get
      // to the ideal distance away; weight this by the spring's weight
      p /= l;
      dir += (springs[tgt][i].wgt * (l - springs[tgt][i].len) / l) * p;
      total_weight += springs[tgt][i].wgt;
    }

    // calculate final update direction
    if (total_weight < 0.1) break;
    dir /= total_weight;

    // dir is not the right amount to move, but it usually overshoots,
    // and it's magnitude tells of something.  So if it's tiny, we should
    // be close to the minimum and give up
    if (len2(dir) < sqr(MIN_RELAX_MOTION)) break;

    // The error metric is quartic, but once we have the gradient direction (dir),
    // it's quadratic in that direction.  So to find the minimum, just calculate
    // three points and fit a parabola.  We use t=0, t=1, and t=0.5.
    float cur_err =  point_error(tgt, springs[tgt], skip_unstable, targets, use_points);
    point cur_pos = targets[tgt];
    targets[tgt] = cur_pos + 0.5f * dir;
    float mid_err = point_error(tgt, springs[tgt], skip_unstable, targets, use_points);
    targets[tgt] = cur_pos + dir;
    new_err = point_error(tgt, springs[tgt], skip_unstable, targets, use_points);

    // fit quadratic to this
    float a = 2 * new_err - 4 * mid_err + 2 * cur_err;
    float b =    -new_err + 4 * mid_err - 3 * cur_err;
    float min_t = -b / (2 * a);

    if (a < 0.01f /* 0.5f */) {
      // very, very flat minimization, so don't bother
      new_pos = cur_pos;
      new_err = cur_err;
      targets[tgt] = cur_pos;
      break;
    }

    // Now we actually do know how far we're
    // moving, so break if it's not very much.
    // Note that this is scaled by len(dir), but
    // we've already checked that len(dir) is not too small.
    if (min_t < MIN_RELAX_MOTION) break;
    // if (len(min_t * dir) > 50) break;

    // ok, now move the point, and update its error
    new_pos = cur_pos + min_t * dir;
    targets[tgt] = new_pos;
    new_err = point_error(tgt, springs[tgt], skip_unstable, targets, use_points);

    if (cur_err < new_err) {
      // we got worse, which is normally due to numerical instability
      // so just restore the old position and go on
      targets[tgt] = cur_pos;
      new_err = cur_err;
      break;
    }
    // assert(errors[tgt] <= old_err);
  }


  // restore the original position of the point and return the new one
  targets[tgt] = old_pos;
  return new_pos;
}

// Iterate through all points, relaxing them in turn.  skip_unstable
// determines whether unstable springs should be considered in the relaxation.
//
// This function is intended to be called many times.
static void move_points(vector<point> &targets, const vector<bool> &use_points,
                        vector<float> &errors, const vector<spring_vector> &springs,
                        bool skip_unstable) {
  set_errors(errors, springs, targets, use_points);
  vector<point> new_targets(targets.size());

  // relax all vertices in turn, starting with those with the highest errors
  for (unsigned int i = 0; i < targets.size(); i++) {
    new_targets[i] = errors[i] < 0 ? targets[i] :
      relax(i, errors, targets, springs, use_points, skip_unstable);
  }

  for (unsigned int i = 0; i < targets.size(); i++) targets[i] = new_targets[i];
}

// Eliminate feature points which move lots more than nearby ones.
// This is done with respect to feature points on a particular mesh,
// the idea being that if a point on a mesh moves a lot more than nearby
// points on the mesh, it will introduce a weird artifact, and is a bad
// correspondence.
//
// This should be called for every range scan in turn.  It has the disadvantage
// that it turns off the feature point entirely.  Better would be to eliminate the
// springs that seem to be bad.
static void prune_big_moving_points(float mtd2, const corr_vector &corrs, const vector<point> &targets, vector<bool> &use_points) {
  vector< pair<float, int> > sorted_dists(corrs.size());
  vector<bool> remove_points(use_points.size());

  // record how much each feature on this mesh will move
  vector<float> dists(corrs.size());
  for (unsigned int i = 0; i < corrs.size(); i++)
    dists[i] = dist2(corrs[i].pnt, targets[corrs[i].tgt]);

  // compute the nearest features to each given one, and find
  // the median distance moved; if the given feature moves
  // more than four times this, it's bad
  for (unsigned int i = 0; i < corrs.size(); i++) {
    // skip springs to disabled features, and unstable springs
    if (!use_points[corrs[i].tgt]) continue;
    if (!corrs[i].stable) continue;

    vector<float> d(corrs.size());
    int k = 0;

    // sort the list of features on this mesh by their distance to feature i
    for (unsigned int j = 0; j < corrs.size(); j++) {
      if (!corrs[j].stable) continue;
      d[k] = dist2(targets[corrs[i].tgt], targets[corrs[j].tgt]);
      sorted_dists[k].first = -d[j];
      sorted_dists[k].second = j;
      k++;
    }
    d.resize(k);
    sorted_dists.resize(k);
    sort(sorted_dists.begin(), sorted_dists.end());

    // select the 9 closest features up to a maximum distance
    d.clear();
    for (unsigned int j = 0; j < sorted_dists.size() && d.size() < 9; j++) {
      int k = sorted_dists[j].second;
      if (dists[k] > 20 * mtd2) break;
      if (use_points[corrs[k].tgt]) d.push_back(dists[k]);
    }
    sort(d.begin(), d.end());

    // If we don't have 9 nearby points to compare with, don't prune
    if (d.size() != 9) continue;

    
    float maxsqdist = 32 * d[d.size() / 2]; // four times median value
    if (dists[i] > maxsqdist) remove_points[corrs[i].tgt] = true;
  }

  // points are actually removed after all features have been pruned, so
  // that it doesn't matter what order we process features in
  for (unsigned int i = 0; i < use_points.size(); i++) {
    if (use_points[i] & remove_points[i]) {
      pruned_points++;
      fprintf(stderr, "big_moving_points pruned point %d\n", (int) i);
    }
    use_points[i] = use_points[i] & !remove_points[i];
  }
}

// This thins out features that are too close together.
// For each feature points, we find all other features that are
// too close to it, then remove all but the one with highest error.
// This is suboptimal, but reasonable.
static void thin_points(float mtd2, vector<point> &targets, vector<bool> &use_points,
                        vector<float> &errors, const vector<spring_vector> &springs) {
  set_errors(errors, springs, targets, use_points);

  for (unsigned int i = 0; i < targets.size(); i++) {
    if (!use_points[i]) continue;
    vector<bool> nearby(targets.size());
    unsigned int best_neighbor = i;
    nearby[i] = true;
    for (unsigned int j = 0; j < targets.size(); j++) {
      if (i == j) continue;
      if (!use_points[j]) continue;

      float dst2 = dist2(targets[i], targets[j]);
      if (dst2 < mtd2) {
        nearby[j] = true;
        if (errors[j] < errors[best_neighbor]) best_neighbor = j;
      }
    }

    for (unsigned int j = 0; j < targets.size(); j++)
      if (nearby[j] && j != best_neighbor) {
        thinned_points++;
        fprintf(stderr, "thin_points pruned point %d\n", (int) j);
        use_points[j] = false;
      }
  }
}

void align_scan(const opts_t &opts, const char *mesh_name, const corr_vector &corrs,
                const vector<point> &targets, const vector<bool> &use_points, const vector<float> &confidence,
                TriMesh *points_mesh, const vector<TriMesh::Face> &all_faces, pthread_mutex_t *mutex) {
  int tps_size = (int) corrs.size();

  // skip meshes that don't have enough features on them
  if (tps_size < 10) {
    fprintf(stderr, "Skipping %s\n", mesh_name);
    return;
  }

  farr x(tps_size, 3), y(tps_size, 3);
  fprintf(stderr, "Opening %s\n", mesh_name);

  xform xfin;
  if (opts.read_xf) {
    // the assumes the mesh has a three-letter extension
    char xfname[1024];
    strcpy(xfname, mesh_name);
    strcpy(xfname + strlen(xfname) - 3, "xf");
    xfin.read(xfname);
    fprintf(stderr, "Read %s:\n", xfname);
    std::cerr << xfin;
  }
  TriMesh *mesh = NULL;

  if (opts.nonrigid_prefix || !opts.write_xf) {
    mesh = TriMesh::read(mesh_name);
    mesh->normals.clear();
    assert(mesh);
  }

  fprintf(stderr, "Original TPS size: %d\n", (int) tps_size);

  // Points_mesh gets written, so we have to lock it for thread safety.
  // It may seem like we're locking a lot because the calls to points_mesh
  // are fairly diffused.  But it doesn't matter because this code executes
  // very fast; the slow part of the loop is actually computing and applying
  // the warps.  That comes later.
  if (mutex) assert(!pthread_mutex_lock(mutex));

  // we have many correspondences listed, but some are unstable, so we
  // need to build a list of usable points; We also keep track of the number
  // of unstable points so we can reject scans if there are too many of them
  int real_tps_size = 0, unstable_corrs = 0;;
  points_mesh->faces.resize(12 * tps_size); // to write cubes around each point

  for (int j = 0; j < tps_size; j++, real_tps_size++) {
    int src = corrs[j].vnum, tgt = corrs[j].tgt;

    if (/* (confidence[tgt] >= CONF_THRESH) && */ (use_points[tgt] == true) && (!corrs[j].stable || confidence[tgt] < CONF_THRESH)) {
      // Throw out features points whose correspondence is unstable, but which
      // would otherwise be fine.  Record the instability.
      real_tps_size--;
      unstable_corrs++;
      continue;
    } else if ((confidence[tgt] < CONF_THRESH) || (use_points[tgt] == false) || !corrs[j].stable) {
      // Throw out any other points we can't use.
      real_tps_size--;
      continue;
    }

    int l = real_tps_size; // for convenience

    // fill in the location of this point on the mesh
    // and its target position in our big matrices
    x[l][0] = corrs[j].pnt[0];
    x[l][1] = corrs[j].pnt[1];
    x[l][2] = corrs[j].pnt[2];

    y[l][0] = targets[tgt][0];
    y[l][1] = targets[tgt][1];
    y[l][2] = targets[tgt][2];

    // copy this cube from the list of all cubes
    for (int k = 0; k < 12; k++) {
      points_mesh->faces[12 * l + k][0] = all_faces[12 * tgt + k][0];
      points_mesh->faces[12 * l + k][1] = all_faces[12 * tgt + k][1];
      points_mesh->faces[12 * l + k][2] = all_faces[12 * tgt + k][2];
    }

    // record the correspondence for posterity
    fprintf(stderr, "%d -> %d: (%f, %f, %f) -> (%f, %f, %f)\n", src, tgt, x[l][0], x[l][1], x[l][2], y[l][0], y[l][1], y[l][2]);

    // Shout if this point is moving a long way; useful for debugging
    float dist = SQ(x[l][0] - y[l][0]) + SQ(x[l][1] - y[l][1]) + SQ(x[l][2] - y[l][2]);
    if (dist > 500) fprintf(stderr, "Big Move! (%s)\n", mesh_name);
  }

  fprintf(stderr, "Real TPS size: %d, unstable_corrs: %d\n", real_tps_size, unstable_corrs);

  if (real_tps_size < 10) {
    fprintf(stderr, "Skipping %s (< 10)\n", mesh_name);
    if (mesh) delete mesh;
    if (mutex) pthread_mutex_unlock(mutex); // unlock mutex before returning
    return;
  }

#if 1
  // if most of the corresondences were unstable, that's a clue
  // that the mesh is inherently unstable, and we should give up
  // on it
  if (real_tps_size < unstable_corrs) {
    fprintf(stderr, "Skipping %s (< 50%%) %d %d\n",
            mesh_name, real_tps_size, unstable_corrs);
    if (mesh) delete mesh;
    if (mutex) pthread_mutex_unlock(mutex); // unlock mutex before returning
    return;
  }
#endif

  // We may have discarded points, so we must shrink the arrays;
  // TNT matrices do not automatically resize.
  if (real_tps_size < tps_size) {
    farr new_x(real_tps_size, 3), new_y(real_tps_size, 3);
    for (int i = 0; i < real_tps_size; i++) {
      new_x[i][0] = x[i][0]; new_x[i][1] = x[i][1]; new_x[i][2] = x[i][2];
      new_y[i][0] = y[i][0]; new_y[i][1] = y[i][1]; new_y[i][2] = y[i][2];
    }
    x = new_x.copy(); y = new_y.copy();
    tps_size = real_tps_size;
  }

  // write out the mesh of all points used for this mesh
  char points_name[1024];
  sprintf(points_name, "points/points_%s", get_mesh_name(mesh_name));
  points_mesh->faces.erase(points_mesh->faces.begin() + 12 * tps_size,
                           points_mesh->faces.end());
  points_mesh->write(points_name);

  // Now we're done with points_mesh, so we're thread-safe from here on out
  if (mutex) pthread_mutex_unlock(mutex);

  char out_name[1024]; // mesh name for output

  // make a backup of the mesh vertices, because the
  // various warps will clobber them
  farr verts;
  if (mesh) {
    verts = farr(mesh->vertices.size(), 3);
    for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
      point p = xfin * mesh->vertices[j];
      for (int k = 0; k < 3; k++) verts[j][k] =  p[k];
    }
  }

  if (opts.rigid_prefix) {
    // do rigid-body alignment, which is mostly in a helper because it's long
    farr R;
    point T;
    compute_rigid_alignment(x, y, R, T);

    xform xfout = xform(R[0][0], R[1][0], R[2][0], 0,
                        R[0][1], R[1][1], R[2][1], 0,
                        R[0][2], R[1][2], R[2][2], 0,
                        T[0], T[1], T[2], 1);

    if (opts.write_xf) {
      sprintf(out_name, "%s%s", opts.rigid_prefix, get_mesh_name(mesh_name));
      strcpy(out_name + strlen(out_name) - 3, "xf");
      fprintf(stderr, "Writing .xf to %s\n", out_name);
      xform xfwrite = xfout * xfin;
      xfwrite.write(out_name);
    } else {
      for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
        point p = point(verts[j][0], verts[j][1], verts[j][2]);
        mesh->vertices[j] = xfout * p;
      }

      sprintf(out_name, "%s%s", opts.rigid_prefix, get_mesh_name(mesh_name));
      mesh->write(out_name);
    }
  }

  if (opts.nonrigid_prefix || opts.affine_prefix) {
    // do TPS warp
    matrix_t lambda = opts.lambda * tps_size;
    farr A, w;
    
    update_transform_mat(x, y, lambda, w, A);

    if (opts.affine_prefix) {
      xform xfout = xform(A[0][0], A[0][1], A[0][2], A[0][3],
                          A[1][0], A[1][1], A[1][2], A[1][3],
                          A[2][0], A[2][1], A[2][2], A[2][3],
                          0, 0, 0, 1);
      if (opts.write_xf) {
        sprintf(out_name, "%s%s", opts.affine_prefix, get_mesh_name(mesh_name));
        strcpy(out_name + strlen(out_name) - 3, "xf");
        xform xfwrite = xfout * xfin;
        xfwrite.write(out_name);
      } else {
        for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
          point p = point(verts[j][0], verts[j][1], verts[j][2]);
          mesh->vertices[j] = xfout * p;
        }

        sprintf(out_name, "%s%s", opts.affine_prefix, get_mesh_name(mesh_name));
        mesh->write(out_name);
      }
    }

    if (opts.nonrigid_prefix) {
      farr tmp_x = warp_points(verts, x, w, A);

      for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
        mesh->vertices[j][0] = tmp_x[j][0];
        mesh->vertices[j][1] = tmp_x[j][1];
        mesh->vertices[j][2] = tmp_x[j][2];
      }

      sprintf(out_name, "%s%s", opts.nonrigid_prefix, get_mesh_name(mesh_name));
      mesh->write(out_name);
    }
  }

  if (mesh) delete mesh;
}

struct point_pair{
	point x;
	point y;
	point_pair(const point &x, const point &y) {
		this->x = x;
		this->y = y;
	}
};

typedef vector<point_pair> point_pair_vector;

//The strategy for ApproxiMultiTPS
void icv2ppvs_ApproxiMultiTPS(const opts_t &opts, const int &num_meshes, const vector<input_corr_vector> &icv, const vector<bool> &use_points,
	vector< vector<point_pair_vector> > &ppvs) {

		std::cout << "original : " << icv.size() << std::endl;

		ppvs.resize(num_meshes, vector<point_pair_vector>(num_meshes));

		int unused = 0, unstable = 0, unallowed = 0;

		//The strategy for ApproxiMultiTPS
		for (int i = 0; i < icv.size(); i++) {
			if(!use_points[i]){
				unused++;
				continue;
			}
			for (int j = 0; j < icv[i].size(); j++) {
				if (icv[i][j].status == input_corr::STABLE) {
					int mesh_index_x = icv[i][j].tgt;
					int mesh_index_y = icv[i][(j+1)%icv[i].size()].tgt;
					point x = icv[i][j].p;
					point y = icv[i][(j+1)%icv[i].size()].p;
					float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
					if(dist < 500) {
						point_pair pp(x, y);
						ppvs[mesh_index_x][mesh_index_y].push_back(pp);
					} else unallowed++;
				}else unstable++;
			}
		}

		std::cout << "unused : " << unused << std::endl;
		std::cout << "unstable : " << unstable << std::endl;
		std::cout << "unallowed : " << unallowed << std::endl;

		for (int i = 0; i < num_meshes; i++) {
			std::cout << i << ": ";
			for (int j = 0; j < num_meshes; j++) {
				std::cout << ppvs[i][j].size() << " ";
			}
			std::cout << std::endl;
		}
}

//The first strategy
void icv2ppvs_1(const opts_t &opts, const int &num_meshes, const vector<input_corr_vector> &icv, const vector<bool> &use_points,
	vector< vector<point_pair_vector> > &ppvs) {

    std::cout << "original : " << icv.size() << std::endl;

	ppvs.resize(num_meshes, vector<point_pair_vector>(num_meshes));

	int unused = 0, unstable = 0, unallowed = 0;

	//The first strategy
	for (int i = 0; i < icv.size(); i++) {
		if(!use_points[i]){
			unused++;
			continue;
		}
		for (int j = 0; j < icv[i].size(); j++) {
			if (icv[i][j].status == input_corr::STABLE) {
				int mesh_index_x = icv[i][j].tgt;
				int mesh_index_y = icv[i][(j+1)%icv[i].size()].tgt;
				point x = icv[i][j].p;
				point y = icv[i][(j+1)%icv[i].size()].p;
				float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
				if(dist < 500) {
					if (mesh_index_x < mesh_index_y) {
						point_pair pp(x, y);
						ppvs[mesh_index_x][mesh_index_y].push_back(pp);
					} 
					else{
						point_pair pp(y, x);
						ppvs[mesh_index_y][mesh_index_x].push_back(pp);
					}
				} else unallowed++;
			}else unstable++;
		}
	}

	std::cout << "unused : " << unused << std::endl;
	std::cout << "unstable : " << unstable << std::endl;
	std::cout << "unallowed : " << unallowed << std::endl;

	for (int i = 0; i < num_meshes; i++) {
		std::cout << i << ": ";
		for (int j = 0; j < num_meshes; j++) {
			std::cout << ppvs[i][j].size() << " ";
		}
		std::cout << std::endl;
	}
}

//The second strategy
void icv2ppvs_2(const opts_t &opts, const int &num_meshes, const vector<int> &icv_num_for_each, const vector<input_corr_vector> &icv, const vector<bool> &use_points,
	vector< vector<point_pair_vector> > &ppvs) {

		std::cout << "original : " << icv.size() << std::endl;

		ppvs.resize(num_meshes, vector<point_pair_vector>(num_meshes));

		int unused = 0, unstable = 0, unallowed = 0;

		//The second strategy
		for (int i = 0, start =0; i < icv_num_for_each.size(); i++) {
			for (int j = start; j - start < icv_num_for_each[i]; j++) {
				if (!use_points[j]) {
					unused++;
					continue;
				}
				int mesh_index_x = i;
				point x;
				for (int k = 0; k < icv[j].size(); k++) {
					if (icv[j][k].tgt == mesh_index_x) {
						x = icv[j][k].p;
						break;
					}
				}
				for(int k = 0; k < icv[j].size(); k++) {
					if(icv[j][k].tgt == mesh_index_x) continue;
					if (icv[j][k].status == input_corr::STABLE) {
						int mesh_index_y = icv[j][k].tgt;
						point y = icv[j][k].p;
						float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
						if(dist < 500 ){
							if (mesh_index_x < mesh_index_y) {
								point_pair pp(x, y);
								ppvs[mesh_index_x][mesh_index_y].push_back(pp);
							} 
							else{
								point_pair pp(y, x);
								ppvs[mesh_index_y][mesh_index_x].push_back(pp);
							}
						}
						else unallowed ++;
					}else unstable++;
				}
			}
			start += icv_num_for_each[i];
		}

		std::cout << "unused : " << unused << std::endl;
		std::cout << "unstable : " << unstable << std::endl;
		std::cout << "unallowed : " << unallowed << std::endl;

		for (int i = 0; i < num_meshes; i++) {
			std::cout << i << ": ";
			for (int j = 0; j < num_meshes; j++) {
				std::cout << ppvs[i][j].size() << " ";
			}
			std::cout << std::endl;
		}
}

//The third strategy
void icv2ppvs_3(const opts_t &opts, const int &num_meshes, const vector<input_corr_vector> &icv, const vector<bool> &use_points,
	vector< vector<point_pair_vector> > &ppvs) {

		std::cout << "original : " << icv.size() << std::endl;

		ppvs.resize(num_meshes, vector<point_pair_vector>(num_meshes));

		int unused = 0, unstable = 0, unallowed = 0;

		//The third strategy
		for (int i = 0; i < icv.size(); i++) {
			if (use_points[i]) {
				vector<int> stable_indices;
				for (int j = 0; j < icv[i].size(); j++) {
					if (icv[i][j].status == input_corr::STABLE) {
						stable_indices.push_back(j);
					} else unstable++;
				}
				srand(time(NULL));
				for (int i = stable_indices.size()-1; i>=1; i--) {
					int j = rand() % (i+1);
					int temp = stable_indices[j];
					stable_indices[j] = stable_indices[i];
					stable_indices[i] = temp;
				}

				if (stable_indices.size() >= 2) {
					for (int j = 0; j < stable_indices.size()-1; j++) {
						int first = stable_indices[j];
						int second = stable_indices[j+1];
						int mesh_index_x = icv[i][first].tgt;
						int mesh_index_y = icv[i][second].tgt;
						point x = icv[i][first].p;
						point y = icv[i][second].p;
						float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
						if(dist < 500 ){
							if (mesh_index_x < mesh_index_y) {
								point_pair pp(x, y);
								ppvs[mesh_index_x][mesh_index_y].push_back(pp);
							} else{
								point_pair pp(y, x);
								ppvs[mesh_index_y][mesh_index_x].push_back(pp);
							}
						} else unallowed ++;
					}	
				}
			} else unused++;
		}

		std::cout << "unused : " << unused << std::endl;
		std::cout << "unstable : " << unstable << std::endl;
		std::cout << "unallowed : " << unallowed << std::endl;

		for (int i = 0; i < num_meshes; i++) {
			std::cout << i << ": ";
			for (int j = 0; j < num_meshes; j++) {
				std::cout << ppvs[i][j].size() << " ";
			}
			std::cout << std::endl;
		}
}

//The fourth strategy
void icv2ppvs_4(const opts_t &opts, const int &num_meshes, const vector<input_corr_vector> &icv, const vector<bool> &use_points,
	vector< vector<point_pair_vector> > &ppvs) {

		std::cout << "original : " << icv.size() << std::endl;

		ppvs.resize(num_meshes, vector<point_pair_vector>(num_meshes));

		int unused = 0, unstable = 0, unallowed = 0;

		//Fourth strategy
		for (int i = 0; i < icv.size(); i++) {
			if (use_points[i]) {
				for (int j = 0; j < icv[i].size()-1; j++) {
					if (icv[i][j].status == input_corr::STABLE) {
						for(int k = j+1; k < icv[i].size(); k++) {
							if (icv[i][k].status == input_corr::STABLE) {
								int mesh_index_x = icv[i][j].tgt;
								int mesh_index_y = icv[i][k].tgt;
								point x = icv[i][j].p;
								point y = icv[i][k].p;
								float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
								if(dist < 500 ) {
									if (mesh_index_x < mesh_index_y) {
										point_pair pp(x, y);
										ppvs[mesh_index_x][mesh_index_y].push_back(pp);
									} else{
										point_pair pp(y, x);
										ppvs[mesh_index_y][mesh_index_x].push_back(pp);
									}
								}
								else unallowed++;
							}
						}
					} else unstable++;
				}
			}else unused++;
		}

		std::cout << "unused : " << unused << std::endl;
		std::cout << "unstable : " << unstable << std::endl;
		std::cout << "unallowed : " << unallowed << std::endl;

		for (int i = 0; i < num_meshes; i++) {
			std::cout << i << ": ";
			for (int j = 0; j < num_meshes; j++) {
				std::cout << ppvs[i][j].size() << " ";
			}
			std::cout << std::endl;
		}
}

void generate_alpha_beta_m_na_nb_total_num_corr_ApproxiMultiTPS(const opts_t &opts, const int &num_meshes, const vector< vector<point_pair_vector> > &ppvs,
	vector<int> &alpha, vector<int> &beta, vector<int> &m, vector<int> &na, vector<int> &nb, int &total_num_corr) {
		
	total_num_corr = 0;

	for (int i = 0; i < num_meshes-1; i++) {
		for (int j = i+1; j < num_meshes; j++) {
			if (!ppvs[i][j].empty() || !ppvs[j][i].empty()) {
				alpha.push_back(i);
				beta.push_back(j);
				na.push_back(ppvs[i][j].size());
				nb.push_back(ppvs[j][i].size());
				m.push_back(ppvs[i][j].size() + ppvs[j][i].size());
				total_num_corr += m.back();
			}
		}
	}
}

void generate_alpha_beta_m_na_nb_total_num_corr_1234(const opts_t &opts, const int &num_meshes, const vector< vector<point_pair_vector> > &ppvs,
	vector<int> &alpha, vector<int> &beta, vector<int> &m, vector<int> &na, vector<int> &nb, int &total_num_corr) {

		total_num_corr = 0;

		int min_each_num_corr = opts.min_each_num_corr;
		int max_each_num_corr = opts.max_each_num_corr;

		for (int i = 0; i < num_meshes-1; i++) {
			for (int j = i+1; j < num_meshes; j++) {
				if (!ppvs[i][j].empty() && ppvs[i][j].size() >= min_each_num_corr) {
					alpha.push_back(i);
					beta.push_back(j);
					if (ppvs[i][j].size() <= max_each_num_corr) {
						m.push_back(ppvs[i][j].size());
						total_num_corr += m.back();
					} else {
						m.push_back(max_each_num_corr);
						total_num_corr += m.back();
					}	
				}
			}
		}
}

void ppvs2XY_ApproxiMultiTPS(const opts_t &opts, const int &num_meshes, const vector< vector<point_pair_vector> > &ppvs, const int &total_num_corr,
	gmnr::PointSet3D &input_X, gmnr::PointSet3D &input_Y) {

	input_X.resize(total_num_corr, 3);
	input_Y.resize(total_num_corr, 3);

	int current_num_corr = 0;

	for (int i = 0; i < num_meshes-1; i++) {
		for (int j = i+1; j < num_meshes; j++) {
			if (!ppvs[i][j].empty()  || !ppvs[j][i].empty()) {
				for (int k = 0; k < ppvs[i][j].size(); k++) {
					input_X(current_num_corr, 0) = ppvs[i][j][k].x[0];
					input_X(current_num_corr, 1) = ppvs[i][j][k].x[1];
					input_X(current_num_corr, 2) = ppvs[i][j][k].x[2];

					input_Y(current_num_corr, 0) = ppvs[i][j][k].y[0];
					input_Y(current_num_corr, 1) = ppvs[i][j][k].y[1];
					input_Y(current_num_corr, 2) = ppvs[i][j][k].y[2];

					current_num_corr++;
				}
				for (int k = 0; k < ppvs[j][i].size(); k++) {
					input_X(current_num_corr, 0) = ppvs[j][i][k].y[0];
					input_X(current_num_corr, 1) = ppvs[j][i][k].y[1];
					input_X(current_num_corr, 2) = ppvs[j][i][k].y[2];

					input_Y(current_num_corr, 0) = ppvs[j][i][k].x[0];
					input_Y(current_num_corr, 1) = ppvs[j][i][k].x[1];
					input_Y(current_num_corr, 2) = ppvs[j][i][k].x[2];

					current_num_corr++;
				}
			}
		}
	}
}

void ppvs2XY_1234(const opts_t &opts, const int &num_meshes, const vector< vector<point_pair_vector> > &ppvs, const int &total_num_corr,
	gmnr::PointSet3D &input_X, gmnr::PointSet3D &input_Y) {

		int min_each_num_corr = opts.min_each_num_corr;
		int max_each_num_corr = opts.max_each_num_corr;

		input_X.resize(total_num_corr, 3);
		input_Y.resize(total_num_corr, 3);

		int current_num_corr = 0;

		for (int i = 0; i < num_meshes-1; i++) {
			for (int j = i+1; j < num_meshes; j++) {
				if (!ppvs[i][j].empty()  && ppvs[i][j].size() >= min_each_num_corr) {
					if (ppvs[i][j].size() <= max_each_num_corr) {
						for (int k = 0; k < ppvs[i][j].size(); k++) {
							input_X(current_num_corr, 0) = ppvs[i][j][k].x[0];
							input_X(current_num_corr, 1) = ppvs[i][j][k].x[1];
							input_X(current_num_corr, 2) = ppvs[i][j][k].x[2];

							input_Y(current_num_corr, 0) = ppvs[i][j][k].y[0];
							input_Y(current_num_corr, 1) = ppvs[i][j][k].y[1];
							input_Y(current_num_corr, 2) = ppvs[i][j][k].y[2];

							current_num_corr++;
						}
					} else {
						vector<int> indices(ppvs[i][j].size());
						for (int i = 0; i < indices.size(); i++) {
							indices[i] = i;
						}
						srand(time(NULL));
						for (int i = indices.size()-1; i>=1; i--) {
							int j = rand() % (i+1);
							int temp = indices[j];
							indices[j] = indices[i];
							indices[i] = temp;
						}

						for (int k = 0; k < max_each_num_corr; k++) {
							input_X(current_num_corr, 0) = ppvs[i][j][indices[k]].x[0];
							input_X(current_num_corr, 1) = ppvs[i][j][indices[k]].x[1];
							input_X(current_num_corr, 2) = ppvs[i][j][indices[k]].x[2];

							input_Y(current_num_corr, 0) = ppvs[i][j][indices[k]].y[0];
							input_Y(current_num_corr, 1) = ppvs[i][j][indices[k]].y[1];
							input_Y(current_num_corr, 2) = ppvs[i][j][indices[k]].y[2];

							current_num_corr++;
						}
					}
				}
			}
		}
}

void write_cube_for_XY(const opts_t &opts, const vector<char *> &mesh_names, TriMesh *points_mesh,
	const gmnr::PointSet3D &input_X, const gmnr::PointSet3D &input_Y, const vector<int> &m, const vector<int> &alpha, const vector<int> &beta) {
	int start = 0;
	for (int i = 0; i < m.size(); i++) {
		points_mesh->vertices.resize(8 * m[i]);
		points_mesh->faces.resize(12 * m[i]);
		for (int j = 0; j  < m[i]; j++) {
			point p;
			p[0] = input_X(j + start, 0);
			p[1] = input_X(j + start, 1);
			p[2] = input_X(j + start, 2);
			write_cube_vertex(&points_mesh->vertices[8 * j], p, opts.cube_size);
			write_cube_face(points_mesh->faces, 12 * j, 8 * j);
		}
		char points_name[1024];
		sprintf(points_name, "points/points_%d-(%d_%d)_(%s_%s)_%d.ply", i, alpha[i], beta[i], get_mesh_name(mesh_names[alpha[i]]), get_mesh_name(mesh_names[beta[i]]), 1);
		points_mesh->write(points_name);
		start += m[i];
	}

	start = 0;
	for (int i = 0; i < m.size(); i++) {
		points_mesh->vertices.resize(8 * m[i]);
		points_mesh->faces.resize(12 * m[i]);
		for (int j = 0; j < m[i]; j++) {
			point p;
			p[0] = input_Y(j + start, 0);
			p[1] = input_Y(j + start, 1);
			p[2] = input_Y(j + start, 2);
			write_cube_vertex(&points_mesh->vertices[8 * j], p, opts.cube_size);
			write_cube_face(points_mesh->faces, 12 * j, 8 * j);
		}
		char points_name[1024];
		sprintf(points_name, "points/points_%d-(%d_%d)_(%s_%s)_%d.ply", i, alpha[i], beta[i], get_mesh_name(mesh_names[alpha[i]]), get_mesh_name(mesh_names[beta[i]]), 2);
		points_mesh->write(points_name);
		start += m[i];
	}
}

void write_view0_target0_xyz(const opts_t &opts, const vector<input_corr_vector> &icv, const vector<point> &targets, const vector<bool> &use_points) {
	vector<point> viewX, viewY;

	int orignial_size = 0;
	int unused = 0;
	int unstable = 0;

	for (int i = 0; i < icv.size(); i++) {
		for (int j = 0; j < icv[i].size(); j++){ 
			if(icv[i][j].tgt == 0) {
				orignial_size++;
				if(icv[i][j].status == input_corr::STABLE) {if(!use_points[i]) unused++;}
				else unstable++;
			}
		}
		if (use_points[i]) {
			for (int j = 0; j < icv[i].size(); j++) {
				if (icv[i][j].status == input_corr::STABLE && icv[i][j].tgt == 0) {
					point x = icv[i][j].p;
					point y = targets[i];
					point_pair pp(x, y);
					float dist = SQ(x[0] - y[0]) + SQ(x[1] - y[1]) + SQ(x[2] - y[2]);
					if(dist < 500 )	{
						viewX.push_back(x);
						viewY.push_back(y);
					}				
				}
			}
		}
	}

	std::cout << "Original size : "<< orignial_size << std::endl;
	std::cout << "Unused size : "<< unused << std::endl;
	std::cout << "Unstable size : "<< unstable << std::endl;

	std::fstream fs_viewX("viewX2.xyz", std::ios::out);
	if(fs_viewX) {
		for(int p = 0; p < viewX.size(); p++) {
			for (int q = 0; q < 3; q++) {
				fs_viewX << viewX[p][q] << " ";
			}
			fs_viewX << std::endl;
		}
	}
	fs_viewX.close();
	std::fstream fs_viewY("viewY2.xyz", std::ios::out);
	if(fs_viewY) {
		for(int p = 0; p < viewY.size(); p++) {
			for (int q = 0; q < 3; q++) {
				fs_viewY << viewY[p][q] << " ";
			}
			fs_viewY << std::endl;
		}
	}
	fs_viewY.close();

	exit(0);
}

void align_scan2(const opts_t &opts, const vector<char *> &mesh_names, const int &num_meshes, const vector<int> &icv_num_for_each,
	const vector<input_corr_vector> &icv, const vector<point> &targets, const vector<bool> &use_points, const vector<float> &confidence,
	TriMesh *points_mesh, const vector<TriMesh::Face> &all_faces, pthread_mutex_t *mutex) {

		//write_view0_target0_xyz(opts, icv, targets, use_points);

		vector< vector<point_pair_vector> > ppvs;
		icv2ppvs_ApproxiMultiTPS(opts, num_meshes, icv, use_points, ppvs);

		vector<int> alpha, beta, m, na, nb;
		int total_num_corr = 0;
		generate_alpha_beta_m_na_nb_total_num_corr_ApproxiMultiTPS(opts, num_meshes, ppvs, alpha, beta, m, na, nb, total_num_corr);

		gmnr::PointSet3D input_X, input_Y;
		ppvs2XY_ApproxiMultiTPS(opts, num_meshes, ppvs, total_num_corr, input_X, input_Y);

		write_cube_for_XY(opts, mesh_names, points_mesh, input_X, input_Y, m, alpha, beta);

		vector<int> corr_for_each_view(num_meshes, 0);
		int total_num = 0;
		for (int i = 0; i < m.size(); i++) {
			corr_for_each_view[alpha[i]] += m[i];
			corr_for_each_view[beta[i]] += m[i];
			total_num += m[i];
		}

		gmnr::Vector kappa(num_meshes);
		for (int i = 0; i < num_meshes; i++) {
			//kappa[i] = opts.kappa * corr_for_each_view[i] * 1.0 / total_num;
			kappa[i] = 0.0001;
		}
		gmnr::Vector lambda(num_meshes);
		for (int i = 0; i < num_meshes; i++) {
			//lambda[i] = opts.lambda * corr_for_each_view[i] * 1.0 / total_num;
			lambda[i] = 0.00001;
		}

		//gmnr::MultiTPS mtps(input_X, input_Y, m, alpha, beta, kappa, lambda, 10);
		gmnr::ApproxiMultiTPS mtps(input_X, input_Y, m, alpha, beta, kappa, lambda, na, nb, 4);

		for (int i = 0; i < num_meshes; i++) {
			const char *mesh_name = mesh_names[i];

			fprintf(stderr, "Opening %s\n", mesh_name);

			xform xfin;
			if (opts.read_xf) {
				// the assumes the mesh has a three-letter extension
				char xfname[1024];
				strcpy(xfname, mesh_name);
				strcpy(xfname + strlen(xfname) - 3, "xf");
				xfin.read(xfname);
				fprintf(stderr, "Read %s:\n", xfname);
				std::cerr << xfin;
			}
			TriMesh *mesh = NULL;

			if (opts.nonrigid_prefix || !opts.write_xf) {
				mesh = TriMesh::read(mesh_name);
				mesh->normals.clear();
				assert(mesh);
			}

			char out_name[1024]; // mesh name for output

			// make a backup of the mesh vertices, because the
			// various warps will clobber them
			farr verts;
			if (mesh) {
				verts = farr(mesh->vertices.size(), 3);
				for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
					point p = xfin * mesh->vertices[j];
					for (int k = 0; k < 3; k++) verts[j][k] =  p[k];
				}
			}

			if (opts.nonrigid_prefix) {
				// do TPS warp

				gmnr::TPSFunction f = mtps.getfs()[i];
				
				gmnr::PointSet3D original_points(mesh->vertices.size(), 3);
				for (int j = 0; j < mesh->vertices.size(); j++) {
					original_points(j, 0) = verts[j][0];
					original_points(j, 1) = verts[j][1];
					original_points(j, 2) = verts[j][2];
				}
				
				gmnr::PointSet3D new_points = f.evaluate(original_points);

				if (opts.nonrigid_prefix) {

					for (unsigned int j = 0; j < mesh->vertices.size(); j++) {
						mesh->vertices[j][0] = new_points(j, 0);
						mesh->vertices[j][1] = new_points(j, 1);
						mesh->vertices[j][2] = new_points(j, 2);
					}

					sprintf(out_name, "%s%s", opts.nonrigid_prefix, get_mesh_name(mesh_name));
					mesh->write(out_name);
				}
			}

			if (mesh) delete mesh;
		}
}

class read_corrs {
public:
  read_corrs(const opts_t &o) : opts(o) {
    points_mesh = new TriMesh;

    // read the input in either old or new style
    read_corr_points();

    errors.resize(num_points);
  }
  
  int num_meshes, num_sources, num_points;
  vector<int> offsets, num_corrs;
  vector<bool> use_points;
  vector<float> confidence, errors;
  vector<point> targets, orig_targets;
  vector<corr_vector> corrs;
  vector<spring_vector> springs;

  vector<char *> mesh_names;
  TriMesh *points_mesh;

  vector<input_corr_vector> icv;
  vector<int> icv_num_for_each;

private:
  const opts_t &opts;

  // read the new-style correspondence input from corr*.txt
  void read_corr_points(void) {
    vector<char *> corr_files;

    // header: number of points and sources
    scanf(" %d %d ", &num_meshes, &num_sources);

    // get list of mesh names
    for (int i = 0; i < num_meshes; i++) {
      char fname[1024];
      scanf(" %s ", fname);
      mesh_names.push_back(strdup(fname));
    }

    // get list of corr files
    for (int i = 0; i < num_meshes; i++) {
      char fname[1024];
      scanf(" %s ", fname);
      corr_files.push_back(strdup(fname));
    }

    offsets.resize(num_meshes + 1);

    // process each source in turn
    num_points = 0;
    confidence.clear();
    for (int i = 0; i < num_sources; i++) {
      char mesh_name[1024];
      int num_scan_points = 0;

      FILE *cf = fopen(corr_files[i], "r");

      // get the number of features on this source, and the mesh name
      fscanf(cf, "%d %s", &num_scan_points, mesh_name);
      strcpy(mesh_name + strlen(mesh_name) - 3, "ply");

      // every source should also be a target, but we don't know which one,
      // so we need to scan the list of targets to figure out
      int offset_num = 0;
      while (offset_num < num_meshes && strcmp(mesh_name, mesh_names[offset_num])) offset_num++;
      assert(offset_num < num_meshes); // true if source is also a target
      offsets[offset_num + 1] = num_scan_points;

	  icv_num_for_each.push_back(num_scan_points);

      // read in the list of correspondeces for this source
      for (int j = 0; j < num_scan_points; j++, num_points++) {
        int v, nc;
        point p;

        // Header tells us where this point was to start, and the number
        // of correspondences.  It also tells us the vertex in the source
        // mesh where it was, but we don't really care about that since
        // we have the coordinates.
        fscanf(cf, " %f %f %f %d %d ", &p[0], &p[1], &p[2], &v, &nc);
        orig_targets.push_back(p);
        icv.push_back(input_corr_vector());

        int num_stable = 0, num_enabled = 0;
        for (int k = 0; k < nc; k++) {
          input_corr ic;
          fscanf(cf, " %d %f %f %f %f %f %f ", &ic.tgt, &ic.p[0], &ic.p[1], &ic.p[2],
                 &ic.stability, &ic.max_stability, &ic.err);
          ic.status = input_corr::STABLE;

          // here we do some pruning of springs based on error and stability thresholds
          if (ic.err > opts.max_allowed_err) ic.status = input_corr::DISABLED;
          // if (ic.err > opts.max_allowed_err) ic.status = input_corr::UNSTABLE;
          else if ((ic.stability < opts.min_allowed_stability) ||
                   (ic.max_stability < opts.min_allowed_max_stability))
            ic.status = input_corr::UNSTABLE;
          assert(num_points < (int) icv.size());
          icv[num_points].push_back(ic);
          num_stable += (ic.status == input_corr::STABLE);
          num_enabled += (ic.status != input_corr::DISABLED);
        }
        confidence.push_back((float) num_stable / (float) num_enabled);
      }

      fclose(cf);
    }

    fprintf(stderr, "Loaded input corerspondences (%d points).\n", (int) num_points);

    //
    // Do any trimming of springs here
    //

    fprintf(stderr, "Done pruning.\n");

    // make the offsets array cumulative
    num_corrs.resize(num_points);
    targets.resize(num_points);
    // confidence.resize(num_points);
    use_points.resize(num_points);

    corrs.resize(num_meshes);
    springs.resize(num_points);
    errors.resize(num_points);

    for (int i = 1; i <= num_meshes; i++) offsets[i] += offsets[i - 1];

    for (int i = 0; i < num_points ; i++) {
      targets[i] = point(0, 0, 0); num_corrs[i] = 0;
      // targets[i] = orig_targets[i]; num_corrs[i] = 1;
      use_points[i] = (confidence[i] >= CONF_THRESH);
      use_points[i] = true;
      if (confidence[i] < CONF_THRESH) {
        conf_points++;
        // fprintf(stderr, "CONF_THRESH pruned point %d\n", i);
      }
    }

    fprintf(stderr, "Total CONF_THRESH pruning, %d\n", conf_points);

    // now write out a mesh with little cubes around each of the points
    points_mesh->vertices.resize(8 * num_points);
    points_mesh->faces.resize(12 * num_points);
    for (int i = 0; i < num_points; i++) {
      if (!use_points[i]) continue;
      write_cube_vertex(&points_mesh->vertices[8 * i], orig_targets[i], opts.cube_size);
      write_cube_face(points_mesh->faces, 12 * i, 8 * i);
    }
    points_mesh->write("points/orig_points.ply");

    std::vector<point> camera_centers;
    if (opts.camera_centers_file) {
      FILE *ccf = fopen(opts.camera_centers_file, "r");
      assert(ccf);
      camera_centers.resize(num_meshes);

      for (int i = 0; i < num_meshes; i++) {
        assert(fscanf(ccf, " %f %f %f ", &camera_centers[i][0],
                      &camera_centers[i][1], &camera_centers[i][2]) == 3);
      }
      fclose(ccf);
      fprintf(stderr, "Camera centers: %d\n", (int) camera_centers.size());
    }

    int unstable_pruning = 0;
    for (int i = 0, src_i = 0; i < num_points; i++) {
      int stable_corrs = 0, unstable_corrs = 0;
      point &ot = orig_targets[i];

      // this incrementally computes the source mesh
      // for point i
      // int old_src_i = src_i;
      while (offsets[src_i + 1] <= i) src_i++;

#if 0
      if (src_i > old_src_i)
        fprintf(stderr, "Adding springs for mesh %d\n", src_i);
#endif

      std::vector<point> rays, centers;
      for (unsigned int j = 0; j < icv[i].size(); j++) {
        if (icv[i][j].status == input_corr::DISABLED) continue; // skipped pruned correspondences
        if (icv[i][j].status == input_corr::UNSTABLE) continue; // also skip unstable ones for spring economy
        point &p  = icv[i][j].p;

        bool stable = (icv[i][j].status == input_corr::STABLE);
#if 1
        float divergence = dist(ot, icv[i][j].p);
        if (divergence > opts.max_allowed_divergence) {
          fprintf(stderr, "Bad correspondence?: %4d -> %4d (%d), %.3f %.3f %.3f -> %.3f %.3f %.3f\n",
                  i, icv[i][j].tgt, src_i, ot[0], ot[1], ot[2], p[0], p[1], p[2]);
		  icv[i][j].status = input_corr::UNSTABLE;
          stable = false;
          unstable_corrs++;
          continue;
        } else {
#if 0
          fprintf(stderr, "Good correspondence: %4d -> %4d (%d), %.3f %.3f %.3f -> %.3f %.3f %.3f\n",
                  i, icv[i][j].tgt, src_i, ot[0], ot[1], ot[2], p[0], p[1], p[2]);
#endif
        }
#endif

        // icv[i] tells us where feature point i lives on every mesh
        // We need to invert this to get corrs, which tells us
        // which feature points live on each mesh.  That means that
        // corrs is indexed by icv[i][j].tgt, which is the mesh associated
        // with correspondence j of feature point i
        corr_t c(src_i, i, p, stable);
        corrs[icv[i][j].tgt].push_back(c);
        stable_corrs += stable; unstable_corrs += !stable;

        // fprintf(stderr, "Created correspondence: %4d -> %4d %f %f %f %d\n", icv[i][j].tgt, i, p[0], p[1], p[2], stable);

        if (!stable) continue;

        // Remember, targets is the array of feature points, which will
        // ultimately be the target locations for the TPS transform.
        if (opts.camera_centers_file) {
          centers.push_back(camera_centers[icv[i][j].tgt]);
          rays.push_back(icv[i][j].p - centers[num_corrs[i]]);
        }
        targets[i] += p;
        num_corrs[i]++;

        // add springs to all features/targets on mesh src_i
        for (int k = offsets[icv[i][j].tgt]; k < offsets[icv[i][j].tgt + 1]; k++) {
          // if (k == icv[i][j].tgt) continue;
          if (!use_points[k]) continue;
          if (k == i) continue; // skip springs to yourself

          float len = dist(p, orig_targets[k]);
          if (len > opts.max_spring_len) continue;

          float wgt = 1 / (opts.wgt_epsilon + len );

          // see if we have this spring already (from the other alignment direction)
          bool found_spring = false;
          bool stable_spring = stable;
          for (unsigned int l = 0; l < springs[i].size(); l++) {
            if (springs[i][l].tgt == k) {
              stable_spring |= springs[i][l].stable;
              bool combine = !(stable ^ springs[i][l].stable);
              len = combine ? (len + springs[i][l].len) / 2.0f : (stable ? len : springs[i][l].len);
              wgt = combine ? (wgt + springs[i][l].wgt) / 2.0f : (stable ? wgt : springs[i][l].wgt);

              springs[i][l].stable = stable_spring;
              springs[i][l].len = len;
              springs[i][l].wgt = wgt;

              found_spring = true;
              break;
            }
          }

          if (found_spring) {
            // if we found the spring in this direction, it exists in the other direction too
            bool found_opp_spring = false;
            for (unsigned int l = 0; l < springs[k].size(); l++) {
              if (springs[k][l].tgt == i) {
                springs[k][l].stable = stable_spring;
                springs[k][l].len = len;
                springs[k][l].wgt = wgt;
                found_opp_spring = true;
                break;
              }
            }
            assert(found_opp_spring);
          } else {
            // the spring didn't exist, so add it in both directions
            springs[k].push_back(spring_t(i, len, wgt, stable));
            springs[i].push_back(spring_t(k, len, wgt, stable));

            // /* if (len < 10) */ fprintf(stderr, "Spring: %d %d %f %f %c\n", k, i, len, wgt, stable ? 't' : 'f');
          }
        }
      }

      if (opts.camera_centers_file && num_corrs[i] >= 2) {
        // compute intersection of camera rays
        farr A(3 * num_corrs[i], num_corrs[i] + 3, 0.0);
        farr B(3 * num_corrs[i], 1);
        for (int eq = 0; eq < num_corrs[i]; eq++) {
          // construct the linear system for solving the ray intersection
          for (int var = 0; var < 3; var++) {
            int row = 3 * eq + var;
            A[row][eq] = rays[eq][var];
            A[row][num_corrs[i] + var] = 1;
            B[row][0] = centers[eq][var];
          }
        }

        // now solve the system
        // fprintf(stderr, "%d %d\n", i, eq);
        // fprintf(stderr, "A =\n");
        // print_matrix(A);
        // fprintf(stderr, "B =\n");
        // print_matrix(B);

        JAMA::LU<matrix_t> lu(A);
        farr res = lu.solve(B);
        assert(res.dim1() == B.dim1());
        assert(res.dim2() == B.dim2());
        targets[i] = point(res[num_corrs[i]][0], res[num_corrs[i] + 1][0], res[num_corrs[i] + 2][0]);
      } else {
        targets[i] /= num_corrs[i];
      }
      if (stable_corrs < max(2, unstable_corrs / 4)) {
        fprintf(stderr, "Point %d pruned for lack of stable correspondences.\n", (int) i);
        unstable_pruning++;
        use_points[i] = false;
      }
    }

    fprintf(stderr, "Total unstable pruning: %d\n", unstable_pruning);

    fprintf(stderr, "Added springs\n");
  }
};

// see opts_t for options
int main(int argc, char *argv[]) {
  vector<char *> global_names;

  opts_t opts(argc, argv);

  read_corrs rc(opts);
  vector<TriMesh::Face> all_faces = rc.points_mesh->faces;

  if (opts.write_individual_points) {
    // write out meshes of locations of each point on its constituent meshes
    for (int n = 0; n < rc.num_points; n++) {
      rc.points_mesh->vertices.clear();
      rc.points_mesh->faces.clear();
      int v = 0, f = 0;

      for (int i = 0; i < rc.num_meshes; i++) {
        for (unsigned int j = 0; j < rc.corrs[i].size(); j++) {
          if (rc.corrs[i][j].tgt != n) continue;
          if (rc.corrs[i][j].stable == false) continue;
          rc.points_mesh->vertices.resize(v + 8);
          rc.points_mesh->faces.resize(f + 12);
          write_cube_vertex(&rc.points_mesh->vertices[v], rc.corrs[i][j].pnt, opts.cube_size);
          write_cube_face(rc.points_mesh->faces, f, v);
          v += 8; f += 12;
        }
      }

      char fname[1024];
      sprintf(fname, "points/p_%08di.ply", n);
      rc.points_mesh->write(fname);
    }
  }

  if (opts.write_individual_points) {
    // write out average point locations individually
    rc.points_mesh->vertices.resize(8);
    rc.points_mesh->faces.resize(12);
    for (int i = 0; i < rc.num_points; i++) {
      if (!rc.use_points[i]) continue;
      write_cube_vertex(&rc.points_mesh->vertices[0], rc.targets[i], opts.cube_size);
      write_cube_face(rc.points_mesh->faces, 0, 0);

      char fname[1024];
      sprintf(fname, "points/p_%08da.ply", i);
      rc.points_mesh->write(fname);
    }
  }

  rc.points_mesh->vertices.resize(8 * rc.num_points);
  rc.points_mesh->faces.resize(12 * rc.num_points);
  for (int i = 0; i < rc.num_points; i++) {
    if (!rc.use_points[i]) continue;
    write_cube_vertex(&rc.points_mesh->vertices[8 * i], rc.targets[i], opts.cube_size);
    write_cube_face(rc.points_mesh->faces, 12 * i, 8 * i);
  }
  rc.points_mesh->write("points/avg_points.ply");

  timestamp opt_start = now();

  if (!opts.suppress_optimization) {
	  for (int i = 0; i < 1; i++) move_points(rc.targets, rc.use_points, rc.errors, rc.springs, false);
	  for (int i = 0; i < 1; i++) move_points(rc.targets, rc.use_points, rc.errors, rc.springs, true);
  }
  int f = 0;

#if 1
  thin_points(opts.min_target_dist2, rc.targets, rc.use_points, rc.errors, rc.springs);

  f = 0;
  rc.points_mesh->faces.resize(12 * rc.num_points);
  for (int i = 0; i < rc.num_points; i++) {
	  if (rc.use_points[i]) continue;
	  write_cube_face(rc.points_mesh->faces, f,  8 * i);
	  f += 12;
  }
  rc.points_mesh->faces.resize(f);
  rc.points_mesh->write("points/thin_points.ply");
#endif

#if 1
  for (int i = 0; i < rc.num_meshes; i++)
	  prune_big_moving_points(opts.min_target_dist2, rc.corrs[i], rc.targets, rc.use_points);

  if (!opts.suppress_optimization) {
    for (int i = 0; i < 1; i++) move_points(rc.targets, rc.use_points, rc.errors, rc.springs, true);
  }
#endif
  float opt_time = now() - opt_start;
  fprintf(stderr, "Optimization time: %f\n", opt_time);

  vector<bool> new_state = rc.use_points;
  if (opts.disable_pts_file) {
    // remove the selected points
    int pt, state;
    while (fscanf(opts.disable_pts_file, " %d %d ", &pt, &state) == 2) {
      new_state[pt] = !state;
    }
  }

  if (opts.pts_file) {
    for (int i = 0; i < rc.num_points; i++) {
      fprintf(opts.pts_file, "%d %f %f %f %u %f ", (int) rc.use_points[i],
              rc.targets[i][0], rc.targets[i][1], rc.targets[i][2],
              (unsigned int) rc.springs[i].size(), rc.errors[i]);
      for (unsigned int j = 0; j < rc.springs[i].size(); j++) {
        int color = 0;
        if (rc.use_points[i] && rc.use_points[rc.springs[i][j].tgt] && rc.springs[i][j].stable) {
          float len2 = sqr(rc.springs[i][j].len);
          float d2 = dist2(rc.targets[i], rc.targets[rc.springs[i][j].tgt]);
          color = (d2 > 1.1f * len2) ? 1 : ((d2 < 0.9f * len2) ? 2 : 3);
        }
        fprintf(opts.pts_file, " %d %g %d", rc.springs[i][j].tgt,
                rc.springs[i][j].len, color);
      }
      putc('\n', opts.pts_file);
    }
  }

  rc.use_points = new_state;

  f = 0;
  rc.points_mesh->faces.resize(12 * rc.num_points);
  for (int i = 0; i < rc.num_points; i++) {
	  if (rc.use_points[i]) continue;
	  write_cube_vertex(&rc.points_mesh->vertices[8 * i], rc.targets[i], opts.cube_size);
	  write_cube_face(rc.points_mesh->faces, f, 8 * i);
	  f += 12;
  }
  rc.points_mesh->faces.resize(f);
  rc.points_mesh->write("points/disabled_points.ply");

  f = 0;
  rc.points_mesh->faces.resize(12 * rc.num_points);
  for (int i = 0; i < rc.num_points; i++) {
    if (!rc.use_points[i]) continue;
    write_cube_vertex(&rc.points_mesh->vertices[8 * i], rc.targets[i], opts.cube_size);
    write_cube_face(rc.points_mesh->faces, f, 8 * i);
    f += 12;
  }
  rc.points_mesh->faces.resize(f);
  rc.points_mesh->write("points/final_points.ply");

  if (opts.write_individual_points) {
    // write out final point locations individually
    rc.points_mesh->vertices.resize(8);
    rc.points_mesh->faces.resize(12);
    for (int i = 0; i < rc.num_points; i++) {
      if (!rc.use_points[i]) continue;
      write_cube_vertex(&rc.points_mesh->vertices[0], rc.targets[i], opts.cube_size);
      write_cube_face(rc.points_mesh->faces, 0, 0);

      char fname[1024];
      sprintf(fname, "points/p_%08df.ply", i);
      rc.points_mesh->write(fname);
    }
  }

  if (opts.correspondences_file) {
    // write out the final correspondences
    FILE *f = fopen(opts.correspondences_file, "w");
    if (f) {
      // write out list of target positions and if the feature is used
      fprintf(f, "%d\n", (int) rc.targets.size());
#if 1
      for (int i = 0; i < rc.num_points; i++) {
        fprintf(f, "%d %f %f %f\n", rc.use_points[i] && rc.confidence[i] >= CONF_THRESH,
                rc.targets[i][0], rc.targets[i][1], rc.targets[i][2]);
      }
#endif

      // write out all correspondences, including stability, mesh number, vertex number, and point
      for (int i = 0; i < rc.num_meshes; i++) {
        for (size_t j = 0; j < rc.corrs[i].size(); j++) {
          const corr_t &c = rc.corrs[i][j];
          if (!c.stable) continue;
          if (!rc.use_points[c.tgt]) continue;
#if 1
          fprintf(f, "%d %d %d %d %f %f %f\n", c.stable, c.tgt, i, c.vnum,
                  c.pnt[0], c.pnt[1], c.pnt[2]);
#else
          fprintf(f, "%d %d %f %f %f\n", c.tgt, i, c.pnt[0], c.pnt[1], c.pnt[2]);
#endif
        }
      }

      fclose(f);
    }
  }

  // return 0; // uncomment this for debugging to see just the point optimizations, and not waste time warping

  // do alignment
  if (opts.nthreads == 1) {
	align_scan2(opts, rc.mesh_names, rc.num_meshes, rc.icv_num_for_each, rc.icv, rc.targets, rc.use_points, rc.confidence, rc.points_mesh, all_faces, NULL);
  } else {
    // multithreaded version
    threaded_alignment ta(opts, rc.mesh_names, rc.corrs, rc.targets, rc.use_points, rc.confidence, rc.points_mesh, all_faces);
  }

  delete rc.points_mesh;

  int num_enabled_points = 0;
  for (unsigned int i = 0; i < rc.use_points.size(); num_enabled_points += rc.use_points[i++]);
  fprintf(stderr, "%u %d %d %d %d\n", (unsigned int) rc.targets.size(), num_enabled_points,
         pruned_points, thinned_points, conf_points);

#ifdef COMPUTE_RMS_ERR
  // calculate all pairs of rms err
  for (unsigned int i = 0; i < global_names.size(); i++) {
    TriMesh *tgt = TriMesh::read(global_names[i]);
    tgt->need_neighbors();
    tgt->need_normals();
    tgt->need_adjacentfaces();
    KDtree *tgt_kd = new KDtree(tgt->vertices);
    for (unsigned int j = i + 1; j < global_names.size(); j++) {
      TriMesh *src = TriMesh::read(global_names[j]);
      src->need_neighbors();
      src->need_normals();
      src->need_adjacentfaces();
      double err = rms_err(src, tgt, tgt_kd);
      printf("\n%s -> %s: Global TPS Err: %f\n", global_names[i], global_names[j], err);
      delete src;
    }
    delete tgt_kd;
    delete tgt;
  }
#endif
}
