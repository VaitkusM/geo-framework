#include "npatch.hh"
#include <functional>

NPatch::NPatch(std::string filename, size_t num_sides) : Object(filename), show_basis_fcn(false), n_(num_sides)
{
  
}

NPatch::~NPatch() {
  
}

void 
NPatch::updateBaseMesh()
{
  mesh = domain_mesh;
  for(auto vv : domain_mesh.vertices()) {
    mesh.point(vv) = evaluateAtParam(vv);
    if (mesh.data(vv).spider_idx[0] == 0 && mesh.data(vv).spider_idx[1] == 0 && mesh.data(vv).spider_idx[2] == 0) {
      std::cerr << std::endl << mesh.point(vv) << std::endl;
    }
  }
  Object::updateBaseMesh(false, false);
}

void 
NPatch::generateSpiderMesh(size_t resolution, BaseMesh &mesh)
{
  mesh.clear();

  std::vector<BaseMesh::VertexHandle> handles;
  size_t meshSize = 1 + n_ * resolution * (resolution + 1) / 2;
  handles.reserve(meshSize);

  // Adding vertices
  Vector center = Vector(0.0, 0.0, 0.0);
  auto vtx = mesh.add_vertex(center);
  handles.push_back(vtx);
  mesh.data(vtx).spider_idx = {0, 0, 0};
  for (size_t j = 1; j <= resolution; ++j) {
    double u = (double)j / (double)resolution;
    for (size_t k = 0; k < n_; ++k) {
      for (size_t i = 0; i < j; ++i) {
        double v = (double)i / (double)j;
        Vector ep = vertices_[prev(k)] * (1.0 - v) + vertices_[k] * v;
        Vector p = center * (1.0 - u) + ep * u;
        vtx = mesh.add_vertex(p);
        handles.push_back(vtx);
        mesh.data(vtx).spider_idx = {i,j,k} ;
      }
    }
  }

  //Adding triangles
  size_t inner_start = 0, outer_vert = 1;
  for (size_t layer = 1; layer <= resolution; ++layer) {
    size_t inner_vert = inner_start, outer_start = outer_vert;
    for (size_t side = 0; side < n_; ++side) {
      size_t vert = 0;
      while (true) {
        size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
        mesh.add_face(handles[inner_vert], handles[outer_vert], handles[next_vert]);
        ++outer_vert;
        if (++vert == layer)
          break;
        size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
        mesh.add_face(handles[inner_vert], handles[next_vert], handles[inner_next]);
        inner_vert = inner_next;
      }
    }
    inner_start = outer_start;
  }
}

double
NPatch::getGBC(double u, double v, size_t i, BarycentricType type) const
{
  std::vector<Vector> vectors; vectors.reserve(n_);
  std::transform(vertices_.begin(), vertices_.end(),
    std::back_inserter(vectors),
    [u, v](const Vector& p)->Vector { return Vector(u, v, 0) - p; });

  DoubleVector areas; areas.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    const Vector& si = vectors[i];
    const Vector& si1 = vectors[next(i)];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n_);

  for (size_t i = 0; i < n_; ++i) {
    size_t i_1 = prev(i), i1 = next(i);
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n_; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector& si_1 = vectors[i_1];
    const Vector& si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    double ri_1 = 1.0, ri = 1.0, ri1 = 1.0;
    switch (type) {
    case BarycentricType::WACHSPRESS:
      break;
    case BarycentricType::MEAN_VALUE:
      ri_1 = vectors[i_1].norm();
      ri = vectors[i].norm();
      ri1 = vectors[i1].norm();
      break;
    case BarycentricType::HARMONIC:
      ri_1 = vectors[i_1].sqrnorm();
      ri = vectors[i].sqrnorm();
      ri1 = vectors[i1].sqrnorm();
      break;
    };
    if (ri < epsilon) {         // at a vertex of the domain (mean/harmonic)
      l.assign(n_, 0.0);
      l[i] = 1.0;
      break;
    }
    l.push_back(ri_1 * Ai_1 + ri1 * Ai - ri * Bi * Ai_1i);
  }

  normalizeValues(l);
  
  if (i >= n_) {
    double prod = 1.0;
    double dl = (vectors[1] - vectors[0]).norm();
    for (auto a : areas) {
      prod *= a/dl;
    }
    return prod;
  }
  else {
    return l[i];
  }
}

void NPatch::normalizeValues(DoubleVector& values)
{
  double sum = std::accumulate(values.begin(), values.end(), 0.0);
  std::transform(values.begin(), values.end(), values.begin(), [sum](double x) { return x / sum; });
}