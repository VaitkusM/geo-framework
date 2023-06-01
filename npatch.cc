#include "npatch.hh"

NPatch::NPatch(std::string filename, size_t num_sides) : Object(filename), n_(num_sides)
{
  
}

NPatch::~NPatch() {
  
}

void 
NPatch::updateBaseMesh()
{
  mesh = domain_mesh;
  for(auto vv : domain_mesh.vertices()) {
    auto pt_uv = mesh.point(vv);
    mesh.point(vv) = evaluateAtParam(pt_uv[0], pt_uv[1]);
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
  handles.push_back(mesh.add_vertex(center));
  for (size_t j = 1; j <= resolution; ++j) {
    double u = (double)j / (double)resolution;
    for (size_t k = 0; k < n_; ++k) {
      for (size_t i = 0; i < j; ++i) {
        double v = (double)i / (double)j;
        Vector ep = vertices_[prev(k)] * (1.0 - v) + vertices_[k] * v;
        Vector p = center * (1.0 - u) + ep * u;
        handles.push_back(mesh.add_vertex(p));
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