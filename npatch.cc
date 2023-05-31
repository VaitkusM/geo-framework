#include "npatch.hh"

NPatch::NPatch(std::string filename) : Object(filename)
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

