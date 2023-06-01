#pragma once

#include "object.hh"

class NPatch : public Object
{
public:
  const double epsilon = 1.0e-8;

  NPatch(std::string filename);
  virtual ~NPatch();
  virtual void updateBaseMesh() override;

  virtual void initDomain() = 0;
  virtual void initDomainMesh(size_t resolution = 30) = 0;
  virtual Vector evaluateAtParam(double u, double v) const = 0;

  BaseMesh domain_mesh;
  std::vector<Vector> vertices_;
};