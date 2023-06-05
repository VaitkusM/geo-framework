#pragma once

#include "object.hh"

class NPatch : public Object
{
public:
  const double epsilon = 1.0e-8;
  using DoubleVector = std::vector<double>;
  enum class BarycentricType { WACHSPRESS, MEAN_VALUE, HARMONIC };

  NPatch(std::string filename, size_t num_sides = 0);
  virtual ~NPatch();
  virtual void updateBaseMesh() override;

  virtual void initDomain() = 0;
  virtual void initDomainMesh(size_t resolution = 60) = 0;
  virtual Vector evaluateAtParam(double u, double v) const = 0;

  void generateSpiderMesh(size_t resolution, BaseMesh &mesh);

  double getGBC(double u, double v, size_t i, BarycentricType type = BarycentricType::WACHSPRESS) const;


  static void normalizeValues(DoubleVector& values);

  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }

  virtual void swapFootpoints() = 0;

  bool show_basis_fcn;
  size_t n_;
  BaseMesh domain_mesh;
  std::vector<Vector> vertices_;
};