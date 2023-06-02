#pragma once

#include "npatch.hh"

class GBPatch : public NPatch {
public:
  using Index = std::array<size_t, 4>;
  using Ribbon = std::vector<std::vector<Vector> >;

  GBPatch(std::string filename);
  GBPatch(size_t num_sides, size_t num_layers);
  virtual ~GBPatch();

  virtual void initDomain() override;
  virtual void initDomainMesh(size_t resolution = 60) override;
  virtual Vector evaluateAtParam(double u, double v) const override;
  virtual bool reload() override;

  virtual void draw(const Visualization& vis) const override;
  virtual void drawWithNames(const Visualization& vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector& pos) override;

  static void bernstein(size_t n, double u, DoubleVector& coeff);

  //virtual void setControlPoint(const Index& index, const Vector& p);

  size_t d_;
  std::vector<Ribbon> ribbons_;
  Vector central_cp;
  bool normalized;
};
