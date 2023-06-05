#pragma once

#include "npatch.hh"

#include <map>
#include <array>

class MPatch : public NPatch {
public:
  using Index = std::array<size_t,4>;

  MPatch(std::string filename);
  MPatch(size_t num_sides, size_t depth);
  virtual ~MPatch();

  virtual void initBasicShape() override;
  virtual void initDomain() override;
  virtual void initDomainMesh(size_t resolution = 60) override;
  virtual Vector evaluateAtParam(double u, double v) const override;
  virtual bool reload() override;

  virtual void draw(const Visualization& vis) const override;
  virtual void drawWithNames(const Visualization& vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector& pos) override;

  static size_t binomial(size_t n, size_t k);
  static size_t multinomial(const Index& index);
  static double multiBernstein(const Index& index, const DoubleVector& bc);

  void getBlendFunctions(double u, double v, DoubleVector &values, bool normalize = true) const;
  double getBlendFunctionUnnormalized(double u, double v, const Index &i) const;

  std::vector<Index> neighbors(const Index& si) const;

  virtual void swapFootpoints() override;

  virtual void setControlPoint(const Index& index, const Vector& p);

  size_t d_;
  std::map<Index, Vector> footpoints_;
  std::map<Index, Vector> net_;
};