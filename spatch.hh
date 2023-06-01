#pragma once

#include "npatch.hh"

#include <map>

class SPatch : public NPatch {
public:
  using DoubleVector = std::vector<double>;
  using Index = std::vector<size_t>;
  enum class BarycentricType { WACHSPRESS, MEAN_VALUE, HARMONIC };

  SPatch(std::string filename);
  SPatch(size_t num_sides, size_t depth);
  virtual ~SPatch();

  virtual void initDomain() override;
  virtual void initDomainMesh(size_t resolution = 30) override;
  virtual Vector evaluateAtParam(double u, double v) const override;
  virtual bool reload() override;

  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;

  virtual void setControlPoint(const Index &index, const Vector& p);

  double getGBC(double u, double v, size_t i, BarycentricType type = BarycentricType::WACHSPRESS) const;

  static size_t binomial(size_t n, size_t k);
  static size_t multinomial(const Index &index);
  static double multiBernstein(const Index &index, const DoubleVector &bc);
  static std::vector<Index> indices(size_t n, size_t d);
  
  std::vector<Index> neighbors(const Index& si) const;
  
  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }
  
  size_t n_;
  size_t d_;
  std::map<Index, Vector> net_;
};