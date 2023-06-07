#pragma once

#include "npatch.hh"

class ToricPatch : public NPatch {
public:
  using Index = std::array<size_t, 3>;
  using Ribbon = std::vector<std::vector<Vector> >;

  ToricPatch(std::string filename);
  ToricPatch(size_t num_sides, size_t num_layers);
  virtual ~ToricPatch();

  virtual void initBasicShape() override;
  virtual void initDomain() override;
  virtual void initDomainMesh(size_t resolution = 60) override;
  virtual Vector evaluateAtParam(double u, double v) const override;
  virtual Vector evaluateAtParam(const BaseMesh::VertexHandle& vtx) const override;
  virtual bool reload() override;

  virtual void draw(const Visualization& vis) const override;
  virtual void drawWithNames(const Visualization& vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector& pos) override;

  static size_t binomial(size_t n, size_t k);
  size_t coefficient(const Index & idx) const;
  static void bernstein(size_t n, double u, DoubleVector& coeff);
  double sideDistance(double u, double v, size_t i) const;
  void getBlendFunctions(double u, double v, std::map<Index, double>& values) const;

  virtual void swapFootpoints() override;

  size_t num_layers() const { return (d_ + 1) / 2; };
  //virtual void setControlPoint(const Index& index, const Vector& p);

  std::vector<Index> neighbors(const Index& si) const;

  size_t d_, tdu_, tdv_;
  std::map<Index, Vector> footpoints_;
  std::map<Index, Vector> net_;
  std::map<Index, std::vector<size_t>> ld_;
  std::vector<std::array<double, 3> > side_eqs_;
  std::vector<std::map<Index, double>> blend_functions;
  bool normalized;
};