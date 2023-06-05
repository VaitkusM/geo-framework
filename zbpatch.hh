#pragma once

#include "npatch.hh"

class ZBPatch : public NPatch {
public:
  using Index = std::vector<size_t>;
  using Parameter = DoubleVector;
  using Ribbon = std::vector<std::vector<Vector> >;

  ZBPatch(std::string filename);
  ZBPatch(size_t num_sides, size_t num_layers);
  virtual ~ZBPatch();

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
  static void bernstein(size_t n, double u, DoubleVector& coeff);
  void getBlendFunctions(const BaseMesh::VertexHandle& vtx, std::map<Index, double>& values) const;

  virtual void swapFootpoints() override;

  size_t num_layers() const { return (d_ + 1) / 2; };
  static size_t numControlPoints(size_t n, size_t m, bool only_inside) {
    size_t T = m % 2 == 0 ? n * (m - 2) * m / 4 + 1 : n * (m - 1) * (m - 1) / 4;
    return only_inside ? T : T + n * m;
  }
  //virtual void setControlPoint(const Index& index, const Vector& p);

  std::vector<Index> neighbors(const Index& si) const;

  static DoubleVector affine(const DoubleVector& a, double x, const DoubleVector& b);
  double coefficient() const;
  static Parameter project(size_t n, const Parameter& p);
  void getParameters(const BaseMesh::VertexHandle& vtx, Parameter& params) const;

  size_t d_;
  std::map<Index, Vector> footpoints_;
  std::map<Index, Vector> net_;
  bool normalized;
};
