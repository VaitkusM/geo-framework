#include "gbpatch.hh"

GBPatch::GBPatch(std::string filename) : NPatch(filename)
{
  reload();
}

GBPatch::GBPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.sp", num_sides), d_(depth), normalized(false)
{
  initDomain();

  net_.clear();
  Vector center = Vector(0.0, 0.0, 0.0);
  net_[{n_, 0, 0}] = center; 
  for (size_t i = 0; i < n_; ++i) {
    auto cim2 = vertices_[prev(prev(i))];
    auto cim1 = vertices_[prev(i)];
    auto ci = vertices_[i];
    auto cip1 = vertices_[next(i)];

    auto ql00 = cim1;
    auto ql01 = 0.5 * (cim1 + ci);
    auto ql10 = 0.5 * (cim2 + cim1);
    auto ql11 = center;
    auto qr00 = ci;
    auto qr01 = 0.5 * (cim1 + ci);
    auto qr10 = 0.5 * (cip1 + ci);
    auto qr11 = center;

    auto cp00 = ql00;
    auto cp01 = 1.0 / 3.0 * ql00 + 2.0 / 3.0 * ql01;
    auto cp02 = 1.0 / 3.0 * qr00 + 2.0 / 3.0 * qr01;
    auto cp03 = qr00;

    auto cp10 = 1.0 / 3.0 * ql00 + 2.0 / 3.0 * ql10;
    auto cp11 = 1. / 9. * ql00 + 2. / 9. * ql01 + 2. / 9. * ql10 + 4. / 9. * ql11;
    auto cp12 = 1. / 9. * qr00 + 2. / 9. * qr01 + 2. / 9. * qr10 + 4. / 9. * qr11;
    auto cp13 = 1.0 / 3.0 * qr00 + 2.0 / 3.0 * qr10;

    net_[{i, 0, 0}] = cp00;
    net_[{i, 0, 1}] = cp01;
    net_[{i, 1, 0}] = cp10;
    net_[{i, 1, 1}] = cp11;
    net_[{next(i), 0, 0}] = cp03;
    net_[{next(i), 0, 1}] = cp13;
    net_[{next(i), 1, 0}] = cp02;
    net_[{next(i), 1, 1}] = cp12;
  }

  footpoints_ = net_;

  initBasicShape();
  initBlendFunctions();

  // std::cerr << "Number of CPs:" << net_.size() << std::endl;

  updateBaseMesh();
}

GBPatch::~GBPatch() {

}

bool
GBPatch::reload()
{
  initDomain();
  updateBaseMesh();
  return true;
}

void
GBPatch::initBasicShape()
{
  for (auto& [id, p] : net_) {
    if (basic_shape_type == BasicShapeType::PARABOLOID) {
      p[2] = 1.0 - p[0] * p[0] - p[1] * p[1];
    }
    else if (basic_shape_type == BasicShapeType::HYPERBOLOID) {
      p[2] = 1.0 - p[0] * p[0] + p[1] * p[1];
    }
    else {
      p[2] = 0.0;
    }
  }
}

void
GBPatch::initDomain()
{
  double alpha = 2.0 * M_PI / n_;
  vertices_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    vertices_[i] = Vector(std::cos(alpha * i), std::sin(alpha * i), 0.0);
  }

  initDomainMesh();
}

void
GBPatch::initDomainMesh(size_t resolution)
{
  generateSpiderMesh(resolution, domain_mesh);
}

Vector
GBPatch::evaluateAtParam(const BaseMesh::VertexHandle& vtx) const
{
  const auto &bf = blend_functions_[vtx.idx()];
  auto uv = domain_mesh.point(vtx);
  auto u = uv[0], v = uv[1];
  size_t nl = num_layers();
  Vector p(0, 0, 0);
  double sum = 0.0;
  size_t idx = 0;
  for (const auto& [id, cp] : net_) {

    auto [i, l, c] = id;

    if (i < n_) {
      double Blc = bf[i][l][c] + bf[prev(i)][c][d_ - l];
      if (!show_basis_fcn) {
        p += Blc * cp;
      }
      else {
        if (idx == (selected_idx % net_.size())) {
          p = Vector(u, v, Blc);
        }
      }
      sum += Blc;
      ++idx;
    }
  }

  if (!normalized) {
    if (!show_basis_fcn) {
      auto ccp = net_.at({ n_,0,0 });
      p += (1.0 - sum) * ccp;
    }
    else {
      if (idx == (selected_idx % net_.size())) {
        p = Vector(u, v, 1.0 - sum);
        ++idx;
      }
    }
  }
  else {
    if (!show_basis_fcn) {
      p /= sum;
    }
    else {
      p[2] /= sum;
    }
  }
  return p;
}

Vector
GBPatch::evaluateAtParam(double u, double v) const
{
  std::vector<std::vector<DoubleVector> > bf;
  getBlendFunctions(u ,v, bf);
  size_t nl = num_layers();
  Vector p(0, 0, 0);
  double sum = 0.0;
  size_t idx = 0;
  for (const auto& [id, cp] : net_) {

    auto [i, l, c] = id;

    if (i < n_) {
      double Blc = bf[i][l][c] + bf[prev(i)][c][d_ - l];
      if (!show_basis_fcn) {
        p += Blc * cp;
      }
      else {
        if (idx == (selected_idx % net_.size())) {
          p = Vector(u, v, Blc);
        }
      }
      sum += Blc;
      ++idx;
    }
  }

  if (!normalized) {
    if (!show_basis_fcn) {
      auto ccp = net_.at({ n_,0,0 });
      p += (1.0 - sum) * ccp;
    }
    else {
      if (idx == (selected_idx % net_.size())) {
        p = Vector(u, v, 1.0 - sum);
        ++idx;
      }
    }
  }
  else {
    if (!show_basis_fcn) {
      p /= sum;
    }
    else {
      p[2] /= sum;
    }
  }
  return p;
}

void
GBPatch::initBlendFunctions()
{
  blend_functions_.clear();
  blend_functions_.resize(domain_mesh.n_vertices());
  for (auto vtx : domain_mesh.vertices()) {
    auto pt = domain_mesh.point(vtx);
    getBlendFunctions(pt[0], pt[1], blend_functions_[vtx.idx()]);
  }
}

void
GBPatch::getBlendFunctions(double u, double v, std::vector<std::vector<DoubleVector> >& values) const
{
  values.clear();
  values.resize(n_, std::vector<DoubleVector>(d_+1, DoubleVector(d_+1)));
  std::vector<double> bc(n_), s(n_, 0.), h(n_, 0.);
  for (size_t i = 0; i < n_; ++i) {
    bc[i] = getGBC(u, v, i);
  }
  //std::rotate(bc.begin(),bc.begin() + 1, bc.end());
  for (size_t i = 0; i < n_; ++i) {
    if (bc[i] + bc[prev(i)] > epsilon) {
      s[i] = bc[i] / (bc[i] + bc[prev(i)]);
    }
    h[i] = 1.0 - bc[i] - bc[prev(i)];
  }

  size_t nl = num_layers();
  for (size_t i = 0; i < n_; ++i) {
    auto a = 1.0, b = 1.0;
    if ((std::pow(h[prev(i)], nl) + std::pow(h[i], nl)) < epsilon) {
      a = 0.5;
    }
    else {
      a = std::pow(h[prev(i)], nl) / (std::pow(h[prev(i)], nl) + std::pow(h[i], nl));
    }

    if ((std::pow(h[next(i)], nl) + std::pow(h[i], nl)) < epsilon) {
      b = 0.5;
    }
    else {
      b = std::pow(h[next(i)], nl) / (std::pow(h[next(i)], nl) + std::pow(h[i], nl));
    }
    DoubleVector Bu, Bv;
    bernstein(d_, s[i], Bu);
    bernstein(d_, h[i], Bv);
    for (size_t l = 0; l < nl; ++l) {
      for (size_t c = 0; c <= d_; ++c) {
        double mu = c < nl ? a : b;
        values[i][l][c] = mu * Bu[c] * Bv[l];
      }
    }
  }
}

void
GBPatch::bernstein(size_t n, double u, DoubleVector& coeff)
{
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

void
GBPatch::swapFootpoints()
{
  std::swap(footpoints_, net_);
}

std::vector<GBPatch::Index>
GBPatch::neighbors(const Index& si) const
{
  auto [i,l,c] = si;
  auto nl = num_layers();
  std::vector<Index> candidates;
  if(i < n_) {
    if(l < nl - 1) {
      candidates.push_back({i, l + 1 , c});
    }
    else {
      candidates.push_back({prev(i), c, nl - 1});
    }
    if(c < nl - 1) {
      candidates.push_back({ i, l , c + 1  });
    }
    else {
      candidates.push_back({ next(i), nl - 1 , l });
    }
  }

  else {
    for (size_t i = 0; i < n_; ++i) {
      candidates.push_back({ i, nl - 1, nl - 1 });
    }
  }

  return candidates;
}


void
GBPatch::draw(const Visualization& vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    for (const auto& [idx,cp] : net_) {
      auto nbs = neighbors(idx);
      for (const auto& nb : nbs) {
        glBegin(GL_LINE_STRIP);
        glVertex3dv(cp.data());
        glVertex3dv(net_.at(nb).data());
        glEnd();
      }
    }
    glLineWidth(1.0);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    size_t idx = 0;
    for(const auto& [id,cp] : net_) {
      if (idx == (selected_idx % net_.size())) {
        glPointSize(16.0);
      }
      else {
        glPointSize(8.0);
      }
      glBegin(GL_POINTS);
      glVertex3dv(cp.data());
      glEnd();
      ++idx;
    }
    glPointSize(1.0);
    glEnable(GL_LIGHTING);
  }
}

void
GBPatch::drawWithNames(const Visualization& vis) const
{
  if (!vis.show_control_points) {
    return;
  }
  size_t i = 0;
  for (const auto& p : net_) {
    glPushName(static_cast<GLuint>(i));
    glRasterPos3dv(p.second.data());
    glPopName();
    ++i;
  }
}

Vector
GBPatch::postSelection(int selected)
{
  size_t i = 0;
  for (const auto& p : net_) {
    if (i == selected) {
      return p.second;
    }
    ++i;
  }
  return Vector(0, 0, 0);
}

void
GBPatch::movement(int selected, const Vector& pos)
{
  size_t i = 0;
  for (auto& p : net_) {
    if (i == selected) {
      p.second = pos;
      break;
    }
    ++i;
  }
}