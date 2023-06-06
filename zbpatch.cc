#include "zbpatch.hh"
#include "powell.hh"
#include <OpenMesh/Core/Utils/PropertyManager.hh>

ZBPatch::ZBPatch(std::string filename) : NPatch(filename)
{
  reload();
}

ZBPatch::ZBPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.zb", num_sides), d_(depth), normalized(true)
{
  initDomain();

  net_.clear();
  Vector center = Vector(0.0, 0.0, 0.0);
  //net_[{n_, 0, 0}] = center;
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

    net_[gb2zb(i, 0, 0)] = cp00;
    net_[gb2zb(i, 0, 1)] = cp01;
    net_[gb2zb(i, 1, 0)] = cp10;
    net_[gb2zb(i, 1, 1)] = cp11;
  }

  footpoints_ = net_;

  initBasicShape();

  // std::cerr << "Number of CPs:" << net_.size() << std::endl;

  updateBaseMesh();
}

ZBPatch::~ZBPatch() {

}

bool
ZBPatch::reload()
{
  initDomain();
  updateBaseMesh();
  return true;
}

void
ZBPatch::initBasicShape()
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
ZBPatch::initDomain()
{
  double alpha = 2.0 * M_PI / n_;
  vertices_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    vertices_[i] = Vector(std::cos(alpha * i), std::sin(alpha * i), 0.0);
  }

  initDomainMesh();
}

void
ZBPatch::initDomainMesh(size_t resolution)
{
  generateSpiderMesh(resolution, domain_mesh);
}

Vector
ZBPatch::evaluateAtParam(const BaseMesh::VertexHandle& vtx) const
{
  //return domain_mesh.point(vtx);
  std::map<Index, double> bf;
  getBlendFunctions(vtx, bf);
  double u = domain_mesh.point(vtx)[0];
  double v = domain_mesh.point(vtx)[1];

  size_t nl = num_layers();
  Vector p(0, 0, 0);
  double sum = 0.0;
  size_t idx = 0;
  for (const auto& [id, cp] : net_) {

    if (true) {
      double Blc = bf.at(id);
      if (!show_basis_fcn) {
        p += Blc * cp;
      }
      else {
        if (idx == (selected_idx % net_.size())) {
          //p = Vector(u, v, mesh.data(vtx).gbc[selected_idx % n_]);
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
    p /= sum;
  }
  return p;
}

Vector 
ZBPatch::evaluateAtParam(double u, double v) const
{
  return Vector(0,0,0);
}

void
ZBPatch::getBlendFunctions(const BaseMesh::VertexHandle& vtx, std::map<Index, double>& values) const
{
  values.clear();
  Parameter params;
  getParameters(vtx, params);
  mesh.data(vtx).gbc = params;

  double B, Bsum = 0;
  for (const auto& [l, p] : net_) {
    auto it = std::find(l.begin(), l.end(), 0);
    if (it == l.end()) {
      // Inside - (4.5b)
      size_t i = std::min_element(l.begin(), l.end()) - l.begin();
      size_t im = prev(i), ip = next(i);
      if (l[im] > l[ip]) {
        im = i;
        i = ip;
      }
      B =  static_cast<double>(binomial(d_, l[im])); 
      B *= static_cast<double>(binomial(d_, l[i]));
      for (size_t j = 0; j < n_; ++j) {
        B *= std::pow(params[j], l[j]);
      }
    }
    else {
      // On side i - (4.5a)
      size_t i = it - l.begin();
      size_t im = prev(i), ip = next(i);
      size_t k = l[im];
      B = binomial(d_, k) * std::pow(params[im], k) * std::pow(params[ip], d_ - k);
      for (size_t j = 0; j < n_; ++j) {
        if (j != im && j != i && j != ip) {
          B *= std::pow(params[j], d_);
        }
      }
      if (n_ == 3) {
        // Special case - (3.2b)
        // Note: the equation there treats side (i-1), and k and m-k are swapped
        double f = d_ - k <= k ? d_ - k + (2 * k - d_) * params[ip] : k + (d_ - 2 * k) * params[im];
        B *= 1 - params[i] * f;
      }
      else {
        double prod = 1;
        for (size_t j = 0; j < n_; ++j) {
          prod *= params[j];
        }
        B *= (1 - d_ * coefficient() * prod);
      }
    }
    values[l] = B;
    Bsum += B;
  }
  // Correction (2.6)
  double S = Bsum - 1;
  size_t T = numControlPoints(n_, d_, true);
  for (const auto& [l, p] : net_) {
    if (std::find(l.begin(), l.end(), 0) == l.end()) {
      values[l] -= (S / T);
    }
  }
}

void 
ZBPatch::getParameters(const BaseMesh::VertexHandle& vtx, Parameter& params) const
{
  // auto p = domain_mesh.point(vtx);
  // auto gbc = getGBCs(p[0], p[1]);
  // Parameter h(n_);
  // for(size_t i = 0; i < n_; ++i) {
  //   h[i] = 1.0 - gbc[i] - gbc[prev(i)];
  // }
  // params = h;//project(n_, h);
  // return;
  std::vector<DoubleVector> vertices_nd;
  vertices_nd.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    DoubleVector v(n_, 1);
    v[prev(i)] = 0;
    v[i] = 0;
    vertices_nd.push_back(v);
  }
  DoubleVector center(n_, n_ == 5 ? (std::sqrt(5) - 1) / 2 : (1.0 / std::sqrt(2)));
  //params.clear();
  //params.reserve(n_);

  size_t resolution = 60;
  const auto [i,j,k] = domain_mesh.data(vtx).spider_idx;
  double u = (double)j / (double)resolution;
  double v = j > 0 ? ((double)i / (double)j) : 0;

  auto ep = affine(vertices_nd[prev(k)], v, vertices_nd[k]);
  auto pt = affine(center, u, ep);
  //params = pt;
  params = project(n_, pt);
}

size_t
ZBPatch::binomial(size_t n, size_t k)
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

void
ZBPatch::bernstein(size_t n, double u, DoubleVector& coeff)
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
ZBPatch::swapFootpoints()
{
  std::swap(footpoints_, net_);
}

std::vector<ZBPatch::Index>
ZBPatch::neighbors(const Index& si) const
{
  // auto [i, l, c] = si;
  // auto nl = num_layers();
  std::vector<Index> candidates;
  // if (i < n_) {
  //   if (l < nl - 1) {
  //     candidates.push_back({ i, l + 1 , c });
  //   }
  //   else {
  //     candidates.push_back({ prev(i), c, nl - 1 });
  //   }
  //   if (c < nl - 1) {
  //     candidates.push_back({ i, l , c + 1 });
  //   }
  //   else {
  //     candidates.push_back({ next(i), nl - 1 , l });
  //   }
  // }

  // else {
  //   for (size_t i = 0; i < n_; ++i) {
  //     candidates.push_back({ i, nl - 1, nl - 1 });
  //   }
  // }

  //return candidates;

  return candidates;
}

ZBPatch::DoubleVector 
ZBPatch::affine(const DoubleVector& a, double x, const DoubleVector& b) {
  size_t n = a.size();
  DoubleVector result(n);
  for (size_t i = 0; i < n; ++i)
    result[i] = a[i] * (1 - x) + b[i] * x;
  return result;
}

double 
ZBPatch::coefficient() const {
  switch (n_) {
  case 3: return 2.0;
  case 5: return 1.0;
  case 6: return 2.0;
  default: assert(false && "coefficient: only implemented for n=3,5,6");
  }
  return 0;
}

ZBPatch::Parameter
ZBPatch::project(size_t n, const Parameter& p)
{
  auto f5 = [&](const Powell::Point& x) {
    for (size_t i = 0; i < 5; ++i)
      if (x[i] < 0)
        return std::numeric_limits<double>::max();
    double err = 0;
    for (size_t i = 0; i < 5; ++i) // 3 should be enough
      err += std::pow(1 - x[i] - x[(i + 2) % n] * x[(i + 3) % n], 2);
    return err;
  };
  auto f6 = [&](const Powell::Point& x) {
    for (size_t i = 0; i < 6; ++i)
      if (x[i] < 0)
        return std::numeric_limits<double>::max();
    double err = 0;
    for (size_t i = 0; i < 6; ++i) { // 4 should be enough
      double xm = x[(i + n - 1) % n], xi = x[i], xp = x[(i + 1) % n];
      err += std::pow(xp * xp * (1 - xm * xi) * (1 - 2 * xm * xi) +
        xp * (2 * xm - 3 * xm * xm * xi + xm * xi * xi) +
        xm * xm - 1, 2);
    }
    return err;
  };
  Powell::Point x = p;
  for (size_t i = 0; i < 15; ++i) // increased iteration #
    if (n == 5)
      Powell::optimize(f5, x, 50, 0.2, 20, 1e-18); // increased tolerance
    else
      Powell::optimize(f6, x, 50, 0.2, 20, 1e-18);
  //std::cerr << "Error: " << f5(x) << std::endl; 
  return x;
}

ZBPatch::Index
ZBPatch::gb2zb(size_t s, size_t l, size_t c) const
{
  Index idx(n_, d_);
  size_t i = c <= d_/2 ? s : next(s);
  size_t im1 = prev(i);
  size_t im2 = prev(im1);
  size_t ip1 = next(i);
  idx[im1] = c;
  idx[i] = l;
  idx[im2] = d_ - idx[i];
  idx[ip1] = d_ - idx[im1];
  auto dj = d_ - std::min(idx[im1], idx[i]);
  for(size_t j = 0; j < n_; ++j) {
    if(j != im2 && j != im1 && j != i && j!= ip1){
      idx[j] = dj;
    }
  }
  return idx;
}

void
ZBPatch::draw(const Visualization& vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    for (const auto& [idx, cp] : net_) {
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
    for (const auto& [id, cp] : net_) {
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
ZBPatch::drawWithNames(const Visualization& vis) const
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
ZBPatch::postSelection(int selected)
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
ZBPatch::movement(int selected, const Vector& pos)
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