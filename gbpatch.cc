#include "gbpatch.hh"

GBPatch::GBPatch(std::string filename) : NPatch(filename)
{
  reload();
}

GBPatch::GBPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.sp", num_sides), d_(depth), normalized(false)
{
  initDomain();

  ribbons_.clear();
  ribbons_.reserve(n_);
  Vector center = Vector(0.0, 0.0, 0.0);
  central_cp = center;
  for (size_t i = 0; i < n_; ++i) {
    auto cim2 = vertices_[prev(prev(i))];
    auto cim1 = vertices_[prev(i)];
    auto ci = vertices_[i];
    auto cip1 = vertices_[next(i)];

    auto ql00 = cim1;
    auto ql01 = 0.5 * (cim1 + ci);
    auto ql10 = 0.5 * (cim2 + cim1);
    auto ql11 = central_cp;
    auto qr00 = ci;
    auto qr01 = 0.5 * (cim1 + ci);
    auto qr10 = 0.5 * (cip1 + ci);
    auto qr11 = central_cp;

    auto cp00 = ql00;
    auto cp01 = 1.0 / 3.0 * ql00 + 2.0 / 3.0 * ql01;
    auto cp02 = 1.0 / 3.0 * qr00 + 2.0 / 3.0 * qr01;
    auto cp03 = qr00;

    auto cp10 = 1.0 / 3.0 * ql00 + 2.0 / 3.0 * ql10;
    auto cp11 = 1. / 9. * ql00 + 2. / 9. * ql01 + 2. / 9. * ql10 + 4. / 9. * ql11;
    auto cp12 = 1. / 9. * qr00 + 2. / 9. * qr01 + 2. / 9. * qr10 + 4. / 9. * qr11;
    auto cp13 = 1.0 / 3.0 * qr00 + 2.0 / 3.0 * qr10;

    ribbons_.push_back({ {cp00,cp01,cp02,cp03},{cp10,cp11,cp12,cp13} });
  }

  for (auto& r : ribbons_) {
    for (auto& l : r) {
      for (auto& cp : l) {
        cp[2] = 1.0 - cp[0] * cp[0] - cp[1] * cp[1];
      }
    }
  }
  central_cp[2] = 1.0 - central_cp[0] * central_cp[0] - central_cp[1] * central_cp[1];

  // for (auto& p : net_) {
  //     p.second[2] = 1.0 - p.second[0] * p.second[0] + p.second[1] * p.second[1];
  // }

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
GBPatch::evaluateAtParam(double u, double v) const
{
  std::vector<double> bc(n_), s(n_, 0.), h(n_, 0.);
  for (size_t i = 0; i < n_; ++i) {
    bc[i] = getGBC(u, v, i);
  }
  //std::rotate(bc.begin(),bc.begin() + 1, bc.end());
  for(size_t i = 0; i < n_; ++i) {
    if(bc[i] + bc[prev(i)] > epsilon) {
      s[i] = bc[i]/(bc[i] + bc[prev(i)]);
    }
    h[i] = 1.0 - bc[i] - bc[prev(i)];
  }

  size_t nl = (d_+1)/2;
  Vector p(0, 0, 0);
  double sum = 0.0;
  for(size_t i = 0; i < n_; ++i) {
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
    DoubleVector Bu,Bv;
    bernstein(d_, s[i], Bu);
    bernstein(d_, h[i], Bv);
    for(size_t l = 0; l < nl; ++l) {
      for(size_t c = 0; c <= d_; ++c) {
        double mu  = c < nl ? a : b;
        double Blc = mu*Bu[c]*Bv[l];
        p   += Blc*ribbons_[i][l][c];
        sum += Blc;
      }
    }
  }
  if(!normalized) {
    p += (1.0 - sum)*central_cp;
  }
  else {
    p /= sum;
  }
  return p;
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
GBPatch::draw(const Visualization& vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    size_t idx = 0;
    for (const auto& r : ribbons_) {
      for(const auto& l : r) {
        for(const auto& cp : l) {
          glBegin(GL_POINTS);
          glVertex3dv(cp.data());
          glEnd();
        }
      }
    }
    glBegin(GL_POINTS);
    glVertex3dv(central_cp.data());
    glEnd();
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
}

Vector
GBPatch::postSelection(int selected)
{
  return Vector(0, 0, 0);
}

void
GBPatch::movement(int selected, const Vector& pos)
{

}