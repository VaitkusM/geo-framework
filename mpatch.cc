#include "mpatch.hh"

#include <fstream>

MPatch::MPatch(std::string filename) : NPatch(filename)
{
  reload();
}

MPatch::MPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.sp", num_sides), d_(depth)
{
  initDomain();

  net_.clear();
  Vector center = Vector(0.0, 0.0, 0.0);
  net_[{0, 0, d_, 0}] = center;
  for (size_t j = 1; j <= d_; ++j) {
    double u = (double)j / (double)d_;
    for (size_t k = 0; k < n_; ++k) {
      for (size_t i = 0; i < j; ++i) {
        double v = (double)i / (double)j;
        Vector ep = vertices_[prev(k)] * (1.0 - v) + vertices_[k] * v;
        Vector p = center * (1.0 - u) + ep * u;
        net_[{j-i, i, d_-j, k}] = p;
      }
    }
  }

  footpoints_ = net_;

  initBasicShape();

  std::cerr << "Number of CPs:" << net_.size() << std::endl;

  updateBaseMesh();
}

MPatch::~MPatch() {

}

bool
MPatch::reload()
{
  initDomain();
  updateBaseMesh();
  return true;
}

void
MPatch::initBasicShape()
{
  for (auto& [id, p] : net_) {
    if (basic_shape_type == BasicShapeType::PARABOLOID) {
      p[2] = 1.0 - p[0] * p[0] - p[1] * p[1];
    }
    else if (basic_shape_type == BasicShapeType::HYPERBOLOID) {
      p[2] = 1.0 - p[0] * p[0] + p[1] * p[1];
    }
  }
}

void
MPatch::initDomain()
{
  double alpha = 2.0 * M_PI / n_;
  vertices_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    vertices_[i] = Vector(std::cos(alpha * i), std::sin(alpha * i), 0.0);
  }

  initDomainMesh();
}

void
MPatch::initDomainMesh(size_t resolution)
{
  generateSpiderMesh(resolution, domain_mesh);
}

Vector
MPatch::evaluateAtParam(double u, double v) const
{
  // return {0,0,0};

  Vector p(0, 0, 0);
  size_t cp_idx = 0;
  DoubleVector bf;
  getBlendFunctions(u, v, bf, true);
  for (const auto& cp : net_) {
    const auto& id = cp.first;
    if (!show_basis_fcn) {
      p += cp.second * bf[cp_idx];
    }
    else {
      if (cp_idx == (selected_idx % net_.size())) {
        p += Vector(u, v, bf[cp_idx]);
      }
    }
    ++cp_idx;
  }
  return p;
}

void 
MPatch::getBlendFunctions(double u, double v, DoubleVector &values, bool normalize) const
{
  values.clear();
  values.reserve(net_.size());

  for(const auto& cp : net_) {
    values.push_back(getBlendFunctionUnnormalized(u, v, cp.first));
  }

  if(normalize){
    normalizeValues(values);
  }
}

double 
MPatch::getBlendFunctionUnnormalized(double u, double v, const Index &i) const
{
  DoubleVector ls(n_);
  const double alpha = 2.0 * M_PI / double(n_);
  const double cs    = std::cos(alpha);
  const double cs2   = std::cos(alpha/2.0);
  const double C     = std::pow(cs2, 2) / std::pow(cs, 2) - std::pow(u, 2) - std::pow(v, 2);
  for (size_t i = 0; i < n_; ++i) {
    ls[i] = -u * std::cos(alpha * i + alpha / 2.0) - v * std::sin(alpha * i + alpha / 2.0) + cs2;
    //ls[i] /= (1 + std::cos(alpha/2.0));
  }
  DoubleVector bc(n_);
  double pr = 1.0;
  for (size_t i = 0; i < n_; ++i) {
    bc[i] = C;
    for (size_t j = 0; j < n_; ++j) {
      if (j != i && j != prev(i)) {
        bc[i] *= ls[j];
      }
    }
    pr *= ls[i];
  }

  normalizeValues(bc);

  double Bi = multiBernstein(
    { i[0], i[1], i[2] },
    { bc[prev(i[3])], bc[i[3]], pr }
  );
  return Bi;
}

size_t
MPatch::binomial(size_t n, size_t k)
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

size_t
MPatch::multinomial(const Index& index)
{
  auto fact = [](size_t n) { return (size_t)std::lround(std::tgamma(n + 1)); };
  size_t numerator = 0, denominator = 1;
  for (auto i : index) {
    numerator += i;
    denominator *= fact(i);
  }
  return fact(numerator) / denominator;
}

double
MPatch::multiBernstein(const Index& index, const DoubleVector& bc)
{
  auto result = static_cast<double>(multinomial(index));
  for (size_t i = 0; i < index.size(); ++i)
    result *= std::pow(bc[i], index[i]);
  return result;
}


void
MPatch::setControlPoint(const Index& i, const Vector& p)
{
    net_[i] = p;
}

std::vector<MPatch::Index> 
MPatch::neighbors(const Index& si) const
{
  auto step = [&](const Index& index, size_t i, int dir, size_t n) {
    size_t j = (i + dir + n) % n;
    Index result(index);
    if (result[i] > 0) {
      result[i]--;
      result[j]++;
    }
    if (result[0] == 0) {
      std::swap(result[0], result[1]);
      result[3] = next(result[3]); 
    }
    return result;
  };

  std::vector<Index> candidates;
  for (size_t i = 0; i < 3; i++) {
    auto sip1 = step(si, i, 1, 3);
    if(net_.find(sip1) != net_.end()) {
      candidates.push_back(sip1);
    }
    if (si[i] > 0) { // add a check here as well to prevent wrap-around
      auto sim1 = step(si, i, -1, 3);
      if (net_.find(sim1) != net_.end()) {
        candidates.push_back(sim1);
      }
    }
  }

  if(si[2] == d_) {
    for(size_t i = 0; i < n_; ++i) {
      candidates.push_back({1, 0, d_-1, i});
    }
  }

  return candidates;
}

void 
MPatch::swapFootpoints()
{
  std::swap(footpoints_, net_);
}

void 
MPatch::draw(const Visualization& vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    for (const auto& cp : net_) {
      auto nbs = neighbors(cp.first);
      for (const auto& nb : nbs) {
        glBegin(GL_LINE_STRIP);
        glVertex3dv(cp.second.data());
        glVertex3dv(net_.at(nb).data());
        glEnd();
      }
    }
    glLineWidth(1.0);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    size_t idx = 0;
    for (const auto& cp : net_) {
      if (idx++ == (selected_idx % net_.size())) {
        glPointSize(16.0);
      }
      else {
        glPointSize(8.0);
      }
      glBegin(GL_POINTS);
      glVertex3dv(cp.second.data());
      glEnd();
    }
    glPointSize(1.0);
    glEnable(GL_LIGHTING);
  }
}

void
MPatch::drawWithNames(const Visualization& vis) const
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
MPatch::postSelection(int selected)
{
  size_t i = 0;
  for (const auto& p : net_) {
    if(i == selected) {
      return p.second;
    }
    ++i;
  }
  return Vector(0,0,0);
}

void 
MPatch::movement(int selected, const Vector& pos)
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