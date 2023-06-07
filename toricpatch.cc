#include "toricpatch.hh"

ToricPatch::ToricPatch(std::string filename) : NPatch(filename)
{
  reload();
}

ToricPatch::ToricPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.trp", num_sides), d_(depth), normalized(false)
{
  initDomain();

  net_.clear();
  double large = std::numeric_limits<double>::max();
  Vector box_min(large, large, large), box_max(-large, -large, -large);
  for (auto v : vertices_) {
    box_min.minimize(v);
    box_max.maximize(v);
  }
  tdu_ = static_cast<size_t>(std::round(box_max[0] - box_min[0]));
  tdv_ = static_cast<size_t>(std::round(box_max[1] - box_min[1]));


  for (size_t i = 0; i <= tdu_; ++i) {
    for(size_t j = 0; j <= tdv_; ++j) {
      bool inside = true;
      for(size_t k = 0; k < n_; ++k) {
        if (sideDistance(static_cast<double>(i), static_cast<double>(j), k) < 0) {
          inside = false;
        }
      }
      if(inside) {
        net_[{i,j,0}] = {i,j,0};
      }
    }
  }


  footpoints_ = net_;

  initBasicShape();
  
  ld_.clear();
  for(auto& [id,pt] : net_) {
    std::vector<size_t> d(n_);
    for(size_t i = 0; i < n_; ++i) {
      d[i] = static_cast<size_t>(std::round(sideDistance(std::round(id[0]), std::round(id[1]), i)));
    }
    ld_[id] = d;
  }



  // std::cerr << "Number of CPs:" << net_.size() << std::endl;

  updateBaseMesh();
}

ToricPatch::~ToricPatch() {

}

bool
ToricPatch::reload()
{
  initDomain();
  updateBaseMesh();
  return true;
}

void
ToricPatch::initBasicShape()
{
  auto c = getDomainCenter();

  for (auto& [id, p] : net_) {
    p = (Vector(id[0], id[1], 0) - c) / d_;
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
ToricPatch::initDomain()
{
  assert(n_ == 5);
  vertices_.clear();
  vertices_.resize(n_);
  vertices_[0] = {0,0,0};
  vertices_[1] = {d_,0,0};
  vertices_[2] = {d_ + d_,d_,0};
  vertices_[3] = {d_, d_ + d_, 0};
  vertices_[4] = {0, d_, 0};

  side_eqs_.clear();
  side_eqs_.reserve(n_);
  for(size_t i = 0; i < n_; ++i) {
    auto v0 = vertices_[i];
    auto v1 = vertices_[next(i)];
    auto d = v1 - v0;
    auto r = Vector(-d[1], d[0], 0);
    auto g = std::gcd(static_cast<size_t>(std::round(std::abs(r[0]))), static_cast<size_t>(std::round(std::abs(r[1]))));
    r /= g;
    side_eqs_.push_back({r[0],r[1],-(r[0]*v0[0] + r[1]*v0[1])});
    std::cerr << "side_eq " << i << " : (" << side_eqs_[i][0] << ", " << side_eqs_[i][1] << ", " << side_eqs_[i][2] << ")" << std::endl;
  }

  initDomainMesh();
}

void
ToricPatch::initDomainMesh(size_t resolution)
{
  generateSpiderMesh(resolution, domain_mesh);
}

Vector
ToricPatch::evaluateAtParam(const BaseMesh::VertexHandle& vtx) const
{
  //return Vector();
  auto p = domain_mesh.point(vtx);
  return evaluateAtParam(p[0], p[1]);
}

Vector
ToricPatch::evaluateAtParam(double u, double v) const
{
  std::map<Index, double> bf;
  getBlendFunctions(u, v, bf);

  size_t nl = num_layers();
  Vector p(0, 0, 0);
  //double sum = 0.0;
  size_t idx = 0;
  for (const auto& [id, cp] : net_) {
    double Bi = bf[id];
    if (!show_basis_fcn) {
      p += Bi * cp;
    }
    else {
      if (idx == (selected_idx % net_.size())) {
        p = Vector(u, v, Bi);
      }
    }
    //sum += Bi;
    ++idx;
  }

  //p /= sum;

  return p;
}

void
ToricPatch::getBlendFunctions(double u, double v, std::map<Index, double>& values) const
{
  values.clear();
  double sum = 0.0;
  for (const auto& [id, cp] : net_) {
    double B = 1.0;
    B *= coefficient(id);
    auto ld = ld_.at(id);
    for (size_t i = 0; i < n_; ++i) {
      auto luv = sideDistance(u, v, i);
      auto lab = ld[i];
      B *= std::pow(luv, lab);
    }
    values[id] = B;
    sum += B;
  }
  for (const auto& [id, cp] : net_) {
    values[id] /= sum;
  }
}

double
ToricPatch::sideDistance(double u, double v, size_t i) const
{
  auto [a,b,c] = side_eqs_[i];
  return a*u + b*v + c;
}

size_t
ToricPatch::binomial(size_t n, size_t k)
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

size_t 
ToricPatch::coefficient(const Index& idx) const
{
  size_t i;
  bool inside = true;
  for(i = 0; i < n_; i++) {
    if(sideDistance(static_cast<double>(idx[0]), static_cast<double>(idx[1]), i) == 0) {
      inside = false;
      break;
    }
  }
  if(!inside) {
    auto v0 = vertices_[i];
    auto v1 = vertices_[next(i)];
    auto dv = (v1- v0)/d_;
    auto li = (Vector(idx[0],idx[1],idx[2]) - v0).norm()/dv.norm(); 
    if(std::round(li) > d_) {
      std::cerr << "Baj van! l" << i << " = " << li << std::endl;
      std::cerr << "v0: " << v0[0] << ", " << v0[1] << ", " << v0[2] << std::endl;
      std::cerr << "v1: " << v1[0] << ", " << v1[1] << ", " << v1[2] << std::endl;
      std::cerr << "cp: " << idx[0] << ", " << idx[1] <<  std::endl;
    }
    return binomial(d_, size_t(std::round(li)));
  }
  else {
    return 1;
  }
}

void
ToricPatch::bernstein(size_t n, double u, DoubleVector& coeff)
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
ToricPatch::swapFootpoints()
{
  std::swap(footpoints_, net_);
}

std::vector<ToricPatch::Index>
ToricPatch::neighbors(const Index& si) const
{
  auto [i, l, c] = si;
  auto nl = num_layers();
  std::vector<Index> candidates;
  if (i < n_) {
    if (l < nl - 1) {
      candidates.push_back({ i, l + 1 , c });
    }
    else {
      candidates.push_back({ prev(i), c, nl - 1 });
    }
    if (c < nl - 1) {
      candidates.push_back({ i, l , c + 1 });
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
ToricPatch::draw(const Visualization& vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    // glLineWidth(3.0);
    // glColor3d(0.3, 0.3, 1.0);
    // for (const auto& [idx, cp] : net_) {
    //   auto nbs = neighbors(idx);
    //   for (const auto& nb : nbs) {
    //     glBegin(GL_LINE_STRIP);
    //     glVertex3dv(cp.data());
    //     glVertex3dv(net_.at(nb).data());
    //     glEnd();
    //   }
    // }
    // glLineWidth(1.0);
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
ToricPatch::drawWithNames(const Visualization& vis) const
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
ToricPatch::postSelection(int selected)
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
ToricPatch::movement(int selected, const Vector& pos)
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