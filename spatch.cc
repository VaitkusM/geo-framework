#include "spatch.hh"

#include <fstream>

SPatch::SPatch(std::string filename) : NPatch(filename)
{
  reload();
}

SPatch::SPatch(size_t num_sides, size_t depth) : NPatch("NOTHING.sp", num_sides), d_(depth)
{
  initDomain();
  auto ids = indices(num_sides, depth);
  for(auto id : ids) {
    Vector pt(0.0, 0.0, 0.0);
    for(size_t i = 0; i < num_sides; ++i) {
      pt += (double(id[i])/double(depth))*vertices_[i];
    }
    net_[id] = pt;
  }
  updateBaseMesh();
}

SPatch::~SPatch() {

}

bool 
SPatch::reload()
{
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  f >> n_ >> d_;
  size_t n_cp = binomial(n_ + d_ - 1, d_);
  Index index(n_);
  for (size_t i = 0; i < n_cp; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      f >> index[j];
    }
    Vector p;
    for (size_t j = 0; j < 3; ++j) {
      f >> p[j];
    }
    setControlPoint(index, p);
  }
  initDomain();
  updateBaseMesh();
  return true;
}

void 
SPatch::initDomain()
{
  double alpha = 2.0 * M_PI / n_;
  vertices_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    vertices_[i] = Vector(std::cos(alpha * i), std::sin(alpha * i), 0.0);
  }

  initDomainMesh();
}

void 
SPatch::initDomainMesh(size_t resolution)
{
  generateSpiderMesh(resolution, domain_mesh);
}

double 
SPatch::getGBC(double u, double v, size_t i, BarycentricType type) const
{
  std::vector<Vector> vectors; vectors.reserve(n_);
  std::transform(vertices_.begin(), vertices_.end(),
                 std::back_inserter(vectors),
                 [u,v](const Vector &p)->Vector { return Vector(u,v,0) - p; });

  DoubleVector areas; areas.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    const Vector &si = vectors[i];
    const Vector &si1 = vectors[next(i)];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n_);

  for (size_t i = 0; i < n_; ++i) {
    size_t i_1 = prev(i), i1 = next(i);
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n_; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector &si_1 = vectors[i_1];
    const Vector &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    double ri_1 = 1.0, ri = 1.0, ri1 = 1.0;
    switch (type) {
    case BarycentricType::WACHSPRESS:
      break;
    case BarycentricType::MEAN_VALUE:
      ri_1 = vectors[i_1].norm();
      ri = vectors[i].norm();
      ri1 = vectors[i1].norm();
      break;
    case BarycentricType::HARMONIC:
      ri_1 = vectors[i_1].sqrnorm();
      ri = vectors[i].sqrnorm();
      ri1 = vectors[i1].sqrnorm();
      break;
    };
    if (ri < epsilon) {         // at a vertex of the domain (mean/harmonic)
      l.assign(n_, 0.0);
      l[i] = 1.0;
      break;
    }
    l.push_back(ri_1 * Ai_1 + ri1 * Ai - ri * Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });
  return l[i];
}

Vector 
SPatch::evaluateAtParam(double u, double v) const
{
  const bool basis_fcn = true;
  
  std::vector<double> bc(n_);
  for(size_t i = 0; i < n_; ++i) {
    bc[i] = getGBC(u, v, i);
  }
  Vector p(0,0,0);
  size_t cp_idx = 0;
  for (const auto &cp : net_) {
    if(!basis_fcn) {
      p += cp.second * multiBernstein(cp.first, bc);
    }
    else {
      if(cp_idx++ == (selected_idx % net_.size())) {
        p += Vector(u, v, multiBernstein(cp.first, bc));
        break;
      }
    }
  }
  return p;
}

size_t
SPatch::binomial(size_t n, size_t k) 
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

size_t 
SPatch::multinomial(const Index &index) 
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
SPatch::multiBernstein(const Index &index, const DoubleVector &bc) 
{
  double result = multinomial(index);
  for (size_t i = 0; i < index.size(); ++i)
    result *= std::pow(bc[i], index[i]);
  return result;
}

std::vector<SPatch::Index>
SPatch::indices(size_t n, size_t d)
{
    if(n == 1) {
      Index temp(1, d);
      std::vector<Index> result(1, temp);
      return result;
    }
    else {
      std::vector<Index> result;
      for(size_t i = 0; i <= d; i++) {
          std::vector<Index> sub_indices = indices(n - 1, d - i);
          for(auto& si : sub_indices) {
              si.insert(si.begin(), i);
              result.push_back(si);
          }
      }
      return result;
    }
}

void
SPatch::setControlPoint(const Index &i, const Vector &p) 
{
  net_[i] = p;
}

// bool
// SPatch::is_neigbour(const Index &ida, const Index &idb) const
// {
//   for(size_t i = 0; i < n_; ++i) {
//     auto sl = shift_idx_l(ida, i);
//     auto sr = shift_idx_r(ida, i);
//     if(std::equal(sl.begin(),sl.end(), idb.begin()) || std::equal(sr.begin(),sr.end(), idb.begin())) {
//       return true;
//     }
//   }
//   return false;
// }

// SPatch::Index 
// SPatch::shift_idx_l(const Index &idx, size_t i) const
// {
//   if(idx[i] > 0) {
//     Index shifted_idx = idx;
//     --shifted_idx[i];
//     ++shifted_idx[prev(i)];
//   }
//   else {
//     return idx;
//   }
// }

// SPatch::Index 
// SPatch::shift_idx_r(const Index &idx, size_t i) const
// {
//   if(idx[i] > 0) {
//     Index shifted_idx = idx;
//     --shifted_idx[i];
//     ++shifted_idx[next(i)];
//   }
//   else {
//     return idx;
//   }
// }

std::vector<SPatch::Index> 
SPatch::neighbors(const Index& si) const
{
  auto step = [&](size_t i, int dir) {
    size_t j = (i + dir + n_) % n_;
    Index result(si);
    if(result[i] > 0) {
      result[i]--;
      result[j]++;
    }
    return result;
  };
  
  std::vector<Index> candidates;
  for(size_t i = 0; i < n_; i++) {
    candidates.push_back(step(i, 1));
    if(si[i] > 0) { // add a check here as well to prevent wrap-around
      candidates.push_back(step(i, -1));
    }
  }
  
  return candidates;
}

void SPatch::draw(const Visualization &vis) const
{
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    for (const auto &cp : net_) {
      auto nbs = neighbors(cp.first);
      for(const auto &nb : nbs) {
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
    for (const auto &cp : net_) {
      if(idx++ == (selected_idx % net_.size())) {
        glPointSize(16.0);
      }
      else{
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
SPatch::drawWithNames(const Visualization &vis) const
{

}

Vector SPatch::postSelection(int selected) 
{
  return {0.0, 0.0, 0.0};
}

void SPatch::movement(int selected, const Vector &pos)
{

}