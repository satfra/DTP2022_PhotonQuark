// Minimal reproducer for an nvc++ OpenACC device-codegen bug that bit the GPU
// port of the K' precalculation (see Kernels_K.hh::k_component).
//
// Symptom: a `#pragma acc routine seq` member function whose body is a large
// switch (>~20 sparse cases) dispatching to *private member methods* silently
// falls through and returns 0 on the GPU. The identical class is correct on the
// host, and the identical arithmetic inlined into the switch cases is correct on
// the GPU. Nothing is reported at compile time — the kernel just yields zeros.
//
// This is why Kernels_K.hh defines the kernel components as a free function with
// the formulas inlined in the switch (k_component) rather than as class methods,
// and why K::get() delegates to it.
//
// Bisection on the real K class put the threshold at 20 -> 21 switch cases.
// A simpler variant (cases calling a one-argument *free* function) does NOT
// reproduce; the member-method form does.
//
// Build & run:
//   nvc++ -O2 -acc=gpu -gpu=cc89 tools/nvc_openacc_switch_bug.cpp -o bug && ./bug
//   nvc++ -O2                    tools/nvc_openacc_switch_bug.cpp -o ok  && ./ok
// Expected: host prints -0.872; GPU prints 0 (bug).
//
// Observed with nvc++ 26.5 on an NVIDIA RTX 4070 Laptop (cc89), CUDA 13.3.
#include <cstdio>

class Kx {
  double y, l2, u, up, V, w, wp, X;
 public:
  #pragma acc routine seq
  Kx(double a,double b,double c,double d,double e,double f,double g,double h)
    : y(a), l2(b), u(c), up(d), V(e), w(f), wp(g), X(h) {}

  #pragma acc routine seq
  double get(unsigned i, unsigned j) const {
    const unsigned key = 100u*(i+1u) + (j+1u);
    switch (key) {
      case 101: return m11(y,l2,u,up,V,w,wp,X);  case 202: return m22(y,l2,u,up,V,w,wp,X);
      case 303: return m33(y,l2,u,up,V,w,wp,X);  case 404: return m44(y,l2,u,up,V,w,wp,X);
      case 505: return m55(y,l2,u,up,V,w,wp,X);  case 606: return m66(y,l2,u,up,V,w,wp,X);
      case 707: return m77(y,l2,u,up,V,w,wp,X);  case 808: return m88(y,l2,u,up,V,w,wp,X);
      case 106: return m16(y,l2,u,up,V,w,wp,X);  case 601: return m61(y,l2,u,up,V,w,wp,X);
      case 107: return m17(y,l2,u,up,V,w,wp,X);  case 701: return m71(y,l2,u,up,V,w,wp,X);
      case 607: return m67(y,l2,u,up,V,w,wp,X);  case 706: return m76(y,l2,u,up,V,w,wp,X);
      case 203: return m23(y,l2,u,up,V,w,wp,X);  case 302: return m32(y,l2,u,up,V,w,wp,X);
      case 208: return m28(y,l2,u,up,V,w,wp,X);  case 802: return m82(y,l2,u,up,V,w,wp,X);
      case 308: return m38(y,l2,u,up,V,w,wp,X);  case 803: return m83(y,l2,u,up,V,w,wp,X);
      case 909: return m99(y,l2,u,up,V,w,wp,X);  case 1010: return mA(y,l2,u,up,V,w,wp,X);
      case 1011: return mB(y,l2,u,up,V,w,wp,X);  case 1110: return mC(y,l2,u,up,V,w,wp,X);
      case 1111: return mD(y,l2,u,up,V,w,wp,X);  case 1212: return mE(y,l2,u,up,V,w,wp,X);
    }
    return 0.0;
  }
 private:
  #pragma acc routine seq
  double m11(double y,double l2,double u,double up,double V,double w,double wp,double X) const { return -(1.+y*y)/2. - y*(1.-y*y)*X; }
  #define M(name,expr) _Pragma("acc routine seq") double name(double y,double l2,double u,double up,double V,double w,double wp,double X) const { return expr; }
  M(m22, -(1.+y*y)/2.*(1.-2.*l2*V*V)+y*(1.-y*y)*X)
  M(m33,  y*(1.-2.*l2*V*V)-(1.-y*y)*X)
  M(m44,  y+(1.-y*y)*X)      M(m55, 3.*y)             M(m66, -y*(1.+2.*l2*V*V))
  M(m77, -y*y*(3.-2.*l2*V*V)+2.*y*(1.-y*y)*X)         M(m88, y*y-2.*y*(1.-y*y)*X)
  M(m16,  1.4142135623730951*(1.-y*y)*up*V)          M(m61, -1.4142135623730951*(1.-y*y)*u*V)
  M(m17, -(1.-y*y)/1.4142135623730951*(1.+2.*wp-2.*y*X))
  M(m71, -(1.-y*y)/1.4142135623730951*(1.+2.*w-2.*y*X))
  M(m67, 2.*y*(up-y*u)*V)    M(m76, -2.*y*(u-y*up)*V)
  M(m23, (2.*y*u-(1.+y*y)*up)*V)   M(m32, -(2.*y*up-(1.+y*y)*u)*V)
  M(m28, -(1.-y*y)/1.4142135623730951*(1.+2.*w-2.*y*X)+1.4142135623730951*(1.-y*y))
  M(m82, -(1.-y*y)/1.4142135623730951*(1.+2.*wp-2.*y*X)+1.4142135623730951*(1.-y*y))
  M(m38, 1.4142135623730951*(1.-y*y)*u*V)   M(m83, -1.4142135623730951*(1.-y*y)*up*V)
  M(m99, 3.)   M(mA, -(1.+2.*l2*V*V))   M(mB, 2.*(up-y*u)*V)   M(mC, -2.*(u-y*up)*V)
  M(mD, -y*(3.-2.*l2*V*V)+2.*(1.-y*y)*X)   M(mE, y-2.*(1.-y*y)*X)
  #undef M
};

int main() {
  const int N = 256;
  double out[N];
  #pragma acc parallel loop copyout(out[0:N])
  for (int t = 0; t < N; ++t) {
    Kx k(0.6, 2.0, 0.3, 0.4, 0.5, 0.1, 0.2, 0.5);
    out[t] = k.get(0u, 0u);   // -> m11 = -(1+0.36)/2 - 0.6*0.64*0.5 = -0.872
  }
  printf("Kx::get(0,0) = %g  (expected -0.872)\n", out[0]);
  const bool bug = (out[0] == 0.0);
  printf("%s\n", bug ? "BUG REPRODUCED on device" : "no bug (correct)");
  return bug ? 1 : 0;
}
