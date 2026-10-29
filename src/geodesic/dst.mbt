// MoonBit port of GeographicLib::DST (Charles Karney, MIT/X11 License).
// The area integral in GeodesicExact is a discrete sine transform of the
// I4 integrand. Transcription of DST.cpp from GeographicLib 2.5. The real
// FFT is a direct DFT: N for double precision is small except at |n| → 1.

///|
let c4_narr : Array[Int32] = [
  22, 21, 19, 18, 17, 17, 16, 15, 15, 15, 14, 14, 14, 13, 13, 13, 13, 13, 13, 12,
  12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10,
  10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5,
  5, 5, 5, 5, 4, 4, 3, 2, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10,
  10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13,
  13, 13, 14, 14, 15, 15, 16, 17, 18, 19,
]

///|
fn c4_nterms(n_ellip : Double) -> Int32 {
  let scaled = 100.0 * n_ellip
  let rounded = if scaled < 0.0 { scaled.floor() } else { scaled.ceil() }
  let mut j = 100 + rounded.to_int()
  if j < 0 {
    j = 0
  }
  if j > 200 {
    j = 200
  }
  let code = c4_narr[j]
  let base = if code % 2 == 0 { 2 } else { 3 }
  let p = code / 2
  base * (1 << p)
}

///|
fn complex_dft(
  re_in : Array[Double],
  im_in : Array[Double],
) -> (Array[Double], Array[Double]) {
  let m = re_in.length()
  let re = Array::make(m, 0.0)
  let im = Array::make(m, 0.0)
  let md = m.to_double()
  for k in 0..<m {
    let mut sr = 0.0
    let mut si = 0.0
    let kd = k.to_double()
    for n in 0..<m {
      let ang = -2.0 * @geomath.pi * kd * n.to_double() / md
      let c = @math.cos(ang)
      let s = @math.sin(ang)
      sr = sr + re_in[n] * c - im_in[n] * s
      si = si + re_in[n] * s + im_in[n] * c
    }
    re[k] = sr
    im[k] = si
  }
  (re, im)
}

///|
fn dst_transform(n : Int32, f : (Double) -> Double) -> Array[Double] {
  if n == 0 {
    return []
  }
  let data = Array::make(4 * n, 0.0)
  let d = @geomath.pi / (2.0 * n.to_double())
  for i in 1..<=n {
    data[i] = f(i.to_double() * d)
  }
  data[0] = 0.0
  for i in 1..<n {
    data[n + i] = data[n - i]
  }
  let two_n = 2 * n
  for i in 0..<two_n {
    data[two_n + i] = -data[i]
  }
  let m = two_n
  let re_in = Array::make(m, 0.0)
  let im_in = Array::make(m, 0.0)
  for k in 0..<m {
    re_in[k] = data[2 * k]
    im_in[k] = data[2 * k + 1]
  }
  let (cre, cim) = complex_dft(re_in, im_in)
  let r0 = cre[0]
  let i0 = cim[0]
  cre[0] = r0 + i0
  cim[0] = r0 - i0
  let mut k = 1
  while 2 * k < m {
    let mk = m - k
    let w_re = 0.5 * (cre[k] + cre[mk])
    let w_im = 0.5 * (cim[k] - cim[mk])
    let z_re = 0.5 * (cim[k] + cim[mk])
    let z_im = 0.5 * (-cre[k] + cre[mk])
    let ang = -@geomath.pi * k.to_double() / m.to_double()
    let t_re = @math.cos(ang)
    let t_im = @math.sin(ang)
    let wz_re = t_re * z_re - t_im * z_im
    let wz_im = t_re * z_im + t_im * z_re
    let new_k_re = w_re + wz_re
    let new_k_im = w_im + wz_im
    let new_mk_re = w_re - wz_re
    let new_mk_im = -(w_im - wz_im)
    cre[k] = new_k_re
    cim[k] = new_k_im
    cre[mk] = new_mk_re
    cim[mk] = new_mk_im
    k = k + 1
  }
  if m % 2 == 0 {
    cim[m / 2] = -cim[m / 2]
  }
  let out = Array::make(n, 0.0)
  let scale = two_n.to_double()
  let mut i = 0
  let mut j = 1
  while i < n {
    out[i] = -cim[j] / scale
    i = i + 1
    j = j + 2
  }
  out
}

///|
fn dst_integral(sinx : Double, cosx : Double, coeff : Array[Double]) -> Double {
  let ar = 2.0 * (cosx - sinx) * (cosx + sinx)
  let mut y0 = 0.0
  let mut y1 = 0.0
  let mut n = coeff.length() - 1
  while n >= 0 {
    let t = ar * y0 - y1 + coeff[n] / (2.0 * n.to_double() + 1.0)
    y1 = y0
    y0 = t
    n = n - 1
  }
  cosx * (y1 - y0)
}

///|
fn dst_integral_diff(
  sinx : Double,
  cosx : Double,
  siny : Double,
  cosy : Double,
  coeff : Array[Double],
) -> Double {
  let ac = 2.0 * (cosy * cosx + siny * sinx) * (cosy * cosx - siny * sinx)
  let as_ = -2.0 * (siny * cosx - cosy * sinx) * (siny * cosx + cosy * sinx)
  let mut y0 = 0.0
  let mut y1 = 0.0
  let mut z0 = 0.0
  let mut z1 = 0.0
  let mut n = coeff.length() - 1
  while n >= 0 {
    let ty = ac * y0 + as_ * z0 - y1 + coeff[n] / (2.0 * n.to_double() + 1.0)
    let tz = as_ * y0 + ac * z0 - z1
    y1 = y0
    y0 = ty
    z1 = z0
    z0 = tz
    n = n - 1
  }
  (y1 - y0) * (cosy - cosx) + (z1 - z0) * (cosy + cosx)
}
