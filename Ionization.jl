#!/usr/bin/env julia

omega = 0.0228
amp = 0.0655
n = 36
ipot = 0.579
ac = 0.001

field(t) = amp*cos(omega*t/(2n))^2*sin(omega*t) + 
  (amp/n)*cos(omega*t)*cos(omega*t/(2*n))*sin(omega*t/(2*n))

fcx(x, y) = -x/(x^2 + y^2 + ac^2)^1.5
fcy(x, y) = -y/(x^2 + y^2 + ac^2)^1.5

t_start = -pi/(2*omega)
t_end = 3*pi/(2*omega)
v0_min = 1.0
v0_max = -1.0

function generate_electron()
  while true
    t0 = (t_end - t_start)*rand() + t_start
    v0 = (v0_max - v0_min)*rand() + v0_min
    xi = rand()
    supremum = 0.01 #00113 #0.000112430432247
    if xi <= tunneling_density(t0, v0)/supremum
      return [t0, v0]
    end
  end
end

function tunneling_density(t0, v0_signed)
  f = abs(field(t0))
  F_a = (2*ipot)^1.5
  W = 8*(2*ipot)^0.5*F_a/f*exp(-2/3.0*F_a/f)
  a = (f/F_a)^0.5*(2*ipot)^0.5
  w_perp = exp(-1*v0_signed^2/a^2)/(pi*a^2)
  W*w_perp*pi
end

function newton_eqs(t, rp)
  px, py, x, y = rp
  [
    fcx(x, y) - field(t),
    fcy(x, y),
    px,
    py,
  ]
end

function momentum(v0, t0)
  times = collect(linspace(t0, pi*n/omega, 1000))
  res = ODE.ode45(newton_eqs, [0.0, v0, -ipot/field(t0), 0.0], times)
  res[2][end]
end

function generate(n)
  particles = Array{Float64}(n, 4)
  for i = 1:n
    t0, v0 = generate_electron()
    particles[i, :] = momentum(v0, t0)
  end
  fname = "results" * string(rand()) * ".h5"
  HDF5.h5write(fname, "particles", particles)
  print("Written $n rows to $fname\n")
end

