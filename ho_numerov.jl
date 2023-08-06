using Plots
using NumericalIntegration

# Potential function.
function V(x)
  return x^2/2
end


function numerov(xmin, xmax, nmax, E)
  # Initialize vector.
  x = range(xmin, xmax, nmax)
  steps = (xmax - xmin)/nmax

  y = zeros(nmax)
  y[1] = 0.
  y[2] = 1.

  # Initialize auxiliary arrays.
  f = zeros(nmax)
  g = zeros(nmax)

  g = @. 2*(E - V(x))
  f = @. 1 + g*(steps^2)/12

  mask = g .> 0

  # Numerov algorithm.
  for i = 3:nmax
    y[i] = ((12 - 10*f[i - 1])*y[i - 1] - f[i - 2]*y[i - 2])/f[i]
  end

  if xmax > xmin
    return x[mask], y[mask]
  else
    return reverse(x), reverse(y)
  end
end

function plot_numerov(xmin, xmax, nmax, E)
  # Solve the oscillator.
  x, y = numerov(xmin, xmax, nmax, E)

  # Plot the solution.
  plot(x, y, xlabel = "x", ylabel = "p", linewidth = 2, title = "Numerical solution of the quantum harmonic oscillator")
end


# Set the parameters.
E = 1.5
nmax = 1000
xmax = 8.
xmin = -8.


# Solve the oscillator.
xl, yl = numerov(xmin, xmax, nmax, E)
xr, yr = numerov(xmax, xl[length(xl)], nmax, E)
yr = yl[length(yl)]/yr[1] * yr



# Plot the solution.
pl = plot(xl, yl)
pr = plot(xr, yr)

x = vcat(xl, xr)
y = vcat(yl, yr)

y2 = y.^2
area = sqrt(integrate(x, y2))
y = y/area
p = plot(x, y)