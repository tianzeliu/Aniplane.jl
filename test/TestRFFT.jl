using FFTW

dt = 0.01
npts = 100
fs = 1/dt

freqax = fftfreq(npts, fs)
freqax_pos = freqax[1:div(npts, 2)+1]

println(length(freqax))
println(length(freqax_pos))