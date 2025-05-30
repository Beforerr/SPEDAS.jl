# Interactive tplot with Speasy
# Visual exploration of OMNI data

## tplot with Speasy product ID strings

```@example interactive_speasy
using Speasy
using SPEDAS
using Dates
using CairoMakie

t0 = DateTime("2008-09-05T10:00:00")
t1 = DateTime("2008-09-05T22:00:00")
tvars = [
    "cda/OMNI_HRO_1MIN/flow_speed",
    "cda/OMNI_HRO_1MIN/E",
    "cda/OMNI_HRO_1MIN/Pressure"
]
f, axes = tplot(tvars, t0, t1)
```

## Interactive tplot

Here we simulate a user interacting with the plot by progressively zooming out in time with `tlims!`.
Note: For real-time interactivity, consider using the `GLMakie` backend instead of `CairoMakie`.

```@example interactive_speasy
dt = Hour(12)

record(f, "interactive_speasy.mp4", 1:5; framerate=1) do n
    tlims!(t0 - n * dt, t1 + n * dt)
    sleep(1)
end
```

```@raw html
<video autoplay loop muted playsinline controls src="../interactive_speasy.mp4" />
```